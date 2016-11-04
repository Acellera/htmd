""""""
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import abc
from abc import ABCMeta
from htmd.molecule.molecule import Molecule
from htmd.metricdata import MetricData
from scipy import stats
from joblib import Parallel, delayed
import joblib.parallel
from collections import defaultdict
from htmd.progress.progress import ProgressBar
import time
import logging
logger = logging.getLogger(__name__)

def progessbar(seq, description=None, total=None):
    p = ProgressBar(total, description=description)
    while True:
        try:
            yield next(seq)
            p.progress() # Had to put progress after yield because last call goes over the total and then I can't decrement in stop()
        except StopIteration:
            p.stop()
            raise

def ParallelExecutor(use_bar='txt', **joblib_args):
    def aprun(bar=use_bar, **tq_args):
        def tmp(op_iter):
            foo = lambda args: lambda x: progessbar(x, **args)
            bar_func = foo(tq_args)
            return Parallel(**joblib_args)(bar_func(op_iter))
        return tmp
    return aprun


class Metric:
    """ Class for calculating projections of a simlist.

    Parameters
    ----------
    simulations : list
        A list of simulations produced by :func:`simlist <htmd.simlist.simlist>`
    skip : int
        Frame skipping. Setting i.e. to 3 will keep only every third frame of each simulation.

    Examples
    --------
    >>> metr = Metric(sims)  # doctest: +SKIP
    >>> metr.projection(MetricSelfDistance('protein and name CA', metric='contacts'))  # doctest: +SKIP
    >>> data = metr.project()  # doctest: +SKIP

    .. currentmodule:: htmd.projections.metric.Metric
    .. rubric:: Methods
    .. autoautosummary:: htmd.projections.metric.Metric
        :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.projections.metric.Metric
        :attributes:
    """
    def __init__(self, simulations, skip=1):
        self.simulations = simulations
        self.skip = skip
        self.projectionlist = []

    def projection(self, metric):
        """ Deprecated
        """
        self.projectionlist.append(metric)

    def set(self, projection):
        """ Sets the projection to be applied to the simulations.

        Parameters
        ----------
        projection : :class:`Projection <htmd.projections.projection.Projection>` object or list of objects
            A projection or a list of projections which to use on the simulations
        """
        self.projectionlist = projection
        if not isinstance(self.projectionlist, list) and not isinstance(self.projectionlist, tuple):
            self.projectionlist = [self.projectionlist]

    def getMapping(self, mol):
        """ Returns the description of each projected dimension.

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object which will be used to calculate the descriptions of the projected dimensions.

        Returns
        -------
        map : :class:`DataFrame <pandas.core.frame.DataFrame>` object
            A DataFrame containing the descriptions of each dimension
        """
        import pandas as pd
        if mol is None:
            return
        pandamap = pd.DataFrame(columns=('type', 'atomIndexes', 'description'))
        for proj in self.projectionlist:
            pandamap = pandamap.append(proj.getMapping(mol), ignore_index=True)
        return pandamap

    def project(self):
        """
        Applies all projections stored in Metric on all simulations.

        Returns
        -------
        data : MetricData object
               Returns a MetricData object containing the projected data.
        """
        if len(self.projectionlist) == 0:
            raise NameError('You need to provide projections using the Metric.projection method.')

        if isinstance(self.simulations, Molecule):
            data = []
            for proj in self.projectionlist:
                data.append(proj.project(self.simulations))
            return data

        numSim = len(self.simulations)

        # Find out if there is a unique molfile. If there is, initialize a single Molecule to speed up calculations
        uqMol = None
        import pandas as pd
        map = pd.DataFrame(columns=('type', 'atomIndexes', 'description'))
        (single, molfile) = _singleMolfile(self.simulations)
        if single:
            uqMol = Molecule(molfile)
            for proj in self.projectionlist:
                proj._precalculate(uqMol)
        else:
            logger.warning('Cannot calculate description of dimensions due to different topology files for each trajectory.')
        map = self.getMapping(uqMol)

        logger.info('Metric: Starting projection of trajectories.')
        metrics = np.empty(numSim, dtype=object)
        ref = np.empty(numSim, dtype=object)
        deletesims = np.zeros(numSim, dtype=bool)
        fstep = np.zeros(numSim)

        from htmd.config import _config
        #results = Parallel(n_jobs=_config['ncpus'], verbose=6)(
        #        delayed(_processSim)(self.simulations[i], self.projectionlist, uqMol, self.skip) for i in range(numSim))

        aprun = ParallelExecutor(n_jobs=_config['ncpus'])
        results = aprun(total=numSim, description='Projecting trajectories')(delayed(_processSim)(self.simulations[i], self.projectionlist, uqMol, self.skip) for i in range(numSim))

        for i in range(len(results)):
            metrics[i] = results[i][0]
            ref[i] = results[i][1]
            fstep[i] = results[i][2]
            deletesims[i] = results[i][3]

        logger.info('Finished projecting the trajectories.')

        # Removing empty trajectories
        emptyM = np.array([True if x is None else False for x in metrics], dtype=bool)
        emptyR = np.array([True if x is None else False for x in ref], dtype=bool)
        assert np.all(deletesims == emptyM) and np.all(emptyR == emptyM)

        metrics = np.delete(metrics, np.where(emptyM)[0])
        ref = np.delete(ref, np.where(emptyM)[0])
        updlist = np.delete(self.simulations, np.where(emptyM)[0])

        if len(metrics) == 0:
            raise NameError('No trajectories were read')

        # Constructing a MetricData object
        data = MetricData(dat=metrics, ref=ref, map=map, simlist=updlist)

        uqfsteps = np.unique(fstep)
        data.fstep = float(stats.mode(fstep).mode)
        if len(uqfsteps) != 1:
            logger.warning('Multiple framesteps were read from the simulations. '
                           'Taking the statistical mode: ' + str(data.fstep) + 'ns. '
                           'If it looks wrong, you can modify it by manually setting the MetricData.fstep property.')
        else:
            logger.info('Frame step {}ns was read from the trajectories. If it looks wrong, redefine it by manually '
                        'setting the MetricData.fstep property.'.format(data.fstep))

        return data

    def _projectSingle(self, index):
        data, ref, fstep, _ = _processSim(self.simulations[index], self.projectionlist, None, self.skip)
        return data, ref, fstep


def _highfreqFilter(mol,steps):
    newframes = int(mol.coords.shape[2]/steps)*steps
    mol.coords =   mol.coords[:,:,:newframes]
    mol.box = mol.box[:,:newframes]
    mol.box = mol.box[:,::steps]
    mol.fstep = mol.fstep*steps
    coo=np.reshape(mol.coords,[mol.coords.shape[0],mol.coords.shape[1],int(mol.coords.shape[2]/steps),steps])
    X=np.mean(coo,axis=3)
#   print("mol: ", mol.coords.shape, " X: ",X.shape)
    mol.coords=X

def _processSim(sim, projectionlist, uqmol, skip):
    pieces = sim.trajectory
    try:
        if uqmol is not None:
            mol = uqmol.copy()
        else:
            mol = Molecule(sim.molfile)
        logger.debug(pieces[0])
       
       
        mol.read(pieces, skip=skip)
        #Gianni testing
        #_highfreqFilter(mol,10)
 
        data = []
        for p in projectionlist:
            pj = p.project(mol)
            if pj.ndim == 1:
                pj = np.atleast_2d(pj).T
            data.append(pj)
        data = np.hstack(data)
        if data.dtype == np.float64:
            data = data.astype(np.float32)
    except Exception as e:
        logger.warning('Error in simulation with id: ' + str(sim.simid) + '. "' + e.__str__() + '"')
        return None, None, None, True

    return data, _calcRef(pieces, mol.fileloc), mol.fstep, False


def _calcRef(pieces, fileloc):
    locs = np.array(list([x[0] for x in fileloc]))
    frames = list([x[1] for x in fileloc])
    ref = np.zeros((len(frames), 2), dtype='u4')
    ref[:, 1] = frames
    for i, p in enumerate(pieces):
        ref[locs == p, 0] = i
    return ref


def _singleMolfile(sims):
    single = False
    molfile = []
    if isinstance(sims, Molecule):
        single = False
    elif isinstance(sims, np.ndarray) and len(set([x.molfile for x in sims])) == 1:
        single = True
        molfile = sims[0].molfile
    return single, molfile


def _projectionGenerator(metric, ncpus):
    for i in range(0, len(metric.simulations), ncpus):
        simrange = range(i, np.min((i+ncpus, len(metric.simulations))))
        results = Parallel(n_jobs=ncpus, verbose=0)(delayed(_projector)(metric, i) for i in simrange)
        yield results


def _projector(metric, i):
    return metric._projectSingle(i)


class MetricPyemma(metaclass=ABCMeta):
    """
    Parent class for all trajectory projecting classes. Implements metrify() and defines abstract functions.
    """

    @staticmethod
    @abc.abstractmethod
    def project(self):
        """ Subclasses need to implement and overload this method """
        return

    @abc.abstractmethod
    def _processTraj(self, traj):
        """ Subclasses need to implement this method """
        return

    def _getMapping(self, mol):
        return []

    def _metrify(self, sims, skip, verbose, update):
        """
        Takes a set of trajectory folders and projects all trajectories within them onto the given space defined by the Metric* class.

        Parameters
        ----------

        simList : numpy list of structs
              A list of structs produced by the simList function.

        skip : int
               Skips every x frames.

        verbose : int
              Verbosity toggle

        update : MetricData object
             Provide a previous MetricData object and only metrify new trajectories.

        Returns
        -------

        data : MetricData object
               Returns a MetricData object containing the projected data and the ref data.

        """

        if isinstance(sims, Molecule):
            return self.processTraj(sims)

        # [updList, oldList] = checkUpdate(simList, update, verbose);
        updList = sims
        numSim = len(updList)

        # Find out if there is a unique molfile. If there is, initialize a single Molecule to speed up calculations
        uniqueMol = 0
        uqMol = []
        map = []
        (single, molfile) = _singleMolfile(updList)
        if single:
            uniqueMol = 1
            uqMol = Molecule(molfile)
            # Calculating the mapping of metric columns to atom pair indeces
            map = self._getMapping(uqMol)

        logger.info('Metric: Starting projection of trajectories.')
        metrics = np.empty(numSim, dtype=object)
        ref = np.empty(numSim, dtype=object)
        deleteSims = np.zeros(numSim, dtype=bool)
        fstep = np.empty(numSim)

        #global parpool
        Parallel(n_jobs=6, backend="threading")(delayed(_processSimPyemma)(self, i, updList, uniqueMol, uqMol, skip, deleteSims, metrics, ref, fstep) for i in range(numSim))

        logger.info('Finished projecting the trajectories.')

        # Removing empty trajectories
        emptyM = [True if np.size(x) == 0 else False for x in metrics]
        emptyR = [True if np.size(x) == 0 else False for x in ref]
        #assert np.all(deleteSims == emptyM)# and np.all(emptyR == emptyM)

        metrics = np.delete(metrics, np.where(emptyM))
        ref = np.delete(ref, np.where(emptyM))
        #updList = np.delete(updList, emptyM)

        if len(metrics) == 0:
            raise NameError('No trajectories were read')

        # Constructing a MetricData object
        if not update:
            data = MetricData(dat=metrics, ref=ref, map=map, simlist=updList)
        else:
            data = update
            data.dat.extend(metrics)
            data.ref.extend(ref)
            data.simList.extend(updList)

        uqfsteps = np.unique(fstep)
        data.fstep = stats.mode(fstep).mode
        if len(uqfsteps) != 1:
            logger.warning('Multiple framesteps were read from the simulations. Taking the statistical mode: ' + str(data.fstep) + 'ns.')
            logger.warning('If it looks wrong, you can modify it by manually setting the MetricData.fstep property.')

        return data

    def _calcRef(self, lengths):
        ref = np.zeros((np.sum(lengths).astype(int), 2), dtype='u4')
        prevL = 0
        for i in range(len(lengths)):
            currL = prevL + lengths[i]
            ref[prevL:currL, 0] = i
            ref[prevL:currL, 1] = range(lengths[i])
            prevL = currL
        return ref


def _processSimPyemma(obj, i, updList, uniqueMol, uqMol, skip, deleteSims, metrics, ref, fstep):
    pieces = updList[i].trajectory
    '''
    try:
        if uniqueMol:
            mol = uqMol.copy()
        else:
            mol = Molecule(updList[i].molfile)
        logging.debug(pieces[0])
        mol.readTraj(pieces, skip=skip)
    except Exception as e:
        logging.warning('Error in simulation with id: ' + str(updList[i].id) + ' ' + e.__str__())
        deleteSims[i] = True
        return
    '''

    #fstep[i] = mol.fstep[1] / 1E6
    (metrics[i], trajLengths) = obj._processTraj(pieces)
    ref[i] = obj._calcRef(trajLengths)


