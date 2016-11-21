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
from htmd.parallelprogress import ParallelExecutor
import logging
logger = logging.getLogger(__name__)



class Metric:
    """ Class for calculating projections of a simlist.

    Parameters
    ----------
    simulations : list
        A list of simulations produced by :func:`simlist <htmd.simlist.simlist>`
    skip : int
        Frame skipping. Setting i.e. to 3 will keep only every third frame of each simulation.
    metricdata : :class:`MetricData <htmd.metricdata.MetricData>` object
        If a MetricData object is passed in the constructor, Metric will try to update it by only adding simulations
        which don't exist in it yet.

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
    def __init__(self, simulations, skip=1, metricdata=None):
        self.simulations = simulations
        self.skip = skip
        self.projectionlist = []
        self.metricdata = metricdata

    def projection(self, metric):
        """ Deprecated
        """
        logger.warning('Projection method is deprecated. Use the .set() method of Metric instead.')
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
        (single, molfile) = _singleMolfile(self.simulations)
        if single:
            uqMol = Molecule(molfile)
            for proj in self.projectionlist:
                proj._precalculate(uqMol)
        else:
            logger.warning('Cannot calculate description of dimensions due to different topology files for each trajectory.')
        mapping = self.getMapping(uqMol)

        logger.debug('Metric: Starting projection of trajectories.')
        from htmd.config import _config
        aprun = ParallelExecutor(n_jobs=_config['ncpus'])
        results = aprun(total=numSim, description='Projecting trajectories')(delayed(_processSim)(self.simulations[i], self.projectionlist, uqMol, self.skip) for i in range(numSim))

        metrics = np.empty(numSim, dtype=object)
        ref = np.empty(numSim, dtype=object)
        deletesims = np.zeros(numSim, dtype=bool)
        fstep = np.zeros(numSim)
        for i in range(len(results)):
            metrics[i] = results[i][0]
            ref[i] = results[i][1]
            fstep[i] = results[i][2]
            deletesims[i] = results[i][3]

        logger.debug('Finished projecting the trajectories.')

        # Removing empty trajectories
        metrics, ref, updlist, fstep = self._removeEmpty(metrics, ref, deletesims, fstep)

        # Constructing a MetricData object
        data = MetricData(dat=metrics, ref=ref, map=mapping, simlist=updlist)

        uqfsteps = np.unique(fstep)
        data.fstep = float(stats.mode(fstep).mode)
        if len(uqfsteps) != 1:
            logger.warning('Multiple framesteps [{}] ns were read from the simulations. '
                           'Taking the statistical mode: {}ns. '
                           'If it looks wrong, you can modify it by manually '
                           'setting the MetricData.fstep property.'.format(', '.join(map(str,uqfsteps)), data.fstep))
        else:
            logger.info('Frame step {}ns was read from the trajectories. If it looks wrong, redefine it by manually '
                        'setting the MetricData.fstep property.'.format(data.fstep))

        return data

    def _projectSingle(self, index):
        data, ref, fstep, _ = _processSim(self.simulations[index], self.projectionlist, None, self.skip)
        return data, ref, fstep

    def _removeEmpty(self, metrics, ref, deletesims, fstep):
        emptyM = np.array([True if x is None else False for x in metrics], dtype=bool)
        emptyR = np.array([True if x is None else False for x in ref], dtype=bool)
        assert np.all(deletesims == emptyM) and np.all(emptyR == emptyM)

        metrics = np.delete(metrics, np.where(emptyM)[0])
        ref = np.delete(ref, np.where(emptyM)[0])
        updlist = np.delete(self.simulations, np.where(emptyM)[0])
        fstep = np.delete(fstep, np.where(emptyM)[0])

        if len(metrics) == 0:
            raise NameError('No trajectories were projected. Check if the simlist is empty or for projection errors.')

        return metrics, ref, updlist, fstep



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

