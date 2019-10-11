""""""
# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from moleculekit.molecule import Molecule
from htmd.metricdata import MetricData
from scipy import stats
from moleculekit.projections.projection import Projection
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
    >>> metr.set(MetricSelfDistance('protein and name CA', metric='contacts'))  # doctest: +SKIP
    >>> data = metr.project()  # doctest: +SKIP
    >>>
    >>> # Or define your own function which accepts as first argument a Molecule object. Further arguments are passed as
    >>> # function/argument tuples
    >>> def foo(mol, ref):
    >>>     from moleculekit.util import molRMSD
    >>>     mol.wrap('protein')
    >>>     mol.align('protein and name CA', refmol=ref)
    >>>     return molRMSD(mol, ref, mol.atomselect('protein and name CA'), ref.atomselect('protein and name CA'))
    >>>
    >>> metr = Metric(sims)
    >>> metr.set( (foo, (ref,)) )
    >>> data2 = metr.project()

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

    def set(self, projection):
        """ Sets the projection to be applied to the simulations.

        Parameters
        ----------
        projection : function or :class:`Projection <moleculekit.projections.projection.Projection>` object or list of objects
            A function or projection or a list of projections/functions which to use on the simulations
        """
        self.projectionlist = projection
        if not (isinstance(self.projectionlist, list) or isinstance(self.projectionlist, tuple)) \
                and (isinstance(self.projectionlist, Projection) or hasattr(self.projectionlist, '__call__')):
            self.projectionlist = [self.projectionlist, ]
        elif isinstance(self.projectionlist, list) or isinstance(self.projectionlist, tuple):
            if hasattr(self.projectionlist[0], '__call__'):
                self.projectionlist = [self.projectionlist, ]
            pass
        else:
            raise AttributeError('Metric.set only accepts Projection objects, functions, function/argument tuples or '
                                 'lists and tuples thereof.')

    def getMapping(self, mol):
        """ Returns the description of each projected dimension.

        Parameters
        ----------
        mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
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
            if isinstance(proj, Projection):
                pandamap = pandamap.append(proj.getMapping(mol), ignore_index=True)
        return pandamap

    def project(self, njobs=None):
        """
        Applies all projections stored in Metric on all simulations.

        Parameters
        ----------
        njobs : int
            Number of parallel jobs to spawn for projection of trajectories. Take care that this can use large amounts 
            of memory as multiple trajectories are loaded at once.  If None it will use the default from htmd.config.

        Returns
        -------
        data : MetricData object
               Returns a MetricData object containing the projected data.
        """
        if len(self.projectionlist) == 0:
            raise RuntimeError('You need to provide projections using the Metric.set method.')

        # Projecting single Molecules
        if isinstance(self.simulations, Molecule):
            data = []
            mol = self.simulations
            for proj in self.projectionlist:
                data.append(_project(proj, mol))
            return data

        numSim = len(self.simulations)

        # Find out if there is a unique molfile. If there is, initialize a single Molecule to speed up calculations
        uqMol = None
        (single, molfile) = _singleMolfile(self.simulations)
        if single:
            uqMol = Molecule(molfile)
            for proj in self.projectionlist:
                if isinstance(proj, Projection):
                    proj._setCache(uqMol)
        else:
            logger.warning('Cannot calculate description of dimensions due to different topology files for each trajectory.')
        mapping = self.getMapping(uqMol)

        logger.debug('Metric: Starting projection of trajectories.')
        from htmd.config import _config
        aprun = ParallelExecutor(n_jobs=njobs if njobs is not None else _config['njobs'])
        results = aprun(total=numSim, desc='Projecting trajectories')(delayed(_processSim)(self.simulations[i], self.projectionlist, uqMol, self.skip) for i in range(numSim))

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
        data = MetricData(dat=metrics, ref=ref, description=mapping, simlist=updlist)

        uqfsteps = np.unique(fstep)
        data.fstep = float(stats.mode(fstep).mode)
        if len(uqfsteps) != 1:
            logger.warning('Multiple framesteps [{}] ns were read from the simulations. '
                           'Taking the statistical mode: {}ns. '
                           'If it looks wrong, you can modify it by manually '
                           'setting the MetricData.fstep property.'.format(', '.join(map(str,uqfsteps)), data.fstep))
        else:
            if data.fstep == 0:
                logger.warning('A fstep of 0 was read from the trajectories. Please manually set the MetricData.fstep'
                               ' property, otherwise calculations in Model and Kinetics classes can fail.')
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


def _project(proj, target):
    if isinstance(proj, Projection):
        res = proj.project(target)
    elif hasattr(proj, '__call__'): # If it's a function
        res = proj(target)
    elif isinstance(proj, tuple) and hasattr(proj[0], '__call__'): # If it's a function with extra arguments
        res = proj[0](target, *proj[1])
    else:
        raise RuntimeError('Invalid projection type {}'.format(type(proj)))

    if isinstance(res, list) or isinstance(res, tuple):
        res = np.array(res)

    _checkProjectionDims(res, target, 'projection')
    return res


def _checkProjectionDims(result, mol, name):
    if (result.ndim == 1 and len(result) != mol.numFrames) or \
       (result.ndim == 2 and result.shape[0] != mol.numFrames):
        raise RuntimeError('The {} produced an incorrect result. It produced a {} shaped array for {} frames. '
                           'The {} should return a numpy array of (nframes, ndim) shape '
                           'where nframes the number of frames in the Molecule it accepts as an '
                           'argument.'.format(name, result.shape, mol.numFrames, name))


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
 
        data = []
        for proj in projectionlist:
            result = _project(proj, mol)
            if result.ndim == 1:
                result = np.atleast_2d(result).T
            if result.size == 0:
                logger.warning(f'No data was produced by projection {proj.__class__} of simulation id: {sim.simid}')
            data.append(result)

        data = np.hstack(data)
        if data.dtype == np.float64:
            data = data.astype(np.float32)

    except Exception as e:
        logger.warning(f'Error while projecting simulation id: {sim.simid}. "{e}"')
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
    from moleculekit.molecule import mol_equal
    from htmd.util import ensurelist
    if isinstance(sims, Molecule):
        return False, []
    elif isinstance(sims, np.ndarray):
        molfiles = []
        for s in sims:
            molfiles.append(tuple(ensurelist(s.molfile)))

        uqmolfiles = list(set(molfiles))

        if len(uqmolfiles) == 0:
            raise RuntimeError('No molfiles found in simlist')
        elif len(uqmolfiles) == 1:
            return True, uqmolfiles[0]
        elif len(uqmolfiles) > 1:  # If more than one molfile load them and see if they are different Molecules
            ref = Molecule(uqmolfiles[0], _logger=False)
            for i in range(1, len(uqmolfiles)):
                mol = Molecule(uqmolfiles[i], _logger=False)
                if not mol_equal(ref, mol, exceptFields=['coords']):
                    return False, []
            return True, uqmolfiles[0]
    return False, []


def _projectionGenerator(metric, njobs):
    for i in range(0, len(metric.simulations), njobs):
        simrange = range(i, np.min((i+njobs, len(metric.simulations))))
        results = Parallel(n_jobs=njobs, verbose=0)(delayed(_projector)(metric, i) for i in simrange)
        yield results


def _projector(metric, i):
    return metric._projectSingle(i)


import unittest
class _TestMetric(unittest.TestCase):
    def test_set(self):
        from moleculekit.molecule import Molecule
        from moleculekit.projections.metricrmsd import MetricRmsd

        # Testing the set method of Metric
        sims = []
        ref = Molecule('3PTB')
        metr = Metric(sims)
        metr.set(MetricRmsd(ref, 'protein and name CA'))
        assert len(metr.projectionlist) == 1

        metr.set([MetricRmsd(ref, 'protein and name CA'), MetricRmsd(ref, 'protein and name CA')])
        assert len(metr.projectionlist) == 2

    def test_function_projections(self):
        def foo(mol, ref):
            from moleculekit.util import molRMSD
            mol.wrap('protein')
            mol.align('protein and name CA', refmol=ref)
            return molRMSD(mol, ref, mol.atomselect('protein and name CA'), ref.atomselect('protein and name CA'))

        ref = Molecule('3PTB')
        sims = []
        metr = Metric(sims)
        metr.set((foo, (ref,)))
        assert len(metr.projectionlist) == 1


if __name__ == '__main__':
    unittest.main(verbosity=2)



