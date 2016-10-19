# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from glob import glob
from os import path
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from htmd.util import _getNcpus
from htmd.adaptive.adaptive import AdaptiveBase
from htmd.simlist import simlist, simfilter
from htmd.model import Model, macroAccumulate
from htmd.projections.tica import TICA
from htmd.projections.metric import Metric, _projectionGenerator
from htmd.protocols.protocolinterface import TYPE_INT, RANGE_0POS, RANGE_POS, TYPE_FLOAT, RANGE_ANY
import logging

logger = logging.getLogger(__name__)


class AdaptiveGoal(AdaptiveBase):
    """ Adaptive class which uses a Markov state model for respawning

    AdaptiveMD uses Markov state models to choose respawning poses for the next epochs. In more detail, it projects all
    currently retrieved simulations according to the specified projection, clusters those and then builds a Markov model using
    the discretized trajectories. From the Markov model it then chooses conformations from the various states based on
    the chosen criteria which will be used for starting new simulations.

    Parameters
    ----------
    app : :class:`App <htmd.apps.app.App>` object, default=None
        An App class object used to retrieve and submit simulations
    project : str, default='adaptive'
        The name of the project
    nmin : int, default=1
        Minimum number of running simulations
    nmax : int, default=1
        Maximum number of running simulations
    nepochs : int, default=100
        Maximum number of epochs
    inputpath : str, default='input'
        The directory used to store input folders
    generatorspath : str, default='generators'
        The directory containing the generators
    dryrun : boolean, default=False
        A dry run means that the adaptive will retrieve and generate a new epoch but not submit the simulations
    updateperiod : float, default=0
        When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds
    coorname : str, default='input.coor'
        Name of the file containing the starting coordinates for the new simulations
    datapath : str, default='data'
        The directory in which the completed simulations are stored
    filter : bool, default=True
        Enable or disable filtering of trajectories.
    filtersel : str, default='not water'
        Filtering atom selection
    filteredpath : str, default='filtered'
        The directory in which the filtered simulations will be stored
    projection : :class:`Projection <htmd.projections.projection.Projection>` object, default=None
        A Projection class object or a list of objects which will be used to project the simulation data before constructing a Markov model
    goalfunction : function, default=None
        This function will be used to convert the goal-projected simulation data to a ranking whichcan be used for the directed component of FAST.
    ucscale : float, default=1
        Scaling factor for undirected component.
    truncation : str, default=None
        Method for truncating the prob distribution (None, 'cumsum', 'statecut'
    statetype : str, default='micro'
        What states (cluster, micro, macro) to use for calculations.
    macronum : int, default=8
        The number of macrostates to produce
    skip : int, default=1
        Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame
    lag : int, default=1
        The lagtime used to create the Markov model
    clustmethod : :class:`ClusterMixin <sklearn.base.ClusterMixin>` object, default=<class 'sklearn.cluster.k_means_.MiniBatchKMeans'>
        Clustering algorithm used to cluster the contacts or distances
    ticalag : int, default=20
        Lagtime to use for TICA in frames. When using `skip` remember to change this accordingly.
    ticadim : int, default=3
        Number of TICA dimensions to use. When set to 0 it disables TICA
    contactsym : str, default=None
        Contact symmetry
    save : bool, default=False
        Save the model generated

    Example
    -------
    >>> crystalSS = MetricSecondaryStructure().project(Molecule('crystal.pdb'))[0]
    >>>
    >>> # First argument of a goal function always has to be a Molecule object
    >>> def ssGoal(mol):
    >>>     proj = MetricSecondaryStructure().project(mol)
    >>>     ss_score = np.sum(proj == crystalSS, axis=1) / proj.shape[1]  # How many predicted SS match
    >>>     return ss_score
    >>>
    >>> adapt = AdaptiveGoal()
    >>> adapt.generatorspath = '../generators/'
    >>> adapt.nmin = 2
    >>> adapt.nmax = 3
    >>> adapt.nepochs = 2
    >>> adapt.ticadim = 3
    >>> adapt.projection = [MetricDistance('name CA', 'name N'), MetricDihedral()]
    >>> adapt.goalfunction = ssGoal
    >>> adapt.app = AcemdLocal()
    >>> adapt.run()
    >>>
    >>> # Or alternatively if we have a multi-argument goal function
    >>> def ssGoalAlt(mol, ss):
    >>>     proj = MetricSecondaryStructure().project(mol)
    >>>     ss_score = np.sum(proj == ss, axis=1) / proj.shape[1]
    >>>     return ss_score
    >>> from joblib import delayed
    >>> adapt.goalfunction = delayed(ssGoalAlt)(crystalSS)
    >>> adapt.app = AcemdLocal()
    >>> adapt.run()
    """

    def __init__(self):
        from sklearn.base import ClusterMixin
        from htmd.projections.projection import Projection
        super().__init__()
        self._cmdString('datapath', 'str', 'The directory in which the completed simulations are stored', 'data')
        self._cmdBoolean('filter', 'bool', 'Enable or disable filtering of trajectories.', True)
        self._cmdString('filtersel', 'str', 'Filtering atom selection', 'not water')
        self._cmdString('filteredpath', 'str', 'The directory in which the filtered simulations will be stored',
                        'filtered')
        self._cmdObject('projection', ':class:`Projection <htmd.projections.projection.Projection>` object',
                        'A Projection class object or a list of objects which will be used to project the simulation '
                        'data before constructing a Markov model', None, Projection)
        self._cmdFunction('goalfunction', 'function',
                          'This function will be used to convert the goal-projected simulation data to a ranking which'
                          'can be used for the directed component of FAST.', None)
        self._cmdValue('ucscale', 'float', 'Scaling factor for undirected component.', 1, TYPE_FLOAT, RANGE_ANY)
        self._cmdString('truncation', 'str', 'Method for truncating the prob distribution (None, \'cumsum\', '
                                             '\'statecut\'', None)
        self._cmdString('statetype', 'str', 'What states (cluster, micro, macro) to use for calculations.', 'micro')
        self._cmdValue('macronum', 'int', 'The number of macrostates to produce', 8, TYPE_INT, RANGE_POS)
        self._cmdValue('skip', 'int',
                       'Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame',
                       1, TYPE_INT, RANGE_POS)
        self._cmdValue('lag', 'int', 'The lagtime used to create the Markov model', 1, TYPE_INT, RANGE_POS)
        self._cmdObject('clustmethod', ':class:`ClusterMixin <sklearn.base.ClusterMixin>` object',
                        'Clustering algorithm used to cluster the contacts or distances', MiniBatchKMeans, ClusterMixin)
        self._cmdValue('ticalag', 'int',
                       'Lagtime to use for TICA in frames. When using `skip` remember to change this accordingly.', 20,
                       TYPE_INT, RANGE_0POS)
        self._cmdValue('ticadim', 'int', 'Number of TICA dimensions to use. When set to 0 it disables TICA', 3,
                       TYPE_INT, RANGE_0POS)
        self._cmdString('contactsym', 'str', 'Contact symmetry', None)
        self._cmdBoolean('save', 'bool', 'Save the model generated', False)

    def _algorithm(self):
        logger.info('Postprocessing new data')
        sims = simlist(glob(path.join(self.datapath, '*', '')), glob(path.join(self.inputpath, '*', 'structure.pdb')),
                       glob(path.join(self.inputpath, '*', '')))
        if self.filter:
            sims = simfilter(sims, self.filteredpath, filtersel=self.filtersel)

        metr = Metric(sims, skip=self.skip)
        metr.set(self.projection)

        # if self.contactsym is not None:
        #    contactSymmetry(data, self.contactsym)

        if self.ticadim > 0:
            # tica = TICA(metr, int(max(2, np.ceil(self.ticalag))))  # gianni: without project it was tooooo slow
            data = metr.project()
            ticalag = int(np.ceil(max(2, min(np.min(data.trajLengths)/2, self.ticalag))))  # 1 < ticalag < (trajLen / 2)
            tica = TICA(data, ticalag)
            datadr = tica.project(self.ticadim)
        else:
            datadr = metr.project()

        datadr.dropTraj()  # Preferably we should do this before any projections. Corrupted sims can affect TICA
        datadr.cluster(self.clustmethod(n_clusters=self._numClusters(datadr.numFrames)))
        model = Model(datadr)
        self._model = model
        self._model.markovModel(self.lag, self._numMacrostates(datadr))
        if self.save:
            self._model.save('adapt_model_e{}.dat'.format(self._getEpoch()))

        # Undirected component
        uc = -model.data.N  # Lower counts should give higher score hence the -
        if self.statetype == 'micro':
            uc = uc[model.cluster_ofmicro]
        if self.statetype == 'macro':
            uc = macroAccumulate(model, uc[model.cluster_ofmicro])

        # Calculating the directed component
        dc = self._calculateDirectedComponent(model.data.simlist, model.data.St, model.data.N)
        if self.statetype == 'micro':
            dc = dc[model.cluster_ofmicro]
        if self.statetype == 'macro':
            dc = macroAccumulate(model, dc[model.cluster_ofmicro])

        uc = self._featScale(uc)
        dc = self._featScale(dc)
        logger.debug('Undirected component: {}'.format(uc))
        logger.debug('Directed component: {}'.format(dc))

        reward = dc + self.ucscale * uc

        relFrames = self._getSpawnFrames(reward, self._model, datadr)
        self._writeInputs(datadr.rel2sim(np.concatenate(relFrames)))

    def _featScale(self, feat):
        return (feat - np.min(feat)) / (np.max(feat) - np.min(feat))

    def _getSpawnFrames(self, reward, model, data):
        (spawncounts, prob) = self._spawn(reward, self.nmax - self._running)
        logger.debug('spawncounts {}'.format(spawncounts))
        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(stateIdx, spawncounts[stateIdx], statetype=self.statetype, replacement=(data.K < 10))
        return relFrames

    def _spawn(self, ranking, N):
        if self.truncation is not None and self.truncation.lower() != 'none':
            if self.truncation == 'cumsum':
                idx = np.argsort(ranking)
                idx = idx[::-1]  # decreasing sort
                errs = ranking[idx]
                H = (N * errs / np.cumsum(errs)) < 1
                ranking[idx[H]] = 0
            if self.truncation == 'statecut':
                idx = np.argsort(ranking)
                idx = idx[::-1]  # decreasing sort
                ranking[idx[N:]] = 0  # Set all states ranked > N to zero.
        prob = ranking / np.sum(ranking)
        logger.debug('Sampling probabilities {}'.format(prob))
        spawnmicro = np.random.multinomial(N, prob)
        return spawnmicro, prob

    def _numClusters(self, numFrames):
        """ Heuristic that calculates number of clusters from number of frames """
        K = int(max(np.round(0.6 * np.log10(numFrames / 1000) * 1000 + 50), 100))  # heuristic
        if K > numFrames / 3:  # Ugly patch for low-data regimes ...
            K = int(numFrames / 3)
        return K

    def _numMacrostates(self, data):
        """ Heuristic for calculating the number of macrostates for the Markov model """
        macronum = self.macronum
        if data.K < macronum:
            macronum = np.ceil(data.K / 2)
            logger.warning(
                'Using less macrostates than requested due to lack of microstates. macronum = ' + str(macronum))

        # Calculating how many timescales are above the lag time to limit number of macrostates
        from pyemma.msm import timescales_msm
        timesc = timescales_msm(data.St.tolist(), lags=self.lag, nits=macronum).get_timescales()
        macronum = min(self.macronum, max(np.sum(timesc > self.lag), 2))
        return macronum

    def _calculateDirectedComponent(self, sims, St, N):
        from joblib import Parallel, delayed
        from htmd.util import _getNcpus
        from htmd.molecule.molecule import Molecule

        if hasattr(self.goalfunction, '__call__'):
            results = Parallel(n_jobs=_getNcpus(), verbose=0)(delayed(self.goalfunction)(Molecule(s)) for s in sims)
        elif isinstance(self.goalfunction, tuple) and hasattr(self.goalfunction[0], '__call__'):
            results = Parallel(n_jobs=_getNcpus(), verbose=0)(delayed(self.goalfunction[0])(Molecule(s), *self.goalfunction[1]) for s in sims)

        goalconcat = np.concatenate(results)
        stconcat = np.concatenate(St)
        clustermeans = np.bincount(stconcat, goalconcat.flatten())
        return clustermeans / N


class _AdaptiveGoalOld(AdaptiveGoal):
    def __init__(self):
        from htmd.projections.projection import Projection
        super().__init__()
        self._cmdObject('goalprojection', ':class:`Projection <htmd.projections.projection.Projection>` object',
                        'A Projection class object or a list of objects which will be used to project the simulation '
                        'data. This data will be used for the directed component of FAST.', None, Projection)

    def _calculateDirectedComponent(self, sims, St, N):
        metr = Metric(sims, skip=self.skip)
        metr.set(self.goalprojection)
        clustermeans = np.zeros(len(N))
        k = 0
        for proj in _projectionGenerator(metr, _getNcpus()):
            for pro in proj:
                clustermeans[:np.max(St[k])+1] += np.bincount(St[k], self.goalfunction(pro[0]).flatten())
                k += 1
        return clustermeans / N


if __name__ == '__main__':
    from htmd import *
    from htmd.adaptive.adaptivegoal import AdaptiveGoal
    import os
    import shutil
    from htmd.util import tempname

    def rmsdgoal(proj):
        return -proj  # Lower RMSDs should give higher score

    tmpdir = tempname()
    shutil.copytree(htmd.home() + '/data/adaptive/', tmpdir)
    os.chdir(tmpdir)
    md = AdaptiveGoal()
    md.dryrun = True
    md.nmin = 1
    md.nmax = 2
    md.nepochs = 3
    md.ticalag = 2
    md.ticadim = 3
    md.updateperiod = 5
    md.projection = MetricDistance('protein and name CA', 'resname BEN and noh')
    # md.goalprojection = MetricRmsd(Molecule(htmd.home() + '/data/adaptive/generators/1/structure.pdb'),
    #                               'protein and name CA')
    md.goalfunction = rmsdgoal
    # md.app = AcemdLocal()
    # md.run()
