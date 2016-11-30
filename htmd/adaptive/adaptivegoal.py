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
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.simlist import simlist, simfilter
from htmd.model import Model, macroAccumulate
from htmd.projections.tica import TICA
from htmd.projections.metric import Metric, _projectionGenerator
from htmd.protocols.protocolinterface import TYPE_INT, RANGE_0POS, RANGE_POS, TYPE_FLOAT, RANGE_ANY
import logging

logger = logging.getLogger(__name__)


class AdaptiveGoal(AdaptiveMD):
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
    nepochs : int, default=1000
        Stop adaptive once we have reached this number of epochs
    nframes : int, default=0
        Stop adaptive once we have simulated this number of aggregate simulation frames.
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
    lock : bool, default=False
        Lock the folder while adaptive is ongoing
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
    method : str, default='1/Mc'
        Criteria used for choosing from which state to respawn from
    ticalag : int, default=20
        Lagtime to use for TICA in frames. When using `skip` remember to change this accordinly.
    ticadim : int, default=3
        Number of TICA dimensions to use. When set to 0 it disables TICA
    contactsym : str, default=None
        Contact symmetry
    save : bool, default=False
        Save the model generated
    goalfunction : function, default=None
        This function will be used to convert the goal-projected simulation data to a ranking whichcan be used for the directed component of FAST.
    ucscale : float, default=1
        Scaling factor for undirected component.
    nosampledc : bool, default=False
        Spawn only from top DC conformations without sampling

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
    >>> ag = AdaptiveGoal()
    >>> ag.generatorspath = '../generators/'
    >>> ag.nmin = 2
    >>> ag.nmax = 3
    >>> ag.projection = [MetricDistance('name CA', 'name N'), MetricDihedral()]
    >>> ag.goalfunction = ssGoal
    >>> ag.app = AcemdLocal()
    >>> ag.run()
    >>>
    >>> # Or alternatively if we have a multi-argument goal function
    >>> def ssGoalAlt(mol, ss):
    >>>     proj = MetricSecondaryStructure().project(mol)
    >>>     ss_score = np.sum(proj == ss, axis=1) / proj.shape[1]
    >>>     return ss_score
    >>> from joblib import delayed
    >>> ag.goalfunction = delayed(ssGoalAlt)(crystalSS)
    >>> ag.app = AcemdLocal()
    >>> ag.run()
    """

    def __init__(self):
        super().__init__()
        self._cmdFunction('goalfunction', 'function',
                          'This function will be used to convert the goal-projected simulation data to a ranking which'
                          'can be used for the directed component of FAST.', None)
        self._cmdValue('ucscale', 'float', 'Scaling factor for undirected component.', 1, TYPE_FLOAT, RANGE_ANY)
        self._cmdBoolean('nosampledc', 'bool', 'Spawn only from top DC conformations without sampling', False)
        self._debug = False

    def _algorithm(self):
        sims = self._getSimlist()

        if self.nosampledc:
            print('Spawning only from top DC conformations without sampling')
            goaldata = self._getGoalData(sims)
            if not self._checkNFrames(goaldata): return False
            datconcat = np.concatenate(goaldata.dat).flatten()
            sortedabs = np.argsort(datconcat)[::-1]
            if self._debug: np.save('debug.npy', sortedabs[:self.nmax - self._running]); return True
            self._writeInputs(goaldata.abs2sim(sortedabs[:self.nmax - self._running]))
            return True

        data = self._getData(sims)
        if not self._checkNFrames(data): return False
        goaldata = self._getGoalData(data.simlist)  # Using the simlist of data in case some trajectories were dropped
        if len(data.simlist) != len(goaldata.simlist):
            logger.warning('The goal function was not able to project all trajectories that the MSM projection could.'
                           'Check for possible errors in the goal function.')
        data.dropTraj(keepsims=goaldata.simlist)  # Ensuring that I use the intersection of projected simulations
        self._createMSM(data)

        model = self._model
        data = self._model.data

        # Undirected component
        uc = -model.data.N  # Lower counts should give higher score hence the -
        if self.statetype == 'micro':
            uc = uc[model.cluster_ofmicro]
        if self.statetype == 'macro':
            uc = macroAccumulate(model, uc[model.cluster_ofmicro])

        # Calculating the directed component
        dc = self._calculateDirectedComponent(goaldata, model.data.St, model.data.N)
        if self.statetype == 'micro':
            dc = dc[model.cluster_ofmicro]
        if self.statetype == 'macro':
            dc = macroAccumulate(model, dc[model.cluster_ofmicro])

        uc = self._featScale(uc)
        dc = self._featScale(dc)
        logger.debug('Undirected component: {}'.format(uc))
        logger.debug('Directed component: {}'.format(dc))

        reward = dc + self.ucscale * uc

        relFrames = self._getSpawnFrames(reward, self._model, data)
        if self._debug: np.save('debug.npy', relFrames); return True
        self._writeInputs(data.rel2sim(np.concatenate(relFrames)))
        return True

    def _getSpawnFrames(self, reward, model, data):
        (spawncounts, prob) = self._spawn(reward, self.nmax - self._running)
        logger.debug('spawncounts {}'.format(spawncounts))
        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(stateIdx, spawncounts[stateIdx], statetype=self.statetype, replacement=True)
        return relFrames

    def _featScale(self, feat):
        denom = np.max(feat) - np.min(feat)
        if denom == 0:  # Handling the trivial case where all states have equal features
            res = np.zeros(len(feat))
            res[:] = 1 / len(feat)
            return res
        return (feat - np.min(feat)) / denom

    def _getGoalData(self, sims):
        logger.debug('Starting projection of directed component')
        metr = Metric(sims)
        metr.set(self.goalfunction)
        data = metr.project()
        logger.debug('Finished calculating directed component')
        return data

    def _calculateDirectedComponent(self, goaldata, St, N):
        goalconcat = np.concatenate(goaldata.dat).flatten()
        stconcat = np.concatenate(St)
        clustermeans = np.bincount(stconcat, goalconcat)
        return clustermeans / N


if __name__ == '__main__':
    from htmd import *
    from htmd.adaptive.adaptivegoal import AdaptiveGoal
    import os
    import shutil
    from htmd.util import tempname
    from joblib import delayed
    from htmd.home import home

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

    # Some real testing now
    from htmd.projections.metricsecondarystructure import MetricSecondaryStructure
    from htmd.projections.metricdistance import MetricSelfDistance
    import numpy as np

    os.chdir(path.join(home(), 'data', 'test-adaptive'))

    goalProjectionDict = {'ss': MetricSecondaryStructure(),
                          'contacts': MetricSelfDistance('protein and name CA', metric='contacts', threshold=10),
                          'ss_contacts': [MetricSecondaryStructure(),
                                          MetricSelfDistance('protein and name CA', metric='contacts', threshold=10)]}

    def getLongContacts(crystal, long=8):
        crystalMap = MetricSelfDistance('protein and name CA', metric='contacts', threshold=10, pbc=False).getMapping(
            crystal)
        indexes = np.vstack(crystalMap.atomIndexes.as_matrix())
        return crystal.resid[indexes[:, 1]] - crystal.resid[indexes[:, 0]] > long

    def getCrystalSS(crystal):
        return MetricSecondaryStructure().project(crystal)[0]

    def getCrystalCO(crystal):
        crystalCO = MetricSelfDistance('protein and name CA', metric='contacts', threshold=10, pbc=False).project(
            crystal)
        longCO = getLongContacts(crystal)
        return crystalCO & longCO

    def ssContactGoal(mol, crystal, project=True, crystalSS=None, crystalCO=None):
        if crystalSS is None:
            crystalSS = getCrystalSS(crystal)
        if crystalCO is None:
            crystalCO = getCrystalCO(crystal)

        if project:
            projss = goalProjectionDict['ss'].copy().project(mol)
            projco = goalProjectionDict['contacts'].copy().project(mol)
        else:
            projss = mol[:, :len(crystalSS)]
            projco = mol[:, len(crystalSS):]

        if len(crystalCO) != projco.shape[1]:
            raise RuntimeError(
                'Different lengths between crystal {} and traj {} contacts for fileloc {}'.format(len(crystalCO),
                                                                                                  projco.shape[1],
                                                                                                  mol.fileloc))
        if len(crystalSS) != projss.shape[1]:
            raise RuntimeError(
                'Different lengths between crystal {} and traj {} SS for fileloc {}'.format(len(crystalSS),
                                                                                            projss.shape[1],
                                                                                            mol.fileloc))

        ss_score = np.sum(projss == crystalSS, axis=1) / projss.shape[1]
        co_score = np.sum(projco[:, crystalCO] == 1, axis=1) / np.sum(crystalCO)  # Predicted conts are True?
        return 0.6 * ss_score + 0.4 * co_score

    refmol = Molecule('ntl9_2hbb.pdb')
    crystalSS = getCrystalSS(refmol)
    crystalCO = getCrystalCO(refmol)

    np.random.seed(10)
    ad = AdaptiveGoal()
    ad.app = LocalGPUQueue()
    ad.nmin = 10
    ad.nmax = 20
    ad.nepochs = 999999
    # ad.nframes = nframes['ntl9'] test that as well
    ad.generatorspath = '../../generators/'
    ad.projection = MetricSelfDistance('protein and name CA')
    ad.goalfunction = delayed(ssContactGoal)(refmol, True, crystalSS, crystalCO)
    ad.statetype = 'micro'
    ad.truncation = 'cumsum'
    ad._debug = True
    ad.nosampledc = True
    ad.run()
    assert np.array_equal(np.load('debug.npy'), np.load('ref_nosampledc.npy'))

    # TODO: Make this test work. Seems to ignore the random seed
    #np.random.seed(10)
    #ad.nosampledc = False
    #ad.run()
    #assert np.array_equal(np.load('debug.npy'), np.load('ref.npy'))
