from os import path, makedirs
import numpy as np
from htmd.adaptive.adaptive import AdaptiveBase
from protocolinterface import val
from htmd.projections.tica import TICA
from htmd.projections.metric import Metric
from htmd.clustering.regular import RegCluster
from htmd.adaptive.util import getParentSimIdxFrame, updatingMean
from sklearn.cluster import MiniBatchKMeans
import logging
logger = logging.getLogger(__name__)


class AdaptiveBandit(AdaptiveBase):
    """

    Parameters
    ----------
    app : :class:`SimQueue <jobqueues.simqueue.SimQueue>` object, default=None
        A SimQueue class object used to retrieve and submit simulations
    project : str, default='adaptive'
        The name of the project
    nmin : int, default=0
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
    projection : :class:`Projection <moleculekit.projections.projection.Projection>` object, default=None
        A Projection class object or a list of objects which will be used to project the simulation data before
        constructing a Markov model
    goalfunction : function, default=None
        This function will be used to convert the goal-projected simulation data to a ranking whichcan be used for the
        directed component of FAST.
    reward_method : str, default='max'
        The reward method
    skip : int, default=1
        Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame
    lag : int, default=1
        The lagtime used to create the Markov model. Units are in frames.
    exploration : float, default=0.5
        Exploration is the coefficient used in UCB algorithm to weight the exploration value
    temperature : int, default=300
        Temperature used to compute the free energy
    ticalag : int, default=20
        Lagtime to use for TICA in frames. When using `skip` remember to change this accordinly.
    ticadim : int, default=3
        Number of TICA dimensions to use. When set to 0 it disables TICA
    clustmethod : :class:`ClusterMixin <sklearn.base.ClusterMixin>` class, default=<class 'sklearn.cluster.k_means_.MiniBatchKMeans'>
        Clustering algorithm used to cluster the contacts or distances
    macronum : int, default=8
        The number of macrostates to produce
    save : bool, default=False
        Save the model generated
    save_qval : bool, default=False
        Save the Q(a) and N values for every epoch
    actionspace : str, default='metric'
        The action space
    recluster : bool, default=False
        If to recluster the action space.
    reclusterMethod : , default=<class 'sklearn.cluster.k_means_.MiniBatchKMeans'>
        Clustering method for reclustering.
    goal_init : float, default=0.3
        The proportional ratio of goal initialization compared to max frames set by nframes
    actionpool : int, default=0
        The number of top scoring actions used to randomly select respawning simulations
    """
    def __init__(self):
        from sklearn.base import ClusterMixin
        from moleculekit.projections.projection import Projection
        super().__init__()
        self._arg('datapath', 'str', 'The directory in which the completed simulations are stored', 'data', val.String())
        self._arg('filter', 'bool', 'Enable or disable filtering of trajectories.', True, val.Boolean())
        self._arg('filtersel', 'str', 'Filtering atom selection', 'not water', val.String())
        self._arg('filteredpath', 'str', 'The directory in which the filtered simulations will be stored', 'filtered', val.String())
        self._arg('projection', ':class:`Projection <moleculekit.projections.projection.Projection>` object',
                  'A Projection class object or a list of objects which will be used to project the simulation '
                   'data before constructing a Markov model', None, val.Object(Projection), nargs='+')
        self._arg('goalfunction', 'function',
                  'This function will be used to convert the goal-projected simulation data to a ranking which'
                  'can be used for the directed component of FAST.', None, val.Function(), nargs='any')
        self._arg('reward_method', 'str', 'The reward method', 'mean', val.String()) # Default should be mean
        self._arg('skip', 'int', 'Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame', 1, val.Number(int, 'POS'))
        self._arg('lag', 'int', 'The lagtime used to create the Markov model. Units are in frames.', 1, val.Number(int, 'POS'))
        self._arg('exploration', 'float', 'Exploration is the coefficient used in UCB algorithm to weight the exploration value', 0.01, val.Number(float, 'OPOS')) # Default changed to 0.01
        self._arg('temperature', 'int', 'Temperature used to compute the free energy', 300, val.Number(int, 'POS'))
        self._arg('ticalag', 'int', 'Lagtime to use for TICA in frames. When using `skip` remember to change this accordinly.', 20, val.Number(int, '0POS'))
        self._arg('ticadim', 'int', 'Number of TICA dimensions to use. When set to 0 it disables TICA', 3, val.Number(int, '0POS'))
        self._arg('clustmethod', ':class:`ClusterMixin <sklearn.base.ClusterMixin>` class', 'Clustering algorithm used to cluster the contacts or distances', MiniBatchKMeans, val.Class(ClusterMixin))
        self._arg('macronum', 'int', 'The number of macrostates to produce', 8, val.Number(int, 'POS'))
        self._arg('save', 'bool', 'Save the model generated', False, val.Boolean())
        self._arg('save_qval', 'bool', 'Save the Q(a) and N values for every epoch', False, val.Boolean())
        self._arg('actionspace', 'str', 'The action space', 'tica', val.String())
        self._arg('recluster', 'bool', 'If to recluster the action space.', False, val.Boolean())
        self._arg('reclusterMethod', '', 'Clustering method for reclustering.', MiniBatchKMeans)
        self._arg('goal_init', 'float', 'The proportional ratio of goal initialization compared to max frames set by nframes', 0.3, val.Number(float, 'POS'))
        self._arg('actionpool', 'int', 'The number of top scoring actions used to randomly select respawning simulations', 0, val.Number(int, 'OPOS'))

    def _algorithm(self):
        from htmd.kinetics import Kinetics
        sims = self._getSimlist()
        metr = Metric(sims, skip=self.skip)
        metr.set(self.projection)

        data = metr.project()
        data.dropTraj()  # Drop before TICA to avoid broken trajectories

        if self.goalfunction is not None:
            goaldata = self._getGoalData(data.simlist)
            if len(data.simlist) != len(goaldata.simlist):
                raise RuntimeError('The goal function was not able to project all trajectories that the MSM projection could. Check for possible errors in the goal function.')
            goaldataconcat = np.concatenate(goaldata.dat)
            if self.save:
                makedirs('saveddata', exist_ok=True)
                goaldata.save(path.join('saveddata', 'e{}_goaldata.dat'.format(self._getEpoch())))

        # tica = TICA(metr, int(max(2, np.ceil(self.ticalag))))  # gianni: without project it was tooooo slow
        if self.ticadim > 0:
            ticalag = int(np.ceil(max(2, min(np.min(data.trajLengths) / 2, self.ticalag))))  # 1 < ticalag < (trajLen / 2)
            tica = TICA(data, ticalag)
            datatica = tica.project(self.ticadim)
            if not self._checkNFrames(datatica): return False
            self._createMSM(datatica)
        else:
            if not self._checkNFrames(data): return False
            self._createMSM(data)

        confstatdist = self.conformationStationaryDistribution(self._model)
        if self.actionspace == 'metric':
            if not data.K:
                data.cluster(self.clustmethod(n_clusters=self._numClusters(data.numFrames)))
            data_q = data.copy()
        elif self.actionspace == 'goal':
            data_q = goaldata.copy()
        elif self.actionspace == 'tica':
            data_q = datatica.copy()
        elif self.actionspace == 'ticapcca':
            data_q = datatica.copy()
            for traj in data_q.trajectories:
                traj.cluster = self._model.macro_ofcluster[traj.cluster]
            data_q.K = self._model.macronum

        if self.recluster:
            print('Reclustering with {}'.format(self.reclusterMethod))
            data_q.cluster(self.reclusterMethod)
        
        numstates = data_q.K
        print('Numstates: {}'.format(numstates))
        currepoch = self._getEpoch()
        q_values = np.zeros(numstates, dtype=np.float32)
        n_values = np.zeros(numstates, dtype=np.int32)

        if self.goalfunction is not None:
            ## For every cluster in data_q, get the max score and initialize
            qstconcat = np.concatenate(data_q.St)
            statemaxes = np.zeros(numstates)
            np.maximum.at(statemaxes, qstconcat, np.squeeze(goaldataconcat))

            goalenergies = -Kinetics._kB * self.temperature * np.log(1-statemaxes)
            q_values = goalenergies
            n_values += int((self.nframes / self._numClusters(self.nframes)) * self.goal_init) ## Needs nframes to be set properly!!!!!!!!

        rewardtraj = np.arange(data_q.numTrajectories) # Recalculate reward for all states
        rewards = self.getRewards(rewardtraj, data_q, confstatdist, numstates, self.reward_method)
        for i in range(numstates):
            if len(rewards[i]) == 0:
                continue
            q_values[i] = updatingMean(q_values[i], n_values[i], rewards[i])
        n_values += np.array([len(x) for x in rewards])


        if self.save_qval:
            makedirs('saveddata', exist_ok=True)
            np.save(path.join('saveddata', 'e{}_qval.npy'.format(currepoch)), q_values)
            np.save(path.join('saveddata', 'e{}_nval.npy'.format(currepoch)), n_values)

        ucb_values = np.array([self.count_ucb(q_values[clust], self.exploration, currepoch + 1, n_values[clust]) for clust in range(numstates)])

        if self.save_qval:
            makedirs('saveddata', exist_ok=True)
            np.save(path.join('saveddata', 'e{}_ucbvals.npy'.format(currepoch)), ucb_values)

        N = self.nmax - self._running
        if self.actionpool <= 0:
            self.actionpool = N
       
        topactions = np.argsort(-ucb_values)[:self.actionpool]
        action = np.random.choice(topactions, N, replace=False)

        action_sel = np.zeros(numstates, dtype=int)
        action_sel[action] += 1
        while np.sum(action_sel) < N:  # When K is lower than N repeat some actions
            for a in action:
                action_sel[a] +=1
                if np.sum(action_sel) == N:
                    break

        if self.save_qval:
            np.save(path.join('saveddata', 'e{}_actions.npy'.format(currepoch)), action_sel)
        relFrames = self._getSpawnFrames_UCB(action_sel, data_q) 
        self._writeInputs(data.rel2sim(np.concatenate(relFrames)))
        return True

    def _getSimlist(self):
        from glob import glob
        from htmd.simlist import simlist, simfilter
        logger.info('Postprocessing new data')

        sims = simlist(glob(path.join(self.datapath, '*', '')), glob(path.join(self.inputpath, '*', '')),
                       glob(path.join(self.inputpath, '*', '')))

        if self.filter:
            sims = simfilter(sims, self.filteredpath, filtersel=self.filtersel)
        return sims


    def count_ucb(self, q_value, exploration, step, n_value):
        return (q_value + (exploration * np.sqrt((np.log(step) / (n_value + 1)))))

    def count_pucb(self, q_value, exploration, predictor, step, n_value):
        return (q_value + (exploration * predictor * np.sqrt((np.log(step) / (n_value + 1)))))

    def getRewards(self, trajidx, data_q, confstatdist, numstates, rewardmethod):
        from htmd.kinetics import Kinetics
        import pandas as pd
        rewards = [[] for _ in range(numstates)]
        for simidx in trajidx:
            # Get the eq distribution of each of the states the sim passed through
            states = data_q.St[simidx]
            statprob = confstatdist[simidx]
            connected = (states != -1) & (statprob != 0)
            if not np.any(connected):
                continue
            states = states[connected]
            statprob = statprob[connected]
            #energies = Kinetics._kB * self.temperature * np.log(statprob)
            energies = -Kinetics._kB * self.temperature * np.log(1-statprob)

            ww = len(energies)

            if rewardmethod == 'mean':
                windowedreward = pd.Series(energies[::-1]).rolling(ww, min_periods=1).mean().values[::-1]
            elif rewardmethod == 'max':
                windowedreward = pd.Series(energies[::-1]).rolling(ww, min_periods=1).max().values[::-1]
            else:
                raise RuntimeError('Reward method {} not available'.format(rewardmethod))

            for st, re in zip(states, windowedreward):
                rewards[st].append(re)

        return rewards

    def conformationStationaryDistribution(self, model):
        statdist = np.zeros(model.data.numFrames) # zero for disconnected set
        dataconcatSt = np.concatenate(model.data.St)
        for i in range(model.micronum):
            microframes = np.where(model.micro_ofcluster[dataconcatSt] == i)[0]
            statdist[microframes] = model.msm.stationary_distribution[i]
        return model.data.deconcatenate(statdist)

    def _checkNFrames(self, data):
        if self.nframes != 0 and data.numFrames >= self.nframes:
            logger.info('Reached maximum number of frames. Stopping adaptive.')
            return False
        return True

    def _getGoalData(self, sims):
        from htmd.projections.metric import Metric
        logger.debug('Starting projection of directed component')
        metr = Metric(sims, skip=self.skip)
        metr.set(self.goalfunction)
        data = metr.project()
        logger.debug('Finished calculating directed component')
        return data

    def _createMSM(self, data):
        from htmd.model import Model
        kmeanserror = True
        while kmeanserror:
            try:
                data.cluster(self.clustmethod(n_clusters=self._numClusters(data.numFrames)))
            except IndexError:
                continue
            kmeanserror = False
            
        self._model = Model(data)
        self._model.markovModel(self.lag, self._numMacrostates(data))
        if self.save:
            makedirs('saveddata', exist_ok=True)
            self._model.save(path.join('saveddata', 'e{}_adapt_model.dat'.format(self._getEpoch())))

    def _getSpawnFrames_UCB(self, reward, data):
        stateIdx = np.where(reward > 0)[0]
        _, relFrames = data.sampleClusters(stateIdx, reward[stateIdx], replacement=True, allframes=False)
        logger.debug('relFrames {}'.format(relFrames))
        return relFrames

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
            logger.warning('Using less macrostates than requested due to lack of microstates. macronum = ' + str(macronum))

        # Calculating how many timescales are above the lag time to limit number of macrostates
        from pyemma.msm import timescales_msm
        timesc = timescales_msm(data.St.tolist(), lags=self.lag, nits=macronum).get_timescales()
        macronum = min(self.macronum, max(np.sum(timesc > self.lag), 2))
        return macronum


import unittest
class _TestAdaptiveBandit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.util import tempname
        from htmd.home import home
        from moleculekit.projections.metricdistance import MetricDistance
        import shutil
        import os

        tmpdir = tempname()
        shutil.copytree(home(dataDir='adaptive'), tmpdir)
        os.chdir(tmpdir)
        print(f"Running AdaptiveBandit test in {tmpdir}")

    def test_adaptive(self):
        from sklearn.cluster import MiniBatchKMeans
        from jobqueues.localqueue import LocalCPUQueue
        from moleculekit.projections.metricdistance import MetricDistance

        import numpy as np
        import random
        np.random.seed(0)  # Needed for the clustering to always give same results
        random.seed(0)

        md = AdaptiveBandit()
        md.app = LocalCPUQueue()
        md.generatorspath = 'generators'
        md.inputpath = 'input'
        md.datapath = 'data'
        md.coorname = 'input.coor'
        md.filter = True
        md.filtersel = 'all'

        md.clustmethod = MiniBatchKMeans
        md.projection = MetricDistance('protein resid 173 and name CA', 'resname BEN and name C1 C2 C3 C7', periodic='selections')
        md.ticadim = 2
        md.nmin=1
        md.nmax=2
        md.nepochs= 9999
        md.nframes = 1000000

        md.reward_method = 'mean'
        md.exploration = 0.01
        md.actionspace = 'tica'
        md.actionpool = 0
        md.recluster = False

        md.save = True
        md.dryrun = True
        md.run()

if __name__ == '__main__':
    unittest.main(verbosity=2)