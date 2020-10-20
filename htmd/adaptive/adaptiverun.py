# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from glob import glob
from os import path, makedirs
import numpy as np
from htmd.adaptive.adaptive import AdaptiveBase
from htmd.simlist import simlist, simfilter
from htmd.model import Model, macroAccumulate
from protocolinterface import val
from htmd.projections.tica import TICA
from htmd.projections.metric import Metric
import logging
logger = logging.getLogger(__name__)


class AdaptiveMD(AdaptiveBase):
    """ Adaptive class which uses a Markov state model for respawning

    AdaptiveMD uses Markov state models to choose respawning poses for the next epochs. In more detail, it projects all
    currently retrieved simulations according to the specified projection, clusters those and then builds a Markov model using
    the discretized trajectories. From the Markov model it then chooses conformations from the various states based on
    the chosen criteria which will be used for starting new simulations.

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
        Atom selection string for filtering.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    filteredpath : str, default='filtered'
        The directory in which the filtered simulations will be stored
    projection : :class:`Projection <moleculekit.projections.projection.Projection>` object, default=None
        A Projection class object or a list of objects which will be used to project the simulation data before constructing a Markov model
    truncation : str, default=None
        Method for truncating the prob distribution (None, 'cumsum', 'statecut'
    statetype : ('micro', 'cluster', 'macro'), str, default='micro'
        What states (cluster, micro, macro) to use for calculations.
    macronum : int, default=8
        The number of macrostates to produce
    skip : int, default=1
        Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame
    lag : int, default=1
        The lagtime used to create the Markov model. The units are in frames.
    clustmethod : :class:`ClusterMixin <sklearn.base.ClusterMixin>` class, default=<class 'htmd.clustering.kcenters.KCenter'>
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

    Example
    -------
    >>> adapt = AdaptiveMD()
    >>> adapt.nmin = 2
    >>> adapt.nmax = 3
    >>> adapt.nepochs = 2
    >>> adapt.ticadim = 3
    >>> adapt.projection = [MetricDistance('name CA', 'resname MOL', periodic='selections'), MetricDihedral()]
    >>> adapt.generatorspath = htmd.home()+'/data/dhfr'
    >>> adapt.app = LocalGPUQueue()
    >>> adapt.run()
    """

    def __init__(self):
        from sklearn.base import ClusterMixin
        from htmd.clustering.kcenters import KCenter
        from moleculekit.projections.projection import Projection
        super().__init__()
        self._arg('datapath', 'str', 'The directory in which the completed simulations are stored', 'data', val.String())
        self._arg('filter', 'bool', 'Enable or disable filtering of trajectories.', True, val.Boolean())
        self._arg('filtersel', 'str', 'Filtering atom selection', 'not water', val.String())
        self._arg('filteredpath', 'str', 'The directory in which the filtered simulations will be stored', 'filtered', val.String())
        self._arg('projection', ':class:`Projection <moleculekit.projections.projection.Projection>` object',
                  'A Projection class object or a list of objects which will be used to project the simulation '
                   'data before constructing a Markov model', None, val.Object(Projection), nargs='+')
        self._arg('truncation', 'str', 'Method for truncating the prob distribution (None, \'cumsum\', \'statecut\'', None, val.String())
        self._arg('statetype', 'str', 'What states (cluster, micro, macro) to use for calculations.', 'micro', val.String(), valid_values=('micro', 'cluster', 'macro'))
        self._arg('macronum', 'int', 'The number of macrostates to produce', 8, val.Number(int, 'POS'))
        self._arg('skip', 'int', 'Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame', 1, val.Number(int, 'POS'))
        self._arg('lag', 'int', 'The lagtime used to create the Markov model', 1, val.Number(int, 'POS'))
        self._arg('clustmethod', ':class:`ClusterMixin <sklearn.base.ClusterMixin>` class', 'Clustering algorithm used to cluster the contacts or distances', KCenter, val.Class(ClusterMixin))
        self._arg('method', 'str', 'Criteria used for choosing from which state to respawn from', '1/Mc', val.String())
        self._arg('ticalag', 'int', 'Lagtime to use for TICA in frames. When using `skip` remember to change this accordinly.', 20, val.Number(int, '0POS'))
        self._arg('ticadim', 'int', 'Number of TICA dimensions to use. When set to 0 it disables TICA', 3, val.Number(int, '0POS'))
        self._arg('contactsym', 'str', 'Contact symmetry', None, val.String())
        self._arg('save', 'bool', 'Save the model generated', False, val.Boolean())

    def _algorithm(self):
        data = self._getData(self._getSimlist())
        if not self._checkNFrames(data): return False
        self._createMSM(data)

        N = self.nmax - self._running
        reward = self._criteria(self._model, self.method)
        reward = self._truncate(reward, N)
        relFrames, _, _ = self._getSpawnFrames(reward, self._model, self._model.data, N)
        self._writeInputs(self._model.data.rel2sim(np.concatenate(relFrames)))
        return True

    def _checkNFrames(self, data):
        if self.nframes != 0 and data.numFrames >= self.nframes:
            logger.info('Reached maximum number of frames. Stopping adaptive.')
            return False
        return True

    def _getSimlist(self):
        logger.info('Postprocessing new data')
        sims = simlist(glob(path.join(self.datapath, '*', '')), glob(path.join(self.inputpath, '*', '')),
                       glob(path.join(self.inputpath, '*', '')))
        if self.filter:
            sims = simfilter(sims, self.filteredpath, filtersel=self.filtersel)
        return sims

    def _getData(self, sims):
        metr = Metric(sims, skip=self.skip)
        metr.set(self.projection)

        # if self.contactsym is not None:
        #    contactSymmetry(data, self.contactsym)

        if self.ticadim > 0:
            # tica = TICA(metr, int(max(2, np.ceil(self.ticalag))))  # gianni: without project it was tooooo slow
            data = metr.project()
            data.dropTraj()  # Drop before TICA to avoid broken trajectories
            ticalag = int(
                np.ceil(max(2, min(np.min(data.trajLengths) / 2, self.ticalag))))  # 1 < ticalag < (trajLen / 2)
            tica = TICA(data, ticalag)
            datadr = tica.project(self.ticadim)
        else:
            datadr = metr.project()
        datadr.dropTraj()  # Preferably we should do this before any projections. Corrupted sims can affect TICA
        return datadr

    def _createMSM(self, data):
        data.cluster(self.clustmethod(n_clusters=self._numClusters(data.numFrames)))
        self._model = Model(data)
        self._model.markovModel(self.lag, self._numMacrostates(data))
        if self.save:
            if not path.exists('saveddata'):
                makedirs('saveddata')
            self._model.save(path.join('saveddata', 'e{}_adapt_model.dat'.format(self._getEpoch())))

    def _getSpawnFrames(self, reward, model, data, N):
        prob = reward / np.sum(reward)
        logger.debug('Sampling probabilities {}'.format(prob))
        spawncounts = np.random.multinomial(N, prob)
        logger.debug('spawncounts {}'.format(spawncounts))

        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(stateIdx, spawncounts[stateIdx], statetype='micro', replacement=True)
        logger.debug('relFrames {}'.format(relFrames))
        return relFrames, spawncounts, prob

    def _criteria(self, model, criteria):
        if criteria == '1/Mc':
            nMicroPerMacro = macroAccumulate(model, np.ones(model.micronum))
            P_I = 1 / macroAccumulate(model, model.data.N[model.cluster_ofmicro])
            P_I = P_I / nMicroPerMacro
            ret = P_I[model.macro_ofmicro]
        elif criteria == 'pi/Mc':
            nMicroPerMacro = macroAccumulate(model, np.ones(model.micronum))
            P_I = 1 / macroAccumulate(model, model.data.N[model.cluster_ofmicro])
            P_I = P_I / nMicroPerMacro
            ret = P_I[model.macro_ofmicro]*model.msm.stationary_distribution
        return ret

    def _truncate(self, ranking, N):
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
        return ranking

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


if __name__ == "__main__":
    import htmd.home
    import os
    import shutil
    from htmd.util import tempname
    from moleculekit.projections.metricdistance import MetricDistance

    tmpdir = tempname()
    shutil.copytree(htmd.home.home()+'/data/adaptive/', tmpdir)
    os.chdir(tmpdir)
    md = AdaptiveMD()
    # md.dryrun = True
    md.nmin = 1
    md.nmax = 2
    md.nepochs = 3
    md.ticalag = 2
    md.ticadim = 3
    md.updateperiod = 5
    md.projection = MetricDistance('protein and name CA', 'resname BEN and noh', periodic='selections')
    md.projection = [MetricDistance('protein and name CA', 'resname BEN and noh', periodic='selections'), MetricDistance('protein and name CA', 'resname BEN and noh', periodic='selections')]


