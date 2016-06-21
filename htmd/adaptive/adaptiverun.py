# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from glob import glob
from os import path
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from htmd.adaptive.adaptive import Adaptive, AdaptiveNew
from htmd.simlist import simlist, simfilter
from htmd.projections.metricdistance import MetricDistance, MetricSelfDistance
from htmd.model import Model, macroAccumulate
from htmd.userinterface import uisetattr
from htmd.projections.tica import TICA
from htmd.projections.metric import Metric
from htmd.protocols.protocolinterface import TYPE_INT, RANGE_0POS, RANGE_POS
import logging
logger = logging.getLogger(__name__)


class AdaptiveRunNew(AdaptiveNew):
    """ Adaptive class which uses a Markov state model for respawning

    AdaptiveRun uses Markov state models to choose respawning poses for the next epochs. In more detail, it projects all
    currently retrieved simulations on either contacts or distances, clusters those and then builds a Markov model using
    currently retrieved simulations on either contacts or distances, clusters those and then builds a Markov model using
    the discretized trajectories. From the Markov model it then chooses conformations from the various states based on
    the chosen criteria which will be used for starting new simulations.

    Parameters
    ----------
    app : :class:`App <htmd.apps.app.App>` object
        An App class object used to retrieve and submit simulations
    project : str
        The name of the project
    nmin : int
        Minimum number of running simulations
    nmax : int
        Maximum number of running simulations
    nepochs : int
        Maximum number of epochs
    inputpath : str, optional
        The directory used to store input folders
    generatorspath : str, optional
        The directory containing the generators
    dryrun : bool, optional
        A dry run means that the adaptive will retrieve and generate a new epoch but not submit the simulations
    updateperiod : float, optional
        When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds
    metricsel1 : str
        The atomselection to be used for projecting coordinates to atom distances
    metricsel2 : str, optional
        Second atomselection for projecting coordinates to distances between the two selections
    metrictype : {'contacts', 'distances'}, optional, default='contacts'
        Metric type to user. Choose between contacts or distances
    datapath : str, optional
        The directory in which the completed simulations are stored
    filteredpath : str, optional
        The directory in which the filtered simulations will be stored
    macronum : int, optional
        The number of macrostates to produce
    skip : int, optional
        Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame
    lag : int, optional
        The lagtime used to create the Markov model
    clustmethod : :class:`ClusterMixin <sklearn.base.ClusterMixin>` object, optional
        Clustering algorithm used to cluster the contacts or distances
    method : str, optional
        Criteria used for choosing from which state to respawn from
    ticadim : int, optional
        Number of TICA dimensions to use. When set to 0 it disables TICA
    filtersel : str, optional
        Filtering atom selection

    Example
    -------
    >>> adapt = AdaptiveRun()
    >>> adapt.nmin = 2
    >>> adapt.nmax = 3
    >>> adapt.nepochs = 2
    >>> adapt.ticadim = 0
    >>> adapt.metricsel1 = 'name CA'
    >>> adapt.generatorspath = htmd.home()+'/data/dhfr'
    >>> adapt.app = AcemdLocal()
    >>> adapt.run()
    """

    def __init__(self):
        from sklearn.base import ClusterMixin
        super().__init__()
        self._cmdString('datapath', 'str', 'The directory in which the completed simulations are stored', 'data')
        self._cmdString('filteredpath', 'str', 'The directory in which the filtered simulations will be stored', 'filtered')
        self._cmdString('metricsel1', 'str', 'The atomselection to be used for projecting coordinates to atom distances', None)
        self._cmdString('metricsel2', 'str', 'Second atomselection for projecting coordinates to distances between the two selections', None)
        self._cmdString('metrictype', 'str', 'Metric type to user. Choose between contacts or distances', 'contacts', valid_values=['contacts', 'distances'])
        self._cmdValue('macronum', 'int', 'The number of macrostates to produce', 8, TYPE_INT, RANGE_POS)
        self._cmdValue('skip', 'int', 'Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame', 1, TYPE_INT, RANGE_POS)
        self._cmdValue('lag', 'int', 'The lagtime used to create the Markov model', 1, TYPE_INT, RANGE_POS)
        self._cmdObject('clustmethod', ':class:`ClusterMixin <sklearn.base.ClusterMixin>` object', 'Clustering algorithm used to cluster the contacts or distances', MiniBatchKMeans, ClusterMixin)
        self._cmdString('method', 'str', 'Criteria used for choosing from which state to respawn from', '1/Mc')
        self._cmdValue('ticadim', 'int', 'Number of TICA dimensions to use. When set to 0 it disables TICA', 0, TYPE_INT, RANGE_0POS)
        self._cmdString('filtersel', 'str', 'Filtering atom selection', 'not water')
        self._cmdString('contactsym', 'str', 'Contact symmetry', None)

    def _algorithm(self):
        logger.info('Postprocessing new data')
        datalist = simlist(glob(path.join(self.datapath, '*', '')), glob(path.join(self.inputpath, '*', 'structure.pdb')),
                           glob(path.join(self.inputpath, '*', '')))
        filtlist = simfilter(datalist, self.filteredpath, filtersel=self.filtersel)

        if hasattr(self, 'metricsel2') and self.metricsel2 is not None:
            proj = MetricDistance(self.metricsel1, self.metricsel2, metric=self.metrictype)
        else:
            proj = MetricSelfDistance(self.metricsel1, metric=self.metrictype)
        metr = Metric(filtlist, skip=self.skip)
        metr.projection(proj)
        data = metr.project()

        #if self.contactsym is not None:
        #    contactSymmetry(data, self.contactsym)

        data.dropTraj()
        if self.ticadim > 0:
            tica = TICA(data, int(max(2, np.ceil(20/self.skip))))
            datadr = tica.project(self.ticadim)
        else:
            datadr = data

        K = int(max(np.round(0.6 * np.log10(datadr.numFrames/1000)*1000+50), 100))  # heuristic
        if K > datadr.numFrames / 3: # Freaking ugly patches ...
            K = int(datadr.numFrames / 3)

        datadr.cluster(self.clustmethod(n_clusters=K), mergesmall=5)
        replacement = False
        if datadr.K < 10:
            datadr.cluster(self.clustmethod(n_clusters=K))
            replacement = True

        model = Model(datadr)
        macronum = self.macronum
        if datadr.K < macronum:
            macronum = np.ceil(datadr.K / 2)
            logger.warning('Using less macrostates than requested due to lack of microstates. macronum = ' + str(macronum))

        from pyemma.msm import timescales_msm
        timesc = timescales_msm(datadr.St.tolist(), lags=self.lag, nits=macronum).get_timescales()
        macronum = min(self.macronum, max(np.sum(timesc > self.lag), 2))

        model.markovModel(self.lag, macronum)
        p_i = self._criteria(model, self.method)
        (spawncounts, prob) = self._spawn(p_i, self.nmax - self._running)
        logger.debug('spawncounts {}'.format(spawncounts))
        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(stateIdx, spawncounts[stateIdx], statetype='micro', replacement=replacement)
        logger.debug('relFrames {}'.format(relFrames))

        self._writeInputs(datadr.rel2sim(np.concatenate(relFrames)))

    def _criteria(self, model, criteria):
        # TODO. REST OF CRITERIA!
        P_I = []
        if criteria == '1/Mc':
            nMicroPerMacro = macroAccumulate(model, np.ones(model.micronum))
            P_I = 1 / macroAccumulate(model, model.data.N[model.cluster_ofmicro])
            P_I = P_I / nMicroPerMacro
        return P_I[model.macro_ofmicro]

    def _spawn(self, ranking, N, truncated=False):
        if truncated:
            idx = np.argsort(ranking)
            idx = idx[::-1]  # decreasing sort
            errs = ranking[idx]
            H = (N * errs / np.cumsum(errs)) < 1
            ranking[idx[H]] = 0
        prob = ranking / np.sum(ranking)
        spawnmicro = np.random.multinomial(N, prob)
        return spawnmicro, prob


class AdaptiveRun(Adaptive):
    """ Adaptive class which uses a Markov state model for respawning

    AdaptiveRun uses Markov state models to choose respawning poses for the next epochs. In more detail, it projects all
    currently retrieved simulations on either contacts or distances, clusters those and then builds a Markov model using
    currently retrieved simulations on either contacts or distances, clusters those and then builds a Markov model using
    the discretized trajectories. From the Markov model it then chooses conformations from the various states based on
    the chosen criteria which will be used for starting new simulations.

    Parameters
    ----------
    app : :class:`App <htmd.apps.app.App>` object
        An App class object used to retrieve and submit simulations
    project : str
        The name of the project
    nmin : int
        Minimum number of running simulations
    nmax : int
        Maximum number of running simulations
    nepochs : int
        Maximum number of epochs
    inputpath : str, optional
        The directory used to store input folders
    generatorspath : str, optional
        The directory containing the generators
    dryrun : bool, optional
        A dry run means that the adaptive will retrieve and generate a new epoch but no submit the simulations
    updateperiod : float, optional
        When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds
    metricsel1 : str
        The atomselection to be used for projecting coordinates to atom distances
    metricsel2 : str, optional
        Second atomselection for projecting coordinates to distances between the two selections
    metrictype : {'contacts', 'distances'}, optional, default='contacts'
        Metric type to user. Choose between contacts or distances
    datapath : str, optional
        The directory in which the completed simulations are stored
    filteredpath : str, optional
        The directory in which the filtered simulations will be stored
    resultspath : str, optional
        The directory in which the results of the MSM analysis of the current epoch will be stored
    macronum : int, optional
        The number of macrostates to produce
    skip : int, optional
        Allows skipping of simulation frames to reduce data. i.e. skip=3 will only keep every third frame
    lag : int, optional
        The lagtime used to create the Markov model
    clustmethod : :class:`ClusterMixin <sklearn.base.ClusterMixin>` object, optional
        Clustering algorithm used to cluster the contacts or distances
    method : str, optional
        Criteria used for choosing from which state to respawn from
    ticadim : int, optional
        Number of TICA dimensions to use. When set to 0 it disables TICA
    filtersel : str, optional
        Filtering atom selection

    Example
    -------
    >>> adapt = AdaptiveRun()
    >>> adapt.nmin = 2
    >>> adapt.nmax = 3
    >>> adapt.nepochs = 2
    >>> adapt.ticadim = 0
    >>> adapt.metricsel1 = 'name CA'
    >>> adapt.generatorspath = htmd.home()+'/data/dhfr'
    >>> adapt.app = AcemdLocal()
    >>> adapt.run()
    """

    _cmds = ['metrictype', 'datapath', 'filteredpath', 'resultspath', 'macronum', 'skip', 'lag', 'clustmethod',
            'method','ticadim','metricsel1','metricsel2','contactsym','filtersel']

    def __init__(self, app=None, project=None, nmin=1, nmax=1, nepochs=1, inputpath='input', generatorspath='generators',
                 dryrun=False, updateperiod=0, metricsel1=None, metricsel2=None, metrictype='contacts',
                 datapath='data', filteredpath='filtered', resultspath='results', macronum=8, skip=1, lag=1,
                 clustmethod=MiniBatchKMeans, method='1/Mc', ticadim=0, filtersel='not water'):

        super().__init__(app, project, nmin, nmax, nepochs, inputpath, generatorspath, dryrun, updateperiod)
        self.metrictype = metrictype
        self.datapath = datapath
        self.filteredpath = filteredpath
        self.resultspath = resultspath
        self.macronum = macronum
        self.skip = skip
        self.lag = lag
        self.clustmethod = clustmethod
        self.method = method
        self.ticadim = ticadim
        self.metricsel1 = metricsel1
        self.metricsel2 = metricsel2
        self.filtersel = filtersel
        self.contactsym = None

    def __setattr__(self, key, value):
        all = self._adaptivecmds + self._cmds
        uisetattr(self, key, value, all)

    def _algorithm(self):

        logger.info('Postprocessing new data')

        datalist = simlist(glob(path.join(self.datapath, '*', '')),
                           glob(path.join(self.inputpath, '*', '*.prmtop'))[0],
                           glob(path.join(self.inputpath, '*', '')))

        filtlist = simfilter(datalist, self.filteredpath, filtersel=self.filtersel)

        # calculate metrics for tica projection
        if hasattr(self, 'metricsel2') and self.metricsel2 is not None:
            proj = MetricDistance(self.metricsel1, self.metricsel2, metric=self.metrictype)
        else:
            proj = MetricSelfDistance(self.metricsel1, metric=self.metrictype)

        metr = Metric(filtlist, skip=self.skip)
        metr.projection(proj)
        data = metr.project()

        #if self.contactsym is not None:
        #    contactSymmetry(data, self.contactsym)

        # projecting data with tICA
        data.dropTraj()
        if self.ticadim > 0:
            tica = TICA(data, int(max(2, np.ceil(20/self.skip))))
            datadr = tica.project(self.ticadim)
        else:
            datadr = data

        # what is this K definition???
        K = int(max(np.round(0.6 * np.log10(datadr.numFrames/1000)*1000+50), 100))  # heuristic
        if K > datadr.numFrames / 3: # Freaking ugly patches ...
            K = int(datadr.numFrames / 3)

        datadr.cluster(self.clustmethod(n_clusters=K), mergesmall=5)
        replacement = False
        if datadr.K < 10:
            datadr.cluster(self.clustmethod(n_clusters=K))
            replacement = True

        model = Model(datadr)
        macronum = self.macronum
        if datadr.K < macronum:
            macronum = np.ceil(datadr.K / 2)
            logger.warning('Using less macrostates than requested due to lack of microstates. macronum = ' + str(macronum))

        from pyemma.msm import timescales_msm
        timesc = timescales_msm(datadr.St.tolist(), lags=self.lag, nits=macronum).get_timescales()
        macronum = min(self.macronum, max(np.sum(timesc > self.lag), 2))

        model.markovModel(self.lag, macronum)
        p_i = self._criteria(model, self.method)
        (spawncounts, prob) = self._spawn(p_i, self.nmax-self.running)
        logger.debug('spawncounts {}'.format(spawncounts))
        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(stateIdx, spawncounts[stateIdx], statetype='micro', replacement=replacement)
        logger.debug('relFrames {}'.format(relFrames))

        self._writeInputs(datadr.rel2sim(np.concatenate(relFrames)))

    def _criteria(self, model, criteria):
        # TODO. REST OF CRITERIA!
        P_I = []
        if criteria == '1/Mc':
            nMicroPerMacro = macroAccumulate(model, np.ones(model.micronum))
            P_I = 1 / macroAccumulate(model, model.data.N[model.cluster_ofmicro])
            P_I = P_I / nMicroPerMacro
        return P_I[model.macro_ofmicro]

    def _spawn(self, ranking, N, truncated=False):
        if truncated:
            idx = np.argsort(ranking)
            idx = idx[::-1]  # decreasing sort
            errs = ranking[idx]
            H = (N * errs / np.cumsum(errs)) < 1
            ranking[idx[H]] = 0
        prob = ranking / np.sum(ranking)
        spawnmicro = np.random.multinomial(N, prob)
        return spawnmicro, prob


if __name__ == "__main__":
    from htmd import AcemdLocal
    import htmd
    from os import chdir
    from tempfile import mkdtemp

    chdir(mkdtemp())
    md = AdaptiveRun()
    md.dryrun = True
    md.nmin = 2
    md.nmax = 3
    md.nepochs = 2
    md.ticadim = 0
    md.metricsel1 = 'name CA'
    md.generatorspath = htmd.home()+'/data/dhfr'
    #md.app = AcemdLocal()
    # md.run()  # Disabled because I can't assume acemd is installed in travis machines

