# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from os import path
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.model import macroAccumulate
from protocolinterface import val
import numpy as np
import os
import logging

logger = logging.getLogger(__name__)


class AdaptiveGoal(AdaptiveMD):
    """Adaptive sampling combining Markov state models with a goal function.

    Extends :class:`AdaptiveMD` by adding a directed component: a user-provided
    goal function scores conformations, and the respawning probability is a
    weighted sum of the undirected (MSM-based) and directed (goal-based) scores.

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
    boxname : str, default='input.xsc'
        Name of the file containing the starting box dimensions for the new simulations. Set to 'none' to disable box writing.
    lock : bool, default=False
        Lock the folder while adaptive is ongoing
    mps : int, default=0
        If mps > 0, it will run simulations using the Multi-Process Service (MPS) with the number of processes specified. If set to 0, mps is disabled
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
        The lagtime used to create the Markov model. Units are in frames.
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
    goalfunction : function, default=None
        This function will be used to convert the goal-projected simulation data to a ranking whichcan be used for the directed component of FAST.
    ucscale : float, default=0.5
        Scaling factor for undirected component. Directed component scaling automatically calculated as (1-uscale)
    nosampledc : bool, default=False
        Spawn only from top DC conformations without sampling
    autoscale : bool, default=False
        Automatically scales exploration and exploitation ratios depending on how stuck the adaptive is at a given goal score.
    autoscalemult : float, default=1
        Multiplier for the scaling factor.
    autoscaletol : float, default=0.2
        Tolerance for the scaling factor.
    autoscalediff : int, default=10
        Diff in epochs to use for scaling factor.
    savegoal : str, default=None
        Save the goal values to the specified file

    Examples
    --------
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
    >>> ag.projection = [MetricDistance('name CA', 'resname MOL', periodic='selections'), MetricDihedral()]
    >>> ag.goalfunction = ssGoal
    >>> ag.app = LocalGPUQueue()
    >>> ag.run()
    >>>
    >>> # Or alternatively if we have a multi-argument goal function
    >>> def ssGoalAlt(mol, ss):
    >>>     proj = MetricSecondaryStructure().project(mol)
    >>>     ss_score = np.sum(proj == ss, axis=1) / proj.shape[1]
    >>>     return ss_score
    >>> from joblib import delayed
    >>> ag.goalfunction = delayed(ssGoalAlt)(crystalSS)
    >>> ag.app = LocalGPUQueue()
    >>> ag.run()
    """

    def __init__(self):
        super().__init__()
        self._arg(
            "goalfunction",
            "function",
            "This function will be used to convert the goal-projected simulation data to a ranking which"
            "can be used for the directed component of FAST.",
            None,
            val.Function(),
            nargs="any",
        )
        self._arg(
            "ucscale",
            "float",
            "Scaling factor for undirected component. Directed component scaling "
            "automatically calculated as (1-uscale)",
            0.5,
            val.Number(float, "ANY"),
        )
        self._arg(
            "nosampledc",
            "bool",
            "Spawn only from top DC conformations without sampling",
            False,
            val.Boolean(),
        )
        self._arg(
            "autoscale",
            "bool",
            "Automatically scales exploration and exploitation ratios depending on "
            "how stuck the adaptive is at a given goal score.",
            False,
            val.Boolean(),
        )
        self._arg(
            "autoscalemult",
            "float",
            "Multiplier for the scaling factor.",
            1,
            val.Number(float, "0POS"),
        )
        self._arg(
            "autoscaletol",
            "float",
            "Tolerance for the scaling factor.",
            0.2,
            val.Number(float, "0POS"),
        )
        self._arg(
            "autoscalediff",
            "int",
            "Diff in epochs to use for scaling factor.",
            10,
            val.Number(int, "POS"),
        )
        self._arg(
            "savegoal",
            "str",
            "Save the goal values to the specified file",
            None,
            val.String(),
        )
        self._debug = False

    def _algorithm(self):
        sims = self._getSimlist()

        if self.nosampledc:
            print("Spawning only from top DC conformations without sampling")
            goaldata = self._getGoalData(sims)
            if not self._checkNFrames(goaldata):
                return False
            datconcat = np.concatenate(goaldata.dat).flatten()
            sortedabs = np.argsort(datconcat)[::-1]
            if self._debug:
                np.save("debug.npy", sortedabs[: self.nmax - self._running])
                return True
            self._writeInputs(goaldata.abs2sim(sortedabs[: self.nmax - self._running]))
            return True

        data = self._getData(sims)
        if self.save:
            os.makedirs("saveddata", exist_ok=True)
            np.savetxt(
                path.join("saveddata", f"e{self._getEpoch()}_report.npy"),
                [self._getEpoch(), data.numFrames, len(data.dat)],
            )

        if not self._checkNFrames(data):
            return False

        # Using the simlist of data in case some trajectories were dropped
        goaldata = self._getGoalData(data.simlist)
        if len(data.simlist) != len(goaldata.simlist):
            logger.warning(
                "The goal function was not able to project all trajectories that the MSM projection could."
                "Check for possible errors in the goal function."
            )

        # Ensuring that I use the intersection of projected simulations
        data.dropTraj(keepsims=goaldata.simlist)
        self._createMSM(data)

        model = self._model
        data = self._model.data

        # Undirected component
        uc = -model.data.N  # Lower counts should give higher score hence the -
        if self.statetype == "micro":
            uc = uc[model.cluster_ofmicro]
        if self.statetype == "macro":
            uc = macroAccumulate(model, uc[model.cluster_ofmicro])

        # Calculating the directed component
        dcmeans = dcstds = None
        if self.statetype == "micro":
            dcmeans, dcstds = self._calculateDirectedComponent(
                goaldata, model.data.St, model.micro_ofcluster
            )
        elif self.statetype == "macro":
            # TODO: Should we weigh by equilibrium population?
            dcmeans, dcstds = self._calculateDirectedComponent(
                goaldata, model.data.St, model.macro_ofcluster
            )

        ucunscaled = uc
        dcunscaled = dcmeans
        uc = self._featScale(uc)
        dc = self._featScale(dcmeans)

        scale = self.ucscale
        if self.autoscale:
            scale = AdaptiveGoal._calculateScale(
                goaldata, self.autoscalediff, self.autoscalemult, self.autoscaletol
            )
        reward = scale * uc + (1 - scale) * dc

        N = self.nmax - self._running
        reward = self._truncate(reward, N)
        relFrames, spawncounts, truncprob = self._getSpawnFrames(
            reward, self._model, data, N
        )

        if self.save:
            os.makedirs("saveddata", exist_ok=True)
            epoch = self._getEpoch()
            tosave = {
                "ucunscaled": -ucunscaled,
                "dcunscaled": dcunscaled,
                "uc": uc,
                "dc": dc,
                "ucscale": scale,
                "spawncounts": spawncounts,
                "truncprob": truncprob,
                "relFrames": relFrames,
                "dcmeans": dcmeans,
                "dcstds": dcstds,
                "reward": reward,
            }
            np.save(path.join("saveddata", f"e{epoch}_goalreport.npy"), tosave)
            np.save(path.join("saveddata", f"e{epoch}_spawnframes.npy"), relFrames)
            goaldata.save(path.join("saveddata", f"e{epoch}_goaldata.dat"))

        if self._debug:
            np.save("debug.npy", relFrames)
            return True
        self._writeInputs(data.rel2sim(np.concatenate(relFrames)))
        return True

    @staticmethod
    def _calculateScale(goaldata, epochdiff, multiplier=1, tolerance=0.2):
        from htmd.adaptive.adaptive import epochSimIndexes

        # Calculate the max goal of each epoch and the total min of the goal to normalize the goals to a [0, 1]
        epochs = epochSimIndexes(goaldata.simlist)
        g = np.zeros(len(epochs))
        totalmin = None

        dat = goaldata.dat
        for i, e in enumerate(sorted(epochs.keys())):
            idx = epochs[e]
            epochgoals = np.concatenate(dat[idx])
            g[i] = epochgoals.max()
            if totalmin is None or epochgoals.min() < totalmin:
                totalmin = epochgoals.min()

        rangeG = g.max() - totalmin

        # Calculate the dG
        # Prepending the first element epochdiff-times to calculate the dG
        g = np.hstack(([g[0]] * epochdiff, g))
        dG = np.abs(g[epochdiff:] - g[:-epochdiff]) / rangeG

        # Tolerance is the range of what we consider a significant change on a scale of [0, 1]
        # Multiplier is a scaling factor in case we want to react stronger to
        grad = -multiplier * (dG - tolerance)
        # Euler integration
        tstep = 1
        y = [0]
        # print(dG, g, grad)
        for t in range(1, len(dG), tstep):
            y.append(max(min(y[-1] + tstep * grad[t], 1), 0))
            # print('time: {} a: {} gradient: {} rangemax: {} rangemin: {}'.format(t, y[-1], grad[t], g.max(), totalmin))
        # print("END")

        # return y, grad, dG, g
        return y[-1]

        # dx = np.abs(np.diff(g))
        # dxn = self._featScale(dx)
        # scale = 0.5
        # for d in dxn[::-1]:
        #     if d < 0.1:
        #         scale = max(1, scale+0.1)
        # return scale

    def _getSpawnFrames(self, reward, model, data, N):
        prob = reward / np.sum(reward)
        logger.debug(f"Sampling probabilities {prob}")
        spawncounts = np.random.multinomial(N, prob)
        logger.debug(f"spawncounts {spawncounts}")

        stateIdx = np.where(spawncounts > 0)[0]
        _, relFrames = model.sampleStates(
            stateIdx, spawncounts[stateIdx], statetype=self.statetype, replacement=True
        )
        logger.debug(f"relFrames {relFrames}")
        return relFrames, spawncounts, prob

    def _featScale(self, feat):
        denom = np.max(feat) - np.min(feat)
        if denom == 0:  # Handling the trivial case where all states have equal features
            res = np.zeros(len(feat))
            res[:] = 1 / len(feat)
            return res
        return (feat - np.min(feat)) / denom

    def _getGoalData(self, sims):
        from htmd.projections.metric import Metric
        from htmd.simlist import _simName
        import pickle

        logger.debug("Starting projection of directed component")
        metr = Metric(sims, skip=self.skip)
        metr.set(self.goalfunction)
        data = metr.project()
        logger.debug("Finished calculating directed component")

        if self.savegoal is not None:
            dirname = os.path.dirname(self.savegoal)
            os.makedirs(dirname, exist_ok=True)
            savedata = {}
            for traj in data.trajectories:
                trajname = _simName(traj.sim.trajectory[0])
                savedata[trajname] = traj.projection
            with open(self.savegoal, "wb") as f:
                pickle.dump(savedata, f)

        return data

    def _calculateDirectedComponent(self, goaldata, St, mapping=None):
        import pandas as pd

        goalconcat = np.concatenate(goaldata.dat).flatten()
        stconcat = np.concatenate(St)
        if mapping is not None:
            stconcat = mapping[stconcat]

        x = pd.DataFrame({"a": stconcat})
        indexes = x.groupby("a").groups

        means = np.zeros(stconcat.max() + 1)
        stds = np.zeros(stconcat.max() + 1)
        for i in indexes:
            if i == -1:
                # Mappings have -1 on disconnected clusters (not used in the MSM)
                continue
            means[i] = np.mean(goalconcat[indexes[i]])
            stds[i] = np.std(goalconcat[indexes[i]])

        return means, stds

