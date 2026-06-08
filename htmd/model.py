"""
Markov state models are a statistical tool for analysing molecular simulations which has met with lots of success.
The Model class here, encapsulates all functionallity for the calculation of Markov models while hiding unnecessary
details under the hood. It uses DeepTime [1]_ internally to calculate Markov models.

References
----------
.. [1] Deeptime: a Python library for machine learning dynamical models from time series data. Moritz Hoffmann et al 2022 Mach. Learn.: Sci. Technol. 3 015009
"""

# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import numpy as np
import random
from moleculekit.molecule import Molecule
from htmd.units import convert as unitconvert
import logging

if TYPE_CHECKING:
    from htmd.metricdata import MetricData

logger = logging.getLogger(__name__)


class Model(object):
    """Constructor for the Model class.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object, optional
        A :class:`MetricData <htmd.metricdata.MetricData>` object containing the discretized trajectories
    file : str, optional
        Path to a previously saved Model file to load. Provide either `data` or `file`.

    Examples
    --------
    >>> model = Model(data)

    .. rubric:: Methods
    .. autoautosummary:: htmd.model.Model
        :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.model.Model
        :attributes:
    """

    def __init__(self, data: "MetricData | None" = None, file: str | None = None):
        if data is None:
            if file is not None:
                self.load(file)
            return
        else:
            self.data = data
        self.hmm = None
        self._modelid = None
        if self.data._clusterid is None:
            raise NameError(
                "You need to cluster your data before making a Markov model"
            )
        if self.data._dataid != self.data._clusterid:
            raise NameError(
                "You have modified the data in data.dat after clustering. Please re-cluster."
            )
        self._clusterid = self.data._clusterid

    @staticmethod
    def _get_model(statelist, lagtime, bayesian_samples=None):
        from deeptime.markov.msm import MaximumLikelihoodMSM, BayesianMSM
        from deeptime.markov import TransitionCountEstimator

        count_mode = "sliding" if bayesian_samples is None else "effective"

        counts = TransitionCountEstimator(
            lagtime=lagtime, count_mode=count_mode
        ).fit_fetch(statelist)

        if bayesian_samples is not None:
            return BayesianMSM(n_samples=bayesian_samples).fit_fetch(
                counts.submodel_largest()
            )
        else:
            return MaximumLikelihoodMSM(
                allow_disconnected=False, use_lcc=False
            ).fit_fetch(counts.submodel_largest())

    def markovModel(
        self,
        lag: float,
        macronum: int,
        units: str = "frames",
        sparse: bool = False,
    ):
        """Build a Markov model at a given lag time and calculate metastable states

        Parameters
        ----------
        lag : float
            The lag time at which to calculate the Markov state model. The units are specified with the `units` argument.
        macronum : int
            The number of macrostates (metastable states) to produce
        units : str
            The units of lag. Can be 'frames' or any time unit given as a string.
        sparse : bool
            Make the transition matrix sparse. Useful if lots (> 4000) states are used for the MSM. Warning: untested.

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(150, 4)  # 150 frames lag, 4 macrostates
        """
        self._integrityCheck(markov=True)

        lag = unitconvert(units, "frames", lag, fstep=self.data.fstep)

        statelist = [traj.cluster for traj in self.data.trajectories]
        self.lag = lag
        self.msm = self._get_model(statelist, lag, bayesian_samples=None)
        modelflag = False
        while not modelflag:
            try:
                self.coarsemsm = self.msm.pcca(macronum)
            except Exception as e:
                macronum -= 1
                if macronum < 2:
                    raise RuntimeError(
                        "Could not create even two macrostates. Please revise your clustering."
                    )
                logger.warning(
                    f"PCCA failed with following error. Reducing the number of macrostates to {macronum}. Error: {e}"
                )
                continue

            if len(np.unique(self.coarsemsm.assignments)) != macronum:
                macronum -= 1
                logger.warning(
                    f"PCCA returned empty macrostates. Reducing the number of macrostates to {macronum}."
                )
            else:
                modelflag = True
            if macronum < 2:
                raise RuntimeError(
                    "Could not create even two macrostates. Please revise your clustering."
                )

        self._modelid = random.random()

        logger.info(f"{self.msm.count_fraction * 100:.1f}% of the data was used")

        _macroTrajectoriesReport(
            self.macronum,
            _macroTrajSt(self.data.St, self.macro_ofcluster),
            self.data.simlist,
        )

    def createState(
        self,
        microstates: list | np.ndarray | None = None,
        indexpairs: list | None = None,
    ):
        """Creates a new state. Works both for new clusters and macrostates.

        If creating a new cluster, it just reassigns the given frames to the new cluster.
        If creating a new macrostate, it removes the given microstates from their previous macrostate, creates a new one
        and assigns them to it.

        Parameters
        ----------
        microstates : list or np.ndarray
            The microstates to split out into a new macrostate.
        indexpairs : list
            List of lists. Each row is a simulation index-frame pair which should be added to a new cluster.
        """
        if microstates is not None and indexpairs is not None:
            raise AttributeError(
                "microstates and indexpairs arguments are mutually exclusive"
            )
        if microstates is not None:
            newmacro = self.macronum

            # Fixing hard assignments
            self.coarsemsm.assignments[microstates] = newmacro

            todelete = np.setdiff1d(np.arange(newmacro), self.coarsemsm.assignments)

            # Fixing memberships. Padding the array with 0s for the new macrostate
            self.coarsemsm._memberships = np.pad(
                self.coarsemsm.memberships,
                ((0, 0), (0, 1)),
                mode="constant",
                constant_values=(0),
            )
            self.coarsemsm._memberships[microstates, :] = 0
            self.coarsemsm._memberships[microstates, -1] = 1

            # Moving probabilities of empty states to new one
            othermicro = np.ones(self.micronum, dtype=bool)
            othermicro[microstates] = False
            othermicro = np.where(othermicro)[0]
            self.coarsemsm._memberships[othermicro, -1] = np.sum(
                self.coarsemsm.memberships[othermicro[:, None], todelete], axis=1
            )

            # Fixing distributions
            self.coarsemsm._metastable_distributions = np.pad(
                self.coarsemsm.metastable_distributions,
                ((0, 1), (0, 0)),
                mode="constant",
                constant_values=(0),
            )
            self.coarsemsm.metastable_distributions[-1, microstates] = 1 / len(
                microstates
            )
        elif indexpairs is not None:
            newcluster = self.data.K
            for ip in indexpairs:
                self.data.trajectories[ip[0]].cluster[ip[1]] = newcluster
            self.data.K += 1
            self.data.N = np.bincount(np.concatenate(self.data.St))

    @property
    def _active_set(self):
        return self.msm.count_model.state_symbols

    @property
    def P(self) -> np.ndarray:
        """The transition probability matrix"""
        return self.msm.transition_matrix

    @property
    def micro_ofcluster(self) -> np.ndarray:
        """Mapping of clusters to microstates

        Numpy array which at index i has the index of the microstate corresponding to cluster i.
        Clusters which were not connected and thus are not in the model have a microstate value of -1.
        """
        self._integrityCheck(postmsm=True)
        micro_ofcluster = -np.ones(self.data.K, dtype=int)
        micro_ofcluster[self._active_set] = np.arange(len(self._active_set))
        return micro_ofcluster

    @property
    def cluster_ofmicro(self) -> np.ndarray:
        """Mapping of microstates to clusters

        Numpy array which at index i has the index of the cluster corresponding to microstate i.
        """
        self._integrityCheck(postmsm=True)
        return self._active_set

    @property
    def micronum(self) -> int:
        """Number of microstates"""
        self._integrityCheck(postmsm=True)
        return len(self._active_set)

    @property
    def macronum(self) -> int:
        """Number of macrostates"""
        self._integrityCheck(postmsm=True)
        return len(set(self.coarsemsm.assignments))

    @property
    def macro_ofmicro(self) -> np.ndarray:
        """Mapping of microstates to macrostates

        Numpy array which at index i has the index of the macrostate corresponding to microstate i.
        """
        self._integrityCheck(postmsm=True)
        # Fixing deeptime macrostate numbering
        mask = np.ones(np.max(self.coarsemsm.assignments) + 1, dtype=int) * -1
        mask[list(set(self.coarsemsm.assignments))] = range(self.macronum)
        return mask[self.coarsemsm.assignments]

    @property
    def macro_ofcluster(self) -> np.ndarray:
        """Mapping of clusters to macrostates

        Numpy array which at index i has the index of the macrostate corresponding to cluster i.
        Clusters which were not connected and thus are not in the model have a macrostate value of -1.
        """
        self._integrityCheck(postmsm=True)
        macro_ofcluster = -np.ones(self.data.K, dtype=int)
        macro_ofcluster[self._active_set] = self.macro_ofmicro
        return macro_ofcluster

    def _plot_implied_timescales(
        self,
        fstep,
        data,
        n_its=None,
        process=None,
        show_mle: bool = True,
        show_samples: bool = True,
        show_sample_mean: bool = True,
        show_sample_confidence: bool = True,
        show_cutoff: bool = True,
        sample_confidence: float = 0.95,
        colors=None,
        ax=None,
        ylog=False,
        **kwargs,
    ):
        r"""Creates an implied timescales plot inside exising matplotlib axes.

        ! Taken from deeptime !
        .. plot:: examples/plot_implied_timescales.py

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The matplotlib axes to use for plotting.
        data : ImpliedTimescales
            A timescales data container object, can be obtained, e.g., via :meth:`ImpliedTimescales.from_models`.
        n_its : int, optional, default=None
            Maximum number of timescales to plot.
        process : int, optional, default=None
            A particular process to plot. This is mutually exclusive with n_its.
        show_mle : bool, default=True
            Whether to show the timescale of the maximum-likelihood estimate.
        show_samples : bool, default=True
            Whether to show sample means and/or confidences.
        show_sample_mean : bool, default=True
            Whether to show the sample mean. Only has an effect if show_samples is True and there are samples in the data.
        show_sample_confidence : bool, default=True
            Whether to show the sample confidence. Only has an effect if show_samples is True and there are samples
            in the data.
        show_cutoff : bool, default=True
            Whether to show the model resolution cutoff as grey filled area.
        sample_confidence : float, default=0.95
            The confidence to plot. The default amounts to a shaded area containing 95% of the sampled values.
        colors : list of colors, optional, default=None
            The colors that should be used for timescales. By default uses the matplotlib default colors as per
            rc-config value "axes.prop_cycle".
        **kwargs
            Keyword arguments which are forwarded into the matplotlib plotting function for timescales.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The matplotlib axes that were used to plot the timescales.

        See Also
        --------
        ImpliedTimescales
        """
        from deeptime.plots.util import default_colors
        from deeptime.util import confidence_interval

        if ax is None:
            import matplotlib.pyplot as plt

            ax = plt.gca()
        if n_its is not None and process is not None:
            raise ValueError("n_its and process are mutually exclusive.")
        if process is not None and process >= data.max_n_processes:
            raise ValueError(
                f"Requested process {process} when only {data.max_n_processes} are available."
            )

        if process is None and n_its is None:
            n_its = data.max_n_processes
        it_indices = [process] if process is not None else np.arange(n_its)
        if colors is None:
            colors = default_colors()
        for it_index in it_indices:
            color = colors[it_index % len(colors)]
            if show_mle:
                ax.plot(
                    data.lagtimes * fstep,
                    data.timescales_for_process(it_index) * fstep,
                    color=color,
                    **kwargs,
                )
            if data.has_samples and show_samples:
                its_samples = data.samples_for_process(it_index)
                if show_sample_mean:
                    sample_mean = np.nanmean(its_samples, axis=1)
                    ax.plot(
                        data.lagtimes * fstep,
                        sample_mean * fstep,
                        marker="o",
                        linestyle="dashed",
                        color=color,
                    )
                if show_sample_confidence:
                    l_conf, r_conf = confidence_interval(
                        its_samples.T, conf=sample_confidence, remove_nans=True
                    )
                    ax.fill_between(
                        data.lagtimes * fstep,
                        l_conf * fstep,
                        r_conf * fstep,
                        alpha=0.2,
                        color=color,
                    )

        if show_cutoff:
            ax.plot(
                data.lagtimes * fstep, data.lagtimes * fstep, linewidth=2, color="black"
            )
            ax.fill_between(
                data.lagtimes * fstep,
                np.full((data.n_lagtimes,), fill_value=ax.get_ylim()[0]) * fstep,
                data.lagtimes * fstep,
                alpha=0.5,
                color="grey",
            )
        if ylog:
            ax.set_yscale("log")

        ax.set_title("Implied timescales")
        ax.set_xlabel("Lag time (ns)")
        ax.set_ylabel("Timescale (ns)")
        return ax

    def plotTimescales(
        self,
        lags: list | range | np.ndarray | None = None,
        minlag: float | None = None,
        maxlag: float | None = None,
        numlags: int = 25,
        units: str = "frames",
        errors: int | None = None,
        nits: int | None = None,
        results: bool = False,
        plot: bool = True,
        save: str | None = None,
        njobs: int = 1,
        ylog: bool = True,
    ):
        """Plot the implied timescales of MSMs of various lag times

        Parameters
        ----------
        lags : list or range or np.ndarray
            Specify explicitly at which lag times to compute the timescales.
        minlag : float
            The minimum lag time for the timescales. Used in combination with `maxlag` and `numlags`.
        maxlag : float
            The maximum lag time for the timescales. If None will default to the mode length of the trajectories.
        numlags : int
            The number of points to place between `minlag` and `maxlag`.
        units : str
            The units of lag. Can be 'frames' or any time unit given as a string.
        errors : int
            Number of Bayesian samples used to calculate errors with Bayesian MSMs. If None, no errors are calculated.
        nits : int
            Number of implied timescales to calculate. If None, all are calculated.
        results : bool
            If the method should return the calculated implied timescales
        plot : bool
            If the method should display the plot of implied timescales
        save : str
            Path of the file in which to save the figure
        njobs : int
            DEPRECATED: Number of parallel jobs to spawn for calculation of timescales. Negative numbers are used for spawning jobs as many as CPU threads.
            -1: for all CPUs -2: for all except one etc.
        ylog : bool
            Set to False to get linear y axis instead of logarithmic

        Returns
        -------
        If given results=True this method will return the following data
        its : np.ndarray
            The calculated implied timescales. 2D array with dimensions (len(`lags`), `nits`)
        lags : np.ndarray
            A list of the lag times that were used to calculate the implied timescales

        Examples
        --------
        >>> model = Model(data)
        >>> model.plotTimescales()
        >>> model.plotTimescales(lags=list(range(1,100,5)))
        >>> model.plotTimescales(minlag=0.1, maxlag=20, numlags=25, units='ns')
        """
        from deeptime.util.validation import implied_timescales
        from tqdm import tqdm

        self._integrityCheck()
        if lags is None:
            lags = self.data._defaultLags(minlag, maxlag, numlags, units)
        else:
            lags = unitconvert(units, "frames", lags, fstep=self.data.fstep).tolist()

        if nits is None:
            nits = np.min((self.data.K - 1, 20))

        statelist = [traj.cluster for traj in self.data.trajectories]
        models = []
        for lagtime in tqdm(lags, desc="Estimating Timescales"):
            models.append(self._get_model(statelist, lagtime, bayesian_samples=errors))

        its_data = implied_timescales(models, n_its=nits)
        if plot or (save is not None):
            from matplotlib import pylab as plt

            if plot:
                plt.ion()
            plt.figure()
            ax = plt.gca()
            try:
                self._plot_implied_timescales(
                    self.data.fstep, data=its_data, n_its=nits, ax=ax, ylog=ylog
                )
            except ValueError as ve:
                plt.close()
                raise ValueError(
                    f"{ve} This is probably caused by badly set fstep in the data ({self.data.fstep}). "
                    + "Please correct the model.data.fstep to correspond to the simulation frame step in nanoseconds."
                )
            if save is not None:
                plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
            if plot:
                plt.show()
        if results:
            return its_data._its, its_data.lagtimes

    def maxConnectedLag(self, lags: list | np.ndarray, njobs: int = 1):
        """Heuristic for getting the lagtime before a timescale drops.

        It calculates the last lagtime before a drop occurs in the first implied timescale due to disconnected states.
        If the top timescale is closer to the second top timescale at the previous lagtime than to itself at the previous
        lagtime it means that a drop occured. The lagtime before the drop is returned.

        Parameters
        ----------
        lags : list or np.ndarray
            A list of lag times for which to calculate the implied timescales
        njobs : int
            Number of parallel jobs (currently unused).

        Returns
        -------
        ml : int
            The maximum lagtime before a drop occurs in the top timescale

        Examples
        --------
        >>> model = Model(data)
        >>> model.maxConnectedLag(list(range(1, 100, 5)))
        """
        from deeptime.util.validation import implied_timescales
        from tqdm import tqdm

        if len(lags) == 1:
            return lags
        if isinstance(lags, np.ndarray):
            lags = lags.astype(int)

        statelist = [traj.cluster for traj in self.data.trajectories]
        models = []
        for lagtime in tqdm(lags, desc="Estimating Timescales"):
            models.append(self._get_model(statelist, lagtime, bayesian_samples=None))

        its_data = implied_timescales(models, n_its=2)
        itime = its_data._its

        for i in range(1, np.size(itime, 0)):
            if abs(itime[i, 0] - itime[i - 1, 1]) < abs(itime[i, 0] - itime[i - 1, 0]):
                lagidx = i - 1
                break
            else:
                lagidx = i
        return lags[lagidx], itime

    def sampleStates(
        self,
        states: list | int | range | np.ndarray | None = None,
        frames: int | list | np.ndarray | None = 20,
        statetype: str = "macro",
        replacement: bool = False,
        samplemode: str = "random",
        allframes: bool = False,
    ):
        """Samples frames from a set of states

        Parameters
        ----------
        states : list or int
            A list of state indexes from which we want to sample
        frames : int or list
            An integer with the number of frames we want to sample per state or a list of same length as
            `states` which contains the number of frames we want from each of the states.
            If set to None it will return all frames of the states.
        statetype : str
            The type of state we want to sample from. Can be 'micro', 'macro' or 'cluster'.
        replacement : bool
            If we want to sample with or without replacement.
        samplemode : str
            What sampling strategy to use. Can be 'random', 'even' or 'weighted'. For `statetype` == 'macro' this can be
            set to 'even' to sample evenly from all microstates in the macrostate or to 'weighted' to sample
            proportional to the equilibrium probability of each microstate in the macrostate.
        allframes : bool
            Deprecated. Use frames=None instead.

        Returns
        -------
        absframes : numpy.ndarray
            An array which contains for each state an array containing absolute trajectory frames
        relframes : numpy.ndarray
            An array which contains for each state a 2D array containing the trajectory ID and frame number for each of
            the sampled frames

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.sampleStates(range(5), [10, 3, 2, 50, 1])  # Sample from all 5 macrostates
        >>> model.sampleStates(range(model.micronum), samplesnum, statetype='micro')  # Sample from all microstates
        """
        if statetype == "cluster":
            if samplemode != "random":
                logger.warning(
                    "'cluster' states incompatible with 'samplemode' other than 'random'. Defaulting to 'random'"
                )
            return self.data.sampleClusters(states, frames, replacement, allframes)

        if states is None:
            if statetype == "macro":
                states = range(self.macronum)
            elif statetype == "micro":
                states = range(self.micronum)
        if isinstance(states, int):
            states = [
                states,
            ]

        if allframes:
            logger.warning(
                "The allframes option will be deprecated. Use frames=None instead."
            )
            frames = None

        self._integrityCheck(postmsm=True)
        if statetype != "macro" and samplemode != "random":
            samplemode = "random"
            logger.warning(
                "'micro' states incompatible with 'samplemode' other than 'random'. Defaulting to 'random'"
            )

        if frames is None or isinstance(frames, int):
            frames = np.repeat(frames, len(states))

        stConcat = np.concatenate(self.data.St)
        absFrames = []
        relFrames = []
        for i in range(len(states)):
            if frames[i] == 0:
                absFrames.append(np.array([], dtype=int))
                relFrames.append(np.array([], dtype=int))
                continue

            st = states[i]
            if statetype == "macro":
                (selFr, selMicro) = _sampleMacro(
                    self, st, stConcat, samplemode, frames[i], replacement
                )
                absFrames.append(selFr)
            elif statetype == "micro":
                absFrames.append(
                    _sampleMicro(self, st, stConcat, frames[i], replacement)
                )
            else:
                raise NameError("No valid state type given (read documentation)")

            if len(absFrames[-1]) == 0:
                raise NameError(
                    "No frames could be sampled from {} state {}. State is empty.".format(
                        statetype, st
                    )
                )

            relFrames.append(self.data.abs2rel(absFrames[-1]))
        return absFrames, relFrames

    def eqDistribution(self, plot: bool = True, save: str | None = None) -> np.ndarray:
        """Obtain and plot the equilibrium probabilities of each macrostate

        Parameters
        ----------
        plot : bool, optional
            Disable plotting of the equilibrium distribution by setting it to False
        save : str
            Path of the file in which to save the figure

        Returns
        -------
        eq : np.ndarray
            An array of equilibrium probabilities of the macrostates

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.eqDistribution()
        """
        # logger.warning('Equilibrium distribution calculations for macrostates are now done using membership '
        #                'probabilities and hence your results might differ from analyses done before this change.')
        self._integrityCheck(postmsm=True)
        macroeq = np.ones(self.macronum) * -1
        macroindexes = list(set(self.coarsemsm.assignments))
        for i in range(self.macronum):
            # macroeq[i] = np.sum(self.msm.stationary_distribution[self.macro_ofmicro == i])
            macroeq[i] = np.sum(
                self.coarsemsm.memberships[:, macroindexes[i]]
                * self.msm.stationary_distribution
            )

        if plot or (save is not None):
            from matplotlib import pylab as plt

            if plot:
                plt.ion()
            plt.figure()
            plt.bar(np.arange(self.macronum) + 0.4, macroeq)
            plt.ylabel("Equilibrium probability")
            plt.xlabel("Macrostates")
            plt.xticks(np.arange(self.macronum) + 0.4, range(self.macronum))
            if save is not None:
                plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
            if plot:
                plt.show()

        return macroeq

    def _coarseP(self):
        M = self.coarsemsm.memberships
        Pcoarse = np.linalg.inv(M.T.dot(M)).dot(M.T).dot(self.P).dot(M)
        if len(np.where(Pcoarse < 0)[0]) != 0:
            raise NameError(
                "Cannot produce coarse P matrix. Ended up with negative probabilities. Try using less macrostates."
            )
        return Pcoarse

    def getStates(
        self,
        states: list | int | range | np.ndarray | None = None,
        statetype: str = "macro",
        wrapsel: str | np.ndarray | None = "protein",
        alignsel: str | np.ndarray | None = "name CA",
        alignmol: "Molecule | None" = None,
        samplemode: str = "weighted",
        numsamples: int = 50,
        simlist: np.ndarray | None = None,
    ) -> list:
        """Get samples of MSM states in Molecule classes

        Parameters
        ----------
        states : list or np.ndarray, optional
            A list of states to visualize
        statetype : str, optional
            The type of state to visualize. Can be 'macro', 'micro' or 'cluster'.
        wrapsel : str or np.ndarray, optional
            Atom selection used for wrapping, as an atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`).
            Set to None to disable wrapping.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignsel : str or np.ndarray, optional
            Atom selection used for aligning all frames, as an atom selection string, a boolean mask, or an integer
            index array. Set to None to disable aligning.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignmol : :class:`Molecule <moleculekit.molecule.Molecule>` object
            A reference molecule onto which to align all others
        samplemode : str, optional
            How to obtain the samples from the states. Can be 'weighted' or 'random'.
        numsamples : int
            Number of samples (conformations) for each state.
        simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
            Optionally pass a different (but matching, i.e. filtered) simlist for creating the Molecules.

        Returns
        -------
        mols : list of :class:`Molecule <moleculekit.molecule.Molecule>` objects
            A list of :class:`Molecule <moleculekit.molecule.Molecule>` objects containing the samples of each state

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> mols = model.getStates()
        >>> for m in mols:
        >>>     m.view()
        """
        from htmd.projections.metric import _singleMolfile

        self._integrityCheck(postmsm=(statetype != "cluster"))
        if simlist is None:
            simlist = self.data.simlist
        else:
            if len(simlist) != len(self.data.simlist):
                raise AttributeError(
                    "Provided simlist has different number of trajectories than the one used by the model."
                )

        (single, molfile) = _singleMolfile(simlist)
        if not single:
            raise NameError(
                "Visualizer does not support yet visualization of systems with different structure files. "
                "The simlist should be created with a single molfile (for example a filtered one)"
            )
        if alignmol is None:
            alignmol = Molecule(molfile)
        if statetype != "macro" and statetype != "micro" and statetype != "cluster":
            raise NameError("'statetype' must be either 'macro', 'micro' or ''cluster'")
        if states is None:
            if statetype == "macro":
                states = range(self.macronum)
            elif statetype == "micro":
                states = range(self.micronum)
            elif statetype == "cluster":
                states = range(self.data.K)
        if len(states) == 0:
            raise NameError("No " + statetype + " states exist in the model")

        (tmp, relframes) = self.sampleStates(
            states, numsamples, statetype=statetype, samplemode=samplemode
        )

        from tqdm import tqdm

        # This loop really iterates over states. sampleStates returns an array of arrays
        mols = []
        for rel in tqdm(relframes, desc="Getting state Molecules"):
            mols.append(
                self._loadMols(rel, molfile, wrapsel, alignsel, alignmol, simlist)
            )

        return mols

    def viewStates(
        self,
        states: list | int | range | np.ndarray | None = None,
        statetype: str = "macro",
        protein: bool | None = None,
        ligand: str | np.ndarray | None = None,
        mols: list | None = None,
        numsamples: int = 50,
        wrapsel: str | np.ndarray | None = "protein",
        alignsel: str | np.ndarray | None = "name CA",
        gui: bool = False,
        simlist: np.ndarray | None = None,
    ):
        """Visualize macro/micro/cluster states in VMD

        Parameters
        ----------
        states : list or np.ndarray, optional
            A list of states to visualize
        statetype : str, optional
            The type of state to visualize. Can be 'macro', 'micro' or 'cluster'.
        protein : bool, optional
            Set to True to enable pure protein system visualization
        ligand : str or np.ndarray, optional
            Atom selection for the ligand, as an atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`).
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        mols : list, optional
            A list of :class:`Molecule <moleculekit.molecule.Molecule>` objects to visualize
        numsamples : int
            Number of samples (conformations) for each state.
        wrapsel : str or np.ndarray, optional
            Atom selection used for wrapping, as an atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`).
            Set to None to disable wrapping.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignsel : str or np.ndarray, optional
            Atom selection used for aligning all frames, as an atom selection string, a boolean mask, or an integer
            index array. Set to None to disable aligning.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        gui : bool
            Set to True to enable the GUI in the NGL viewer.
        simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
            Optionally pass a different (but matching, i.e. filtered) simlist for visualizing the states.

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.viewStates(protein=True)

        >>> model.viewStates(ligand='resname MOL')
        """
        from moleculekit.config import _config

        self._integrityCheck(postmsm=(statetype != "cluster"))

        if _config["viewer"] is not None and _config["viewer"].lower() in (
            "ngl",
            "webgl",
        ):
            if protein is None and ligand is None:
                raise RuntimeError(
                    "You need to specify either the `protein` of `ligand` arguments for state visualization."
                )

            return self._viewStatesNGL(
                states, statetype, protein, ligand, mols, numsamples, gui=gui
            )

        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(
                states,
                statetype,
                numsamples=numsamples,
                wrapsel=wrapsel,
                alignsel=alignsel,
                simlist=simlist,
            )
        colors = [0, 1, 3, 4, 5, 6, 7, 9]
        for i, s in enumerate(states):
            mols[i].viewname = f"{statetype} {s}"
            mols[i].reps.add(
                sel="protein or nucleic",
                style="NewCartoon",
                frames=range(mols[i].numFrames) if ligand is None else None,
            )
            if ligand is not None:
                mols[i].reps.add(
                    sel=ligand,
                    style="Lines",
                    color=colors[np.mod(i, len(colors))],
                    frames=range(mols[i].numFrames),
                )
            mols[i].view()

    def _viewStatesNGL(
        self, states, statetype, protein, ligand, mols, numsamples, gui=False
    ):
        from moleculekit.util import sequenceID

        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(states, statetype, numsamples=min(numsamples, 15))
        hexcolors = {
            0: "#0000ff",
            1: "#ff0000",
            2: "#333333",
            3: "#ff6600",
            4: "#ffff00",
            5: "#4c4d00",
            6: "#b2b2cc",
            7: "#33cc33",
            8: "#ffffff",
            9: "#ff3399",
            10: "#33ccff",
        }
        if protein is None and ligand is None:
            raise NameError(
                'Please provide either the "protein" or "ligand" parameter for viewStates.'
            )
        k = 0
        from nglview import NGLWidget, HTMDTrajectory

        view = NGLWidget(gui=gui)
        ref = mols[0].copy()
        for i, s in enumerate(states):
            if protein:
                mol = Molecule()
            if ligand:
                mol = ref.copy()
                mol.remove(ligand, _logger=False)
                mol.dropFrames(keep=0)
                mols[i].filter(ligand, _logger=False)
            mols[i].set("chain", "{}".format(s))
            tmpcoo = mols[i].coords
            for j in range(mols[i].numFrames):
                mols[i].coords = np.atleast_3d(tmpcoo[:, :, j])
                if ligand:
                    mols[i].set("segid", sequenceID(mols[i].resid) + k)
                    k = int(mols[i].segid[-1])
                mol.append(mols[i])
            view.add_trajectory(HTMDTrajectory(mol))
            # Setting up representations
            if ligand:
                view[i].add_cartoon("protein", color="sstruc")
                view[i].add_hyperball(
                    ":{}".format(s), color=hexcolors[np.mod(i, len(hexcolors))]
                )
            if protein:
                view[i].add_cartoon("protein", color="residueindex")

        self._nglButtons(view, statetype, states)
        return view

    def _nglButtons(self, ngl_widget, statetype, states):
        # Adds buttons for enabling and disabling macrostate visualizations
        import ipywidgets
        from IPython.display import display

        container = []
        for s in states:
            w = ipywidgets.Checkbox(description="{} {}".format(statetype, s))
            w.value = True
            container.append(w)

        def updateReps(name):
            on = []
            for i, w in enumerate(container):
                if w.value:
                    on.append(i)
            ngl_widget.show_only(on)

        for w in container:
            w.on_trait_change(updateReps, "value")

        # container.append(ipywidgets.Checkbox(description="all"))

        hb = ipywidgets.HBox(container)
        display(hb)

    def save(self, filename: str):
        """Save a :class:`Model <htmd.model.Model>` object to disk

        Parameters
        ----------
        filename : str
            Path of the file in which to save the object

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.save('./model.dat')
        """
        import pickle

        # Temporarily store data object and replace with dicts
        tmpdata = self.data
        if self.data.parent is not None:
            tmpparentdata = self.data.parent
            self.data.parent = self.data.parent.__dict__
        self.data = self.data.__dict__

        # Dump the dict
        f = open(filename, "wb")
        pickle.dump(self.__dict__, f)
        f.close()

        # Restore data to classes
        self.data = tmpdata
        if self.data.parent is not None:
            self.data.parent = tmpparentdata

    def load(self, filename: str):
        """Load a :class:`MetricData <htmd.metricdata.MetricData>` object from disk

        Parameters
        ----------
        filename : str
            Path to the saved MetricData object

        Examples
        --------
        >>> model = Model()
        >>> model.load('./model.dat')
        """
        import sys
        import pickle
        from htmd.metricdata import MetricData

        try:
            import pandas.indexes
        except ImportError:
            import pandas.core.indexes

            sys.modules["pandas.indexes"] = (
                pandas.core.indexes
            )  # Hacky fix for new pandas version

        f = open(filename, "rb")
        z = pickle.load(f)
        f.close()
        for k in z:
            if k == "data":
                m = MetricData()
                m.load(z[k])
                self.__dict__[k] = m
            else:
                self.__dict__[k] = z[k]

    def copy(self) -> "Model":
        """Produces a deep copy of the object

        Returns
        -------
        data : :class:`Model` object
            A copy of the current object

        Examples
        --------
        >>> model2 = model.copy()
        """
        from copy import deepcopy

        return deepcopy(self)

    def _cktest_lags(
        self,
        lags: list | range | np.ndarray | None = None,
        maxlag: float | None = None,
        numlags: int = 25,
        units: str = "frames",
    ) -> list:
        """Compute the lag times (in frames) for the Chapman-Kolmogorov test.

        When `lags` is None the test is anchored on the production lag
        (``self.lag``) and the returned lags are the integer multiples
        ``tau, 2*tau, 3*tau, ...`` bounded by `maxlag` and `numlags`. deeptime
        uses the smallest-lag model as the test reference, so the production-lag
        model is the one being validated.
        """
        if lags is not None:
            return unitconvert(units, "frames", lags, fstep=self.data.fstep).tolist()

        if maxlag is None:
            from scipy import stats

            maxlag = int(stats.mode(self.data.trajLengths, keepdims=False).mode) - 1
        else:
            maxlag = unitconvert(units, "frames", maxlag, fstep=self.data.fstep)

        tau = int(self.lag)
        nmultiples = min(numlags, max(2, int(maxlag // tau)))
        return [tau * k for k in range(1, nmultiples + 1)]

    def cktest(
        self,
        lags: list | range | np.ndarray | None = None,
        maxlag: float | None = None,
        numlags: int = 25,
        units: str = "frames",
        plot: bool = True,
        save: str | None = None,
        errors: int | None = None,
    ):
        """Conducts a Chapman-Kolmogorov test.

        The Chapman-Kolmogorov test validates the model estimated at the
        production lag time tau (the lag passed to :meth:`markovModel`) by
        checking the identity ``T(k*tau) = T(tau)^k``. When ``lags`` is not
        given, the test is anchored on the production lag and the test lag times
        are the integer multiples ``tau, 2*tau, 3*tau, ...``. The production-lag
        model is used as the reference, so the model you actually built is the
        one being validated.

        Parameters
        ----------
        lags : list or range or np.ndarray, optional
            Explicit list of lag times to test. If given, the smallest lag is
            used as the reference model, so it should be the production lag.
        maxlag : float, optional
            Largest lag time to test. Used if `lags` is None. By default it will
            use the mode of the trajectory lengths. The number of multiples is
            chosen so that ``k*tau`` does not exceed `maxlag`.
        numlags : int, optional
            Upper bound on the number of multiples of the production lag to test.
            Used if `lags` is None.
        units : str, optional
            Units of the lag times. By default it will use frames.
        plot : bool
            If the method should display the plot of the CK test
        save : str
            Path of the file in which to save the figure
        errors : int, optional
            Number of Bayesian samples to use for the error bars. If set to None,
            it will not plot error bars and use a Maximum Likelihood MSM.
        """
        from deeptime.plots.chapman_kolmogorov import plot_ck_test
        from matplotlib import pylab as plt

        self._integrityCheck(postmsm=True)

        lags = self._cktest_lags(
            lags=lags, maxlag=maxlag, numlags=numlags, units=units
        )

        statelist = [traj.cluster for traj in self.data.trajectories]
        models = []
        for lag in lags:
            models.append(self._get_model(statelist, lag, bayesian_samples=errors))

        res = models[0].ck_test(models, n_metastable_sets=self.macronum)
        plot_ck_test(res, legend=True)

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
        if plot:
            plt.show()

    def createCoreSetModel(self, threshold: float = 0.5) -> tuple:
        """Given an MSM this function detects the states belonging to a core set and returns a new model consisting
        only of these states.

        Parameters
        ----------
        threshold : float
            Membership probability threshold over which microstates belong to the core of a macrostate

        Returns
        -------
        newmodel : :class:`Model <htmd.model.Model>` object
            A new model object
        frames : list
            A list of the frames that were kept in the new model
        """
        if (threshold >= 1) or (threshold <= 0):
            raise AttributeError(
                "threshold argument only accepts values between (0, 1)"
            )

        def calcCoreSet(distr, assign, threshold):
            coreset = []
            for i, md in enumerate(distr):
                microofmacro = np.where(assign == i)[0]
                prob = md[microofmacro]
                tt = threshold * (prob.max() - prob.min())
                coreset.append(microofmacro[np.where(prob > tt)[0]])
            return coreset

        def coreDtraj(data, micro_ofcluster, coreset):
            corestates = np.concatenate(coreset)
            newmapping = np.ones(corestates.max() + 1, dtype=int) * -1
            newmapping[corestates] = np.arange(len(corestates))
            discretetraj = [st.copy() for st in data.St]
            frames = []
            newdiscretetraj = []
            newcounts = np.zeros(len(corestates))
            for t, st in enumerate(discretetraj):
                oldmicro = None
                newtraj = []
                tframes = []
                for f, cl in enumerate(st):
                    newmicro = None
                    micro = micro_ofcluster[cl]
                    if micro == -1:  # If we are in an dropped cluster, keep old index
                        newmicro = oldmicro
                    else:
                        for co in coreset:
                            if micro in co:
                                newmicro = micro
                                oldmicro = micro
                                break
                    if newmicro is None and oldmicro is not None:
                        newtraj.append(oldmicro)
                        tframes.append(f)
                    elif newmicro is not None:
                        newtraj.append(newmicro)
                        tframes.append(f)
                mappedtraj = newmapping[np.array(newtraj, dtype=int)]
                newdiscretetraj.append(mappedtraj)
                if len(mappedtraj):
                    newcounts[: mappedtraj.max() + 1] += np.bincount(mappedtraj)
                frames.append(np.array(tframes))
            # kept = np.array([i for i, x in enumerate(newdiscretetraj) if len(x) != 0])
            return (
                newdiscretetraj,
                len(corestates),
                newcounts,
                frames,
            )

        coreset = calcCoreSet(
            self.coarsemsm.metastable_distributions,
            self.coarsemsm.assignments,
            threshold,
        )
        newdata = self.data.copy()
        newSt, newdata.K, newdata.N, frames = coreDtraj(
            self.data, self.micro_ofcluster, coreset
        )
        for i, (s, fr) in enumerate(zip(newSt, frames)):
            if len(s):
                newdata.trajectories[i].cluster[fr] = s

        logger.info(
            f"Kept {[len(x) for x in coreset]} microstates from each macrostate."
        )

        dataobjects = [newdata]
        if newdata.parent is not None:
            dataobjects.append(newdata.parent)
        for data in dataobjects:
            dat = []
            ref = []
            simstmp = []
            cluster = []
            for i, fr in enumerate(frames):
                if len(fr):
                    dat.append(data.trajectories[i].projection[fr, :])
                    ref.append(data.trajectories[i].reference[fr, :])
                    if len(data.trajectories[i].cluster):
                        cluster.append(data.trajectories[i].cluster[fr])
                    simstmp.append(data.trajectories[i].sim)

            data._loadTrajectories(dat, ref, simstmp, cluster if len(cluster) else None)

        return Model(newdata), frames

    def plotFES(
        self,
        dimX: int,
        dimY: int,
        temperature: float,
        states: bool = False,
        s: float = 10,
        cmap=None,
        fescmap=None,
        statescmap=None,
        plot: bool = True,
        save: str | None = None,
        data: "MetricData | None" = None,
        levels: int = 7,
    ):
        """Plots the free energy surface on any given two dimensions. Can also plot positions of states on top.

        Parameters
        ----------
        dimX : int
            Index of projected dimension to use for the X axis.
        dimY : int
            Index of projected dimension to use for the Y axis.
        temperature : float
            Simulation temperature.
        states : bool
            If True, will plot scatter plot of microstates coloured by macro state on top of FES.
        s : float
            Marker size for states.
        cmap : str or matplotlib colormap
            Sets the Matplotlib colormap for both `fescmap` and `statescmap`
        fescmap : str or matplotlib colormap
            Matplotlib colormap for the free energy surface
        statescmap : str or matplotlib colormap
            Matplotlib colormap for the states
        plot : bool
            If the method should display the plot of the FES. If both plot=False and save=None, the method will return the figure and axes.
        save : str
            Path of the file in which to save the figure. If both plot=False and save=None, the method will return the figure and axes.
        data : :class:`MetricData` object
            Optionally you can pass a different MetricData object than the one used to build the model. For example
            if the user wants to build a model on distances but wants to plot the FES on top of RMSD values. The
            MetricData object needs to have the same simlist as the Model.
        levels : int
            Number of contour levels to use for the free energy surface.

        Returns
        -------
        f : matplotlib.figure.Figure
            The matplotlib figure object
        ax : matplotlib.axes.Axes
            The matplotlib axes object

        Examples
        --------
        >>> import matplotlib as plt
        >>> model.plotFES(0, 1, 300)
        >>> model.plotFES(2, 3, 300, data=otherdata, states=True, fescmap=plt.cm.gray)
        """
        self._integrityCheck(postmsm=True)
        from matplotlib import pylab as plt
        from htmd.kinetics import Kinetics

        if data is None:
            data = self.data
            microcenters = self.data.Centers[self.cluster_ofmicro, :]
        else:
            if self.data.numFrames != data.numFrames or ~np.all(
                [s1 == s2 for s1, s2 in zip(self.data.simlist, data.simlist)]
            ):
                raise RuntimeError(
                    "The data argument you provided uses a different simlist than the Model."
                )
            microcenters = np.vstack(
                getStateStatistic(self, data, range(self.micronum), statetype="micro")
            )

        if fescmap is None:
            fescmap = "viridis"
        if statescmap is None:
            statescmap = plt.cm.jet
        if cmap is not None:
            fescmap = cmap
            statescmap = cmap

        xlabel = f"Dimension {dimX}"
        if data.description is not None:
            xlabel = data.description.description[dimX]

        ylabel = f"Dimension {dimY}"
        if data.description is not None:
            ylabel = data.description.description[dimY]

        title = "Free energy surface"

        counts, xbins, ybins = data._getFEShistogramCounts(
            dimX,
            dimY,
            nbins=80,
            pad=0.5,
            micro_ofcluster=self.micro_ofcluster,
            stationary_distribution=self.msm.stationary_distribution,
            St=self.data.St,
        )
        counts /= counts.sum()  # Normalize probabilites
        nonzero = counts != 0
        energy = counts.copy()
        energy[nonzero] = -Kinetics._kB * temperature * np.log(counts[nonzero])
        energy[~nonzero] = energy[nonzero].max() + 100
        f, ax, cf = data._contourPlot(
            energy,
            xbins,
            ybins,
            levels=levels,
            nonzero=nonzero,
            cmap=fescmap,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )

        data._setColorbar(f, cf, "kcal/mol", scientific=False)
        if states:
            colors = statescmap(np.linspace(0, 1, self.macronum))
            for m in range(self.macronum):
                macromicro = microcenters[self.macro_ofmicro == m, :]
                _ = ax.scatter(
                    macromicro[:, dimX],
                    macromicro[:, dimY],
                    s=s,
                    color=colors[m],
                    label=f"Macro {m}",
                    edgecolors="none",
                )
            ax.legend(prop={"size": 8})

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
            plt.close()
        elif plot:
            plt.show()
        else:
            return f, ax

    def _integrityCheck(self, postmsm=False, markov=False):
        if postmsm and self._modelid is None:
            raise NameError("You need to call markovModel before calling this command")
        if (
            not markov
            and self.data._dataid == self.data._clusterid
            and self.data._dataid != self._clusterid
        ):
            raise NameError(
                "After updating the MetricData object you need to call the markovModel command anew."
            )
        if self.data._dataid != self.data._clusterid:
            raise NameError(
                "After modifying the data in the MetricData object you need to recluster and reconstruct the markov model."
            )

    def _loadMols(self, rel, molfile, wrapsel, alignsel, refmol, simlist):
        frames = self.data.rel2sim(rel, simlist=simlist)
        mol = Molecule(molfile)
        trajs = []
        frs = []
        for f in frames:
            trajs.append(f.sim.trajectory[f.piece])
            frs.append(f.frame)
        mol.read(np.array(trajs), frames=np.array(frs))
        if wrapsel is not None and len(wrapsel):
            mol.wrap(wrapsel)
        if alignsel is not None and len(alignsel):
            mol.align(alignsel, refmol=refmol)
        return mol


def getStateStatistic(
    reference: "Model | MetricData",
    data: "MetricData",
    states: list | range | np.ndarray,
    statetype: str = "macro",
    weighted: bool = False,
    method=np.mean,
    axis: int | None = 0,
) -> list:
    """Calculates properties of the states.

    Calculates properties of data corresponding to states. Can calculate for example the mean distances of atoms in a
    state, or the standard deviation of the RMSD in a state.

    Parameters
    ----------
    reference : :class:`Model` object or :class:`MetricData` object
        A model containing the state definitions or a MetricData object containing cluster definitions
    data : :class:`MetricData` object
        A projection corresponding to the conformations in the states of the model. The data and the model need to share
        the same `simlist`
    states : list or np.ndarray
        A list of the states for which we want to calculate the properties
    statetype : str
        The state type. Can be 'macro', 'micro' or 'cluster'.
    weighted : bool
        If the properties of the macrostates should be calculated as the weighted average of their microstates, where the
        weights are the equilibrium probabilities of the microstates
    method :
        A function pointer for the method that calculates our desired property. e.g. np.mean, np.std
    axis : int
        On which axis of the data the method should operate on

    Returns
    -------
    statistic : list
        A list which contains in each element the desired statistic of a specific state specifies in `states`

    Examples
    --------
    >>> data = MetricDistance.project(sims, 'protein and name CA', 'resname MOL')
    >>> model = Model(data)
    >>> model.markovModel(100, 5)
    >>> # Get the standard deviation of distances in all macrostates
    >>> getStateStatistic(model, data, list(range(5)), method=np.std)
    """
    from htmd.metricdata import MetricData

    if axis != 0:
        logger.warning("Axis different than 0 might not work correctly yet")

    if isinstance(reference, Model):
        refdata = reference.data
    elif isinstance(reference, MetricData):
        refdata = reference
        if statetype != "cluster":
            raise RuntimeError(
                "You can only use statetype cluster with reference MetricData object. To use other statetypes build a Model."
            )
    else:
        raise RuntimeError("Invalid argument type")

    if refdata.numTrajectories > 0 and np.any(refdata.trajLengths != data.trajLengths):
        raise NameError(
            "Data trajectories need to match in size and number to the trajectories in the model"
        )
    stconcat = np.concatenate(refdata.St)
    datconcat = np.concatenate(data.dat)

    statistic = []
    for i, st in enumerate(states):
        if statetype == "macro":
            frames = reference.macro_ofcluster[stconcat] == st
        elif statetype == "micro":
            frames = reference.micro_ofcluster[stconcat] == st
        elif statetype == "cluster":
            frames = stconcat == st
        else:
            raise NameError("No valid state type given (read documentation)")

        if statetype == "macro" and weighted:
            statistic.append(
                _weightedMethod(reference, method, stconcat, datconcat, st, axis)
            )
        else:
            if axis is None:
                statistic.append(method(datconcat[frames, ...]))
            else:
                statistic.append(method(datconcat[frames, ...], axis=axis))
    return statistic


def _weightedMethod(model, method, stconcat, datconcat, st, axis):
    microsofmacro = np.where(model.macro_ofmicro == st)[0]
    eq = model.msm.stationary_distribution
    weights = eq / np.sum(eq[microsofmacro])
    avgstatistic = np.zeros(np.size(datconcat, 1))
    for m in microsofmacro:
        frames = model.micro_ofcluster[stconcat] == m
        if axis is None:
            stat = method(datconcat[frames, :])
        else:
            stat = method(datconcat[frames, :], axis=axis)
        avgstatistic = avgstatistic + stat * weights[m]
    return avgstatistic


def macroAccumulate(model: "Model", microvalue: np.ndarray) -> np.ndarray:
    """Accumulate values of macrostates from a microstate array

    Parameters
    ----------
    model : :class:`Model <htmd.model.Model>` object
        The model which to use to accumulate
    microvalue : np.ndarray
        An array of values corresponding to the microstates of the model

    Returns
    -------
    macrovalue : np.ndarray
        An array of values corresponding to the macrostates of the model
    """
    res = np.zeros(model.macronum)
    for i in range(len(microvalue)):
        macro = model.macro_ofmicro[i]
        res[macro] = res[macro] + microvalue[i]
    return res


def _sampleMacro(obj, macro, stConcat, mode, numFrames, replacement):
    from htmd.metricdata import _randomSample

    if mode == "random":
        frames = np.where(obj.macro_ofcluster[stConcat] == macro)[0]
        selFrames = _randomSample(frames, numFrames, replacement)
        selMicro = obj.micro_ofcluster[stConcat[selFrames]]
    elif mode == "even":
        micros = np.where(obj.macro_ofmicro == macro)[0]
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames.append(
                _sampleMicro(obj, micros[i], stConcat, numFrames, replacement)
            )
            selMicro.append(np.ones(numFrames) * micros[i])
    elif mode == "weighted":
        micros = np.where(obj.macro_ofmicro == macro)[0]
        eq = obj.msm.stationary_distribution
        weights = eq[micros] / np.sum(eq[micros])
        framespermicro = np.random.multinomial(numFrames, weights)
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames = np.append(
                selFrames,
                _sampleMicro(obj, micros[i], stConcat, framespermicro[i], replacement),
            )
            selMicro = np.append(selMicro, [micros[i]] * framespermicro[i])
    elif mode == "weightedTrunc":
        micros = np.where(obj.macro_ofmicro == macro)[0]
        eq = obj.msm.stationary_distribution
        weights = eq[micros] / np.sum(eq[micros])
        idx = np.argsort(weights)
        cs = np.cumsum(weights[idx])
        under50 = cs < 0.5
        weights[idx[under50]] = 0
        weights /= np.sum(weights)
        framespermicro = np.random.multinomial(numFrames, weights)
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames = np.append(
                selFrames,
                _sampleMicro(obj, micros[i], stConcat, framespermicro[i], replacement),
            )
            selMicro = np.append(selMicro, [micros[i]] * framespermicro[i])
    else:
        raise NameError("No valid mode given (read documentation)")
    return selFrames, selMicro


def _sampleMicro(obj, micro, stConcat, numFrames, replacement):
    from htmd.metricdata import _randomSample

    frames = np.where(obj.micro_ofcluster[stConcat] == micro)[0]
    return _randomSample(frames, numFrames, replacement)


def _macroTrajectoriesReport(macronum, macrost, simlist=None):
    """
    Prints out the number of trajectories that visited each macrostate and prints warnings

    Parameters
    ----------
    macronum
    macrost
    simlist
    """
    macrotraj = dict()
    macrotrajnum = np.zeros(macronum, dtype=int)
    for i, st in enumerate(macrost):
        macros = np.setdiff1d(np.unique(st), [-1])
        for m in macros:
            if m not in macrotraj:
                macrotraj[m] = []
            macrotraj[m].append(i)
            macrotrajnum[m] += 1
    logger.info("Number of trajectories that visited each macrostate:")
    logger.info(macrotrajnum)

    for m in range(macronum):
        ratio = macrotrajnum[m] / len(macrost)
        if macrotrajnum[m] <= 3 and ratio <= 0.2:
            logger.info(
                "Take care! Macro {} has been visited only in {} trajectories"
                " ({:.1f}% of total):".format(m, macrotrajnum[m], ratio * 100)
            )
            if simlist is not None:
                for s in macrotraj[m]:
                    logger.info(simlist[s])


def _macroTrajSt(St, macro_ofcluster):
    mst = []
    for i in range(len(St)):
        mst.append(macro_ofcluster[St[i]])
    return mst


"""def _macroP(C, macro_ofmicro):
    macronum = np.max(macro_ofmicro) + 1
    macroC = np.zeros((macronum, macronum))

    for m1 in range(macronum):
        sourcemicros = np.where(macro_ofmicro == m1)[0]
        for m2 in range(macronum):
            sinkmicros = np.where(macro_ofmicro == m2)[0]
            ixgrid = np.ix_(sourcemicros, sinkmicros)
            macroC[m1, m2] = np.sum(C[ixgrid].flatten())

    from deeptime.estimation import transition_matrix
    return transition_matrix(macroC, reversible=True)"""


