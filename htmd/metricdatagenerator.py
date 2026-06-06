from typing import TYPE_CHECKING
import numpy as np
from htmd.metricdata import MetricData

if TYPE_CHECKING:
    from htmd.model import Model


class MetricDataGenerator:
    """Generate synthetic trajectories from existing projected data.

    Given a clustered :class:`MetricData <htmd.metricdata.MetricData>` object (and
    optionally a Markov state model), this class produces new synthetic trajectories
    by resampling frames from the clusters of the original data. The various
    ``newTrajectories*`` methods implement different resampling strategies.

    Parameters
    ----------
    fulldata : :class:`MetricData <htmd.metricdata.MetricData>` object
        The clustered MetricData object from which to sample frames.
    model : :class:`Model <htmd.model.Model>` object, optional
        A Markov state model built on `fulldata`. Required for the MSM-based
        sampling strategies.
    is_adaptive : bool
        If True, starting frames are drawn from the first epoch of an adaptive
        sampling run.
    """

    def __init__(
        self,
        fulldata: MetricData,
        model: "Model | None" = None,
        is_adaptive: bool = False,
    ):
        self.fulldata = fulldata
        if model:
            self.micronum = model.micronum
            self.cluster2micro = model.micro_ofcluster
            self.micro2cluster = model.cluster_ofmicro
            self.P = model.P
        self.reference = []
        self.stconcat = np.concatenate(fulldata.St)
        self.lengths = fulldata.trajLengths
        self.is_adaptive = is_adaptive

    @staticmethod
    def _collectTrajectory(datasource, traj):
        dat = []
        ref = []
        for t in traj:
            dat.append(datasource.trajectories[t[0]].projection[t[1]])
            ref.append(datasource.trajectories[t[0]].reference[t[1]])
        dat = np.array(dat, dtype=datasource.trajectories[0].projection.dtype)
        ref = np.array(ref, dtype=datasource.trajectories[0].reference.dtype)
        return dat, ref, datasource.trajectories[t[0]].sim

    def _convertRelFrames(self, relFrames):
        rel = []
        for r in relFrames:
            rel.append(self.reference[r[0]][r[1], :])
        return rel

    def _pickFromCluster(self, relFrames, simlen, maintainlen=True):
        rel = []
        for r in relFrames:
            clu = self.fulldata.trajectories[r[0]].cluster[r[1]]
            clufr = self.fulldata.abs2rel(np.where(self.stconcat == clu)[0])
            if maintainlen:
                # Basically selects those frames from the cluster where there is still trajectory enough to fill simlen
                longenough = np.where(
                    (self.lengths[clufr[:, 0]] - clufr[:, 1]) > simlen
                )[0]
                if len(longenough) < len(relFrames):
                    pass

                if len(longenough) == 0:
                    raise RuntimeError(
                        f"No traj found in cluster {clu} of length {simlen}"
                    )
                clufr = clufr[longenough, :]
                rel.append(clufr[np.random.randint(clufr.shape[0]), :])
            else:
                sel_clu = clufr[np.random.randint(clufr.shape[0]), :]
                len_clu = self.lengths[sel_clu[0]] - sel_clu[1]
                rel.append((sel_clu, len_clu))
        return rel

    def _startingFrames(
        self,
        ntraj,
        startFrames=None,
        simlen=None,
        maintainlen=True,
    ):
        if startFrames is None:
            if self.is_adaptive:
                from htmd.adaptive.adaptive import epochSimIndexes

                idx = epochSimIndexes(self.fulldata.simlist)
                startFrames = []
                # When starting, it picks random starting nmax trajs from epoch 1's trajectories. ntraj == nmax
                for i in np.random.choice(idx[1], ntraj, replace=False):
                    startFrames.append([i, 0])
            else:
                traj = np.random.randint(len(self.fulldata.trajectories))
                nframes = self.fulldata.trajectories[traj].reference.shape[0]
                frame = np.random.randint(nframes)
                startFrames = [[traj, frame]]
        else:
            startFrames = self._convertRelFrames(startFrames)
            startFrames = self._pickFromCluster(
                startFrames, simlen, maintainlen=maintainlen
            )
        return startFrames

    def newTrajectoriesSimple(
        self,
        simlen: int,
        ntraj: int,
        startFrames: list | np.ndarray | None = None,
    ) -> list:
        """Generate trajectories by sampling whole pieces from a cluster.

        For each respawning conformation, selects a random trajectory from the
        conformations in its cluster.

        Parameters
        ----------
        simlen : int
            The length (in frames) of each new trajectory.
        ntraj : int
            The number of new trajectories to generate.
        startFrames : list or np.ndarray, optional
            Starting frames as trajectory index-frame pairs. If None, starting frames
            are picked automatically.

        Returns
        -------
        ret : list of np.ndarray
            One array per new trajectory, each row a trajectory index-frame pair.
        """
        startFrames = self._startingFrames(ntraj, startFrames, simlen, maintainlen=True)
        ret = []
        for r in startFrames:
            traj = np.tile(r, (simlen, 1))
            traj[:, 1] = np.arange(r[1], r[1] + simlen)
            # traj ends being a list of simidx,frame ; starting from [idx,0] until [idx, simlen]
            ret.append(traj)
        self.reference += ret
        return ret

    def newTrajectoriesFiller(
        self,
        simlen: int,
        ntraj: int,
        startFrames: list | np.ndarray | None = None,
    ) -> list:
        """Generate trajectories by chaining cluster pieces until the length is reached.

        Starting from a frame, it appends pieces sampled from the corresponding clusters
        until each trajectory reaches `simlen` frames.

        Parameters
        ----------
        simlen : int
            The length (in frames) of each new trajectory.
        ntraj : int
            The number of new trajectories to generate.
        startFrames : list or np.ndarray, optional
            Starting frames as trajectory index-frame pairs. If None, starting frames
            are picked automatically.

        Returns
        -------
        ret : list of np.ndarray
            One array per new trajectory, each row a trajectory index-frame pair.
        """
        startFrames = self._startingFrames(
            ntraj, startFrames, simlen, maintainlen=False
        )
        ret = []
        for rr, ll in startFrames:
            traj = np.tile(rr, (simlen, 1))
            tref = 0
            while simlen > tref:
                if tref + ll > simlen:
                    traj[tref:, 0] = rr[0]
                    traj[tref:, 1] = np.arange(rr[1], rr[1] + (simlen - tref))
                else:
                    traj[tref : tref + ll, 0] = rr[0]
                    traj[tref : tref + ll, 1] = np.arange(rr[1], rr[1] + ll)
                rr, ll = self._pickFromCluster(
                    [traj[tref]],
                    simlen,
                    maintainlen=False,
                )[0]
                tref += ll
            # traj ends being a list of simidx,frame ; starting from [idx,0] until [idx, simlen]
            ret.append(traj)
        self.reference += ret
        return ret

    def newTrajectoriesClusterJumping(
        self,
        simlen: int,
        ntraj: int,
        startFrames: list | np.ndarray | None = None,
        jumpprob: float = 0.1,
    ) -> list:
        """Generate trajectories by advancing one frame at a time with random cluster jumps.

        At each step the trajectory advances one frame ahead. With probability
        `jumpprob` (or when the end of the current trajectory is reached) it jumps to a
        random frame of the same cluster to continue sampling.

        Parameters
        ----------
        simlen : int
            The length (in frames) of each new trajectory.
        ntraj : int
            The number of new trajectories to generate.
        startFrames : list or np.ndarray, optional
            Starting frames as trajectory index-frame pairs. If None, starting frames
            are picked automatically.
        jumpprob : float
            The per-frame probability of jumping to a new frame within the same cluster.

        Returns
        -------
        ret : list of np.ndarray
            One array per new trajectory, each row a trajectory index-frame pair.
        """
        startFrames = self._startingFrames(ntraj, startFrames, 2)
        ret = []
        for r in startFrames:
            traj = np.zeros((simlen, 2), dtype=int)
            traj[0, :] = r
            for i in range(1, simlen):
                traj[i, :] = traj[i - 1, :] + [0, 1]  # Go to next frame
                # change by random chance or if we are at the last frame of the current trajectory
                if (np.random.rand() < jumpprob) or (
                    traj[i, 1] >= (self.lengths[traj[i, 0]] - 1)
                ):
                    traj[i, :] = self._pickFromCluster(
                        [traj[i, :]], 2, maintainlen=True
                    )[0]
            ret.append(traj)
        self.reference += ret
        return ret

    def newTrajectoriesMSM(
        self,
        simlen: int,
        ntraj: int,
        startFrames: list | np.ndarray | None = None,
    ) -> list:
        """Generate new synthetic trajectories sampled from the Markov state model.

        At each step the next microstate is drawn from the transition probability matrix
        of the model, and a random frame of the corresponding cluster is selected.

        Parameters
        ----------
        simlen : int
            The length (in frames) of each new trajectory.
        ntraj : int
            The number of new trajectories to generate.
        startFrames : list or np.ndarray, optional
            Starting frames as trajectory index-frame pairs. If None, starting frames
            are picked automatically.

        Returns
        -------
        ret : list of np.ndarray
            One array per new trajectory, each row a trajectory index-frame pair.
        """
        if startFrames is None:
            startFrames = self._startingFrames(ntraj, startFrames, simlen)
        else:
            startFrames = self._convertRelFrames(startFrames)

        ret = []
        for r in startFrames:
            traj = np.zeros((simlen, 2), dtype=int)
            traj[0, :] = r
            for i in range(1, simlen):
                micro = self.cluster2micro[
                    self.fulldata.trajectories[traj[i - 1, 0]].cluster[traj[i - 1, 1]]
                ]
                micro_tp = self.P[micro, :]
                if micro_tp.getformat() == "csr":
                    newmicro = np.random.choice(
                        self.micronum, 1, p=np.array(micro_tp.todense()).flatten()
                    )
                    while newmicro not in self.stconcat:
                        newmicro = np.random.choice(
                            self.micronum, 1, p=np.array(micro_tp.todense()).flatten()
                        )
                else:
                    newmicro = np.random.choice(self.micronum, 1, p=micro_tp)
                newclu = self.micro2cluster[newmicro]
                clufr = self.fulldata.abs2rel(np.where(self.stconcat == newclu)[0])
                # TODO: Use sampleStates
                traj[i, :] = clufr[np.random.randint(clufr.shape[0]), :]
            ret.append(traj)
        self.reference += ret
        return ret

    def newMetricData(
        self,
        datasource: MetricData,
        trajectories: list | None = None,
        olddata: "MetricData | None" = None,
    ) -> MetricData:
        """Convert generated trajectory indexes into a new MetricData object.

        Parameters
        ----------
        datasource : :class:`MetricData <htmd.metricdata.MetricData>` object
            The MetricData object from which to collect projections and references for
            the sampled frames.
        trajectories : list, optional
            A list of trajectories, each given as trajectory index-frame pairs (such as
            those returned by the ``newTrajectories*`` methods).
        olddata : :class:`MetricData <htmd.metricdata.MetricData>` object, optional
            If given, the new data is appended to this object and the merged object is
            returned.

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            A MetricData object containing the generated trajectories.
        """
        dat = []
        ref = []
        sim = []
        for traj in trajectories:
            d, r, s = self._collectTrajectory(datasource, traj)
            dat.append(d)
            ref.append(r)
            sim.append(s)

        newdata = MetricData(
            dat=dat,
            ref=ref,
            simlist=sim,
            description=datasource.description,
            fstep=datasource.fstep,
        )
        if olddata is not None:
            olddata.append(newdata)  # Merge with old data
            return olddata
        else:
            return newdata

    def parallelTest(
        self,
        simlen: int,
        ntraj: int,
        startFrames: list | np.ndarray | None = None,
    ) -> list:
        """Generate MSM-sampled trajectories in parallel across multiple processes.

        Parameters
        ----------
        simlen : int
            The length (in frames) of each new trajectory.
        ntraj : int
            The number of new trajectories to generate.
        startFrames : list or np.ndarray, optional
            Starting frames as trajectory index-frame pairs. If None, starting frames
            are picked automatically.

        Returns
        -------
        ret : list of np.ndarray
            One array per new trajectory, each row a trajectory index-frame pair.
        """
        if startFrames is None:
            startFrames = self._startingFrames(ntraj, startFrames, simlen)
        else:
            startFrames = self._convertRelFrames(startFrames)

        from joblib import delayed
        from htmd.parallelprogress import ParallelExecutor

        aprun = ParallelExecutor(n_jobs=-2)
        ret = aprun(total=len(startFrames))(
            delayed(_pickFromMicro)(
                relFrame,
                simlen,
                np.cumsum(self.fulldata.trajLengths),
                self.fulldata.trajectories,
                self.cluster2micro,
                self.micro2cluster,
                self.micronum,
                self.P,
                self.stconcat,
            )
            for relFrame in startFrames
        )
        self.reference += ret
        return ret


def _pickFromMicro(
    r,
    simlen,
    trajLengths,
    trajectories,
    cluster2micro,
    micro2cluster,
    micronum,
    P,
    stconcat,
):

    traj = np.zeros((simlen, 2), dtype=int)
    traj[0, :] = r
    for i in range(1, simlen):
        micro = cluster2micro[trajectories[traj[i - 1, 0]].cluster[traj[i - 1, 1]]]
        micro_tp = P[micro, :]
        newmicro = np.random.choice(micronum, 1, p=micro_tp)
        newclu = micro2cluster[newmicro]
        clufr = abs2rel(np.where(stconcat == newclu)[0], trajLengths)
        traj[i, :] = clufr[np.random.randint(clufr.shape[0]), :]
    return traj


def abs2rel(
    absFrames: int | list | np.ndarray, trajLengths: np.ndarray
) -> np.ndarray:  # trajLengths need to be summed up (np.cumsum())
    """Convert absolute frame indexes into trajectory index-frame pairs.

    Parameters
    ----------
    absFrames : int or list or np.ndarray
        An absolute frame index or a list of absolute frame indexes.
    trajLengths : np.ndarray
        The cumulative sum of the trajectory lengths (i.e. ``np.cumsum(trajLengths)``).

    Returns
    -------
    relframe : np.ndarray
        An array where each row is a trajectory index-frame pair.
    """
    if not hasattr(absFrames, "__len__"):
        absFrames = [absFrames]
    endFrames = np.append(0, trajLengths)

    relframe = np.zeros((len(absFrames), 2), dtype=int)
    for i in range(len(absFrames)):
        trajIdx = np.where(absFrames[i] < endFrames)[0][0] - 1
        trajFr = absFrames[i] - endFrames[trajIdx]
        relframe[i, :] = [trajIdx, trajFr]
    return relframe
