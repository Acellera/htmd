import numpy as np
from htmd.metricdata import MetricData


class MetricDataGenerator:
    def __init__(self, fulldata, model=None):
        self.fulldata = fulldata
        if model:
            self.micronum = model.micronum
            self.cluster2micro = model.micro_ofcluster
            self.micro2cluster = model.cluster_ofmicro
            self.P = model.P
        self.reference = []
        self.stconcat = np.concatenate(fulldata.St)
        self.lengths = fulldata.trajLengths

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

    def _startingFrames(self, ntraj, startFrames=None, simlen=None, maintainlen=True):
        if startFrames is None:
            from htmd.adaptive.adaptive import epochSimIndexes

            idx = epochSimIndexes(self.fulldata.simlist)
            startFrames = []
            # When starting, it picks random starting nmax trajs from epoch 1's trajectories. ntraj == nmax
            for i in np.random.choice(idx[1], ntraj, replace=False):
                startFrames.append([i, 0])
        else:
            startFrames = self._convertRelFrames(startFrames)
            startFrames = self._pickFromCluster(
                startFrames, simlen, maintainlen=maintainlen
            )
        return startFrames

    def newTrajectoriesSimple(self, simlen, ntraj, startFrames=None):
        """TrajectoriesSimple selects a random trajectory from the conformations in the cluster of the respawning conformations"""
        startFrames = self._startingFrames(ntraj, startFrames, simlen, maintainlen=True)
        ret = []
        for r in startFrames:
            traj = np.tile(r, (simlen, 1))
            traj[:, 1] = np.arange(r[1], r[1] + simlen)
            # traj ends being a list of simidx,frame ; starting from [idx,0] until [idx, simlen]
            ret.append(traj)
        self.reference += ret
        return ret

    def newTrajectoriesFiller(self, simlen, ntraj, startFrames=None):
        """"""
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
        self, simlen, ntraj, startFrames=None, jumpprob=0.1
    ):
        """clusterJumping only jumps one frame ahead and uses random chance to change the cluster from where to obtain new frames"""
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

    def newTrajectoriesMSM(self, simlen, ntraj, startFrames=None):
        """Generates new synthetic (fake) trajectories sampled from the Markov State Model"""
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

    def newMetricData(self, datasource, trajectories=None, olddata=None):
        """Converts trajectory indexes to a new MetricData object"""
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

    def parallelTest(self, simlen, ntraj, startFrames=None):
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


def abs2rel(absFrames, trajLengths):  # trajLengths need to be summed up (np.cumsum())
    if not hasattr(absFrames, "__len__"):
        absFrames = [absFrames]
    endFrames = np.append(0, trajLengths)

    relframe = np.zeros((len(absFrames), 2), dtype=int)
    for i in range(len(absFrames)):
        trajIdx = np.where(absFrames[i] < endFrames)[0][0] - 1
        trajFr = absFrames[i] - endFrames[trajIdx]
        relframe[i, :] = [trajIdx, trajFr]
    return relframe
