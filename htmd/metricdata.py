# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from copy import deepcopy
import uuid
import h5py
import pickle
import unittest
import logging

logger = logging.getLogger(__name__)


def _getsizes(x):
    if x is not None:
        if isinstance(x, np.ndarray):
            return x.shape[0]
        else:
            return len(x)


class Trajectory(object):
    """Class used for storing trajectory projections, their clustering and state assignments.

    Parameters
    ----------
    projection : numpy.ndarray
        The projected metrics
    reference : numpy.ndarray
        Reference indices to the simulations and frames that generated the metrics
    sim : :class:`Sim <htmd.simlist.Sim>` object
        A simulation object
    cluster : numpy.ndarray
        A list of cluster indexes for each frame in the simulation
    """

    def __init__(self, projection=None, reference=None, sim=None, cluster=None):
        self._projection = projection
        self._reference = reference
        self._cluster = cluster
        self.sim = sim
        self._checkframes((self.projection, self.reference, self.cluster))

    def toHDF5(self, h5group: h5py.Group):
        h5group.create_dataset("_projection", data=self._projection)
        h5group.create_dataset("_reference", data=self._reference)
        if self._cluster is not None:
            h5group.create_dataset("_cluster", data=self._cluster)
        simgroup = h5group.create_group("sim")
        self.sim.toHDF5(simgroup)

    @staticmethod
    def fromHDF5(h5group: h5py.Group):
        from htmd.simlist import Sim

        self = Trajectory()
        self._projection = np.array(h5group["_projection"])
        self._reference = np.array(h5group["_reference"])
        if "_cluster" in h5group:
            self._cluster = np.array(h5group["_cluster"])
        if "sim" in h5group:
            self.sim = Sim.fromHDF5(h5group["sim"])
        return self

    @property
    def projection(self):
        """The projected metrics of this simulation trajectory"""
        return self._projection

    @projection.setter
    def projection(self, value):
        self._checkframes((value, self.reference, self.cluster))
        self._projection = value

    @property
    def reference(self):
        """The reference indices to the simulations and frames that generated the projections"""
        return self._reference

    @reference.setter
    def reference(self, value):
        self._checkframes((self.projection, value, self.cluster))
        self._reference = value

    @property
    def cluster(self):
        """The cluster indexes for each frame in the simulation"""
        return self._cluster

    @cluster.setter
    def cluster(self, value):
        self._checkframes((self.projection, self.reference, value))
        self._cluster = value

    @property
    def numFrames(self):
        """The number of frames in this trajectory"""
        return self.projection.shape[0]

    @property
    def numDimensions(self):
        """The number of dimensions in the projection of this trajectory"""
        return self.projection.shape[1]

    def _numframes(self, args):
        return np.unique([x for x in list(map(_getsizes, args)) if x is not None])

    def _checkframes(self, args):
        if len(self._numframes(args)) > 1:
            raise RuntimeError(
                "projection, reference and cluster must have same lengths"
            )

    def dropFrames(self, frames):
        """Drop frames from the trajectory

        Parameters
        ----------
        frames : list
            A list of frame indexes to drop
        """
        if np.min(frames) < 0 or np.max(frames) >= self.numFrames:
            raise RuntimeError(
                "Frames to drop must be > 0 and < {}".format(self.numFrames)
            )
        if self._projection is not None:
            self._projection = np.delete(self._projection, frames, axis=0)
        if self._reference is not None:
            self._reference = np.delete(self._reference, frames, axis=0)
        if self._cluster is not None:
            self._cluster = np.delete(self._cluster, frames, axis=0)
        self._checkframes((self.projection, self.reference, self.cluster))

    def copy(self):
        """Produces a deep copy of the object"""
        return deepcopy(self)

    def __repr__(self):
        return (
            "<{}.{} object at {}>\n".format(
                self.__class__.__module__, self.__class__.__name__, hex(id(self))
            )
            + self.__str__()
        )

    def __str__(self):
        return "sim: {}\nprojection: {}\nreference: {}\ncluster: {}\n".format(
            "simid = {}".format(self.sim.simid) if self.sim is not None else None,
            (
                "np.array(shape={})".format(np.shape(self.projection))
                if self.projection is not None
                else None
            ),
            (
                "np.array(shape={})".format(np.shape(self.reference))
                if self.reference is not None
                else None
            ),
            (
                "np.array(shape={})".format(np.shape(self.cluster))
                if self.cluster is not None
                else None
            ),
        )


class MetricData(object):
    """Class used for storing projected trajectories, their clustering and state assignments. Objects of this class
    are constructed by the `project` methods of the other projection classes. Only construct this class if you want to
    load saved data.

    Parameters
    ----------
    dat : numpy.ndarray
        The projected metrics
    ref : numpy.ndarray
        Reference indices to the simulations and frames that generated the metrics
    description : pandas.DataFrame
        A description of the metrics
    simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function
    fstep : float
        Size of simulation step in ns
    parent : :class:`MetricData` object
        The MetricData object that was used to generate this object
    file : str
        Path to a saved MetricData object
    trajectories : list of :class:`Trajectory` objects
        A list of Trajectory objects
    cluster : numpy.ndarray
        A list of cluster indexes for each frame in the data

    Attributes
    ----------
    trajectories : list of :class:`Trajectory` objects
        A list of Trajectory objects
    fstep : float
        Size of simulation step in ns
    description : pandas.DataFrame
        A description of the metrics
    K : int
        Number of clusters
    N : numpy.ndarray
        Number of frames in each cluster
    Centers : numpy.ndarray
        Cluster centers (if available)
    """

    def __init__(
        self,
        dat=None,
        ref=None,
        description=None,
        simlist=None,
        fstep=0,
        parent=None,
        file=None,
        trajectories=None,
        cluster=None,
    ):
        if trajectories is None:
            self._loadTrajectories(dat, ref, simlist, cluster)
        else:
            self.trajectories = trajectories

        self.fstep = fstep
        self.description = description
        self.parent = parent
        self.K = None
        self.N = None
        self.Centers = None

        if file is not None:
            self.load(file)

        self._dataid = str(uuid.uuid4())
        self._clusterid = None
        return

    def _loadTrajectories(
        self, projection=None, reference=None, simlist=None, cluster=None
    ):
        size = np.unique(
            [
                x
                for x in list(map(_getsizes, (projection, reference, simlist, cluster)))
                if x is not None
            ]
        )
        if len(size) == 0:
            size = 0
        elif len(size) > 1:
            raise RuntimeError("dat, ref and simlist must have same lengths")
        if projection is None:
            projection = np.empty(size, dtype=object)
        if reference is None:
            reference = np.empty(size, dtype=object)
        if simlist is None:
            simlist = np.empty(size, dtype=object)
        if cluster is None:
            cluster = np.empty(size, dtype=object)
        self.trajectories = [
            Trajectory(d, r, s, c)
            for d, r, s, c in zip(projection, reference, simlist, cluster)
        ]

    @property
    def dat(self):
        return [t.projection for t in self.trajectories]

    @property
    def ref(self):
        return [t.reference for t in self.trajectories]

    @property
    def St(self):
        return [t.cluster for t in self.trajectories]

    @property
    def simlist(self):
        return [x.sim for x in self.trajectories]

    @property
    def map(self):
        return self.description

    @property
    def trajLengths(self):
        """Get the lengths of all trajectories

        Returns
        -------
        lens : list
            The lengths of all trajectories in the object

        Examples
        --------
        >>> data.trajLengths
        """
        return np.array([x.numFrames for x in self.trajectories])

    @property
    def numFrames(self):
        """Get the total number of frames in all trajectories

        Returns
        -------
        nframes : int
            Total number of frames in all trajectories

        Examples
        --------
        >>> data.numFrames
        """
        return sum(self.trajLengths)

    @property
    def numTrajectories(self):
        """The number of trajectories

        Examples
        --------
        >>> data.numTrajectories
        """
        return len(self.trajectories)

    @property
    def numDimensions(self):
        """The number of dimensions

        Examples
        --------
        >>> data.numDimensions
        """
        return self.trajectories[0].numDimensions

    @property
    def aggregateTime(self):
        """The total aggregate simulation time

        Examples
        --------
        >>> data.aggTime
        """
        if self.fstep > 0:
            return self.numFrames * self.fstep

    def cluster(self, clusterobj, mergesmall=None, batchsize=False):
        """Cluster the metrics

        Parameters
        ----------
        clusterobj : :class:`ClusterMixin <sklearn.cluster.ClusterMixin>` object
            The object of a clustering class from sklearn or with the same interface
        mergesmall : int
            Clusters containing less than `mergesmall` conformations will be joined into their closest well-populated
            neighbour.
        batchsize : int
            Batch sizes bigger than 0 will enable batching.

        Examples
        --------
        >>> from sklearn.cluster import MiniBatchKMeans
        >>> data = MetricDistance.project(sims, 'protein and name CA', 'resname MOL')
        >>> data.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=5)
        """
        # cluster_obj = coor.cluster_kmeans(self.dat, k=20, stride=1)
        if batchsize > 0:
            lengths = self.trajLengths
            currsum = 0
            starts = [0]
            for i, l in enumerate(lengths):
                currsum += l
                if currsum > batchsize:
                    starts.append(i + 1)
                    currsum = 0
            starts.append(self.numTrajectories)
            for i in range(len(starts) - 1):
                clusterobj.partial_fit(
                    np.concatenate(self.dat[starts[i] : starts[i + 1]])
                )
            labels = []
            for i in range(self.numTrajectories):
                labels.append(clusterobj.predict(self.dat[i].astype(np.float32)))
            # This is retarded
            labels = np.concatenate(labels)
        else:
            datconcat = np.concatenate(self.dat)
            if np.ndim(datconcat) == 1:
                datconcat = np.transpose(np.atleast_2d(datconcat))
            import warnings  # Following 3 lines are BS because sklearn refuse to make releases more often than 1 per year...

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                clusterobj.fit(datconcat)
            labels = clusterobj.labels_

        uqclu = np.unique(labels)
        self.Centers = clusterobj.cluster_centers_[uqclu, :]
        self.K = len(uqclu)
        # ---- Fixing missing clusters crap...
        map = np.zeros(np.max(labels) + 1, dtype=int) * -1
        map[uqclu] = range(self.K)
        labels = map[labels]
        # ------------------------------------
        for i, s in enumerate(self.deconcatenate(labels)):
            self.trajectories[i].cluster = s
        self.N = np.bincount(labels)

        if mergesmall is not None:
            oldK = self.K
            self.K, St, self.Centers, self.N, xxx = _mergeSmallClusters(
                mergesmall, datconcat, labels, self.Centers, self.N
            )
            for i, s in enumerate(self.deconcatenate(St)):
                self.trajectories[i].cluster = s
            logger.info(
                "Mergesmall removed {} clusters. Original ncluster {}, new ncluster {}.".format(
                    oldK - self.K, oldK, self.K
                )
            )

        self._dataid = str(uuid.uuid4())
        self._clusterid = self._dataid

    def combine(self, otherdata):
        """Combines two different metrics into one by concatenating them.

        Parameters
        ----------
        otherdata : :class:`MetricData` object
            Concatenates the metrics of otherdata to the current objects metrics

        Examples
        --------
        >>> dataRMSD = MetricRmsd.project(sims)
        >>> dataDist = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> dataRMSD.combine(dataDist)
        """
        import pandas as pd

        if not np.array_equal(self.trajLengths, otherdata.trajLengths):
            raise NameError(
                "Trying to combine MetricData objects with different number/lengths of trajectories. Check the trajLengths property."
            )
        for i in range(self.numTrajectories):
            if self.simlist[i].simid != otherdata.simlist[i].simid:
                raise NameError(
                    "Simulation ids do not match. Cannot combine. Please generate both data from the same simlist"
                )
        for t1, t2 in zip(self.trajectories, otherdata.trajectories):
            t1.projection = np.concatenate((t1.projection, t2.projection), axis=1)
        self.description = pd.concat(
            [self.description, otherdata.description], ignore_index=True
        )
        self._dataid = str(uuid.uuid4())

    def dropDimensions(self, drop=None, keep=None):
        """Drop some dimensions of the data given their indexes

        Parameters
        ----------
        drop : list
            A list of integer indexes of the dimensions to drop
        keep : list
            A list of integer indexes of the dimensions to keep

        Examples
        --------
        >>> data.dropDimensions([1, 24, 3])
        >>> data.dropDimensions(keep=[2, 10])
        """
        if drop is not None and not isinstance(drop, np.ndarray):
            drop = np.array(drop)
        if keep is not None and not isinstance(keep, np.ndarray):
            keep = np.array(keep)
        if drop is not None and keep is not None:
            raise AttributeError(
                "drop and keep arguments for dropDimensions are mutually exclusive. Pass only one."
            )
        if keep is not None:
            keepidx = keep
            dropidx = np.arange(self.numDimensions)
            dropidx = np.setdiff1d(dropidx, keepidx)
        else:
            dropidx = drop
            keepidx = np.arange(self.numDimensions)
            keepidx = np.setdiff1d(keepidx, dropidx)

        for t in self.trajectories:
            t.projection = t.projection[:, keepidx]
        self.description = self.description.drop(self.description.index[dropidx])
        self.description = self.description.reset_index(drop=True)

    def dropTraj(
        self, limits=None, multiple=None, partial=None, idx=None, keepsims=None
    ):
        """Drops trajectories based on their lengths

        By default, drops all trajectories which are not of statistical mode (most common) length.

        Parameters
        ----------
        limits : list, optional
            Lower and upper limits of trajectory lengths we want to keep. e.g. [100, 500]
        multiple : list, optional
            Drops trajectories whose length is not a multiple of lengths in the list. e.g. [50, 80]
        partial : bool
            Not implemented yet
        idx : list, optional
            A list of trajectory indexes to drop
        keepsims : list of :class:`Sim <htmd.simlist.Sim>` objects
            A list of sims which we want to keep

        Examples
        --------
        >>> data = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> data.dropTraj()
        >>> data.dropTraj(multiple=[100])
        """
        trajLengths = self.trajLengths
        orgNum = self.numTrajectories

        if limits is not None:
            drop = (trajLengths < limits[0]) | (trajLengths > limits[1])
        elif multiple is not None:
            if partial is not None:
                raise NameError("TODO")
                pass
            idx = range(orgNum)
            for i in range(len(multiple)):
                idxNew = np.where(np.mod(trajLengths, multiple[i]) != 0)
                idx = np.intersect1d(idxNew, idx)
            drop = np.zeros(orgNum, dtype=bool)
            drop[idx] = True
        elif idx is not None:
            drop = np.zeros(orgNum, dtype=bool)
            drop[idx] = True
        elif keepsims is not None:
            # Fast check to see if simlists are identical
            if len(keepsims) == self.numTrajectories:
                allsame = True
                for i, t in enumerate(self.trajectories):
                    if keepsims[i] != t.sim:
                        allsame = False
                        break
                if allsame:
                    return

            # Slow check where all sims are checked against each other
            drop = np.ones(orgNum, dtype=bool)
            for i, s in enumerate(self.simlist):
                for k in keepsims:
                    if s == k:
                        drop[i] = False
                        break
        else:
            from scipy import stats

            drop = trajLengths != np.array(stats.mode(trajLengths, keepdims=False).mode)

        keep = np.invert(drop)
        dropIdx = np.where(drop)[0]

        self.trajectories = [self.trajectories[x] for x in np.where(keep)[0]]
        self._dataid = str(uuid.uuid4())
        if self.parent:
            self.parent.trajectories = [
                self.parent.trajectories[x] for x in np.where(keep)[0]
            ]
            self.parent._dataid = str(uuid.uuid4())

        logger.info(
            "Dropped "
            + str(np.sum(drop))
            + " trajectories from "
            + str(orgNum)
            + " resulting in "
            + str(self.numTrajectories)
        )
        return dropIdx

    def dropFrames(self, idx, frames):
        self.trajectories[idx].dropFrames(frames)
        self._dataid = str(uuid.uuid4())
        if self.parent:
            self.parent.trajectories[idx].dropFrames(frames)
            self.parent._dataid = str(uuid.uuid4())

    def sampleClusters(
        self, clusters=None, frames=20, replacement=False, allframes=False
    ):
        """Samples frames from a set of clusters

        Parameters
        ----------
        clusters : Union[None, list]
            A list of cluster indexes from which we want to sample
        frames : Union[None, int, list]
            An integer with the number of frames we want to sample from each state or a list of same length as
            `clusters` which contains the number of frames we want from each of the clusters.
            If None is given it will return all frames.
        replacement : bool
            If we want to sample with or without replacement.
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
        >>> data.sampleClusters(range(5), [10, 3, 2, 50, 1])  # Sample from first 5 clusters, 10, 3, etc frames respectively
        """
        if clusters is None:
            clusters = range(self.K)
        if isinstance(clusters, int):
            clusters = [
                clusters,
            ]
        if allframes:
            logger.warning("allframes option is deprecated. Please use frames=None")
            frames = None
        if frames is None or isinstance(frames, int):
            frames = np.repeat(frames, len(clusters))

        stConcat = np.concatenate(self.St)
        absFrames = []
        relFrames = []
        for i in range(len(clusters)):
            if frames[i] == 0 and not allframes:
                continue
            st = clusters[i]
            absFrames.append(_sampleCluster(st, stConcat, frames[i], replacement))
            if len(absFrames[-1]) == 0:
                raise NameError(
                    "No frames could be sampled from cluster {}. Cluster is empty.".format(
                        st
                    )
                )

            relFrames.append(self.abs2rel(absFrames[-1]))
        return absFrames, relFrames

    def bootstrap(self, ratio, replacement=False):
        """Randomly sample a set of trajectories

        Parameters
        ----------
        ratio : float
            What ratio of trajectories to keep. e.g. 0.8
        replacement : bool
            If we should sample with replacement

        Returns
        -------
        bootdata : :class:`MetricData` object
            A new :class:`MetricData` object containing only the sampled trajectories

        Examples
        --------
        >>> data = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> databoot = data.bootstrap(0.8)
        """
        numtraj = self.numTrajectories
        numtokeep = int(np.floor(numtraj * ratio))
        if replacement:
            rndtraj = np.random.randint(numtraj, size=numtokeep)
        else:
            rndtraj = np.random.permutation(numtraj)[0:numtokeep]
        rndtraj = sorted(
            rndtraj
        )  # Important to keep the sorting! i.e. for data.dropTraj(keepsims=sims)

        pp = None
        if self.parent is not None:
            pp = self.parent.copy()
            pp.trajectories = [self.parent.trajectories[x] for x in rndtraj]
            pp._dataid = str(uuid.uuid4())
        bootdata = MetricData(
            trajectories=[self.trajectories[x].copy() for x in rndtraj],
            description=self.description,
            parent=pp,
            fstep=self.fstep,
        )
        return bootdata

    def plotTrajSizes(self):
        """Plot the lengths of all trajectories in a sorted bar plot

        Examples
        --------
        >>> data = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> data.plotTrajSizes()
        """
        trajLengths = self.trajLengths * self.fstep
        import matplotlib.pyplot as plt

        plt.ion()
        _ = plt.hist(trajLengths)
        # plt.bar(range(len(trajLengths)), np.sort(trajLengths), color='b', edgecolor='b')
        plt.ylabel("Num trajectories")
        plt.xlabel("Length of trajectories (in ns)")
        plt.show()
        return

    def splitCols(self):
        raise NameError("Not implemented yet")

    def deconcatenate(self, array):
        indeces = np.cumsum(self.trajLengths)
        if np.ndim(array) == 1:
            return np.split(array, indeces[:-1])
        else:
            return np.vsplit(array, indeces[:-1])

    def abs2rel(self, absFrames):
        """Convert absolute frame indexes into trajectory index-frame pairs

        Useful when doing calculations on a concatenated data array of all trajectories. When you find a frame of
        interest you can `deconcatenate` the frame index to the corresponding trajectory index-frame pair.

        Parameters
        ----------
        absFrames : list of int
            A list of absolute index frames

        Returns
        -------
        pairs : np.ndarray
            A array where each row is a trajectory index-frame pair

        Examples
        --------
        >>> relidx = data.abs2rel(536)
        """
        if not hasattr(absFrames, "__len__"):
            absFrames = [absFrames]
        endFrames = np.append(0, np.cumsum(self.trajLengths))

        relframe = np.zeros((len(absFrames), 2), dtype=int)
        for i in range(len(absFrames)):
            trajIdx = np.where(absFrames[i] < endFrames)[0][0] - 1
            trajFr = absFrames[i] - endFrames[trajIdx]
            relframe[i, :] = [trajIdx, trajFr]
        return relframe

    def rel2sim(self, relFrames, simlist=None):
        """Converts trajectory index-frame pairs into Sim-frame pairs

        Parameters
        ----------
        relFrames : 2D np.ndarray
            An array containing in each row trajectory index and frame pairs
        simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
            Optionally pass a different (but matching, i.e. filtered) simlist for creating the Frames.

        Returns
        -------
        frames : np.ndarray
            An array of :class:`Frame <htmd.simlist.Frame>` objects containing the simulation object, the trajectory
            piece ID and the frame index.

        Examples
        --------
        >>> simframes = data.rel2sim([100, 56])  # 100th simulation frame 56
        """
        from htmd.simlist import Frame

        if simlist is None:
            simlist = self.simlist
        else:
            if len(simlist) != len(self.simlist):
                raise AttributeError(
                    "Provided simlist has different number of trajectories than the one used by this object."
                )

        relFrames = np.array(relFrames)
        if relFrames.ndim == 1:
            relFrames = relFrames[np.newaxis, :]

        frames = []
        for trajID, trajFrame in relFrames:
            ref = self.trajectories[trajID].reference
            frames.append(Frame(simlist[trajID], ref[trajFrame, 0], ref[trajFrame, 1]))
        return np.array(frames)

    def abs2sim(self, absFrames):
        """Converts absolute frame indexes into Sim-frame pairs

        Parameters
        ----------
        absFrames : list of int
            A list of absolute index frames

        Returns
        -------
        frames : np.ndarray
            An array of :class:`Frame <htmd.simlist.Frame>` objects containing the simulation object, the trajectory
            piece ID and the frame index.

        Examples
        --------
        >>> simframes = data.abs2sim(563)  # 563rd frame to simulation/frame pairs
        """
        return self.rel2sim(self.abs2rel(absFrames))

    def copy(self):
        """Produces a deep copy of the object

        Returns
        -------
        data : :class:`MetricData` object
            A copy of the current object

        Examples
        --------
        >>> data = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> data2 = data.copy()
        """
        return deepcopy(self)

    def toHDF5(self, filename=None, h5group=None):
        import h5py

        if h5group is None and filename is None:
            raise RuntimeError("Either h5group or filename has to be set")

        h5f = None
        if h5group is None:
            h5f = h5py.File(filename, "w")
            h5group = h5f.create_group("MetricData")

        def _set_attrs(grp, keys, attr=True):
            for key in keys:
                if getattr(self, key):
                    if attr:
                        grp.attrs[key] = getattr(self, key)
                    else:
                        grp.create_dataset(key, data=getattr(self, key))

        _set_attrs(h5group, ["K", "_dataid", "_clusterid", "fstep"])
        h5group.attrs["description"] = self.description.to_json()
        _set_attrs(h5group, ["N", "Centers"], False)

        trajgrp = h5group.create_group("trajectories")
        for i in range(len(self.trajectories)):
            trajg = trajgrp.create_group(str(i))
            self.trajectories[i].toHDF5(trajg)

        if self.parent is not None:
            parentgroup = h5group.create_group("parent")
            self.parent.toHDF5(h5group=parentgroup)

        if h5f is not None:
            h5f.close()

    @staticmethod
    def fromHDF5(filename=None, h5group=None):
        import pandas as pd
        from natsort import natsorted

        if h5group is None and filename is None:
            raise RuntimeError("Either h5group or filename has to be set")

        data = MetricData()

        h5f = None
        if h5group is None:
            h5f = h5py.File(filename, "r")
            h5group = h5f["MetricData"]

        def _get_attrs(grp, keys, attr=True):
            for key in keys:
                if attr and key in grp.attrs:
                    setattr(data, key, grp.attrs[key])
                elif key in grp:
                    setattr(data, key, np.array(grp[key]))

        _get_attrs(h5group, ["K", "_dataid", "_clusterid", "fstep"])
        data.description = pd.read_json(h5group.attrs["description"])
        _get_attrs(h5group, ["N", "Centers"], False)

        data.trajectories = []
        if "trajectories" in h5group:
            for key in natsorted(h5group["trajectories"]):
                data.trajectories.append(
                    Trajectory.fromHDF5(h5group["trajectories"][key])
                )

        if "parent" in h5group:
            data.parent = MetricData.fromHDF5(h5group=h5group["parent"])

        if h5f is not None:
            h5f.close()

        return data

    def save(self, filename):
        """Save a :class:`MetricData` object to disk

        Parameters
        ----------
        filename : str
            Path of the file in which to save the object

        Examples
        --------
        >>> data = MetricSelfDistance.project(sims, 'protein and name CA')
        >>> data.save('./data.dat')
        """
        # np.save(filename, [self.__dict__[k] for k in self.__dict__])
        parentpointer = self.parent
        if self.parent is not None:
            self.parent = self.parent.__dict__

        f = open(filename, "wb")
        pickle.dump(self.__dict__, f)
        f.close()

        if self.parent is not None:
            self.parent = parentpointer

    def load(self, filename):
        """Load a :class:`MetricData` object from disk

        Parameters
        ----------
        filename : str
            Path to the saved MetricData object

        Examples
        --------
        >>> data = MetricData()
        >>> data.load('./data.dat')
        """
        import sys

        try:
            import pandas.indexes
        except ImportError:
            import pandas.core.indexes

            sys.modules["pandas.indexes"] = (
                pandas.core.indexes
            )  # Hacky fix for new pandas version

        # Patch for old HTMD versions
        if type(filename).__name__ == "MetricData":
            filename = filename.__dict__

        if isinstance(filename, str):
            f = open(filename, "rb")
            vardict = pickle.load(f)
            f.close()
        elif isinstance(filename, dict):
            vardict = filename

        for k in self.__dict__:
            if k == "description" and "map" in vardict:  # Patch for loading old data
                self.description = vardict["map"]
            elif k == "trajectories" and "dat" in vardict:  # Patch for loading old data
                self._loadTrajectories(
                    vardict["dat"], vardict["ref"], vardict["simlist"], vardict["St"]
                )
            elif k != "parent":
                try:
                    self.__dict__[k] = vardict[k]
                except Exception:
                    logger.warning(
                        "Could not find class property {} in file {}".format(
                            k, filename
                        )
                    )

        if "parent" in vardict and vardict["parent"] is not None:
            self.parent = MetricData()
            self.parent.load(vardict["parent"])

    def _defaultLags(self, minlag=None, maxlag=None, numlags=None, units="frames"):
        from htmd.units import convert as unitconvert

        if maxlag is None:
            from scipy import stats

            maxlag = (
                stats.mode(self.trajLengths, keepdims=False).mode - 1
            )  # -1 to avoid warnings in timescales calc
        else:
            maxlag = unitconvert(units, "frames", maxlag, fstep=self.fstep)

        if minlag is None:
            if maxlag > 20:
                minlag = 10
            else:
                minlag = 2
        else:
            minlag = unitconvert(units, "frames", minlag, fstep=self.fstep)

        return np.append(1, np.round(np.linspace(minlag, maxlag, numlags))).astype(int)

    def _getFEShistogramCounts(
        self,
        dimx,
        dimy,
        nbins=80,
        pad=0.5,
        micro_ofcluster=None,
        stationary_distribution=None,
        St=None,
    ):
        from tqdm import tqdm

        if St is None:
            St = self.St

        clusters = np.hstack(St)
        data_dim = []
        for traj in self.trajectories:
            data_dim.append(traj.projection[:, [dimx, dimy]])
        data_dim = np.vstack(data_dim)

        xmin, ymin = data_dim.min(axis=0)
        xmax, ymax = data_dim.max(axis=0)
        hist_range = np.array([[xmin - pad, xmax + pad], [ymin - pad, ymax + pad]])

        if micro_ofcluster is None:
            counts, xbins, ybins = np.histogram2d(
                data_dim[:, 0], data_dim[:, 1], bins=nbins, range=hist_range
            )
            return counts.T, xbins, ybins

        counts = np.zeros((nbins, nbins))
        for m in tqdm(range(len(stationary_distribution))):
            frames = micro_ofcluster[clusters] == m
            if frames.sum():
                state_mask, xbins, ybins = np.histogram2d(
                    data_dim[frames, 0],
                    data_dim[frames, 1],
                    bins=nbins,
                    range=hist_range,
                )
                state_mask = state_mask.T > 0
                counts += state_mask * stationary_distribution[m]
        return counts, xbins, ybins

    def _contourPlot(
        self,
        values,
        xbins,
        ybins,
        levels=7,
        nonzero=None,
        cmap="viridis",
        title=None,
        xlabel=None,
        ylabel=None,
    ):
        import matplotlib.pylab as plt

        if nonzero is None:
            nonzero = np.ones_like(values).astype(bool)

        def getcmap(cmap, bounds):
            import matplotlib.colors as clr

            cmap = plt.get_cmap(cmap)
            cmaplist = [cmap(i) for i in range(cmap.N)]
            cmap = clr.LinearSegmentedColormap.from_list(
                "Custom cmap", cmaplist, cmap.N
            )
            norm = clr.BoundaryNorm(bounds, cmap.N)
            cmap.set_under(color="white")
            cmap.set_over(color="white")
            return cmap, norm

        xcenters = (xbins[:-1] + xbins[1:]) / 2
        ycenters = (ybins[:-1] + ybins[1:]) / 2
        meshx, meshy = np.meshgrid(xcenters, ycenters)

        levels = np.linspace(values[nonzero].min(), values[nonzero].max(), levels)
        cmap, norm = getcmap(cmap, levels)
        f = plt.figure()
        plt.contour(
            meshx, meshy, values, levels=levels, colors="black", vmin=0, vmax=levels[-1]
        )
        cf = plt.contourf(
            meshx,
            meshy,
            values,
            levels=levels,
            cmap=cmap,
            vmin=0,
            vmax=levels[-1],
            norm=norm,
        )

        ax = f.gca()
        _ = ax.axis("equal")
        if title is not None:
            ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        return f, ax, cf

    def _setColorbar(self, f, mappable, label=None, scientific=True):
        import matplotlib.ticker as ticker

        def fmt(x, pos):
            a, b = "{:.2e}".format(x).split("e")
            b = int(b)
            return r"${} \times 10^{{{}}}$".format(a, b)

        if scientific:
            f.colorbar(mappable, format=ticker.FuncFormatter(fmt), label=label)
        else:
            f.colorbar(mappable, label=label)

    def _plotCounts(
        self,
        dimX,
        dimY,
        resolution=100,
        logplot=False,
        levels=7,
        cmap="viridis",
        title=None,
        xlabel=None,
        ylabel=None,
    ):
        counts, xbins, ybins = self._getFEShistogramCounts(dimX, dimY, nbins=resolution)
        nonzero = counts != 0
        if logplot:
            counts[nonzero] = np.log(counts[nonzero])
            nonzero = counts != 0
        f, ax, cf = self._contourPlot(
            counts,
            xbins,
            ybins,
            levels=levels,
            nonzero=nonzero,
            cmap=cmap,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )

        self._setColorbar(f, cf, "Counts")
        return f, ax, cf

    def plotCounts(
        self,
        dimX,
        dimY,
        resolution=100,
        logplot=False,
        plot=True,
        save=None,
        levels=7,
        cmap="viridis",
    ):
        """Plots a histogram of counts on any two given dimensions.

        Parameters
        ----------
        dimX : int
            Index of projected dimension to use for the X axis.
        dimY : int
            Index of projected dimension to use for the Y axis.
        resolution : int
            Resolution of bincount grid.
        logplot : bool
            Set True to plot the logarithm of counts.
        plot : bool
            If the method should display the plot
        save : str
            Path of the file in which to save the figure
        """
        from matplotlib import pylab as plt

        if self.description is not None:
            xlabel = self.description.description[dimX]
        else:
            xlabel = "Dimension {}".format(dimX)
        if self.description is not None:
            ylabel = self.description.description[dimY]
        else:
            ylabel = "Dimension {}".format(dimY)
        title = "Counts histogram"
        if logplot:
            title = "Logarithmic counts histogram"

        self._plotCounts(
            dimX,
            dimY,
            resolution=resolution,
            logplot=logplot,
            levels=levels,
            cmap=cmap,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
        if plot:
            plt.show()

    def plotClusters(
        self,
        dimX,
        dimY,
        resolution=100,
        s=4,
        c=None,
        cmap="Greys",
        logplot=False,
        plot=True,
        save=None,
        data=None,
        levels=7,
    ):
        """Plot a scatter-plot of the locations of the clusters on top of the count histogram.

        Parameters
        ----------
        dimX : int
            Index of projected dimension to use for the X axis.
        dimY : int
            Index of projected dimension to use for the Y axis.
        resolution : int
            Resolution of bincount grid.
        s : float
            Marker size for clusters.
        c : list
            Colors or indexes for each cluster.
        cmap : matplotlib.colors.Colormap
            Matplotlib colormap for the scatter plot.
        logplot : bool
            Set True to plot the logarithm of counts.
        plot : bool
            If the method should display the plot
        save : str
            Path of the file in which to save the figure
        data : :class:`MetricData` object
            Optionally you can pass a different MetricData object than the one used for clustering. For example
            if the user wants to cluster on distances but wants to plot the centers on top of RMSD values. The new
            MetricData object needs to have the same simlist as this object.
        """
        if self.Centers is None:
            raise RuntimeError("Data has not been clustered yet. Cannot plot clusters.")
        from matplotlib import pylab as plt

        if data is None:
            data = self
            centers = self.Centers
        else:
            from htmd.model import getStateStatistic

            if self.numFrames != data.numFrames or ~np.all(
                [s1 == s2 for s1, s2 in zip(self.simlist, data.simlist)]
            ):
                raise RuntimeError(
                    "The data argument you provided uses a different simlist than this object."
                )
            centers = np.vstack(
                getStateStatistic(self, data, range(self.K), statetype="cluster")
            )

        if data.description is not None:
            xlabel = data.description.description[dimX]
        else:
            xlabel = "Dimension {}".format(dimX)

        if data.description is not None:
            ylabel = data.description.description[dimY]
        else:
            ylabel = "Dimension {}".format(dimY)

        title = "Clusters plotted onto counts histogram"
        if logplot:
            title = "Clusters plotted onto logarithmic counts histogram"
        f, ax, cf = self._plotCounts(
            dimX,
            dimY,
            resolution=resolution,
            logplot=logplot,
            levels=levels,
            cmap=cmap,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
        )

        y = ax.scatter(
            centers[:, dimX],
            centers[:, dimY],
            s=s,
            c=c if c is not None else "r",
            cmap=cmap,
            linewidths=0,
            marker="o",
        )
        if c is not None:
            self._setColorbar(f, y, "Cluster groups")

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
        if plot:
            plt.show()

    def __repr__(self):
        return (
            "<{}.{} object at {}>\n".format(
                self.__class__.__module__, self.__class__.__name__, hex(id(self))
            )
            + self.__str__()
        )

    def __str__(self):
        def formatstr(name, field):
            if isinstance(field, np.ndarray) or isinstance(field, list):
                rep = "{} shape: {}".format(name, np.shape(field))
            else:
                rep = "{}: {}".format(name, field)
            return rep

        rep = "MetricData object with {} trajectories".format(self.numTrajectories)
        if self.fstep > 0:
            rep += " of {:.1f}ns aggregate simulation time".format(self.aggregateTime)
        for j in sorted(self.__dict__.keys()):
            if j[0] == "_":
                continue
            if j == "parent":
                rep += "\nparent: {} at {}".format(
                    type(self.parent), hex(id(self.parent))
                )
            elif j == "description":
                rep += "\ndescription: {} at {}".format(
                    type(self.description), hex(id(self.description))
                )
            else:
                rep += "\n"
                rep += formatstr(j, self.__dict__[j])

        return rep

    def append(self, other):
        self.trajectories += other.trajectories
        for i, t in enumerate(self.trajectories):
            t.sim.simid = i

        if self.parent and other.parent:
            self.parent.trajectories += other.parent.trajectories
            for i, t in enumerate(self.parent.trajectories):
                t.sim.simid = i

        self._dataid = str(uuid.uuid4())
        self._resetClustering()
        return self

    def _resetClustering(self):
        for t in self.trajectories:
            t.cluster = None
        self.K = None
        self.N = None
        self.Centers = None

    def sampleRegion(
        self, point=None, radius=None, limits=None, nsamples=20, singlemol=False
    ):
        """Samples conformations from a region in the projected space.

        Parameters
        ----------
        point : list or np.ndarray
            A point in the projected space. Undefined dimensions should have None value.
        radius : float
            The radius in around the point in which to sample conformations.
        limits : np.ndarray
            A (2, ndim) dimensional array containing the min (1st row) and max (2nd row) limits for each dimension.
            None values will be interpreted as no limit in that dimension, or min/max value.
        nsamples : int
            The number of conformations to sample.
        singlemol : bool
            If True it will return all samples within a single Molecule instead of a list of Molecules.

        Returns
        -------
        absFrames : list
            A list of the absolute frame indexes sampled
        relFrames : list of tuples
            A list of (trajNum, frameNum) tuples sampled
        mols : Molecule or list of Molecules
            The conformations stored in a Molecule or a list of Molecules

        Examples
        --------
        >>> # Working with 4 dimensional data for example
        >>> abs, rel, mols = data.sampleRegion(point=(0.5, 3, None, None), radius=0.1)  # Point undefined in dim 3, 4
        >>> minlims = [-1, None, None, 4]  # No min limit for 2, 3 dim
        >>> maxlims = [2,     3, None, 7]  # No max limit for 3 dim
        >>> abs, rel, mols = data.sampleRegion(limits=np.array([minlims, maxlims]))
        """
        from scipy.spatial.distance import cdist

        datconcat = np.concatenate(self.dat)
        numdim = datconcat.shape[1]
        if point is not None:
            if radius is None:
                raise RuntimeError("You must define a radius with a point.")
            point = np.array(point)
            if len(point) != numdim:
                raise RuntimeError(
                    "Argument `point` should be same dimensionality as your data ({} dimensions)".format(
                        numdim
                    )
                )
            keepdim = np.array([p is not None for p in point])
            dists = cdist(datconcat[:, keepdim], [point[keepdim]])
            confs = np.where(dists < radius)[0]
        elif limits is not None:
            if limits.shape != (2, numdim):
                raise RuntimeError(
                    "Argument `limits` should be of shape (2, {})".format(numdim)
                )
            mask = np.ones(datconcat.shape[0], dtype=bool)
            for i in range(numdim):
                if limits[0, i] is not None:
                    mask &= datconcat[:, i] > limits[0, i]
                if limits[1, i] is not None:
                    mask &= datconcat[:, i] < limits[1, i]
            confs = np.where(mask)[0]

        if len(confs) > nsamples:
            confs = np.random.choice(confs, nsamples, replace=False)
        sims = self.abs2sim(confs)

        from moleculekit.molecule import Molecule

        if singlemol:
            mol = Molecule(sims[0])
            for i in range(1, len(sims)):
                m = Molecule(sims[i])
                mol.appendFrames(m)
        else:
            mol = []
            for s in sims:
                mol.append(Molecule(s))
        return confs, self.abs2rel(confs), mol


def _sampleCluster(cluster, stConcat, numFrames, replacement):
    frames = np.where(stConcat == cluster)[0]
    return _randomSample(frames, numFrames, replacement)


def _randomSample(frames, numFr, replacement):
    if numFr == 0:
        return []
    if numFr is None or (numFr >= len(frames) and not replacement):
        rnd = list(range(len(frames)))
    else:
        rnd = np.random.randint(len(frames), size=numFr)
    return frames[rnd]


def _mergeSmallClusters(mergesmall, data, stconcat, centers, N, metric=None):
    if data.dtype == "bool":
        metric = "hamming"
    else:
        metric = "euclidean"
    badclusters = N < mergesmall
    goodclusters = np.invert(badclusters)
    N[badclusters] = 0
    badcluidx = np.where(badclusters)[0]
    goodcluidx = np.where(goodclusters)[0]

    if len(badcluidx) == 0:
        return len(N), stconcat, centers, N, badclusters

    # Keep only good centers
    centers = centers[goodclusters, :]

    # Creating a mapping for new cluster numbers as we will have to remove some
    newidx = np.zeros(len(N), dtype=int)
    newidx[goodcluidx] = range(len(goodcluidx))

    # Find all frames which belong to bad clusters
    frames = _ismember(stconcat, badcluidx)
    badframeidx = np.where(frames >= 0)[0]

    # Calculate distance of all frames belonging to bad clusters to the good cluster centers
    from scipy.spatial import distance

    dists = distance.cdist(
        np.atleast_2d(data[badframeidx, :]), np.atleast_2d(centers), metric
    )
    minidx = np.argmin(
        dists, axis=1
    )  # Find closest center. Indexes are relative to goodidx
    newclu = goodcluidx[minidx]  # Back to absolute cluster indexes

    # Reassign bad frames to good clusters
    stconcat[badframeidx] = newclu  # Assign them to new clusters
    stconcat = newidx[stconcat]  # Convert all cluster indexes to the new indexes

    N = np.bincount(stconcat)
    K = len(N)
    return K, stconcat, centers, N, badclusters


def _ismember(a, b):
    bind = {}
    for i, elt in enumerate(list(set(b))):
        bind[elt] = i
    return np.array(
        [bind.get(itm, -1) for itm in a]
    )  # None can be replaced by any other "not in b" value


def _generate_toy_data(trans_prob, n_traj=1000, n_frames=1000, seed=None, cluster=True):
    from htmd.metricdata import MetricData
    from pandas import DataFrame

    assert np.all((np.sum(trans_prob, axis=1)) - 1 < 1e-15)

    # Create real transition probability matrix
    n_states = trans_prob.shape[0]

    # Generate fake trajectories from it
    dat = []
    ref = []
    if seed is not None:
        np.random.seed(seed)

    for i in range(n_traj):  # ntraj
        trajref = np.zeros((n_frames, 2), dtype=np.int32)
        trajref[:, 0] = i
        trajref[:, 1] = np.arange(n_frames)
        trajdat = np.zeros((n_frames, n_states), dtype=np.int32)
        curr_state = np.random.randint(0, n_states)
        trajdat[0, curr_state] = 1
        for j in range(1, n_frames):  # nsteps
            curr_state = np.random.choice(n_states, p=trans_prob[curr_state])
            trajdat[j, curr_state] = 1
        dat.append(trajdat)
        ref.append(trajref)

    # Create fake test data
    data = MetricData(dat=dat, ref=ref)
    data._dataid = "fake"
    data.fstep = 1

    if cluster:
        for traj in data.trajectories:
            traj.cluster = np.where(traj.projection)[1].copy()
        data._clusterid = "fake"
        data.K = n_states
        data.N = np.bincount(np.concatenate(data.St))
        data.Centers = np.zeros((n_states, n_states), dtype=np.int32)

    types = ["coordinate"] * n_states
    indexes = np.arange(n_states)
    description = ["coordinate"] * n_states
    data.description = DataFrame(
        {"type": types, "atomIndexes": indexes, "description": description}
    )

    return data


class _TestMetricData(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.simlist import simlist, simfilter
        from glob import glob
        from htmd.projections.metric import Metric
        from moleculekit.projections.metricdistance import MetricDistance
        from moleculekit.projections.metricdihedral import MetricDihedral
        from moleculekit.util import tempname
        from htmd.home import home
        from os.path import join

        sims = simlist(
            glob(join(home(dataDir="adaptive"), "data", "*", "")),
            glob(join(home(dataDir="adaptive"), "input", "*")),
        )
        fsims = simfilter(sims, tempname(), "not water")

        metr = Metric(fsims)
        metr.set(
            MetricDistance(
                "protein and resid 10 and name CA",
                "resname BEN and noh",
                periodic="selections",
                metric="contacts",
                groupsel1="residue",
                threshold=4,
            )
        )
        self.data1 = metr.project()

        metr.set(MetricDihedral())
        self.data2 = metr.project()

    def test_combine(self):
        # Testing combining of metrics
        data1 = self.data1.copy()
        data1.combine(self.data2)

        # Testing dimensions
        assert np.array_equal(
            data1.description.shape, (897, 3)
        ), "combine not working correct"
        assert np.array_equal(
            data1.trajectories[0].projection.shape, (6, 897)
        ), "combine not working correct"
        assert np.array_equal(
            np.where(data1.description.type == "contact")[0],
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
        ), "combine not working correct"

    def test_dropping(self):
        # Testing dimension dropping / keeping
        data1 = self.data1.copy()
        assert np.array_equal(data1.description.shape, (9, 3))
        data1.dropDimensions(range(9))
        assert np.array_equal(
            data1.description.shape, (0, 3)
        ), "dropDimensions not working correct"
        assert np.array_equal(
            data1.trajectories[0].projection.shape, (6, 0)
        ), "dropDimensions not working correct"
        assert (
            len(np.where(data1.description.type == "contact")[0]) == 0
        ), "dropDimensions not working correct"

        data2 = self.data2.copy()
        assert np.array_equal(data2.description.shape, (888, 3))
        data2.dropDimensions(keep=range(9))
        assert np.array_equal(
            data2.description.shape, (9, 3)
        ), "dropDimensions not working correct"
        assert np.array_equal(
            data2.trajectories[0].projection.shape, (6, 9)
        ), "dropDimensions not working correct"
        assert (
            len(np.where(data2.description.type == "dihedral")[0]) == 9
        ), "dropDimensions not working correct"

    def test_saving_loading(self):
        from moleculekit.util import tempname

        def checkCorrectness(newdata):
            assert newdata.numTrajectories == 2, "Failed to load trajectories"
            assert newdata.description.shape == (9, 3), "Failed to load pandas data"
            assert newdata.trajectories[0].projection.shape == (6, 9), "Wrong data size"

        savefile = tempname(suffix=".dat")
        self.data1.save(savefile)

        newdata = MetricData(file=savefile)
        checkCorrectness(newdata)

        newdata = MetricData()
        newdata.load(savefile)
        checkCorrectness(newdata)

        # Saving with a parent
        data1 = self.data1.copy()
        data1.parent = self.data2.copy()
        data1.save(savefile)
        newdata = MetricData(file=savefile)
        checkCorrectness(newdata)


if __name__ == "__main__":
    import unittest

    unittest.main(verbosity=2)
