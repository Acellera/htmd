from .kcenters import KCenter
import numpy as np
import random as rd
from scipy.spatial.distance import cdist
import logging
logger = logging.getLogger(__name__)


class RegCluster:
    """ Class to perform regular clustering of a given data set

    RegCluster can be passed a radius or an approximate number of clusters. If a number of clusters is passed, KCenter
    clustering is used to estimate the necessary radius. RegCluster randomly chooses a point and assigns all points
    within the radius of this point to the same cluster. Then it proceeds with the nearest point, which is not yet
    assigned to a cluster and puts all unassigned points within the radius of this point in the next cluster and so on.

    Parameters
    ----------
    radius: float
        radius of clusters
    n_clusters: int
        desired number of clusters

    Examples
    --------
    >>> cluster = RegCluster(radius=5.1)
    >>> cluster.fit(data)

    Properties
    ----------
    cluster_centers: list
        list with the points, which are the centers of the clusters
    centerFrames: list
        list of indices of center points in data array
    labels_: list
        list with number of cluster of each frame
    clusterSize_: list
        list with number of frames in each cluster
    """
    def __init__(self, radius=None, n_clusters=None):
        if radius is None and n_clusters is None:
            raise(NameError, "radius or n_clusters needs to be set")

        self.radius = radius
        self.cluster_centers_ = []
        self.centerFrames = []
        self.labels_ = []
        self.clusterSize = []
        self.n_clusters = n_clusters

    def fit(self, data):
        """ performs clustering of data

        Parameters
        ----------
        data: np.ndarray
                array of data points to cluster
        merge: int
                minimal number of frames within each cluster. Smaller clusters are merged into next big one
        """

        if len(self.cluster_centers_) != 0:
            logger.warning('Clustering already exists. Reclustering data!')

        # if n_clusters is given and no r, estimate n_clusters
        if self.radius is None:
            estClust = KCenter(n_clusters=self.n_clusters)
            estClust.fit(data)
            self.radius = estClust.distance.max()
            logger.log("Estimated radius = {}".format(self.radius))

        self.cluster_centers_ = []
        self.centerFrames = []
        self.labels_ = np.ones(data.shape[0], dtype=int) * -1

        # find initial center
        unassigned = np.where(self.labels_ == -1)[0]
        idxCenter = rd.choice(unassigned)
        self.cluster_centers_.append(data[idxCenter, :])
        self.centerFrames.append(idxCenter)

        # find all points within radius from initial center

        dist = self._dist(self.cluster_centers_, data[unassigned, :])
        self.labels_[dist <= self.radius] = 0

        countCluster = 1

        # Progress bar will not work here
        while len(np.where(self.labels_ == -1)[0]) > 0:
            # find closest point to existing clusters
            unassigned = np.where(self.labels_ == -1)[0]
            dist = cdist(self.cluster_centers_, data[unassigned, :])
            newCenterIdx = np.argmin(dist)
            newCenter = data[unassigned[newCenterIdx], :]

            # make it a center
            self.cluster_centers_.append(newCenter)
            self.centerFrames.append(unassigned[newCenterIdx])

            # find unassigned points within radius from new cluster
            dist = self._dist(self.cluster_centers_[countCluster], data[unassigned, :])
            self.labels_[dist <= self.radius] = countCluster
            countCluster += 1

        self.cluster_centers_ = np.array(self.cluster_centers_)
        self.clusterSize = np.bincount(self.labels_)

    @staticmethod
    def _dist(centers, data):
        dist = np.squeeze(cdist(np.atleast_2d(data), np.atleast_2d(centers)))
        return dist