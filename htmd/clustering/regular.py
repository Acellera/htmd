# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import random as rd
from scipy.spatial.distance import cdist
import logging
from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
logger = logging.getLogger(__name__)


class RegCluster(BaseEstimator, ClusterMixin, TransformerMixin):
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

    Attributes
    ----------
    cluster_centers:  list
        list with the points, which are the centers of the clusters
    centerFrames : list
        list of indices of center points in data array
    labels_ : list
        list with number of cluster of each frame
    clusterSize_ : list
        list with number of frames in each cluster
    """
    def __init__(self, radius=None, n_clusters=None):
        if radius is None and n_clusters is None:
            raise RuntimeError("radius or n_clusters needs to be set")

        self.radius = radius
        self.n_clusters = n_clusters
        self.labels_ = []

    def fit(self, data):
        """ performs clustering of data

        Parameters
        ----------
        data: np.ndarray
                array of data points to cluster
        merge: int
                minimal number of frames within each cluster. Smaller clusters are merged into next big one
        """
        # if n_clusters is given and no r, estimate n_clusters
        if self.radius is None:
            from htmd.clustering.kcenters import KCenter
            estClust = KCenter(n_clusters=self.n_clusters)
            estClust.fit(data)
            self.radius = estClust.distance.max()
            logger.info("Estimated radius = {}".format(self.radius))

        from pyemma.coordinates.clustering.regspace import RegularSpaceClustering
        self._reg = RegularSpaceClustering(dmin=self.radius)
        self.labels_ = self._reg.fit_transform(data).flatten()

    @property
    def cluster_centers_(self):
        return self._reg.clustercenters

    @property
    def clusterSize(self):
        return np.bincount(self.labels_)

#
# class RegCluster:
#     """ Class to perform regular clustering of a given data set
#
#     RegCluster can be passed a radius or an approximate number of clusters. If a number of clusters is passed, KCenter
#     clustering is used to estimate the necessary radius. RegCluster randomly chooses a point and assigns all points
#     within the radius of this point to the same cluster. Then it proceeds with the nearest point, which is not yet
#     assigned to a cluster and puts all unassigned points within the radius of this point in the next cluster and so on.
#
#     Parameters
#     ----------
#     radius: float
#         radius of clusters
#     n_clusters: int
#         desired number of clusters
#
#     Examples
#     --------
#     >>> cluster = RegCluster(radius=5.1)
#     >>> cluster.fit(data)
#
#     Attributes
#     ----------
#     cluster_centers:  list
#         list with the points, which are the centers of the clusters
#     centerFrames : list
#         list of indices of center points in data array
#     labels_ : list
#         list with number of cluster of each frame
#     clusterSize_ : list
#         list with number of frames in each cluster
#     """
#     def __init__(self, radius=None, n_clusters=None):
#         if radius is None and n_clusters is None:
#             raise RuntimeError("radius or n_clusters needs to be set")
#
#         self.radius = radius
#         self.cluster_centers_ = []
#         self.centerFrames = []
#         self.labels_ = []
#         self.clusterSize = []
#         self.n_clusters = n_clusters
#
#     def fit(self, data):
#         """ performs clustering of data
#
#         Parameters
#         ----------
#         data: np.ndarray
#                 array of data points to cluster
#         merge: int
#                 minimal number of frames within each cluster. Smaller clusters are merged into next big one
#         """
#
#         if len(self.cluster_centers_) != 0:
#             logger.warning('Clustering already exists. Reclustering data!')
#
#         # if n_clusters is given and no r, estimate n_clusters
#         if self.radius is None:
#             estClust = KCenter(n_clusters=self.n_clusters)
#             estClust.fit(data)
#             self.radius = estClust.distance.max()
#             logger.info("Estimated radius = {}".format(self.radius))
#
#         self.cluster_centers_ = []
#         self.centerFrames = []
#         self.labels_ = np.ones(data.shape[0], dtype=int) * -1
#
#         # find initial center
#         unassigned = np.where(self.labels_ == -1)[0]
#         idxCenter = rd.choice(unassigned)
#         self.cluster_centers_.append(data[idxCenter, :])
#         self.centerFrames.append(idxCenter)
#
#         # find all points within radius from initial center
#         dist = self._dist(self.cluster_centers_, data)
#         self.labels_[dist <= self.radius] = 0
#
#         countCluster = 1
#
#         # Progress bar will not work here
#         while len(np.where(self.labels_ == -1)[0]) > 0:
#             if np.mod(countCluster, 100) == 0:
#                 logger.info('{} clusters created'.format(countCluster))
#             # find closest unassigned  point to existing clusters
#             unassigned = np.atleast_1d(np.where(self.labels_ == -1)[0])
#             newCenterIdx = unassigned[np.argmin(dist[unassigned])]
#             newCenter = data[newCenterIdx, :]
#
#             # make it a center
#             self.cluster_centers_.append(newCenter)
#             self.centerFrames.append(newCenterIdx)
#
#             # find unassigned points within radius from new cluster
#             newdist = self._dist(self.cluster_centers_[-1], data[unassigned, :])
#             self.labels_[unassigned[newdist <= self.radius]] = countCluster
#             countCluster += 1
#
#             # Update minimum distances
#             smallerdist = newdist < dist[unassigned]
#             dist[unassigned[smallerdist]] = newdist[smallerdist]
#
#         self.cluster_centers_ = np.array(self.cluster_centers_)
#         self.clusterSize = np.bincount(self.labels_)
#
#     @staticmethod
#     def _dist(centers, data):
#         dist = np.squeeze(cdist(np.atleast_2d(data), np.atleast_2d(centers)))
#         return np.atleast_1d(dist)

if __name__ == '__main__':
    import numpy as np
    X = np.random.random((100, 2))
    reg = RegCluster(radius=0.4)
    reg.fit(X)

    # from matplotlib import pylab as plt
    # plt.figure()
    # plt.scatter(X[:, 0], X[:, 1], c=reg.labels_)
    # centers = np.array(reg.cluster_centers_)
    # plt.scatter(centers[:, 0], centers[:, 1], s=100, c=range(centers.shape[0]))
    # plt.show()
