# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
import logging
logger = logging.getLogger(__name__)


class KCenter(BaseEstimator, ClusterMixin, TransformerMixin):
    """ Class to perform KCenter clustering of a given data set

    KCenter randomly picks one point from the data, which is now the center of the first cluster. All points are put
    into the first cluster. In general the furthest point from its center is chosen to be the new center. All points,
    which are closer to the new center than the old one are assigned to the new cluster. This goes on, until K clusters
    have been created.

    Parameters
    ----------
    n_clusters: int
        desired number of clusters

    Examples
    --------
    >>> cluster = KCenter(n_cluster=200)
    >>> cluster.fit(data)

    Attributes
    ----------
    cluster_centers : list
        list with the points, which are the centers of the clusters
    centerFrames : list
        list of indices of center points in data array
    labels_ : list
        list with number of cluster of each frame
    clusterSize_ : list
        list with number of frames in each cluster
    distance : list
        list with the distance of each frame from the nearest center
    """

    def __init__(self, n_clusters):
        self.n_clusters = n_clusters
        self.cluster_centers_ = []
        self.centerFrames = []
        self.labels_ = []
        self.clusterSize = []
        self.distance = []

    def fit(self, data):
        """ Compute the centroids of data.

        Parameters
        ----------
        data : np.ndarray
            A 2D array of data. Columns are features and rows are data examples.
        """
        if len(self.cluster_centers_) != 0:
            logger.warning('Clustering already exists. Reclustering data!')
            self.cluster_centers_ = []
            self.centerFrames = []
            self.clusterSize = []

        # Initialization
        # select random point and assign all points to cluster 0
        numpoints = np.size(data, 0)

        idxCenter = np.random.randint(numpoints)
        self.cluster_centers_.append(data[idxCenter, :])
        self.centerFrames.append(idxCenter)
        self.labels_ = np.zeros(numpoints, dtype=int)

        dist = self._dist(self.cluster_centers_, data)
        countCluster = 1

        while len(self.cluster_centers_) < self.n_clusters:
            if np.max(dist) == 0:
                break

            # find point furthest away from all center
            newCenterIdx = np.argmax(dist)
            newCenter = data[newCenterIdx, :]
            self.centerFrames.append(newCenterIdx)
            self.cluster_centers_.append(newCenter)

            # find all points closer to new center than old center
            newdist = self._dist(newCenter, data)
            switchIdx = dist > newdist

            # assign them to new cluster
            self.labels_[switchIdx] = countCluster
            dist[switchIdx] = newdist[switchIdx]

            countCluster += 1

        # update clusterSize
        self.clusterSize = np.bincount(self.labels_)
        self.distance = dist
        self.cluster_centers_ = np.array(self.cluster_centers_)

    @staticmethod
    def _dist(centers, data):
        dist = np.squeeze(cdist(np.atleast_2d(data), np.atleast_2d(centers)))
        return dist


if __name__ == '__main__':
    """
    infile = open("../../clusterdata/R15.txt")
    data = []
    xdata = []
    ydata =[]
    properLabels = []

    for line in infile:
        stuff = line.split()
        data.append((float(stuff[0]), float(stuff[1])))
        xdata.append(float(stuff[0]))
        ydata.append(float(stuff[1]))
        properLabels.append(stuff[2])

    K = len(set(properLabels))
    """

    data = np.random.rand(100, 10)

    cluster1 = KCenter(n_clusters=20)
    cluster1.fit(data)

    """
    from matplotlib import pylab as plt
    plt.figure(0)
    plt.scatter(xdata, ydata, c=cluster1.labels_)

    xcenters1 = []
    ycenters1 = []
    xcenters2 = []
    ycenters2 = []
    for center in cluster1.cluster_centers_:
        xcenters1.append(center[0])
        ycenters1.append(center[1])

    plt.figure(0)
    plt.title("KCenter K="+str(K))
    plt.scatter(xcenters1, ycenters1, marker='H', s=100, c='red')
    plt.show()
    """
