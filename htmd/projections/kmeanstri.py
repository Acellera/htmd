# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import random
import logging
logger = logging.getLogger(__name__)


class KMeansTri(object):
    """ Class for calculating the Kmeans-triangle projections of a MetricData  object

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        The object whose data we wish to project.

    Example
    -------
    >>> tri = KMeansTri(data)
    >>> dataproj = tri.project(50)
    """

    def __init__(self, data):
        self.data = data

    def project(self, ndim=None):
        """ Projects the data object given to the constructor onto `ndim` dimensions

        Parameters
        ----------
        ndim : int
            The number of dimensions we want to project the data on.

        Returns
        -------
        dataTica : :class:`MetricData <htmd.metricdata.MetricData>` object
            A new :class:`MetricData <htmd.metricdata.MetricData>` object containing the projected data

        Example
        -------
        >>> tri = KMeansTri(data)
        >>> datatri = tri.project(5)
        """
        import scipy.spatial.distance as scidist
        from sklearn.cluster import MiniBatchKMeans
        from htmd.metricdata import MetricData

        datconcat = np.concatenate(self.data.dat)
        mb = MiniBatchKMeans(n_clusters=ndim)
        mb.fit(datconcat)

        # TODO: Could make it into a loop to waste less memory
        dist = scidist.cdist(datconcat, mb.cluster_centers_)
        dist = np.mean(dist, axis=1)[:, np.newaxis] - dist
        dist[dist < 0] = 0

        return MetricData(dat=self.data.deconcatenate(dist), ref=self.data.ref, simlist=self.data.simlist,
                          fstep=self.data.fstep, parent=self.data)
