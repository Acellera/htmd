# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import numpy as np
import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from htmd.metricdata import MetricData


class KMeansTri(object):
    """Class for calculating the KMeans-triangle projections of a MetricData object.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        The object whose data we wish to project.

    Examples
    --------
    >>> tri = KMeansTri(data)
    >>> dataproj = tri.project(50)
    """

    def __init__(self, data: "MetricData"):
        self.data = data

    def project(self, ndim: int | None = None) -> "MetricData":
        """Project the data object given to the constructor onto ``ndim`` dimensions.

        Parameters
        ----------
        ndim : int, optional
            The number of dimensions (cluster centers) to project the data on.

        Returns
        -------
        projdata : :class:`MetricData <htmd.metricdata.MetricData>` object
            A new MetricData object containing the projected data.

        Examples
        --------
        >>> tri = KMeansTri(data)
        >>> datatri = tri.project(50)
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

        return MetricData(
            dat=self.data.deconcatenate(dist),
            ref=self.data.ref,
            simlist=self.data.simlist,
            fstep=self.data.fstep,
            parent=self.data,
        )
