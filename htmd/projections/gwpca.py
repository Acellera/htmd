# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import numpy as np
from htmd.units import convert as unitconvert
import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from htmd.metricdata import MetricData


class GWPCA(object):
    """Class for calculating the globally-weighted PCA projections of a MetricData object.

    References
    ----------

    * N. Blöchliger, A. Caflisch, and A. Vitalis. Weighted Distance Functions Improve Analysis of High-Dimensional
      Data: Application to Molecular Dynamics Simulations, J. Chem. Theory Comput., 2015, 11 (11), pp 5481-5492. doi:
      `10.1021/acs.jctc.5b00618`_

    .. _10.1021/acs.jctc.5b00618: http://dx.doi.org/10.1021/acs.jctc.5b00618

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        The object whose data we wish to project.
    lag : float
        The correlation lagtime used to compute autocorrelation weights.
    units : str, optional
        The units of ``lag``. Can be ``'frames'`` or any time unit given as a string.

    Examples
    --------
    >>> gw = GWPCA(data, lag=10)
    >>> dataproj = gw.project(5)
    """

    def __init__(self, data: "MetricData", lag: float, units: str = "frames"):
        lag = unitconvert(units, "frames", lag, data.fstep)
        if lag == 0:
            raise RuntimeError(
                "Lag time conversion resulted in 0 frames. Please use a larger lag-time for TICA."
            )
        self.data = data

        datconcat = np.concatenate(self.data.dat)
        self.weights = self._autocorrelation(datconcat, lag)

    def _autocorrelation(self, data, lag):
        # Calculate the autocorrelation of each feature
        autocorr = []
        for i in range(data.shape[1]):
            autocorr.append(np.correlate(data[0:-lag, i], data[lag:, i]))
        return np.squeeze(autocorr)

    def project(self, ndim: int | None = None) -> "MetricData":
        """Project the data object given to the constructor onto ``ndim`` dimensions.

        Parameters
        ----------
        ndim : int, optional
            The number of dimensions to project the data on.

        Returns
        -------
        projdata : :class:`MetricData <htmd.metricdata.MetricData>` object
            A new MetricData object containing the projected data.

        Examples
        --------
        >>> gw = GWPCA(data, lag=10)
        >>> dataproj = gw.project(5)
        """
        from sklearn.decomposition import IncrementalPCA
        from tqdm import tqdm

        pca = IncrementalPCA(n_components=ndim, batch_size=10000)
        for t in tqdm(self.data.trajectories):
            pca.partial_fit(t.projection * self.weights)

        projdata = self.data.copy()
        for i, t in enumerate(tqdm(self.data.trajectories)):
            projdata.trajectories[i].projection = pca.transform(
                t.projection * self.weights
            )

        # projdataconc = pca.fit_transform(self.weighedconcat)
        # projdata.dat = projdata.deconcatenate(projdataconc)
        return projdata
