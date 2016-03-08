# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import warnings
import numpy as np
import random
import logging
logger = logging.getLogger(__name__)


class TICA(object):
    """ Class for calculating the TICA projections of a MetricData  object

    Time-based Independent Component Analysis
    Projects your data on the slowest coordinates identified for a
    given lagtime.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        The object whose data we wish to project onto the top TICA dimensions
    lag : int
        The correlation lagtime to use for TICA

    Example
    -------
    >>> from htmd.projections.tica import TICA
    >>> tica = TICA(data,20)

    References
    ----------
    Perez-Hernandez, G. and Paul, F. and Giogino, T. and de Fabritiis, G.
    and Noe, F. (2013) Identification of slow molecular order parameters
    for Markov model construction. J. Chem. Phys., 139 . 015102.
    """

    def __init__(self, data, lag):
        from pyemma.coordinates import tica
        # data.dat.tolist() might be better?
        self.data = data
        self.tic = tica(data.dat.tolist(), lag=lag)

    def project(self, ndim=None):
        """ Projects the data object given to the constructor onto the top `ndim` TICA dimensions

        Parameters
        ----------
        ndim : int
            The number of TICA dimensions we want to project the data on. If None is given it will use choose a number
            of dimensions to cover 95% of the kinetic variance.

        Returns
        -------
        dataTica : :class:`MetricData <htmd.metricdata.MetricData>` object
            A new :class:`MetricData <htmd.metricdata.MetricData>` object containing the TICA projected data

        Example
        -------
        >>> from htmd.projections.tica import TICA
        >>> tica = TICA(data,20)
        >>> dataTica = tica.project(5)
        """
        if ndim is not None:
            self.tic._dim = ndim
        proj = self.tic.get_output()
        if ndim is None:
            logger.info('Kept {} dimension(s) to cover 95% of kinetic variance.'.format(self.tic.dimension()))
        #print(np.shape(proj))
        datatica = self.data.copy()
        #datatica.dat = self.data.deconcatenate(np.squeeze(proj))
        datatica.dat = np.array(proj, dtype=object)
        datatica.parent = self.data
        datatica.St = None
        datatica.Centers = None
        datatica.N = None
        datatica.K = None
        datatica._dataid = random.random()
        datatica._clusterid = None
        return datatica
