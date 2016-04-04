# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import warnings
import numpy as np
import random
from htmd.projections.metric import Metric
from htmd.progress.progress import ProgressBar
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
    Perez-Hernandez, G. and Paul, F. and Giorgino, T. and de Fabritiis, G.
    and Noe, F. (2013) Identification of slow molecular order parameters
    for Markov model construction. J. Chem. Phys., 139 . 015102.
    """

    def __init__(self, data, lag):
        from pyemma.coordinates import tica
        # data.dat.tolist() might be better?
        self.data = data
        if isinstance(data, Metric):
            from pyemma.coordinates.transform.tica import TICA
            self.tic = TICA(lag)

            p = ProgressBar(len(data.simulations))
            for i in range(len(data.simulations)):
                # Fix for pyemma bug. Remove eventually:
                d, _, _ = data._projectSingle(i)
                if d is None or d.shape[0] < lag:
                    continue
                self.tic.partial_fit(d)
                p.progress()
            p.stop()
        else:
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
            self.tic.set_params(dim=ndim)

        if isinstance(self.data, Metric):  # Doesn't project on correct number of dimensions
            proj = []
            refs = []
            fstep = None

            '''from htmd.config import _config
            from joblib import Parallel, delayed
            results = Parallel(n_jobs=_config['ncpus'], verbose=11)(
                delayed(_test)(self.data, self.tic, i) for i in range(len(self.data.simulations)))

            for i in range(len(results)):
                proj.append(results[i][0])
                refs.append(results[i][1])
                fstep.append(results[i][2])'''

            droppedsims = []
            p = ProgressBar(len(self.data.simulations))
            for i in range(len(self.data.simulations)):
                d, r, f = self.data._projectSingle(i)
                if d is None:
                    droppedsims.append(i)
                    continue
                if fstep is None:
                    fstep = f
                refs.append(r)
                proj.append(self.tic.transform(d))
                p.progress()
            p.stop()
            simlist = self.data.simulations
            simlist = np.delete(simlist, droppedsims)
            ref = np.array(refs, dtype=object)
            #fstep = 0
            parent = None
        else:
            proj = self.tic.get_output()
            simlist = self.data.simlist
            ref = self.data.ref
            fstep = self.data.fstep
            parent = self.data

        if ndim is None:
            logger.info('Kept {} dimension(s) to cover 95% of kinetic variance.'.format(self.tic.dimension()))
        #print(np.shape(proj))


        from htmd.metricdata import MetricData
        datatica = MetricData(dat=np.array(proj, dtype=object), simlist=simlist, ref=ref, fstep=fstep, parent=parent)

        '''datatica = self.data.copy()
        #datatica.dat = self.data.deconcatenate(np.squeeze(proj))
        datatica.dat = np.array(proj, dtype=object)
        datatica.parent = self.data
        datatica.St = None
        datatica.Centers = None
        datatica.N = None
        datatica.K = None
        datatica._dataid = random.random()
        datatica._clusterid = None'''
        return datatica


def _test(data, tic, i):
    d, r, f = data._projectSingle(i)
    if d is None:
        return
    else:
        return tic.transform(d), r, f
