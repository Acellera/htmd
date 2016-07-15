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
from htmd.units import convert as unitconvert
from joblib import Parallel, delayed
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
    units : str
        The units of lag. Can be 'frames' or any time unit given as a string.

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

    def __init__(self, data, lag, units='frames'):
        from pyemma.coordinates import tica
        # data.dat.tolist() might be better?
        self.data = data
        if isinstance(data, Metric):
            if units != 'frames':
                raise RuntimeError('Cannot use delayed projection TICA with units other than frames for now. Report this to HTMD issues.')
            metr = data
            from pyemma.coordinates.transform.tica import TICA
            self.tic = TICA(lag)

            p = ProgressBar(len(metr.simulations))
            for proj in _projectionGenerator(metr, _getNcpus()):
                for pro in proj:
                    self.tic.partial_fit(pro[0])
                p.progress(len(proj))
            p.stop()
        else:
            lag = unitconvert(units, 'frames', lag, data.fstep)
            if lag == 0:
                raise RuntimeError('Lag time conversion resulted in 0 frames. Please use a larger lag-time for TICA.')
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
            # self.tic._dim = ndim  # Old way of doing it. Deprecated since pyEMMA 2.1
            self.tic.set_params(dim=ndim)  # Change to this in 2.1 pyEMMA version

        if isinstance(self.data, Metric):  # Doesn't project on correct number of dimensions
            proj = []
            refs = []
            fstep = None

            metr = self.data
            p = ProgressBar(len(metr.simulations))
            k = -1
            droppedsims = []
            for projecteddata in _projectionGenerator(metr, _getNcpus()):
                for pro in projecteddata:
                    k += 1
                    if pro is None:
                        droppedsims.append(k)
                        continue
                    proj.append(self.tic.transform(pro[0]))
                    refs.append(pro[1])
                    if fstep is None:
                        fstep = pro[2]
                p.progress(len(projecteddata))
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

        from htmd.metricdata import MetricData
        datatica = MetricData(dat=np.array(proj, dtype=object), simlist=simlist, ref=ref, fstep=fstep, parent=parent)

        return datatica


def _projectionGenerator(metric, ncpus):
    for i in range(0, len(metric.simulations), ncpus):
        simrange = range(i, np.min((i+ncpus, len(metric.simulations))))
        results = Parallel(n_jobs=ncpus, verbose=0)(delayed(_projector)(metric, i) for i in simrange)
        yield results


def _projector(metric, i):
    return metric._projectSingle(i)


def _getNcpus():
    from htmd.config import _config
    ncpus = _config['ncpus']
    if ncpus < 0:
        import multiprocessing
        ncpus = multiprocessing.cpu_count() + ncpus + 1
    return ncpus