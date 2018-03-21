# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from htmd.kinetics import Rates
import logging
logger = logging.getLogger(__name__)


class KineticsHMM(object):
    _kB = 0.0019872041  # Boltzman constant kcal/(mol K)

    def __init__(self, model, temperature, concentration=1, source=None, sink=None):
        self.model = model
        self.concentration = concentration
        self.temperature = temperature
        self.source = source
        self.sink = sink

        # if model._modelid is None:
        #     raise RuntimeError('You need to call model.markovModel() before calculating kinetics')
        #
        # self._modelid = model._modelid

        if self.source is None:
            logger.info('Detecting source state...')
            self._detectSource()
            logger.info('Source macro = ' + str(model.hmm.active_set[self.source]))
        if self.sink is None:
            logger.info('Detecting sink state...')
            self._detectSink()
            logger.info('Sink macro = ' + str(model.hmm.active_set[self.sink]))
        if model.data.fstep is None:
            raise RuntimeError('Please define a framestep (fstep) for the data object')

    def _detectSource(self):
        if self.model.data.parent is None:
            dataobj = self.model.data
        else:
            dataobj = self.model.data.parent
        St = self.model.hmm.discrete_trajectories_full
        K = self.model.data.K

        distcols = []
        contcols = []
        if not dataobj.description.empty:  # Search for distances or contacts in the data
            distcols = np.where(dataobj.description.type == 'distance')[0]
            contcols = np.where(dataobj.description.type == 'contact')[0]
        if len(distcols) == 0 and len(contcols) == 0:
            raise RuntimeError('Could not detect source state. Please specify it in the Kinetics constructor.')

        if len(distcols) != 0:
            cols = distcols
        else:
            cols = contcols

        averages = np.zeros(self.model.hmm.nstates)
        N = np.zeros(self.model.hmm.nstates)
        dat = dataobj.dat
        for i in range(len(St)):
            if dataobj.trajectories[i].projection.ndim == 2:
                rowsums = np.sum(dat[i][:, cols], axis=1)
            else:
                rowsums = dat[i][:, cols]
            macroSt = self.model.hmm.metastable_assignments[St[i]]
            totalsums = np.bincount(macroSt, weights=rowsums)
            N[0:np.max(macroSt)+1] += np.bincount(macroSt)
            idx = list(set(macroSt))
            averages[idx] += totalsums[idx]
        avg = averages / N
        avg = avg[self.model.hmm.active_set]

        if len(contcols) != 0:
            logger.info('Guessing the source state as the state with minimum contacts.')
            self.source = np.argmin(avg)
        elif len(distcols) != 0:
            logger.info('Guessing the source state as the state with maximum distances.')
            self.source = np.argmax(avg)

    def _detectSink(self):
        if self.source is None:
            raise NameError('Source undefined. Cannot define sink.')
        peq = self.model.eqDistribution(plot=False)
        idx = np.argsort(peq)
        idx = np.delete(idx, np.where(idx == self.source))
        self.sink = idx[-1]

    def getRates(self, source=None, sink=None):
        """ Get the rates between two states

        Parameters
        ----------
        source : int, optional
            The state index to use as source
        sink : int, optional
            The state index to use as sink

        Returns
        -------
        rates : :class:`Rates` object
            A Rates object containing the rates

        Example
        -------
        >>> kin = Kinetics(model, temperature=300, concentration=0.015)
        >>> r = kin.getRates()
        >>> print(r)
        >>> dg = r.g0eq
        """
        actset = self.model.hmm.active_set
        if source is None:
            source = self.source
        else:
            source = np.where(actset == source)[0]

        if sink is None:
            sink = self.sink
        else:
            sink = np.where(actset == sink)[0]
        logger.info('Calculating rates between source: {} and sink: {} states.'.format(actset[source], actset[sink]))
        if source == sink:
            logger.info('Calculating rates between state and itself gives 0')
            r = Rates(); r.mfpton = 0; r.mfptoff=0; r.koff=0; r.kon=0; r.g0eq=0; r.kdeq=0;
            return r

        if source == self.source:  # Apply concentration only on the bulk state
            conc = self.concentration
        elif sink == self.source:
            conc = 1 / self.concentration
        else:
            conc = 1

        model = self.model
        from msmtools.analysis import mfpt
        r = Rates()
        r.mfpton = model.data.fstep * model.hmm.lag * mfpt(self.model.hmm.transition_matrix, origin=source, target=sink)
        r.mfptoff = model.data.fstep * model.hmm.lag * mfpt(self.model.hmm.transition_matrix, origin=sink, target=source)
        r.koff = 1E9 / r.mfptoff
        r.kon = 1E9 / (r.mfpton * conc)
        eq = model.hmm.stationary_distribution
        sinkeq = np.sum(eq[sink])
        sourceeq = np.sum(eq[source])
        if conc != 1:
            logger.info('Concentration correction of {:.2f} kcal/mol.'.format(-self._kBT * np.log(1 / conc)))
        r.g0eq = -self._kBT * np.log(sinkeq / (conc * sourceeq))
        r.kdeq = np.exp(r.g0eq / self._kBT)
        return r

    @property
    def _kBT(self):
        return self._kB * self.temperature



