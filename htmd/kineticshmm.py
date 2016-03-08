# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import warnings
import matplotlib.pyplot as plt
import logging
logger = logging.getLogger(__name__)


class KineticsHMM(object):
    """ Constructor for the Kinetics class

    Parameters
    ----------
    model : :class:`Model <htmd.model.Model>`
        The model from which this class will calculate the kinetics
    concentration : float, optional
        If a ligand is contained in the simulation, give the concentration of the ligand.
    temperature : float, optional
        The simulation temperature
    source : int, optional
        The source macrostate. Default it will be detected as the most disassociated state or the most extended conformation.
    sink : int, optional
        The sink macrostate. Default it will be calculated as the macrostate with the highest equilibrium population.

    Examples
    --------
    >>> kin = KineticsHMM(modelHMM, temperature=300, concentration=0.015)
    """
    _kB = 0.0019872041  # Boltzman constant

    def __init__(self, model, concentration=1, temperature=0, source=None, sink=None):
        self.concentration = concentration
        self.temperature = temperature
        self.model = model
        self.source = source
        self.sink = sink

        if model._modelid is None:
            raise NameError('You need to call model.markovModel() before calculating kinetics')

        self._modelid = model._modelid

        if self.source is None:
            logger.info('Detecting source state...')
            self._detectSource()
            logger.info('Source macro = ' + str(self.source))
        if self.sink is None:
            logger.info('Detecting sink state...')
            self._detectSink()
            logger.info('Sink macro = ' + str(self.sink))
        if model.data.fstep is None:
            raise NameError('Please define a framestep (fstep) for the data object')

    def _detectSource(self):
        dataobj = self.model.data
        if dataobj.parent is None:
            data = dataobj.dat
        else:
            data = dataobj.parent.dat

        averages = np.zeros(dataobj.K)
        for i in range(len(dataobj.St)):
            if data[i].ndim == 2:
                rowsums = np.sum(data[i], axis=1)
            else:
                rowsums = data[i]
            totalsums = np.bincount(dataobj.St[i], weights=rowsums)
            idx = list(set(dataobj.St[i]))
            averages[idx] += totalsums[idx]
        avg = averages / dataobj.N
        avg = avg[self.model.cluster_ofmicro]

        if data[0].dtype == bool:
            logger.info('Guessing the source state as the state with minimum contacts.')
            sourcemicro = np.argmin(avg)
        else:
            logger.info('Guessing the source state as the state with maximum distances.')
            sourcemicro = np.argmax(avg)
        self.source = self.model.macro_ofmicro[sourcemicro]

    def _detectSink(self):
        if self.source is None:
            raise NameError('Source undefined. Cannot define sink.')
        peq = self.model.eqDistribution(plot=False)
        idx = np.argsort(peq)
        idx = np.delete(idx, np.where(idx == self.source))
        self.sink = idx[-1]

    def getRates(self, source=None, sink=None, states='macro'):
        """ Get the rates between two states

        Parameters
        ----------
        source : int, optional
            The state index to use as source
        sink : int, optional
            The state index to use as sink
        states : ['macro','micro'], optional
            The state type of the states given before

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
        self._intergrityCheck()
        if source is None:
            source = self.source
        if sink is None:
            sink = self.sink
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
        r.mfpton = model.data.fstep * model.lag * mfpt(self.model.P, origin=source, target=sink)
        r.mfptoff = model.data.fstep * model.lag * mfpt(self.model.P, origin=sink, target=source)
        r.koff = 1E9 / r.mfptoff
        r.kon = 1E9 / (r.mfpton * conc)
        eq = model.hmm.stationary_distribution
        sinkeq = eq[sink]
        sourceeq = eq[source]
        r.g0eq = -self._kBT * np.log(sinkeq / (conc * sourceeq))
        r.kdeq = np.exp(r.g0eq / self._kBT)
        return r

    def plotRates(self, rates=('mfptoff', 'mfpton', 'g0eq')):
        """ Plot the MFPT on, off and DG of all the macrostates to the sink state

        Parameters
        ----------
        rates : tuple
            Specify which rates you want to plot. Options are ('mfptoff','mfpton','g0eq','kdeq','kon','koff')

        Examples
        --------
        >>> kin = Kinetics(model, temperature=300, concentration=0.015)
        >>> kin.plotRates()
        """
        self._intergrityCheck()
        macronum = self.model.macronum
        mfptoff = np.zeros(macronum)
        mfpton = np.zeros(macronum)
        dg = np.zeros(macronum)
        kon = np.zeros(macronum)
        koff = np.zeros(macronum)
        kdeq = np.zeros(macronum)

        for m in range(macronum):
            if m == self.source:
                continue
            r = self.getRates(source=self.source, sink=m)
            mfptoff[m] = r.mfptoff
            mfpton[m] = r.mfpton
            dg[m] = r.g0eq
            kon[m] = r.kon
            koff[m] = r.koff
            kdeq[m] = r.kdeq

        plt.ion()  # Interactive figure mode on
        if 'mfptoff' in rates:
            self._plotRate(mfptoff, 'MFPT off (ns)', log=True, lim1=True)
        if 'mfpton' in rates:
            self._plotRate(mfpton, 'MFPT on (ns)', log=True, lim1=True)
        if 'g0eq' in rates:
            self._plotRate(dg, 'Standard free energy (kcal/mol)')
        if 'kon' in rates:
            self._plotRate(kon, 'Kon (1/M 1/s)', log=True, lim1=True)
        if 'koff' in rates:
            self._plotRate(koff, 'Koff (1/s)', log=True, lim1=True)
        if 'kdeq' in rates:
            self._plotRate(kdeq, 'Kd (M)')

    def _plotRate(self, rate, ylabel, log=False, lim1=False):
            plt.figure()
            plt.bar(range(self.model.macronum), rate)
            plt.ylabel(ylabel)
            plt.xlabel('Macrostates')
            plt.xticks(np.arange(0.4, self.model.macronum+0.4, 1), range(self.model.macronum))
            if log:
                plt.yscale('log')
            if lim1:
                plt.ylim(ymin=1)
            plt.show()

    def mfptGraph(self):
        """ Not yet implemented
        """
        import pyemma.plots as mplt
        self._intergrityCheck()
        #pos = np.array([[3,3],[4.25,0],[0,1],[1.75,0],[6,1.0]])
        #state_colors = ['green', 'blue', 'yellow', 'cyan', 'purple']
        # Color states by committor! on gradient from blue to red
        hmm = self.model.hmm
        state_sizes = hmm.stationary_distribution ** 0.25  # Scale eq prob down to make states visible
        fig, pos = mplt.plot_markov_model(hmm, state_sizes=state_sizes)

    def fluxPathways(self):
        import pyemma.plots as mplt
        import pyemma.msm as msm
        import matplotlib.pylab as plt
        self._intergrityCheck()

        plt.ion()
        tpt = msm.tpt(self.model.coarsemsm, [self.source], [self.sink])
        mplt.plot_flux(tpt[1], show_committor=False)
        plt.show()


    @property
    def _kBT(self):
        return self._kB * self.temperature

    def _intergrityCheck(self):
        if self.model._modelid != self._modelid:
            raise NameError('After updating the Model object you need to call the Kinetics constructor again')


class Rates(object):
    """ A class containing a set of rates

    Attributes
    ----------
    mfpton : float
        The mean first passage time of going from source to sink
    mfptoff : float
        The mean first passage time of going from sink to source
    kon : float
        The Kon rate (association constant) from source to sink
    koff : float
        The Koff rate (dissociation constant) from sink to source
    kdeq : float
        The Kd, calculated from the equilibrium probability
    g0eq : float
        The free energy between source and sink, calculated from the equilibrium probability
    """
    def __str__(self):
        s = ''
        s += 'mfpton = {:.2E} (ns)\n'.format(self.mfpton)
        s += 'mfptoff = {:.2E} (ns)\n'.format(self.mfptoff)
        s += 'kon = {:.2E} (1/M 1/s)\n'.format(self.kon)
        s += 'koff = {:.2E} (1/s)\n'.format(self.koff)
        s += 'kdeq = {:.2E} (M)\n'.format(self.kdeq)
        s += 'g0eq = {:.2f} (kcal/mol)\n'.format(self.g0eq)
        return s
