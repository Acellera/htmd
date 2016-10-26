# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import warnings
import logging
logger = logging.getLogger(__name__)


class Kinetics(object):
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
    >>> kin = Kinetics(model, temperature=300, concentration=0.015)

    .. currentmodule:: htmd.kinetics.Kinetics
    .. rubric:: Methods
    .. autoautosummary:: htmd.kinetics.Kinetics
        :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.kinetics.Kinetics
        :attributes:
    """
    _kB = 0.0019872041  # Boltzman constant kcal/(mol K)

    def __init__(self, model, temperature, concentration=1, source=None, sink=None):
        self.concentration = concentration
        self.temperature = temperature
        self.model = model
        self.source = source
        self.sink = sink
        if source is not None:
            self.sourcemicro = np.where(model.macro_ofmicro == self.source)[0]
        if sink is not None:
            self.sinkmicro = np.where(model.macro_ofmicro == self.sink)[0]

        if model._modelid is None:
            raise RuntimeError('You need to call model.markovModel() before calculating kinetics')

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
            raise RuntimeError('Please define a framestep (fstep) for the data object')

    def _detectSource(self):
        dataobj = self.model.data
        distcols = []
        contcols = []

        if not dataobj.map.empty:  # Search for distances or contacts in the data
            if dataobj.parent is None:
                distcols = np.where(dataobj.map.type == 'distance')[0]
                contcols = np.where(dataobj.map.type == 'contact')[0]
            else:
                distcols = np.where(dataobj.parent.map.type == 'distance')[0]
                contcols = np.where(dataobj.parent.map.type == 'contact')[0]

        if len(distcols) == 0 and len(contcols) == 0:
            raise RuntimeError('Could not detect source state. Please specify it in the Kinetics constructor.')

        if len(distcols) != 0:
            cols = distcols
        else:
            cols = contcols

        if dataobj.parent is None:
            data = dataobj.dat
        else:
            data = dataobj.parent.dat

        averages = np.zeros(dataobj.K)
        for i in range(len(dataobj.St)):
            if data[i].ndim == 2:
                rowsums = np.sum(data[i][:, cols], axis=1)
            else:
                rowsums = data[i][:, cols]
            totalsums = np.bincount(dataobj.St[i], weights=rowsums)
            idx = list(set(dataobj.St[i]))
            averages[idx] += totalsums[idx]
        avg = averages / dataobj.N
        avg = avg[self.model.cluster_ofmicro]

        if len(contcols) != 0:
            logger.info('Guessing the source state as the state with minimum contacts.')
            sourcemicro = np.argmin(avg)
        elif len(distcols) != 0:
            logger.info('Guessing the source state as the state with maximum distances.')
            sourcemicro = np.argmax(avg)
        self.source = self.model.macro_ofmicro[sourcemicro]
        self.sourcemicro = sourcemicro

    def _detectSink(self):
        if self.source is None:
            raise NameError('Source undefined. Cannot define sink.')
        peq = self.model.eqDistribution(plot=False)
        idx = np.argsort(peq)
        idx = np.delete(idx, np.where(idx == self.source))
        self.sink = idx[-1]
        nonsource = np.where(idx != self.source)[0]
        self.sinkmicro = nonsource[np.argmax(self.model.msm.stationary_distribution[nonsource])]

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
        logger.info('Calculating rates between source: {} and sink: {} states.'.format(source, sink))
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
        if states == 'macro':
            source = np.where(model.macro_ofmicro == source)[0]
            sink = np.where(model.macro_ofmicro == sink)[0]

        from msmtools.analysis import mfpt
        r = Rates()
        r.mfpton = model.data.fstep * model.lag * mfpt(self.model.P, origin=source, target=sink)
        r.mfptoff = model.data.fstep * model.lag * mfpt(self.model.P, origin=sink, target=source)
        r.koff = 1E9 / r.mfptoff
        r.kon = 1E9 / (r.mfpton * conc)
        eq = model.msm.stationary_distribution
        sinkeq = np.sum(eq[sink])
        sourceeq = np.sum(eq[source])
        if conc != 1:
            logger.info('Concentration correction of {:.2f} kcal/mol.'.format(-self._kBT * np.log(1 / conc)))
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

        import matplotlib.pyplot as plt
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
        import matplotlib.pyplot as plt
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
        msm = self.model.msm
        state_sizes = msm.stationary_distribution ** 0.25  # Scale eq prob down to make states visible
        fig, pos = mplt.plot_markov_model(msm, state_sizes=state_sizes)

    def plotFluxPathways(self, statetype='macro', mode='net_flux', fraction=1.0):
        """ Plot flux pathways between source and sink state.

        Parameters
        ----------
        statetype : {'macro','coarse','micro'}
            What type of states to plot
        mode : {'net_flux', 'gross_flux'}
            Type of fluxes to plot
        fraction : float
            Fraction of fluxes for which to report pathways. Doesn't change the plot, only the text output.
        """
        # Make mode a radio button with interactive plot
        from pyemma import msm
        from pyemma.plots import plot_flux
        self._intergrityCheck()

        if statetype == 'micro':
            tpt = msm.tpt(self.model.msm, [self.sourcemicro], [self.sinkmicro])
            plot_flux(tpt, attribute_to_plot=mode)
        elif statetype == 'macro' or statetype == 'coarse':
            metastable_sets = []
            for i in range(self.model.macronum):
                metastable_sets.append(np.where(self.model.macro_ofmicro == i)[0])
            tpt = msm.tpt(self.model.msm, metastable_sets[self.source], metastable_sets[self.sink])
            #from IPython.core.debugger import Tracer
            #Tracer()()
            newsets, tpt = tpt.coarse_grain(metastable_sets)
            setmap = []
            # getting the mapping of new sets to old sets
            for set1 in newsets:
                for idx2, set2 in enumerate(metastable_sets):
                    if set(set1) == set(set2):
                        setmap.append(idx2)
                        continue
            setmap = np.array(setmap)
            plot_flux(tpt, attribute_to_plot=mode, state_labels=setmap)

        paths, pathfluxes = tpt.pathways(fraction=fraction)
        cumflux = 0
        print("Path flux\t\t%path\t%of total\tpath")
        for i in range(len(paths)):
            cumflux += pathfluxes[i]
            print('{}\t{:3.1f}%\t{:3.1f}%\t\t{}'.format(
                    pathfluxes[i], 100.0*pathfluxes[i]/tpt.total_flux,
                    100.0*cumflux/tpt.total_flux, setmap[paths[i]]))
            #print(pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t\t',paths[i])

        # transition rates
        # plot_markov_model(P)


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
    def __repr__(self):
        s = ''
        s += 'mfpton = {:.2E} (ns)\n'.format(self.mfpton)
        s += 'mfptoff = {:.2E} (ns)\n'.format(self.mfptoff)
        s += 'kon = {:.2E} (1/M 1/s)\n'.format(self.kon)
        s += 'koff = {:.2E} (1/s)\n'.format(self.koff)
        s += 'koff/kon = {:.2E} (M)\n'.format(self.koff/self.kon)
        s += 'kdeq = {:.2E} (M)\n'.format(self.kdeq)
        s += 'g0eq = {:.2f} (kcal/mol)\n'.format(self.g0eq)
        return s
