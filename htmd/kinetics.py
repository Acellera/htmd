# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import logging
import unittest

logger = logging.getLogger(__name__)


class Kinetics(object):
    """Constructor for the Kinetics class

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
            raise RuntimeError(
                "You need to call model.markovModel() before calculating kinetics"
            )

        self._modelid = model._modelid

        if self.source is None:
            logger.info("Detecting source state...")
            self._detectSource()
            logger.info("Source macro = " + str(self.source))
        if self.sink is None:
            logger.info("Detecting sink state...")
            self._detectSink()
            logger.info("Sink macro = " + str(self.sink))
        if model.data.fstep is None:
            raise RuntimeError("Please define a framestep (fstep) for the data object")

    def _detectSource(self):
        dataobj = self.model.data
        # This is an ordered tuple based on priority
        supportedProj = ("contact", "distance", "rmsd", "secondary structure")
        cols = []
        detectedProj = None

        if dataobj.parent is None:
            mapping = dataobj.description
            data = dataobj.dat
        else:
            mapping = dataobj.parent.description
            data = dataobj.parent.dat

        if (
            not dataobj.description.empty
        ):  # Search for supported projections in the data
            for proj in supportedProj:
                cols = np.where(mapping.type == proj)[0]
                if len(cols) != 0:
                    detectedProj = proj
                    break

        if len(cols) == 0:
            raise RuntimeError(
                "Could not detect source state. Please specify it in the Kinetics constructor. "
                f"This class only detects automatic detection for the following metrics. {supportedProj}"
            )

        averages = np.zeros(dataobj.K)
        St = dataobj.St
        for i in range(dataobj.numTrajectories):
            if data[i].ndim == 2:
                rowsums = np.sum(data[i][:, cols], axis=1)
            else:
                rowsums = data[i][:, cols]
            totalsums = np.bincount(St[i], weights=rowsums)
            idx = list(set(St[i]))
            averages[idx] += totalsums[idx]
        avg = averages / dataobj.N
        avg = avg[self.model.cluster_ofmicro]

        if detectedProj == "contact" or detectedProj == "secondary structure":
            logger.info(
                "Guessing the source state as the state with minimum {}s.".format(
                    detectedProj
                )
            )
            sourcemicro = np.argmin(avg)
        elif detectedProj == "distance" or detectedProj == "rmsd":
            logger.info(
                "Guessing the source state as the state with maximum {}s.".format(
                    detectedProj
                )
            )
            sourcemicro = np.argmax(avg)
        # TODO: I could also detect the source by largest variance
        self.source = self.model.macro_ofmicro[sourcemicro]
        self.sourcemicro = np.array([sourcemicro])

    def _detectSink(self):
        if self.source is None:
            raise NameError("Source undefined. Cannot define sink.")
        peq = self.model.eqDistribution(plot=False)
        idx = np.argsort(peq)
        idx = np.delete(idx, np.where(idx == self.source))
        self.sink = idx[-1]
        nonsource = np.where(idx != self.source)[0]
        self.sinkmicro = np.array(
            [nonsource[np.argmax(self.model.msm.stationary_distribution[nonsource])]]
        )

    def getRates(self, source=None, sink=None, states="macro", _logger=True):
        """Get the rates between two (sets of) states

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
        >>> r = kin.getRates(source=4, sink=[0,1,2,3])
        """
        import numbers
        from deeptime.markov.tools.analysis import mfpt

        self._intergrityCheck()
        if source is None:
            source = self.source
        if sink is None:
            sink = self.sink
        if isinstance(source, numbers.Integral):
            source = [
                source,
            ]
        if isinstance(sink, numbers.Integral):
            sink = [
                sink,
            ]
        if _logger:
            logger.info(
                "Calculating rates between source: {} and sink: {} states.".format(
                    source, sink
                )
            )

        if len(np.intersect1d(source, sink)) != 0:
            if _logger:
                logger.info("Calculating rates between state and itself gives 0")
            r = Rates()
            return r

        if self.source in source:  # Apply concentration only on the bulk state
            conc = self.concentration
        elif self.source in sink:  # Invert concentration is bulk state is in sink
            if _logger:
                logger.info(
                    "Bulk state detected in sink. Applying concentration correction to sink instead of source."
                )
            conc = 1 / self.concentration
        else:
            conc = 1

        model = self.model
        if states == "macro":  # Finding the microstates of the macrostates
            eq = model.eqDistribution(
                plot=False
            )  # If macro, use the membership probs to calculate eq prob
            sinkeq = np.sum(
                eq[sink]
            )  # sink and source might be multiple macros so we have to sum them
            sourceeq = np.sum(eq[source])
            sourcemicros = []
            for s in source:
                sourcemicros += list(np.where(model.macro_ofmicro == s)[0])
            sinkmicros = []
            for s in sink:
                sinkmicros += list(np.where(model.macro_ofmicro == s)[0])
        elif states == "micro":
            eq = model.msm.stationary_distribution
            sinkeq = np.sum(eq[sink])
            sourceeq = np.sum(eq[source])
            sourcemicros = source
            sinkmicros = sink

        r = Rates()
        r.mfpton = (
            model.data.fstep
            * model.lag
            * mfpt(self.model.P, origin=sourcemicros, target=sinkmicros)
        )
        r.mfptoff = (
            model.data.fstep
            * model.lag
            * mfpt(self.model.P, origin=sinkmicros, target=sourcemicros)
        )
        r.koff = 1e9 / r.mfptoff
        r.kon = 1e9 / (r.mfpton * conc)
        if conc != 1:
            if _logger:
                logger.info(
                    "Concentration correction of {:.2f} kcal/mol.".format(
                        -self._kBT * np.log(1 / conc)
                    )
                )
        r.g0eq = -self._kBT * np.log(sinkeq / (conc * sourceeq))
        r.kdeq = np.exp(r.g0eq / self._kBT)
        return r

    def plotRates(self, rates=("mfptoff", "mfpton", "g0eq")):
        """Plot the MFPT on, off and DG of all the macrostates to the sink state

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
            r = self.getRates(source=self.source, sink=m, _logger=False)
            mfptoff[m] = r.mfptoff
            mfpton[m] = r.mfpton
            dg[m] = r.g0eq
            kon[m] = r.kon
            koff[m] = r.koff
            kdeq[m] = r.kdeq

        import matplotlib.pyplot as plt

        plt.ion()  # Interactive figure mode on
        if "mfptoff" in rates:
            self._plotRate(mfptoff, "MFPT off (ns)", log=True, lim1=True)
        if "mfpton" in rates:
            self._plotRate(mfpton, "MFPT on (ns)", log=True, lim1=True)
        if "g0eq" in rates:
            self._plotRate(dg, "Standard free energy (kcal/mol)")
        if "kon" in rates:
            self._plotRate(kon, "Kon (1/M 1/s)", log=True, lim1=True)
        if "koff" in rates:
            self._plotRate(koff, "Koff (1/s)", log=True, lim1=True)
        if "kdeq" in rates:
            self._plotRate(kdeq, "Kd (M)")

    def _plotRate(self, rate, ylabel, log=False, lim1=False):
        import matplotlib.pyplot as plt

        plt.figure()
        plt.bar(np.arange(self.model.macronum) + 0.4, rate)
        plt.ylabel(ylabel)
        plt.xlabel("Macrostates")
        plt.xticks(np.arange(self.model.macronum) + 0.4, range(self.model.macronum))
        if log:
            plt.yscale("log")
        if lim1:
            plt.ylim(ymin=1)
        plt.show()

    def plotMarkovModel(self, plot=True, save=None):
        """Plot graph of transition probabilities

        Parameters
        ----------
        plot : bool
            If set it False the plot will not show up in a figure
        save : str
            If a path is passed to save, the plot will be saved to the specified file
        """
        from deeptime.plots import plot_markov_model
        from matplotlib import pylab as plt

        self._intergrityCheck()
        # pos = np.array([[3,3],[4.25,0],[0,1],[1.75,0],[6,1.0]])
        # state_colors = ['green', 'blue', 'yellow', 'cyan', 'purple']
        # Color states by committor! on gradient from blue to red
        msm = self.model.msm
        state_sizes = (
            msm.stationary_distribution**0.25
        )  # Scale eq prob down to make states visible
        plot_markov_model(msm, state_sizes=state_sizes)

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
        if plot:
            plt.show()

    def plotFluxPathways(
        self,
        statetype="macro",
        mode="net_flux",
        fraction=1.0,
        plot=True,
        save=None,
        results=False,
    ):
        """Plot flux pathways between source and sink state.

        The flux is in units of transition events per lag time used to construct the Markov model.

        Parameters
        ----------
        statetype : {'macro','coarse','micro'}
            What type of states to plot
        mode : {'net_flux', 'gross_flux'}
            Type of fluxes to plot
        fraction : float
            Fraction of fluxes for which to report pathways. Doesn't change the plot, only the text output.
        plot : bool
            If set it False the plot will not show up in a figure
        save : str
            If a path is passed to save, the plot will be saved to the specified file
        results : bool
            Set to True to return fluxes, paths and path_fluxes
        """
        # Make mode a radio button with interactive plot
        from deeptime.plots import plot_flux
        from matplotlib import pylab as plt

        self._intergrityCheck()

        plt.figure()

        if statetype == "micro":
            setmap = np.arange(self.model.micronum)
            flux = self.model.msm.reactive_flux(
                self.sourcemicro.tolist(), self.sinkmicro.tolist()
            )
        elif statetype == "macro" or statetype == "coarse":
            metastable_sets = []
            for i in range(self.model.macronum):
                metastable_sets.append(
                    np.where(self.model.macro_ofmicro == i)[0].tolist()
                )

            flux = self.model.msm.reactive_flux(
                metastable_sets[self.source], metastable_sets[self.sink]
            )
            newsets, flux = flux.coarse_grain(metastable_sets)
            setmap = []
            # getting the mapping of new sets to old sets
            for set1 in newsets:
                for idx2, set2 in enumerate(metastable_sets):
                    if set(set1) == set(set2):
                        setmap.append(idx2)
                        continue
            setmap = np.array(setmap)

        ax, _ = plot_flux(flux, attribute_to_plot=mode, state_labels=setmap)
        ax.set_aspect("equal")
        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight", pad_inches=0.2)
        if plot:
            plt.show()

        paths, pathfluxes = flux.pathways(fraction=fraction)
        cumflux = 0
        print("Path flux\t\t%path\t%of total\tpath")
        for i in range(len(paths)):
            cumflux += pathfluxes[i]
            print(
                "{}\t{:3.1f}%\t{:3.1f}%\t\t{}".format(
                    pathfluxes[i],
                    100.0 * pathfluxes[i] / flux.total_flux,
                    100.0 * cumflux / flux.total_flux,
                    setmap[paths[i]],
                )
            )
        if results:
            return flux, paths, pathfluxes

    @property
    def _kBT(self):
        return self._kB * self.temperature

    def _intergrityCheck(self):
        if self.model._modelid != self._modelid:
            raise NameError(
                "After updating the Model object you need to call the Kinetics constructor again"
            )


class Rates(object):
    """A class containing a set of rates

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

    def __init__(
        self, mfpton=None, mfptoff=None, kon=None, koff=None, kdeq=None, g0eq=None
    ):
        if mfpton is None:
            self.mfpton = 0
        if mfptoff is None:
            self.mfptoff = 0
        if kon is None:
            self.kon = 0
        if koff is None:
            self.koff = 0
        if kdeq is None:
            self.kdeq = 0
        if g0eq is None:
            self.g0eq = 0

    def __add__(self, other):
        if isinstance(other, Rates):
            for k in self.__dict__.keys():
                self.__dict__[k] += other.__dict__[k]
        else:
            for k in self.__dict__.keys():
                self.__dict__[k] += other

    def __truediv__(self, other):
        if isinstance(other, Rates):
            for k in self.__dict__.keys():
                self.__dict__[k] /= other.__dict__[k]
        else:
            for k in self.__dict__.keys():
                self.__dict__[k] /= other

    def __repr__(self):
        s = ""
        s += "mfpton = {:.2E} (ns)\n".format(self.mfpton)
        s += "mfptoff = {:.2E} (ns)\n".format(self.mfptoff)
        s += "kon = {:.2E} (1/M 1/s)\n".format(self.kon)
        s += "koff = {:.2E} (1/s)\n".format(self.koff)
        if self.kon != 0:
            s += "koff/kon = {:.2E} (M)\n".format(self.koff / self.kon)
        s += "kdeq = {:.2E} (M)\n".format(self.kdeq)
        s += "g0eq = {:.2f} (kcal/M)\n".format(self.g0eq)
        return s


class _TestKinetics(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.trans_prob = np.array([[0.6, 0.2, 0.2], [0.3, 0.4, 0.3], [0.2, 0.3, 0.5]])

    def test_kinetics(self):
        from htmd.metricdata import _generate_toy_data
        from htmd.model import Model
        import tempfile
        import os

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)

        kin = Kinetics(model, 300, source=0)
        assert kin.sink == 2
        assert np.array_equal(kin.sinkmicro, [0])
        assert kin.source == 0
        assert np.array_equal(kin.sourcemicro, [1])

        rates = kin.getRates(0, 2, states="macro")
        assert np.allclose(rates.mfpton, 3.8699044076299556)
        assert np.allclose(rates.mfptoff, 4.312080466715925)
        assert np.allclose(rates.kon, 258404315.6281552)
        assert np.allclose(rates.koff, 231906618.56122518)
        assert np.allclose(rates.kdeq, 0.7542502489709302)
        assert np.allclose(rates.g0eq, -0.16813599009752134)

        rates = kin.getRates(0, 2, states="micro")
        assert np.allclose(rates.mfpton, 4.5232360505688645)
        assert np.allclose(rates.mfptoff, 4.17783500152719)
        assert np.allclose(rates.kon, 221080657.48066255)
        assert np.allclose(rates.koff, 239358423.59366855)
        assert np.allclose(rates.kdeq, 1.1703043390973993)
        assert np.allclose(rates.g0eq, 0.09375460063521723)

        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            kin.plotFluxPathways(
                statetype="macro",
                mode="net_flux",
                fraction=0.1,
                plot=False,
                save=outplot,
            )
            assert os.path.exists(outplot)

            outplot = os.path.join(tmpdir, "mfpt.png")
            kin.plotMarkovModel(plot=False, save=outplot)
            assert os.path.exists(outplot)


if __name__ == "__main__":
    unittest.main(verbosity=2)
