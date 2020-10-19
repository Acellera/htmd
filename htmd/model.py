"""
Markov state models are a statistical tool for analysing molecular simulations which has met with lots of success.
The Model class here, encapsulates all functionallity for the calculation of Markov models while hiding unnecessary
details under the hood. It uses PyEMMA [1] internally to calculate Markov models.

References
----------
.. [1] PyEMMA 2: A Software Package for Estimation, Validation, and Analysis of Markov Models. Martin K. Scherer et al. JCTC 2015.
"""
# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy import stats
import warnings
import random
from moleculekit.molecule import Molecule
from moleculekit.vmdviewer import getCurrentViewer
from htmd.units import convert as unitconvert
import logging
logger = logging.getLogger(__name__)


class Model(object):
    """ Constructor for the Model class.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        A :class:`MetricData <htmd.metricdata.MetricData>` object containing the discretized trajectories

    Example
    -------
    >>> model = Model(data)

    .. rubric:: Methods
    .. autoautosummary:: htmd.model.Model
        :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.model.Model
        :attributes:
    """

    def __init__(self, data=None, file=None):
        if data is None:
            if file is not None:
                self.load(file)
            return
        else:
            self.data = data
        self.hmm = None
        self._modelid = None
        if self.data._clusterid is None:
            raise NameError('You need to cluster your data before making a Markov model')
        if self.data._dataid != self.data._clusterid:
            raise NameError('You have modified the data in data.dat after clustering. Please re-cluster.')
        self._clusterid = self.data._clusterid

    def markovModel(self, lag, macronum, units='frames', sparse=False, hmm=False):
        """ Build a Markov model at a given lag time and calculate metastable states

        Parameters
        ----------
        lag : int
            The lag time at which to calculate the Markov state model. The units are specified with the `units` argument.
        macronum : int
            The number of macrostates (metastable states) to produce
        units : str
            The units of lag. Can be 'frames' or any time unit given as a string.
        sparse : bool
            Make the transition matrix sparse. Useful if lots (> 4000) states are used for the MSM. Warning: untested.

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(150, 4)  # 150 frames lag, 4 macrostates
        """
        import pyemma.msm as msm
        self._integrityCheck(markov=True)

        lag = unitconvert(units, 'frames', lag, fstep=self.data.fstep)

        self.lag = lag
        self.msm = msm.estimate_markov_model(self.data.St.tolist(), self.lag, sparse=sparse)
        modelflag = False
        while not modelflag:
            self.coarsemsm = self.msm.pcca(macronum)
            if len(np.unique(self.msm.metastable_assignments)) != macronum:
                macronum -= 1
                logger.warning('PCCA returned empty macrostates. Reducing the number of macrostates to {}.'.format(macronum))
            else:
                modelflag = True
            if macronum < 2:
                raise RuntimeError('Could not create even two macrostates. Please revise your clustering.')

        self._modelid = random.random()

        if hmm:  # Still in development
            self.hmm = self.msm.coarse_grain(self.macronum)

        logger.info('{:.1f}% of the data was used'.format(self.msm.active_count_fraction * 100))

        _macroTrajectoriesReport(self.macronum, _macroTrajSt(self.data.St, self.macro_ofcluster), self.data.simlist)

    def createState(self, microstates=None, indexpairs=None):
        """ Creates a new state. Works both for new clusters and macrostates.

        If creating a new cluster, it just reassigns the given frames to the new cluster.
        If creating a new macrostate, it removes the given microstates from their previous macrostate, creates a new one
        and assigns them to it.

        Parameters
        ----------
        microstates : list
            The microstates to split out into a new macrostate.
        indexpairs : list
            List of lists. Each row is a simulation index-frame pair which should be added to a new cluster.
        """
        if microstates is not None and indexpairs is not None:
            raise AttributeError('microstates and indexpairs arguments are mutually exclusive')
        if microstates is not None:
            newmacro = self.macronum

            # Fixing sets. Remove microstates from previous macrostates and add new set
            for i, metset in enumerate(self.msm.metastable_sets):
                self.msm.metastable_sets[i] = np.array(np.sort(list(set(metset) - set(microstates))))
            self.msm.metastable_sets.append(np.array(microstates, dtype=np.int64))

            todelete = np.where([len(x) == 0 for x in self.msm.metastable_sets])[0]

            # Fixing hard assignments
            self.msm.metastable_assignments[microstates] = newmacro

            # Fixing memberships. Padding the array with 0s for the new macrostate
            self.msm._metastable_memberships = np.pad(self.msm.metastable_memberships, ((0, 0), (0, 1)), mode='constant', constant_values=(0))
            self.msm._metastable_memberships[microstates, :] = 0
            self.msm._metastable_memberships[microstates, -1] = 1

            # Moving probabilities of empty states to new one
            othermicro = np.ones(self.micronum, dtype=bool)
            othermicro[microstates] = False
            othermicro = np.where(othermicro)[0]
            self.msm._metastable_memberships[othermicro, -1] = np.sum(self.msm._metastable_memberships[othermicro[:, None], todelete], axis=1)

            # Fixing distributions
            self.msm._metastable_distributions = np.pad(self.msm.metastable_distributions, ((0, 1), (0, 0)), mode='constant', constant_values=(0))
            self.msm._metastable_distributions[-1, microstates] = 1 / len(microstates)
        if indexpairs is not None:
            newcluster = self.data.K
            for ip in indexpairs:
                self.data.trajectories[ip[0]].cluster[ip[1]] = newcluster
            self.data.K += 1
            self.data.N = np.bincount(np.concatenate(self.data.St))

    @property
    def P(self):
        """ The transition probability matrix """
        return self.msm.transition_matrix

    @property
    def micro_ofcluster(self):
        """ Mapping of clusters to microstates

        Numpy array which at index i has the index of the microstate corresponding to cluster i.
        Clusters which were not connected and thus are not in the model have a microstate value of -1.
        """
        self._integrityCheck(postmsm=True)
        micro_ofcluster = -np.ones(self.data.K, dtype=int)
        micro_ofcluster[self.msm.active_set] = np.arange(len(self.msm.active_set))
        return micro_ofcluster

    @property
    def cluster_ofmicro(self):
        """ Mapping of microstates to clusters

        Numpy array which at index i has the index of the cluster corresponding to microstate i.
        """
        self._integrityCheck(postmsm=True)
        return self.msm.active_set

    @property
    def micronum(self):
        """ Number of microstates """
        self._integrityCheck(postmsm=True)
        return len(self.msm.active_set)

    @property
    def macronum(self):
        """ Number of macrostates """
        self._integrityCheck(postmsm=True)
        return len(set(self.msm.metastable_assignments))

    @property
    def macro_ofmicro(self):
        """ Mapping of microstates to macrostates

        Numpy array which at index i has the index of the macrostate corresponding to microstate i.
        """
        self._integrityCheck(postmsm=True)
        # Fixing pyemma macrostate numbering
        mask = np.ones(np.max(self.msm.metastable_assignments) + 1, dtype=int) * -1
        mask[list(set(self.msm.metastable_assignments))] = range(self.macronum)
        return mask[self.msm.metastable_assignments]

    @property
    def macro_ofcluster(self):
        """ Mapping of clusters to microstates

        Numpy array which at index i has the index of the macrostate corresponding to cluster i.
        Clusters which were not connected and thus are not in the model have a macrostate value of -1.
        """
        self._integrityCheck(postmsm=True)
        macro_ofcluster = -np.ones(self.data.K, dtype=int)
        macro_ofcluster[self.msm.active_set] = self.macro_ofmicro
        return macro_ofcluster

    def plotTimescales(self, lags=None, minlag=None, maxlag=None, numlags=25, units='frames', errors=None, nits=None,
                       results=False, plot=True, save=None, njobs=-2):
        """ Plot the implied timescales of MSMs of various lag times

        Parameters
        ----------
        lags : list
            Specify explicitly at which lag times to compute the timescales.
        minlag: float
            The minimum lag time for the timescales. Used in combination with `maxlag` and `numlags`.
        maxlag: float
            The maximum lag time for the timescales. If None will default to the mode length of the trajectories.
        numlags: int
            The number of points to place between `minlag` and `maxlag`.
        units : str
            The units of lag. Can be 'frames' or any time unit given as a string.
        errors : errors
            Calculate errors using Bayes (Refer to pyEMMA documentation)
        nits : int
            Number of implied timescales to calculate. Default: all
        results : bool
            If the method should return the calculated implied timescales
        plot : bool
            If the method should display the plot of implied timescales
        save : str
            Path of the file in which to save the figure
        njobs : int
            Number of parallel jobs to spawn for calculation of timescales. Negative numbers are used for spawning jobs as many as CPU threads. 
            -1: for all CPUs -2: for all except one etc.

        Returns
        -------
        If given results=True this method will return the following data
        its : np.ndarray
            The calculated implied timescales. 2D array with dimensions (len(`lags`), `nits`)
        lags : np.ndarray
            A list of the lag times that were used to calculate the implied timescales

        Examples
        --------
        >>> model = Model(data)
        >>> model.plotTimescales()
        >>> model.plotTimescales(lags=list(range(1,100,5)))
        >>> model.plotTimescales(minlag=0.1, maxlag=20, numlags=25, units='ns')
        """
        import pyemma.plots as mplt
        import pyemma.msm as msm
        self._integrityCheck()
        if lags is None:
            lags = self.data._defaultLags(minlag, maxlag, numlags, units)
        else:
            lags = unitconvert(units, 'frames', lags, fstep=self.data.fstep).tolist()

        if nits is None:
            nits = np.min((self.data.K, 20))

        from htmd.config import _config
        its = msm.its(self.data.St.tolist(), lags=lags, errors=errors, nits=nits, n_jobs=njobs) # Use all CPUs minus one
        if plot or (save is not None):
            from matplotlib import pylab as plt
            plt.ion()
            plt.figure()
            try:
                mplt.plot_implied_timescales(its, dt=self.data.fstep, units='ns')
            except ValueError as ve:
                plt.close()
                raise ValueError('{} This is probably caused by badly set fstep in the data ({}). '.format(ve, self.data.fstep) +
                                 'Please correct the model.data.fstep to correspond to the simulation frame step in nanoseconds.')
            if save is not None:
                plt.savefig(save, dpi=300, bbox_inches='tight', pad_inches=0.2)
            if plot:
                plt.show()
        if results:
            return its.get_timescales(), its.lags

    def maxConnectedLag(self, lags):
        """ Heuristic for getting the lagtime before a timescale drops.

        It calculates the last lagtime before a drop occurs in the first implied timescale due to disconnected states.
        If the top timescale is closer to the second top timescale at the previous lagtime than to itself at the previous
        lagtime it means that a drop occured. The lagtime before the drop is returned.

        Parameters
        ----------
        lags : np.ndarray or list
            A list of lag times for which to calculate the implied timescales

        Returns
        -------
        ml : int
            The maximum lagtime before a drop occurs in the top timescale

        Examples
        --------
        >>> model = Model(data)
        >>> model.maxConnectedLag(list(range(1, 100, 5)))
        """
        if len(lags) == 1:
            return lags
        if isinstance(lags, np.ndarray):
            lags = lags.astype(int)

        import pyemma.msm as msm
        itime = msm.its(self.data.St.tolist(), lags=lags, nits=2).get_timescales()

        for i in range(1, np.size(itime, 0)):
            if abs(itime[i, 0] - itime[i-1, 1]) < abs(itime[i, 0] - itime[i-1, 0]):
                lagidx = i-1
                break
            else:
                lagidx = i
        return lags[lagidx], itime

    def sampleStates(self, states=None, frames=20, statetype='macro', replacement=False, samplemode='random', allframes=False):
        """ Samples frames from a set of states

        Parameters
        ----------
        states : Union[list, int]
            A list of state indexes from which we want to sample
        frames : Union[None, int, list]
            An integer with the number of frames we want to sample per state or a list of same length as
            `states` which contains the number of frames we want from each of the states.
            If set to None it will return all frames of the states.
        statetype : ['micro','macro','cluster'], optional
            The type of state we want to sample from.
        replacement : bool
            If we want to sample with or without replacement.
        samplemode : ['random','even','weighted']
            What sampling strategy to use. For `statetype` == 'macro' this can be set to 'even' to sample evenly from
            all microstates in the macrostate or to 'weighted' to sample proportional to the equilibium probability of
            each microstate in the macrostate.
        allframes : bool
            Deprecated. Use frames=None instead.

        Returns
        -------
        absframes : numpy.ndarray
            An array which contains for each state an array containing absolute trajectory frames
        relframes : numpy.ndarray
            An array which contains for each state a 2D array containing the trajectory ID and frame number for each of
            the sampled frames

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.sampleStates(range(5), [10, 3, 2, 50, 1])  # Sample from all 5 macrostates
        >>> model.sampleStates(range(model.micronum), samplesnum, statetype='micro')  # Sample from all microstates
        """
        if statetype == 'cluster':
            if samplemode != 'random':
                logger.warning("'cluster' states incompatible with 'samplemode' other than 'random'. Defaulting to 'random'")
            return self.data.sampleClusters(states, frames, replacement, allframes)

        if states is None:
            if statetype == 'macro':
                states = range(self.macronum)
            elif statetype == 'micro':
                states = range(self.micronum)
        if isinstance(states, int):
            states = [states, ]

        if allframes:
            logger.warning('The allframes option will be deprecated. Use frames=None instead.')
            frames = None

        self._integrityCheck(postmsm=True)
        if statetype != 'macro' and samplemode != 'random':
            samplemode = 'random'
            logger.warning("'micro' states incompatible with 'samplemode' other than 'random'. Defaulting to 'random'")

        if frames is None or isinstance(frames, int):
            frames = np.repeat(frames, len(states))

        stConcat = np.concatenate(self.data.St)
        absFrames = []
        relFrames = []
        for i in range(len(states)):
            if frames[i] == 0:
                absFrames.append(np.array([], dtype=int))
                relFrames.append(np.array([], dtype=int))
                continue

            st = states[i]
            if statetype == 'macro':
                (selFr, selMicro) = _sampleMacro(self, st, stConcat, samplemode, frames[i], replacement)
                absFrames.append(selFr)
            elif statetype == 'micro':
                absFrames.append(_sampleMicro(self, st, stConcat, frames[i], replacement))
            else:
                raise NameError('No valid state type given (read documentation)')

            if len(absFrames[-1]) == 0:
                raise NameError('No frames could be sampled from {} state {}. State is empty.'.format(statetype, st))

            relFrames.append(self.data.abs2rel(absFrames[-1]))
        return absFrames, relFrames

    def eqDistribution(self, plot=True, save=None):
        """ Obtain and plot the equilibrium probabilities of each macrostate

        Parameters
        ----------
        plot : bool, optional, default=True
            Disable plotting of the equilibrium distribution by setting it to False
        save : str
            Path of the file in which to save the figure

        Returns
        -------
        eq : ndarray
            An array of equilibrium probabilities of the macrostates

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.eqDistribution()
        """
        # logger.warning('Equilibrium distribution calculations for macrostates are now done using membership '
        #                'probabilities and hence your results might differ from analyses done before this change.')
        self._integrityCheck(postmsm=True)
        macroeq = np.ones(self.macronum) * -1
        macroindexes = list(set(self.msm.metastable_assignments))
        for i in range(self.macronum):
            # macroeq[i] = np.sum(self.msm.stationary_distribution[self.macro_ofmicro == i])
            macroeq[i] = np.sum(self.msm.metastable_memberships[:, macroindexes[i]] * self.msm.stationary_distribution)

        if plot or (save is not None):
            from matplotlib import pylab as plt
            plt.ion()
            plt.figure()
            plt.bar(np.arange(self.macronum)+0.4, macroeq)
            plt.ylabel('Equilibrium probability')
            plt.xlabel('Macrostates')
            plt.xticks(np.arange(self.macronum)+0.4, range(self.macronum))
            if save is not None:
                plt.savefig(save, dpi=300, bbox_inches='tight', pad_inches=0.2)
            if plot:
                plt.show()

        return macroeq

    def _coarseP(self):
        M = self.msm.metastable_memberships
        Pcoarse = np.linalg.inv(M.T.dot(M)).dot(M.T).dot(self.P).dot(M)
        if len(np.where(Pcoarse < 0)[0]) != 0:
            raise NameError('Cannot produce coarse P matrix. Ended up with negative probabilities. Try using less macrostates.')
        return Pcoarse

    def getStates(self, states=None, statetype='macro', wrapsel='protein', alignsel='name CA', alignmol=None, samplemode='weighted', numsamples=50, simlist=None):
        """ Get samples of MSM states in Molecule classes

        Parameters
        ----------
        states : ndarray, optional
            A list of states to visualize
        statetype : ['macro','micro','cluster'], optional
            The type of state to visualize
        wrapsel : str, optional, default='protein'
            Atom selection string to use for wrapping.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignsel : str, optional, default='name CA'
            Atom selection string used for aligning all frames. Set to None to disable aligning.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignmol : :class:`Molecule <moleculekit.molecule.Molecule>` object
            A reference molecule onto which to align all others
        samplemode : ['weighted','random'], optional, default='weighted'
            How to obtain the samples from the states
        numsamples : int
            Number of samples (conformations) for each state.
        simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
            Optionally pass a different (but matching, i.e. filtered) simlist for creating the Molecules.

        Returns
        -------
        mols : ndarray of :class:`Molecule <moleculekit.molecule.Molecule>` objects
            A list of :class:`Molecule <moleculekit.molecule.Molecule>` objects containing the samples of each state

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> mols = model.getStates()
        >>> for m in mols:
        >>>     m.view()
        """
        from htmd.projections.metric import _singleMolfile
        self._integrityCheck(postmsm=(statetype != 'cluster'))
        if simlist is None:
            simlist = self.data.simlist
        else:
            if len(simlist) != len(self.data.simlist):
                raise AttributeError('Provided simlist has different number of trajectories than the one used by the model.')

        (single, molfile) = _singleMolfile(simlist)
        if not single:
            raise NameError('Visualizer does not support yet visualization of systems with different structure files. '
                            'The simlist should be created with a single molfile (for example a filtered one)')
        if alignmol is None:
            alignmol = Molecule(molfile)
        if statetype != 'macro' and statetype != 'micro' and statetype != 'cluster':
            raise NameError("'statetype' must be either 'macro', 'micro' or ''cluster'")
        if states is None:
            if statetype == 'macro':
                states = range(self.macronum)
            elif statetype == 'micro':
                states = range(self.micronum)
            elif statetype == 'cluster':
                states = range(self.data.K)
        if len(states) == 0:
            raise NameError('No ' + statetype + ' states exist in the model')

        (tmp, relframes) = self.sampleStates(states, numsamples, statetype=statetype, samplemode=samplemode)

        from htmd.config import _config
        from htmd.parallelprogress import ParallelExecutor, delayed
        # This loop really iterates over states. sampleStates returns an array of arrays
        # Don't increase njobs because it was giving errors on some systems.
        aprun = ParallelExecutor(n_jobs=1)
        mols = aprun(total=len(relframes), desc='Getting state Molecules')\
            (delayed(_loadMols)(self, rel, molfile, wrapsel, alignsel, alignmol, simlist) for rel in relframes)
        return np.array(mols, dtype=object)

    def viewStates(self, states=None, statetype='macro', protein=None, ligand=None, viewer=None, mols=None,
                   numsamples=50, wrapsel='protein', alignsel='name CA', gui=False, simlist=None):
        """ Visualize macro/micro/cluster states in VMD

        Parameters
        ----------
        states : ndarray, optional
            A list of states to visualize
        statetype : ['macro','micro','cluster'], optional
            The type of state to visualize
        protein : bool, optional
            Set to True to enable pure protein system visualization
        ligand : str, optional
            Atom selection string for the ligand.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        viewer : :class:`VMD <moleculekit.vmdviewer.VMD>` object, optional
            A viewer in which to visualize the states
        mols : ndarray, optional
            An array of :class:`Molecule <moleculekit.molecule.Molecule>` objects to visualize
        numsamples : int
            Number of samples (conformations) for each state.
        wrapsel : str, optional, default='protein'
            Atom selection string to use for wrapping.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        alignsel : str, optional, default='name CA'
            Atom selection string used for aligning all frames. See to None to disable aligning.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
            Optionally pass a different (but matching, i.e. filtered) simlist for visualizing the states.

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.viewStates(protein=True)

        >>> model.viewStates(ligand='resname MOL')
        """
        from htmd.config import _config
        self._integrityCheck(postmsm=(statetype != 'cluster'))

        if protein is None and ligand is None:
            raise RuntimeError('You need to specify either the `protein` of `ligand` arguments for state visualization.')

        if _config['viewer'].lower() == 'ngl' or _config['viewer'].lower() == 'webgl':
            return self._viewStatesNGL(states, statetype, protein, ligand, mols, numsamples, gui=gui)

        if viewer is None:
            viewer = getCurrentViewer()
        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(states, statetype, numsamples=numsamples, wrapsel=wrapsel, alignsel=alignsel, simlist=simlist)
        colors = [0, 1, 3, 4, 5, 6, 7, 9]
        for i, s in enumerate(states):
            viewer.loadMol(mols[i], name=statetype+' '+str(states[i]))
            if ligand is not None:
                viewer.rep('ligand', sel=ligand, color=colors[np.mod(i, len(colors))])
            if protein is not None:
                viewer.rep('protein')
            viewer.send('start_sscache')

    def _viewStatesNGL(self, states, statetype, protein, ligand, mols, numsamples, gui=False):
        from moleculekit.util import sequenceID
        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(states, statetype, numsamples=min(numsamples, 15))
        colors = [0, 1, 3, 4, 5, 6, 7, 9]
        hexcolors = {0: '#0000ff', 1: '#ff0000', 2: '#333333', 3: '#ff6600', 4: '#ffff00', 5: '#4c4d00', 6: '#b2b2cc',
                     7: '#33cc33', 8: '#ffffff', 9: '#ff3399', 10: '#33ccff'}
        if protein is None and ligand is None:
            raise NameError('Please provide either the "protein" or "ligand" parameter for viewStates.')
        k = 0
        from nglview import NGLWidget, HTMDTrajectory
        view = NGLWidget(gui=gui)
        ref = mols[0].copy()
        for i, s in enumerate(states):
            if protein:
                mol = Molecule()
            if ligand:
                mol = ref.copy()
                mol.remove(ligand, _logger=False)
                mol.dropFrames(keep=0)
                mols[i].filter(ligand, _logger=False)
            mols[i].set('chain', '{}'.format(s))
            tmpcoo = mols[i].coords
            for j in range(mols[i].numFrames):
                mols[i].coords = np.atleast_3d(tmpcoo[:, :, j])
                if ligand:
                    mols[i].set('segid', sequenceID(mols[i].resid)+k)
                    k = int(mols[i].segid[-1])
                mol.append(mols[i])
            view.add_trajectory(HTMDTrajectory(mol))
            # Setting up representations
            if ligand:
                view[i].add_cartoon('protein', color='sstruc')
                view[i].add_hyperball(':{}'.format(s), color=hexcolors[np.mod(i, len(hexcolors))])
            if protein:
                view[i].add_cartoon('protein', color='residueindex')

        self._nglButtons(view, statetype, states)
        return view

    def _nglButtons(self, ngl_widget, statetype, states):
        # Adds buttons for enabling and disabling macrostate visualizations
        import ipywidgets
        from IPython.display import display

        container = []
        for s in states:
            w = ipywidgets.Checkbox(description="{} {}".format(statetype, s))
            w.value = True
            container.append(w)

        def updateReps(name):
            on = []
            for i, w in enumerate(container):
                if w.value:
                    on.append(i)
            ngl_widget.show_only(on)

        for w in container:
            w.on_trait_change(updateReps, "value")

        #container.append(ipywidgets.Checkbox(description="all"))

        hb = ipywidgets.HBox(container)
        display(hb)

    def save(self, filename):
        """ Save a :class:`Model <htmd.model.Model>` object to disk

        Parameters
        ----------
        filename : str
            Path of the file in which to save the object

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.save('./model.dat')
        """
        import pickle

        # Temporarily store data object and replace with dicts
        tmpdata = self.data
        if self.data.parent is not None:
            tmpparentdata = self.data.parent
            self.data.parent = self.data.parent.__dict__
        self.data = self.data.__dict__

        # Dump the dict
        f = open(filename, 'wb')
        pickle.dump(self.__dict__, f)
        f.close()

        # Restore data to classes
        self.data = tmpdata
        if self.data.parent is not None:
            self.data.parent = tmpparentdata
        

    def load(self, filename):
        """ Load a :class:`MetricData <htmd.metricdata.MetricData>` object from disk

        Parameters
        ----------
        filename : str
            Path to the saved MetricData object

        Examples
        --------
        >>> model = Model()
        >>> model.load('./model.dat')
        """
        import sys
        import pickle
        from htmd.metricdata import MetricData
        try:
            import pandas.indexes
        except ImportError:
            import pandas.core.indexes
            sys.modules['pandas.indexes'] = pandas.core.indexes  # Hacky fix for new pandas version

        f = open(filename, 'rb')
        z = pickle.load(f)
        f.close()
        for k in z:
            if k == 'data':
                m = MetricData()
                m.load(z[k])
                self.__dict__[k] = m
            else:
                self.__dict__[k] = z[k]

    def copy(self):
        """ Produces a deep copy of the object

        Returns
        -------
        data : :class:`Model` object
            A copy of the current object

        Examples
        --------
        >>> model2 = model.copy()
        """
        from copy import deepcopy
        return deepcopy(self)

    def cktest(self):
        """ Conducts a Chapman-Kolmogorow test.

        Returns
        -------

        """
        from copy import deepcopy
        from pyemma.plots import plot_cktest
        msm = deepcopy(self.msm)
        ck = msm.cktest(self.macronum)
        plot_cktest(ck)

    def createCoreSetModel(self, threshold=0.5):
        """ Given an MSM this function detects the states belonging to a core set and returns a new model consisting
        only of these states.

        Parameters
        ----------
        threshold : float
            Membership probability threshold over which microstates belong to the core of a macrostate

        Returns
        -------
        newmodel :
            A new model object
        frames : list
            A list of the frames that were kept in the new model
        """
        if (threshold >= 1) or (threshold <= 0):
            raise AttributeError('threshold argument only accepts values between (0, 1)')

        def calcCoreSet(distr, assign, threshold):
            coreset = []
            for i, md in enumerate(distr):
                microofmacro = np.where(assign == i)[0]
                prob = md[microofmacro]
                tt = threshold * (prob.max() - prob.min())
                coreset.append(microofmacro[np.where(prob > tt)[0]])
            return coreset

        def coreDtraj(data, micro_ofcluster, coreset):
            corestates = np.concatenate(coreset)
            newmapping = np.ones(corestates.max()+1, dtype=int) * -1
            newmapping[corestates] = np.arange(len(corestates))
            discretetraj = [st.copy() for st in data.St]
            frames = []
            newdiscretetraj = []
            newcounts = np.zeros(len(corestates))
            for t, st in enumerate(discretetraj):
                oldmicro = None
                newtraj = []
                tframes = []
                for f, cl in enumerate(st):
                    newmicro = None
                    micro = micro_ofcluster[cl]
                    if micro == -1:  # If we are in an dropped cluster, keep old index
                        newmicro = oldmicro
                    else:
                        for co in coreset:
                            if micro in co:
                                newmicro = micro
                                oldmicro = micro
                                break
                    if newmicro is None and oldmicro is not None:
                        newtraj.append(oldmicro)
                        tframes.append(f)
                    elif newmicro is not None:
                        newtraj.append(newmicro)
                        tframes.append(f)
                mappedtraj = newmapping[np.array(newtraj, dtype=int)]
                newdiscretetraj.append(mappedtraj)
                if len(mappedtraj):
                    newcounts[:mappedtraj.max()+1] += np.bincount(mappedtraj)
                frames.append(np.array(tframes))
            #kept = np.array([i for i, x in enumerate(newdiscretetraj) if len(x) != 0])
            return np.array(newdiscretetraj, dtype=object), len(corestates), newcounts, frames

        coreset = calcCoreSet(self.msm.metastable_distributions, self.msm.metastable_assignments, threshold)
        newdata = self.data.copy()
        newSt, newdata.K, newdata.N, frames = coreDtraj(self.data, self.micro_ofcluster, coreset)
        for i, (s, fr) in enumerate(zip(newSt, frames)):
            if len(s):
                newdata.trajectories[i].cluster[fr] = s

        logger.info('Kept {} microstates from each macrostate.'.format([len(x) for x in coreset]))

        dataobjects = [newdata]
        if newdata.parent is not None:
            dataobjects.append(newdata.parent)
        for data in dataobjects:
            dat = []
            ref = []
            simstmp = []
            cluster = []
            for i, fr in enumerate(frames):
                if len(fr):
                    dat.append(data.trajectories[i].projection[fr, :])
                    ref.append(data.trajectories[i].reference[fr, :])
                    if len(data.trajectories[i].cluster):
                        cluster.append(data.trajectories[i].cluster[fr])
                    simstmp.append(data.trajectories[i].sim)
                    
            data._loadTrajectories(dat, ref, simstmp, cluster if len(cluster) else None)

        return Model(newdata), frames

    def plotFES(self, dimX, dimY, temperature, states=False, s=10, cmap=None, fescmap=None, statescmap=None, plot=True, save=None, data=None):
        """ Plots the free energy surface on any given two dimensions. Can also plot positions of states on top.

        Parameters
        ----------
        dimX : int
            Index of projected dimension to use for the X axis.
        dimY : int
            Index of projected dimension to use for the Y axis.
        temperature : float
            Simulation temperature.
        states : bool
            If True, will plot scatter plot of microstates coloured by macro state on top of FES.
        s : float
            Marker size for states.
        cmap :
            Sets the Matplotlib colormap for both `fescmap` and `statescmap`
        fescmap :
            Matplotlib colormap for the free energy surface
        statescmap:
            Matplotlib colormap for the states
        plot : bool
            If the method should display the plot of the FES
        save : str
            Path of the file in which to save the figure
        data : :class:`MetricData` object
            Optionally you can pass a different MetricData object than the one used to build the model. For example
            if the user wants to build a model on distances but wants to plot the FES on top of RMSD values. The
            MetricData object needs to have the same simlist as the Model.

        Examples
        --------
        >>> import matplotlib as plt
        >>> model.plotFES(0, 1, 300)
        >>> model.plotFES(2, 3, 300, data=otherdata, states=True, fescmap=plt.cm.gray)
        """
        self._integrityCheck(postmsm=True)
        from matplotlib import pylab as plt
        from htmd.kinetics import Kinetics

        if data is None:
            data = self.data
            microcenters = self.data.Centers[self.cluster_ofmicro, :]
        else:
            if self.data.numFrames != data.numFrames or ~np.all([s1 == s2 for s1, s2 in zip(self.data.simlist, data.simlist)]):
                raise RuntimeError('The data argument you provided uses a different simlist than the Model.')
            microcenters = np.vstack(getStateStatistic(self, data, range(self.micronum), statetype='micro'))

        if fescmap is None:
            fescmap = plt.cm.jet
        if statescmap is None:
            statescmap = plt.cm.jet
        if cmap is not None:
            fescmap = cmap
            statescmap = cmap

        if data.description is not None:
            xlabel = data.description.description[dimX]
        else:
            xlabel = 'Dimension {}'.format(dimX)
        if data.description is not None:
            ylabel = data.description.description[dimY]
        else:
            ylabel = 'Dimension {}'.format(dimY)
        title = 'Free energy surface'

        energy = -Kinetics._kB * temperature * np.log(self.msm.stationary_distribution)
        f, ax, cf = data._contourPlot(microcenters[:, dimX], microcenters[:, dimY], energy, cmap=fescmap, xlabel=xlabel, ylabel=ylabel, title=title)
        data._setColorbar(f, cf, 'kcal/mol', scientific=False)
        if states:
            colors = statescmap(np.linspace(0, 1, self.macronum))
            for m in range(self.macronum):
                macromicro = microcenters[self.macro_ofmicro == m, :]
                _ = ax.scatter(macromicro[:, dimX], macromicro[:, dimY], s=s, c=colors[m], label='Macro {}'.format(m), edgecolors='none')
            ax.legend(prop={'size': 8})

        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches='tight', pad_inches=0.2)
        if plot:
            plt.show()

    def _integrityCheck(self, postmsm=False, markov=False):
        if postmsm and self._modelid is None:
            raise NameError('You need to call markovModel before calling this command')
        if not markov and self.data._dataid == self.data._clusterid and self.data._dataid != self._clusterid:
            raise NameError('After updating the MetricData object you need to call the markovModel command anew.')
        if self.data._dataid != self.data._clusterid:
            raise NameError('After modifying the data in the MetricData object you need to recluster and reconstruct the markov model.')


def _loadMols(self, rel, molfile, wrapsel, alignsel, refmol, simlist):
    frames = self.data.rel2sim(rel, simlist=simlist)
    mol = Molecule(molfile)
    trajs = []
    frs = []
    for f in frames:
        trajs.append(f.sim.trajectory[f.piece])
        frs.append(f.frame)
    mol.read(np.array(trajs), frames=np.array(frs))
    if len(wrapsel) > 0:
        mol.wrap(wrapsel)
    if (refmol is not None) and (alignsel is not None):
        mol.align(alignsel, refmol=refmol)
    return mol


def getStateStatistic(reference, data, states, statetype='macro', weighted=False, method=np.mean, axis=0, existing=False, target=None):
    """ Calculates properties of the states.

    Calculates properties of data corresponding to states. Can calculate for example the mean distances of atoms in a
    state, or the standard deviation of the RMSD in a state.

    Parameters
    ----------
    reference : :class:`Model` object or :class:`MetricData` object
        A model containing the state definitions or a MetricData object containing cluster definitions
    data : :class:`MetricData` object
        A projection corresponding to the conformations in the states of the model. The data and the model need to share
        the same `simlist`
    states : list
        A list of the states for which we want to calculate the properties
    statetype : ['macro','micro','cluster']
        The state type
    weighted : bool
        If the properties of the macrostates should be calculated as the weighted average of their microstates, where the
        weights are the equilibrium probabilities of the microstates
    method :
        A function pointer for the method that calculates our desired property. e.g. np.mean, np.std
    axis : int
        On which axis of the data the method should operate on
    existing : bool
        Not used
    target :
        Not used

    Returns
    -------
    statistic : list
        A list which contains in each element the desired statistic of a specific state specifies in `states`

    Examples
    --------
    >>> data = MetricDistance.project(sims, 'protein and name CA', 'resname MOL')
    >>> model = Model(data)
    >>> model.markovModel(100, 5)
    >>> # Get the standard deviation of distances in all macrostates
    >>> getStateStatistic(model, data, list(range(5)), method=np.std)
    """
    from htmd.metricdata import MetricData
    if axis != 0:
        logger.warning('Axis different than 0 might not work correctly yet')

    if isinstance(reference, Model):
        refdata = reference.data
    elif isinstance(reference, MetricData):
        refdata = reference
        if statetype != 'cluster':
            raise RuntimeError('You can only use statetype cluster with reference MetricData object. To use other statetypes build a Model.')
    else:
        raise RuntimeError('Invalid argument type')

    if refdata.numTrajectories > 0 and np.any(refdata.trajLengths != data.trajLengths):
        raise NameError('Data trajectories need to match in size and number to the trajectories in the model')
    stconcat = np.concatenate(refdata.St)
    datconcat = np.concatenate(data.dat)

    statistic = []
    for i, st in enumerate(states):
        if statetype == 'macro':
            frames = reference.macro_ofcluster[stconcat] == st
        elif statetype == 'micro':
            frames = reference.micro_ofcluster[stconcat] == st
        elif statetype == 'cluster':
            frames = stconcat == st
        else:
            raise NameError('No valid state type given (read documentation)')

        if statetype == 'macro' and weighted:
            statistic.append(_weightedMethod(reference, method, stconcat, datconcat, st, axis))
        else:
            if axis is None:
                statistic.append(method(datconcat[frames, ...]))
            else:
                statistic.append(method(datconcat[frames, ...], axis=axis))
    return statistic


def _weightedMethod(model, method, stconcat, datconcat, st, axis):
    microsofmacro = np.where(model.macro_ofmicro == st)[0]
    eq = model.msm.stationary_distribution
    weights = eq / np.sum(eq[microsofmacro])
    avgstatistic = np.zeros(np.size(datconcat, 1))
    for m in microsofmacro:
        frames = model.micro_ofcluster[stconcat] == m
        if axis is None:
            stat = method(datconcat[frames, :])
        else:
            stat = method(datconcat[frames, :], axis=axis)
        avgstatistic = avgstatistic + stat * weights[m]
    return avgstatistic


def macroAccumulate(model, microvalue):
    """ Accumulate values of macrostates from a microstate array

    Parameters
    ----------
    model : :class:`Model <htmd.model.Model>` object
        The model which to use to accumulate
    microvalue : ndarray
        An array of values corresponding to the microstates of the model

    Returns
    -------
    macrovalue : ndarray
        An array of values corresponding to the macrostates of the model
    """
    res = np.zeros(model.macronum)
    for i in range(len(microvalue)):
        macro = model.macro_ofmicro[i]
        res[macro] = res[macro] + microvalue[i]
    return res


def _sampleMacro(obj, macro, stConcat, mode, numFrames, replacement):
    from htmd.metricdata import _randomSample
    if mode == 'random':
        frames = np.where(obj.macro_ofcluster[stConcat] == macro)[0]
        selFrames = _randomSample(frames, numFrames, replacement)
        selMicro = obj.micro_ofcluster[stConcat[selFrames]]
    elif mode == 'even':
        micros = np.where(obj.macro_ofmicro == macro)[0]
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames.append(_sampleMicro(obj, micros[i], stConcat, numFrames, replacement))
            selMicro.append(np.ones(numFrames) * micros[i])
    elif mode == 'weighted':
        micros = np.where(obj.macro_ofmicro == macro)[0]
        eq = obj.msm.stationary_distribution
        weights = eq[micros] / np.sum(eq[micros])
        framespermicro = np.random.multinomial(numFrames, weights)
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames = np.append(selFrames, _sampleMicro(obj, micros[i], stConcat, framespermicro[i], replacement))
            selMicro = np.append(selMicro, [micros[i]] * framespermicro[i])
    elif mode == 'weightedTrunc':
        micros = np.where(obj.macro_ofmicro == macro)[0]
        eq = obj.msm.stationary_distribution
        weights = eq[micros] / np.sum(eq[micros])
        idx = np.argsort(weights)
        cs = np.cumsum(weights[idx])
        under50 = cs < 0.5
        weights[idx[under50]] = 0
        weights /= np.sum(weights)
        framespermicro = np.random.multinomial(numFrames, weights)
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames = np.append(selFrames, _sampleMicro(obj, micros[i], stConcat, framespermicro[i], replacement))
            selMicro = np.append(selMicro, [micros[i]] * framespermicro[i])
    else:
        raise NameError('No valid mode given (read documentation)')
    return selFrames, selMicro


def _sampleMicro(obj, micro, stConcat, numFrames, replacement):
    from htmd.metricdata import _randomSample
    frames = np.where(obj.micro_ofcluster[stConcat] == micro)[0]
    return _randomSample(frames, numFrames, replacement)


def _macroTrajectoriesReport(macronum, macrost, simlist=None):
    """
    Prints out the number of trajectories that visited each macrostate and prints warnings

    Parameters
    ----------
    macronum
    macrost
    simlist
    """
    macrotraj = dict()
    macrotrajnum = np.zeros(macronum, dtype=int)
    for i, st in enumerate(macrost):
        macros = np.setdiff1d(np.unique(st), [-1])
        for m in macros:
            if m not in macrotraj:
                macrotraj[m] = []
            macrotraj[m].append(i)
            macrotrajnum[m] += 1
    logger.info('Number of trajectories that visited each macrostate:')
    logger.info(macrotrajnum)

    for m in range(macronum):
        ratio = macrotrajnum[m] / len(macrost)
        if macrotrajnum[m] <= 3 and ratio <= 0.2:
            logger.info('Take care! Macro {} has been visited only in {} trajectories'
                        ' ({:.1f}% of total):'.format(m, macrotrajnum[m], ratio*100))
            if simlist is not None:
                for s in macrotraj[m]:
                    logger.info(simlist[s])


def _macroTrajSt(St, macro_ofcluster):
    mst = np.empty(np.shape(St), dtype=object)
    for i in range(len(St)):
        mst[i] = macro_ofcluster[St[i]]
    return mst


'''def _macroP(C, macro_ofmicro):
    macronum = np.max(macro_ofmicro) + 1
    macroC = np.zeros((macronum, macronum))

    for m1 in range(macronum):
        sourcemicros = np.where(macro_ofmicro == m1)[0]
        for m2 in range(macronum):
            sinkmicros = np.where(macro_ofmicro == m2)[0]
            ixgrid = np.ix_(sourcemicros, sinkmicros)
            macroC[m1, m2] = np.sum(C[ixgrid].flatten())

    from msmtools.estimation import transition_matrix
    return transition_matrix(macroC, reversible=True)'''

import unittest
class _TestModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.simlist import simlist, simfilter
        from glob import glob
        from htmd.projections.metric import Metric
        from moleculekit.projections.metricdistance import MetricDistance
        from moleculekit.projections.metricdihedral import MetricDihedral
        from moleculekit.util import tempname
        from sklearn.cluster import MiniBatchKMeans
        from htmd.home import home
        from os.path import join

        sims = simlist(glob(join(home(dataDir='adaptive'), 'data', '*', '')), glob(join(home(dataDir='adaptive'), 'input', '*')))
        fsims = simfilter(sims, tempname(), 'not water')

        metr = Metric(fsims)
        metr.set(MetricDistance('protein and resid 10 and name CA', 'resname BEN and noh', periodic='selections', metric='contacts', groupsel1='residue', threshold=4))
        data = metr.project()
        data.cluster(MiniBatchKMeans(n_clusters=4))

        self.model = Model(data)

    def test_model_saving_loading(self):
        from moleculekit.util import tempname

        modelfile = tempname(suffix='.dat')
        self.model.save(modelfile)

        newmodel = Model(file=modelfile)
        assert newmodel.data.numTrajectories == 2

        # Testing model saving when data has parent
        self.model.data.parent = self.model.data.copy()
        self.model.save(modelfile)
        newmodel = Model(file=modelfile)
        assert newmodel.data.numTrajectories == 2
        assert newmodel.data.parent.numTrajectories == 2
        self.model.data.parent = None

if __name__ == '__main__':
    unittest.main(verbosity=2)



