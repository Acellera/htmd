# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy import stats
import warnings
import random
from htmd.projections.metric import _singleMolfile
from htmd.molecule.molecule import Molecule
from htmd.vmdviewer import getCurrentViewer
from htmd.units import convert as unitconvert
import logging
logger = logging.getLogger(__name__)


class Model(object):
    """ Constructor for the Model class. Model uses PyEMMA [1] internally to calculate Markov models.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        A :class:`MetricData <htmd.metricdata.MetricData>` object containing the discretized trajectories

    Example
    -------
    >>> model = Model(data)

    References
    ----------
    .. [1] PyEMMA 2: A Software Package for Estimation, Validation, and Analysis of Markov Models. Martin K. Scherer et al. JCTC 2015.
    """

    def __init__(self, data=None):
        if data is None:
            return
        self.data = data
        self.hmm = None
        self._modelid = None
        if data._clusterid is None:
            raise NameError('You need to cluster your data before making a Markov model')
        if data._dataid != data._clusterid:
            raise NameError('You have modified the data in data.dat after clustering. Please re-cluster.')
        self._clusterid = data._clusterid

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
        self.coarsemsm = self.msm.pcca(macronum)

        self._modelid = random.random()

        if hmm:  # Still in development
            self.hmm = self.msm.coarse_grain(self.macronum)

        logger.info('{:.1f}% of the data was used'.format(self.msm.active_count_fraction * 100))

        _macroTrajectoriesReport(self.macronum, _macroTrajSt(self.data.St, self.macro_ofcluster), self.data.simlist)

    @property
    def P(self):
        return self.msm.transition_matrix

    @property
    def micro_ofcluster(self):
        self._integrityCheck(postmsm=True)
        micro_ofcluster = -np.ones(self.data.K+1, dtype=int)
        micro_ofcluster[self.msm.active_set] = np.arange(len(self.msm.active_set))
        return micro_ofcluster

    @property
    def cluster_ofmicro(self):
        self._integrityCheck(postmsm=True)
        return self.msm.active_set

    @property
    def micronum(self):
        self._integrityCheck(postmsm=True)
        return len(self.msm.active_set)

    @property
    def macronum(self):
        self._integrityCheck(postmsm=True)
        return len(set(self.msm.metastable_assignments))

    @property
    def macro_ofmicro(self):
        self._integrityCheck(postmsm=True)
        # Fixing pyemma macrostate numbering
        mask = np.ones(np.max(self.msm.metastable_assignments) + 1, dtype=int) * -1
        mask[list(set(self.msm.metastable_assignments))] = range(self.macronum)
        return mask[self.msm.metastable_assignments]

    @property
    def macro_ofcluster(self):
        self._integrityCheck(postmsm=True)
        macro_ofcluster = -np.ones(self.data.K+1, dtype=int)
        macro_ofcluster[self.msm.active_set] = self.macro_ofmicro
        return macro_ofcluster

    def plotTimescales(self, lags=None, units='frames', errors=None, nits=None, results=False, plot=True):
        """ Plot the implied timescales of MSMs of various lag times

        Parameters
        ----------
        lags : list
            The lag times at which to compute the timescales. By default it spreads out 25 lag times linearly from lag
            10 until the mode length of the trajectories.
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

        Returns
        -------
        If given `results`=True this method will return the following data
        its : np.ndarray
            The calculated implied timescales. 2D array with dimensions (len(`lags`), `nits`)
        lags : np.ndarray
            A list of the lag times that were used to calculate the implied timescales

        Examples
        --------
        >>> model = Model(data)
        >>> model.plotTimescales()
        >>> model.plotTimescales(lags=list(range(1,100,5)))
        """
        import pyemma.plots as mplt
        import pyemma.msm as msm
        self._integrityCheck()
        if lags is None:
            lags = self._defaultLags()
        else:
            lags = unitconvert(units, 'frames', lags, fstep=self.data.fstep).tolist()

        if nits is None:
            nits = np.min((self.data.K, 20))

        from htmd.config import _config
        its = msm.its(self.data.St.tolist(), lags=lags, errors=errors, nits=nits, n_jobs=_config['ncpus'])
        if plot:
            from matplotlib import pylab as plt
            plt.ion()
            plt.figure()
            mplt.plot_implied_timescales(its, dt=self.data.fstep, units='ns')
            plt.show()
        if results:
            return its.get_timescales(), its.lags

    def maxConnectedLag(self, lags):
        """ Calculates the last lagtime before a drop occurs in the first implied timescale

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

    def _defaultLags(self):
        return self.data._defaultLags()

    def sampleStates(self, states, frames, statetype='micro', replacement=False, samplemode='random', allframes=False):
        """ Samples frames from a set of states

        Parameters
        ----------
        states : list
            A list of state indexes from which we want to sample
        frames : list
            A list of same length as `states` which contains the number of frames we want from each of the states
        statetype : ['micro','macro','cluster'], optional
            The type of state we want to sample from.
        replacement : bool
            If we want to sample with or without replacement.
        samplemode : ['random','even','weighted']
            What sampling strategy to use. For `statetype` == 'macro' this can be set to 'even' to sample evenly from
            all microstates in the macrostate or to 'weighted' to sample proportional to the equilibium probability of
            each microstate in the macrostate.
        allframes : bool
            If we want to simply retrieve all frames of the states. Ignores the `frames` argument.

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
        self._integrityCheck(postmsm=(statetype != 'cluster'))
        if statetype != 'macro' and samplemode != 'random':
            samplemode = 'random'
            logger.warning("'micro' and 'cluster' states incompatible with 'samplemode' other than 'random'. Defaulting to 'random'")

        stConcat = np.concatenate(self.data.St)
        absFrames = []
        relFrames = []
        for i in range(len(states)):
            if frames[i] == 0 and not all:
                continue
            st = states[i]
            if statetype == 'macro':
                self._integrityCheck(postmsm=True)
                (selFr, selMicro) = _sampleMacro(self, st, stConcat, samplemode, frames[i], allframes, replacement)
                absFrames.append(selFr)
            elif statetype == 'micro':
                absFrames.append(_sampleMicro(self, st, stConcat, frames[i], allframes, replacement))
            elif statetype == 'cluster':
                absFrames.append(_sampleCluster(st, stConcat, frames[i], allframes, replacement))
            else:
                raise NameError('No valid state type given (read documentation)')

            if len(absFrames[-1]) == 0:
                raise NameError('No frames could be sampled from {} state {}. State is empty.'.format(statetype, st))

            relFrames.append(self.data.abs2rel(absFrames[-1]))
        return absFrames, relFrames

    def eqDistribution(self, plot=True):
        """ Obtain and plot the equilibrium probabilities of each macrostate

        Parameters
        ----------
        plot : bool, optional, default=True
            Disable plotting of the probabilities by setting it to False

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
        self._integrityCheck(postmsm=True)
        macroeq = np.ones(self.macronum) * -1
        for i in range(self.macronum):
            macroeq[i] = np.sum(self.msm.stationary_distribution[self.macro_ofmicro == i])

        if plot:
            from matplotlib import pylab as plt
            plt.ion()
            plt.figure()
            plt.bar(range(self.macronum), macroeq)
            plt.ylabel('Equilibrium probability')
            plt.xlabel('Macrostates')
            plt.xticks(np.arange(0.4, self.macronum+0.4, 1), range(self.macronum))
            plt.show()
        return macroeq

    def coarseP(self):
        M = self.msm.metastable_memberships
        Pcoarse = np.linalg.inv(M.T.dot(M)).dot(M.T).dot(self.P).dot(M)
        if len(np.where(Pcoarse < 0)[0]) != 0:
            raise NameError('Cannot produce coarse P matrix. Ended up with negative probabilities. Try using less macrostates.')
        return Pcoarse

    def getStates(self, states=None, statetype='macro', wrapsel='protein', alignsel='name CA', alignmol=None, samplemode='weighted', numsamples=50):
        """ Get samples of MSM states in Molecule classes

        Parameters
        ----------
        states : ndarray, optional
            A list of states to visualize
        statetype : ['macro','micro','cluster'], optional
            The type of state to visualize
        wrapsel : str, optional, default='protein'
            A selection to use for wrapping
        alignsel : str, optional, default='name CA'
            A selection used for aligning all frames
        alignmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A reference molecule onto which to align all others
        samplemode : ['weighted','random'], optional, default='weighted'
            How to obtain the samples from the states
        numsamples : int
            Number of samples (conformations) for each state.

        Returns
        -------
        mols : ndarray of :class:`Molecule <htmd.molecule.molecule.Molecule>` objects
            A list of :class:`Molecule <htmd.molecule.molecule.Molecule>` objects containing the samples of each state

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> mols = model.getStates()
        >>> for m in mols:
        >>>     m.view()
        """
        self._integrityCheck(postmsm=(statetype != 'cluster'))
        (single, molfile) = _singleMolfile(self.data.simlist)
        refmol = None
        if not single:
            raise NameError('Visualizer does not support yet visualization of systems with different number of atoms')
        if alignmol is None:
            alignmol = molfile
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
        if len(alignsel) > 0 and len(alignmol) > 0:
            refmol = Molecule(alignmol)

        (tmp, relframes) = self.sampleStates(states, [numsamples]*len(states), statetype=statetype, samplemode=samplemode)

        from joblib import Parallel, delayed
        from htmd.config import _config
        # This loop really iterates over states. sampleStates returns an array of arrays
        # Removed ncpus because it was giving errors on some systems.
        mols = Parallel(n_jobs=1, verbose=11)(delayed(_loadMols)(self, i, rel, molfile, wrapsel, alignsel, refmol)
                                                  for i, rel in enumerate(relframes))
        return np.array(mols, dtype=object)

    def viewStates(self, states=None, statetype='macro', protein=None, ligand=None, viewer=None, mols=None,
                   numsamples=50, wrapsel='protein', alignsel='name CA'):
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
            Atomselection string for the ligand
        viewer : :class:`VMD <htmd.vmdviewer.VMD>` object, optional
            A viewer in which to visualize the states
        mols : ndarray, optional
            An array of :class:`Molecule <htmd.molecule.molecule.Molecule>` objects to visualize
        numsamples : int
            Number of samples (conformations) for each state.
        wrapsel : str, optional, default='protein'
            A selection to use for wrapping
        alignsel : str, optional, default='name CA'
            A selection used for aligning all frames

        Examples
        --------
        >>> model = Model(data)
        >>> model.markovModel(100, 5)
        >>> model.viewStates(protein=True)

        >>> model.viewStates(ligand='resname MOL')
        """
        from htmd.config import _config
        self._integrityCheck(postmsm=(statetype != 'cluster'))

        if _config['viewer'].lower() == 'ngl':
            return self._viewStatesNGL(states, statetype, protein, ligand, mols, numsamples)

        if viewer is None:
            viewer = getCurrentViewer()
        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(states, statetype, numsamples=numsamples, wrapsel=wrapsel, alignsel=alignsel)
        colors = [0, 1, 3, 4, 5, 6, 7, 9]
        for i, s in enumerate(states):
            viewer.loadMol(mols[i], name=statetype+' '+str(states[i]))
            if ligand is not None:
                viewer.rep('ligand', sel=ligand, color=colors[np.mod(i, len(colors))])
            if protein is not None:
                viewer.rep('protein')
            viewer.send('start_sscache')

    def _viewStatesNGL(self, states, statetype, protein, ligand, mols, numsamples):
        if states is None:
            states = range(self.macronum)
        if isinstance(states, int):
            states = [states]
        if mols is None:
            mols = self.getStates(states, statetype, numsamples=min(numsamples, 15))
        colors = [0, 1, 3, 4, 5, 6, 7, 9]
        if protein is None and ligand is None:
            raise NameError('Please provide either the "protein" or "ligand" parameter for viewStates.')
        if protein:
            mol = Molecule()
        if ligand:
            mol = mols[0].copy()
            mol.remove(ligand, _logger=False)
            mol.coords = np.atleast_3d(mol.coords[:, :, 0])
            mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
        for i, s in enumerate(states):
            if protein:
                mol.reps.add(sel='segid ST{}'.format(s), style='NewCartoon', color='Index')
            if ligand:
                mol.reps.add(sel='segid ST{}'.format(s), style='Licorice', color=colors[np.mod(i, len(colors))])
                mols[i].filter(ligand, _logger=False)

            mols[i].set('segid', 'ST{}'.format(s))
            tmpcoo = mols[i].coords
            for j in range(mols[i].numFrames):
                mols[i].coords = np.atleast_3d(tmpcoo[:, :, j])
                mol.append(mols[i])

        w = mol.view(viewer='ngl')
        self._nglButtons(w, statetype, states)
        return w

    def _nglButtons(self, ngl_widget, statetype, states):
        # Adds buttons for enabling and disabling macrostate visualizations
        import IPython.html.widgets as widgets
        from IPython.display import display
        originalreps = ngl_widget.representations.copy()
        otherreps = originalreps[:-len(states)]
        originalreps = originalreps[-len(states):]

        container = []
        for s in states:
            w = widgets.Checkbox(description="{} {}".format(statetype, s))
            w.value = True
            container.append(w)

        def updateReps(name):
            ngl_widget.isClick = True
            reps = otherreps.copy()
            for i, w in enumerate(container):
                if w.value:
                    reps.append(originalreps[i])
            ngl_widget.representations = reps

        for w in container:
            w.on_trait_change(updateReps, "value")
        hb = widgets.HBox(container)
        display(hb)

    def save(self, filename):
        """ Save a :class:`Model` object to disk

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
        f = open(filename, 'wb')
        pickle.dump(self.__dict__, f)
        f.close()

    def load(self, filename):
        """ Load a :class:`MetricData` object from disk

        Parameters
        ----------
        filename : str
            Path to the saved MetricData object

        Examples
        --------
        >>> model = Model()
        >>> model.load('./model.dat')
        """
        import pickle
        f = open(filename, 'rb')
        z = pickle.load(f)
        f.close()
        for k in z:
            self.__dict__[k] = z[k]

    def _integrityCheck(self, postmsm=False, markov=False):
        if postmsm and self._modelid is None:
            raise NameError('You need to call markovModel before calling this command')
        if not markov and self.data._dataid == self.data._clusterid and self.data._dataid != self._clusterid:
            raise NameError('After updating the MetricData object you need to call the markovModel command anew.')
        if self.data._dataid != self.data._clusterid:
            raise NameError('After modifying the data in the MetricData object you need to recluster and reconstruct the markov model.')


def _loadMols(self, i, rel, molfile, wrapsel, alignsel, refmol):
    frames = self.data.rel2sim(rel)
    mol = Molecule(molfile)
    trajs = np.empty(0, dtype=str)
    frs = np.empty(0, dtype=int)
    for f in frames:
        trajs = np.append(trajs, f.sim.trajectory[f.piece])
        frs = np.append(frs, f.frame)
    mol.read(trajs, frames=frs)
    if len(wrapsel) > 0:
        mol.wrap(wrapsel)
    if refmol is not None:
        mol.align(alignsel, refmol=refmol)
    return mol


def getStateStatistic(model, data, states, statetype='macro', weighted=False, method=np.mean, axis=0, existing=False, target=None):
    """ Calculates properties of the states.

    Calculates properties of data corresponding to states. Can calculate for example the mean distances of atoms in a
    state, or the standard deviation of the RMSD in a state.

    Parameters
    ----------
    model : :class:`Model` object
        A model containing the state definitions
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
    statistic : np.ndarray
        A array which contains in each element the desired statistic of a specific state specifies in `states`

    Examples
    --------
    >>> data = MetricDistance.project(sims, 'protein and name CA', 'resname MOL')
    >>> model = Model(data)
    >>> model.markovModel(100, 5)
    >>> # Get the standard deviation of distances in all macrostates
    >>> getStateStatistic(model, data, method=np.std)
    """
    if axis != 0:
        logger.warning('Axis different than 0 might not work correctly yet')
    if len(model.data.dat) > 0 and np.any(model.data.trajLengths != data.trajLengths):
        raise NameError('Data trajectories need to match in size and number to the trajectories in the model')
    stconcat = np.concatenate(model.data.St)
    datconcat = np.concatenate(data.dat)

    statistic = []
    for i, st in enumerate(states):
        if statetype == 'macro':
            frames = model.macro_ofcluster[stconcat] == st
        elif statetype == 'micro':
            frames = model.micro_ofcluster[stconcat] == st
        elif statetype == 'cluster':
            frames = stconcat == st
        else:
            raise NameError('No valid state type given (read documentation)')

        if statetype == 'macro' and weighted:
            statistic.append(_weightedMethod(model, method, stconcat, datconcat, st, axis))
        else:
            if axis is None:
                statistic.append(method(datconcat[frames, ...]))
            else:
                statistic.append(method(datconcat[frames, ...], axis=axis))
    return statistic


def reconstructContactMap(map, datavec):
    """ Plots a given vector as a contact map

    Parameters
    ----------
    map : np.ndarray 2D
        The map from a MetricData object
    datavec : np.ndarray
        The data we want to plot in a 2D map
    """
    map = np.array(map, dtype=int)
    atomidx = np.unique(map.flatten()).astype(int)
    mask = np.zeros(max(atomidx)+1, dtype=int)
    mask[atomidx] = range(len(atomidx))

    # Create a new map which maps from vector indexes to matrix indexes
    newmap = np.zeros(np.shape(map), dtype=int)
    newmap[:, 0] = mask[map[:, 0]]
    newmap[:, 1] = mask[map[:, 1]]

    contactmap = np.zeros((len(atomidx), len(atomidx)))
    for i in range(len(datavec)):
        contactmap[newmap[i, 0], newmap[i, 1]] = datavec[i]
        contactmap[newmap[i, 1], newmap[i, 0]] = datavec[i]

    from matplotlib import pylab as plt
    plt.imshow(contactmap, interpolation='nearest', aspect='equal')
    plt.colorbar()
    #plt.axis('off')
    #plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    #plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    plt.show()


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


def _sampleMacro(obj, macro, stConcat, mode, numFrames, allFrames, replacement):
    if mode == 'random':
        frames = np.where(obj.macro_ofcluster[stConcat] == macro)[0]
        selFrames = _randomSample(frames, numFrames, allFrames, replacement)
        selMicro = obj.micro_ofcluster[stConcat[selFrames]]
    elif mode == 'even':
        micros = np.where(obj.macro_ofmicro == macro)[0]
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames.append(_sampleMicro(obj, micros[i], stConcat, numFrames, allFrames, replacement))
            selMicro.append(np.ones(numFrames) * micros[i])
    elif mode == 'weighted':
        micros = np.where(obj.macro_ofmicro == macro)[0]
        eq = obj.msm.stationary_distribution
        weights = eq[micros] / np.sum(eq[micros])
        framespermicro = np.random.multinomial(numFrames, weights)
        selFrames = []
        selMicro = []
        for i in range(len(micros)):
            selFrames = np.append(selFrames, _sampleMicro(obj, micros[i], stConcat, framespermicro[i], allFrames, replacement))
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
            selFrames = np.append(selFrames, _sampleMicro(obj, micros[i], stConcat, framespermicro[i], allFrames, replacement))
            selMicro = np.append(selMicro, [micros[i]] * framespermicro[i])
    else:
        raise NameError('No valid mode given (read documentation)')
    return selFrames, selMicro


def _sampleMicro(obj, micro, stConcat, numFrames, allFrames, replacement):
    frames = np.where(obj.micro_ofcluster[stConcat] == micro)[0]
    return _randomSample(frames, numFrames, allFrames, replacement)


def _sampleCluster(cluster, stConcat, numFrames, allFrames, replacement):
    frames = np.where(stConcat == cluster)[0]
    return _randomSample(frames, numFrames, allFrames, replacement)


def _randomSample(frames, numFr, allFrames, replacement):
    if numFr == 0:
        return []
    if allFrames or (numFr >= len(frames) and not replacement):
        rnd = list(range(len(frames)))
    else:
        rnd = np.random.randint(len(frames), size=numFr)
    return frames[rnd]


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
        if macrotrajnum[m] <= 3:
            logger.info('Take care! Macro {} has been visited only in {} trajectories:'.format(m, macrotrajnum[m]))
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


