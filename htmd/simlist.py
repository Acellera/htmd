"""
HTMD can handle a large amount of simulations.
Simulation lists allow to create a simple list containing all relevant information about the simulations to later
perform any type of analysis.
"""
# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import os.path as path
from os import makedirs
from glob import glob as glob
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Frame(object):
    """ Class used for referencing a specific frame of a :class:`Sim <htmd.simlist.Sim>` object object.

    Parameters
    ----------
    sim : :class:`Sim <htmd.simlist.Sim>` object
        The simulation.
    piece : int
        Which trajectory piece the frame corresponds to.
    frame : int
        The frame of the specific trajectory piece.
    """
    __slots__ = ['sim', 'piece', 'frame']

    def __init__(self, sim, piece, frame):
        self.sim = sim
        self.piece = piece
        self.frame = frame

    def __repr__(self):
        return 'sim = {}\npiece = {}\nframe = {}'.format(self.sim, self.piece, self.frame)


class Sim(object):
    """ Information class for a single simulation.

    Do not use directly. Objects of this class are constructed by the :func:`simlist` and :func:`simfilter` functions.
    This class is used for storing information on simulations. This information includes the path to the simulation,
    the path to a structure file (pdb) which corresponds to the simulation, the folder containing the input files used
    to generate the simulation (useful for adaptive), the parent of the simulation (if it was filtered it will point to
    the original simulation) and a unique simulation id.

    Attributes
    ----------
    simid : int
        A unique simulation ID
    parent : :class:`Sim <htmd.simlist.Sim>` object
        The parent of the simulations
    input : str
        The path to the input folder which generated this simulation
    trajectory : list
        A list of trajectory files
    molfile : str
        The path to the structural information about the simulation. Usually a PDB file
    numframes : list
        Number of frames in trajectories
    """
    # __slots__ = ['simid', 'parent', 'input', 'trajectory', 'molfile']

    def __init__(self, simid, parent, input, trajectory, molfile, numframes=None):
        self.simid = simid
        self.parent = parent
        self.input = input
        self.trajectory = trajectory
        self.molfile = molfile
        self.numframes = numframes

    def __repr__(self):
        parent = self.parent
        if self.parent is not None:
            parent = self.parent.simid
        return '\nsimid = {}\nparent = {}\ninput = {}\ntrajectory = {}\nmolfile = {}\nnumframes = {}\n'.format(self.simid,
                                                                                                         parent,
                                                                                                         self.input,
                                                                                                         self.trajectory,
                                                                                                         self.molfile,
                                                                                                         self.numframes)

    def __eq__(self, other):
        iseq = True
        iseq &= self.simid == other.simid
        iseq &= np.all([x == y for x, y in zip(self.molfile, other.molfile)])
        iseq &= self.input == other.input
        iseq &= self.parent == other.parent
        iseq &= len(self.trajectory) == len(other.trajectory)
        if not iseq:
            return False
        for i in range(len(self.trajectory)):
            iseq &= self.trajectory[i] == other.trajectory[i]
        return iseq

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)


class _simlist2(np.ndarray):
    def __new__(cls, datafolders, topologies, inputfolders=None):
        obj = simlist(datafolders, topologies, inputfolders).view(cls)
        return obj

    # def __array_finalize__(self, obj):
    #     # see InfoArray.__array_finalize__ for comments
    #     if obj is None: return
    #     self.info = getattr(obj, 'info', None)

    def numFrames(self):
        return [x.numframes for x in self]

    def append(self):
        pass


def simlist(datafolders, topologies, inputfolders=None):
    """Creates a list of simulations

    Parameters
    ----------
    datafolders : str list
        A list of directories, each containing a single trajectory
    topologies : str list
        A list of topology files or folders containing a topology file corresponding to the trajectories in dataFolders.
        Can also be a single string to a single structure which corresponds to all trajectories.
    inputfolders : optional, str list
        A list of directories, each containing the input files used to produce the trajectories in dataFolders

    Return
    ------
    sims : np.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A list of simulations

    Examples
    --------
    >>> simlist(glob('./test/data/*/'), glob('./test/input/*/'), glob('./test/input/*/'))
    >>> simlist(glob('./test/data/*/'), glob('./test/input/*/*.pdb'), glob('./test/input/*/'))
    """
    from htmd.util import ensurelist
    import natsort

    if not datafolders:
        raise FileNotFoundError('No data folders were given, check your arguments.')
    if not topologies:
        raise FileNotFoundError('No molecule files were given, check your arguments.')
    topologies = ensurelist(topologies)
    datafolders = ensurelist(datafolders)
    for folder in datafolders:
        if not os.path.isdir(folder):
            raise NotADirectoryError('{}'.format(folder))
    if inputfolders:
        inputfolders = ensurelist(inputfolders)
        for folder in inputfolders:
            if not os.path.isdir(folder):
                raise NotADirectoryError('{}'.format(folder))

    # I need to match the simulation names inside the globs given. The
    # reason is that there can be more input folders in the glob than in
    # the data glob due to not having been retrieved. Hence I need to match
    # the folder names.

    # Create a hash map of data folder names
    datanames = dict()
    for folder in datafolders:
        if _simName(folder) in datanames:
            raise RuntimeError('Duplicate simulation name detected. Cannot name-match directories.')
        datanames[_simName(folder)] = folder

    molnames = dict()
    for mol in topologies:
        if not os.path.exists(mol):
            raise FileNotFoundError('File {} does not exist'.format(mol))
        molnames[_simName(mol)] = mol

    if inputfolders:
        inputnames = dict()
        for inputf in inputfolders:
            inputnames[_simName(inputf)] = inputf

    logger.debug('Starting listing of simulations.')
    sims = []

    keys = natsort.natsorted(datanames.keys())
    i = 0
    from tqdm import tqdm
    for k in tqdm(keys, desc='Creating simlist'):
        trajectories = _autoDetectTrajectories(datanames[k])

        if not trajectories:
            continue

        if len(topologies) > 1:
            if k not in molnames:
                raise FileNotFoundError('Did not find molfile with folder name ' + k + ' in the given glob')
            molfile = molnames[k]
        else:
            molfile = topologies[0]

        if os.path.isdir(molfile):
            molfile = _autoDetectTopology(molfile)

        inputf = []
        if inputfolders:
            if k not in inputnames:
                raise FileNotFoundError('Did not find input with folder name ' + k + ' in the given glob')
            inputf = inputnames[k]

        numframes = [_readNumFrames(f) for f in trajectories]
        sims.append(Sim(simid=i, parent=None, input=inputf, trajectory=trajectories, molfile=molfile, numframes=numframes))
        i += 1
    logger.debug('Finished listing of simulations.')
    return np.array(sims, dtype=object)


def simfilter(sims, outfolder, filtersel, njobs=None):
    """ Filters a list of simulations generated by :func:`simlist`

    This function takes as input a list of simulations produced by `simList` and writes new trajectories containing only
    the desired atoms in a new directory.

    Parameters
    ----------
    sims : list
        A simulation list produced by the `simList` function
    outfolder : str
        The folder in which to write the modified trajectories
    filtersel : str
        Atom selection string describing the atoms we want to keep.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    njobs : int
        Number of parallel jobs to spawn for filtering of trajectories. If None it will use the default from htmd.config.

    Returns
    -------
    fsims : np.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A list of filtered simulations

    Example
    -------
    >>> sims  = simlist(glob('data/*/'), glob('input/*/structure.pdb'))
    >>> fsims = simfilter(sims, 'filtered', filtersel='not water')
    """
    if not path.exists(outfolder):
        makedirs(outfolder)

    if len(sims) > 0:
        _filterTopology(sims[0], outfolder, filtersel)

    logger.debug('Starting filtering of simulations.')

    from htmd.config import _config
    from htmd.parallelprogress import ParallelExecutor, delayed
    aprun = ParallelExecutor(n_jobs=njobs if njobs is not None else _config['njobs'])
    filtsims = aprun(total=len(sims), desc='Filtering trajectories')(delayed(_filtSim)(i, sims, outfolder, filtersel) for i in range(len(sims)))

    logger.debug('Finished filtering of simulations')
    return np.array(filtsims)


def simmerge(simlist1, simlist2):
    """ Merges two simlists by updating their `simid` fields

    Parameters
    ----------
    simlist1 : numpy array of :class:`Sim <htmd.simlist.Sim>` objects
        First list
    simlist2 : numpy array of :class:`Sim <htmd.simlist.Sim>` objects
        Second list

    Returns
    -------
    newlist : np.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A new list containing all simulations
    """
    if len(simlist1) == 0:
        return simlist2
    if len(simlist2) == 0:
        return simlist1

    newsimlist = np.append(simlist1, simlist2)
    for i, s in enumerate(newsimlist):
        s.simid = i
    return newsimlist


def _filtSim(i, sims, outFolder, filterSel):
    name = _simName(sims[i].trajectory[0])
    directory = path.join(outFolder, name)
    if not path.exists(directory):
        makedirs(directory)

    logger.debug('Processing trajectory ' + name)

    fmolfile = [path.join(outFolder, 'filtered.psf'), path.join(outFolder, 'filtered.pdb')]
    (traj, outtraj) = _renameSims(sims[i].trajectory, name, outFolder)
    if not traj:
        ftrajectory = _autoDetectTrajectories(path.join(outFolder, name))
        numframes = _getNumFrames(sims[i], ftrajectory)
        return Sim(simid=sims[i].simid, parent=sims[i], input=None, trajectory=ftrajectory, molfile=fmolfile, numframes=numframes)

    try:
        from moleculekit.molecule import Molecule
        mol = Molecule(sims[i].molfile)
    except Exception as e:
        logger.warning('Error! Skipping simulation ' + name + ' due to error: ' + str(e))
        return

    sel = mol.atomselect(filterSel)

    for j in range(0, len(traj)):
        try:
            mol.read(traj[j])
        except IOError as e:
            logger.warning('{}, skipping trajectory'.format(e))
            break

        mol.write(outtraj[j], sel)

    ftrajectory = _autoDetectTrajectories(path.join(outFolder, name))
    numframes = _getNumFrames(sims[i], ftrajectory)
    return Sim(simid=sims[i].simid, parent=sims[i], input=None, trajectory=ftrajectory, molfile=fmolfile, numframes=numframes)


def _getNumFrames(sim, trajectories):
    if trajectories is None:
        return None
    numframes = sim.numframes
    if numframes is None or np.any([f is None for f in sim.numframes]):
        numframes = [_readNumFrames(f) for f in trajectories]
    return numframes


def _readNumFrames(filepath):
    filepath = os.path.abspath(filepath)
    filedir = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    numframefile = os.path.join(filedir, '.{}.numframes'.format(basename))
    numframes = None
    if os.path.exists(numframefile):
        with open(numframefile, 'r') as f:
            try:
                numframes = int(f.readline())
            except:
                raise RuntimeError('{} does not contain an integer number of frames. '
                                   'Please delete this file.'.format(numframefile))
    return numframes


def _renameSims(trajectory, simname, outfolder):
    traj = list()
    outtraj = list()

    for t in range(0, len(trajectory)):
        (tmp, fname) = path.split(trajectory[t])
        (fname, ext) = path.splitext(fname)
        outname = path.join(outfolder, simname, fname + '.filtered{}'.format(ext))

        if not path.isfile(outname) or (path.getmtime(outname) < path.getmtime(trajectory[t])):
            traj.append(trajectory[t])
            outtraj.append(outname)

    return traj, outtraj


def _filterTopology(sim, outfolder, filtsel):
    from htmd.util import ensurelist
    try:
        from moleculekit.molecule import Molecule
        mol = Molecule(sim.molfile)
    except IOError as e:
        raise RuntimeError('simFilter: {}. Cannot read topology file {}'.format(e, sim.molfile))

    if mol.coords.size == 0:  # If we read for example psf or prmtop which have no coords, just add 0s everywhere
        mol.coords = np.zeros((mol.numAtoms, 3, 1), dtype=np.float32)

    extensions = ['pdb', 'psf']  # Adding pdb and psf to make sure they are always written
    for m in ensurelist(sim.molfile):
        extensions.append(os.path.splitext(m)[1][1:])

    for ext in list(set(extensions)):
        filttopo = path.join(outfolder, 'filtered.{}'.format(ext))
        if not path.isfile(filttopo):
            try:
                mol.write(filttopo, filtsel)
            except Exception as e:
                logger.warning('Filtering was not able to write {} due to error: {}'.format(filttopo, e))


def _autoDetectTrajectories(folder):
    from moleculekit.readers import _TRAJECTORY_READERS
    import natsort
    for tt in _TRAJECTORY_READERS:
        if tt in ("xsc",):  # Some trajectory readers don't really load trajectories like xsc
            continue
        trajectories = glob(path.join(folder, '*.{}'.format(tt)))
        if len(trajectories) > 0:
            return natsort.natsorted(trajectories)


from moleculekit.readers import _TOPOLOGY_READERS
__readers = list(_TOPOLOGY_READERS.keys())
__defaultReaders = ['pdb', 'prmtop', 'psf']
__otherReaders = list(np.setdiff1d(__readers, __defaultReaders))
__topotypes =  __defaultReaders + __otherReaders  # Prepending PDB, PSF, PRMTOP so that they are the default

def _autoDetectTopology(folder):
    topo = {}
    for tt in __topotypes:
        files = glob(path.join(folder, '*.{}'.format(tt)))
        if len(files) > 0:
            if len(files) > 1:
                logger.warning('Multiple "{}" files were found in folder {}. '
                               'Picking {} as the topology'.format(tt, folder, files[0]))
            topo[tt] = files[0]
    if len(topo) == 0:
        raise RuntimeError('No topology file found in folder {}. '
                           'Supported extensions are {}'.format(folder, list(_TOPOLOGY_READERS.keys())))
    return list(topo.values())


def _simName(foldername):
    # Uses the name of the last folder as the simulation name
    if os.path.isdir(foldername):
        name = os.path.basename(os.path.normpath(foldername))
    else:
        name = os.path.basename(os.path.dirname(foldername))
    return name


import unittest
class _TestSimlist(unittest.TestCase):
    def test_simlist_auto_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(glob(path.join(home(dataDir='adaptive'), 'data', '*', '')), glob(path.join(home(dataDir='adaptive'), 'input', '*')))
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert len(sims[0].molfile) == 2
        assert _singleMolfile(sims)[0]

    def test_simlist_many_structures(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(glob(path.join(home(dataDir='adaptive'), 'data', '*', '')), glob(path.join(home(dataDir='adaptive'), 'input', '*', 'structure.pdb')))
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]

    def test_simlist_single_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(glob(path.join(home(dataDir='adaptive'), 'data', '*', '')), path.join(home(dataDir='adaptive'), 'input', 'e1s1_1', 'structure.pdb'))
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]


if __name__ == '__main__':
    unittest.main(verbosity=2)


