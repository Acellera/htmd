"""
HTMD can handle a large amount of simulations.
Simulation lists allow to create a simple list containing all relevant information about the simulations to later
perform any type of analysis.
"""
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import os.path as path
from os import makedirs
from glob import glob as glob
import natsort
import numpy as np
from htmd.molecule.molecule import Molecule
from joblib import Parallel, delayed
from htmd.progress.progress import ProgressBar
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
        iseq &= self.molfile == other.molfile
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

    if not datafolders:
        raise FileNotFoundError('No data folders were given, check your arguments.')
    if not topologies:
        raise FileNotFoundError('No molecule files were given, check your arguments.')
    topologies = ensurelist(topologies)
    datafolders = ensurelist(datafolders)
    if inputfolders:
        inputfolders = ensurelist(inputfolders)

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
    bar = ProgressBar(len(keys), description='Creating simlist')
    for k in keys:
        trajectories = _autoDetectTrajectories(datanames[k])

        if not trajectories:
            bar.progress()
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
        bar.progress()
    bar.stop()
    logger.debug('Finished listing of simulations.')
    return np.array(sims, dtype=object)


def simfilter(sims, outfolder, filtersel):
    """ Filters a list of simulations generated by :func:`simlist`

    This function takes as input a list of simulations produced by `simList` and writes new trajectories containing only
    the desired atoms in a new directory.

    Parameters
    ----------
    sims : list
        A simulation list produced by the `simList` function
    outFolder : str
        The folder in which to write the modified trajectories
    filterSel : str
        An atomselection string describing the atoms we want to keep

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
        _filterPDBPSF(sims[0], outfolder, filtersel)

    logger.debug('Starting filtering of simulations.')

    from htmd.config import _config
    from htmd.parallelprogress import ParallelExecutor
    aprun = ParallelExecutor(n_jobs=_config['ncpus'])
    filtsims = aprun(total=len(sims), description='Filtering trajectories')(delayed(_filtSim)(i, sims, outfolder, filtersel) for i in range(len(sims)))

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

    maxid = simlist1[-1].simid + 1
    for s in simlist2:
        s.simid = maxid
        maxid += 1
    return np.append(simlist1, simlist2)


def _filtSim(i, sims, outFolder, filterSel):
    name = _simName(sims[i].trajectory[0])
    directory = path.join(outFolder, name)
    if not path.exists(directory):
        makedirs(directory)

    logger.debug('Processing trajectory ' + name)

    fmolfile = path.join(outFolder, 'filtered.pdb')
    (traj, outtraj) = _renameSims(sims[i].trajectory, name, outFolder)
    if not traj:
        ftrajectory = _autoDetectTrajectories(path.join(outFolder, name))
        numframes = _getNumFrames(sims[i], ftrajectory)
        return Sim(simid=sims[i].simid, parent=sims[i], input=None, trajectory=ftrajectory, molfile=fmolfile, numframes=numframes)

    try:
        mol = Molecule(sims[i].molfile)
    except:
        logger.warning('Error! Skipping simulation ' + name)
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
            numframes = int(f.readline())
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


def _filterPDBPSF(sim, outfolder, filtsel):
    try:
        mol = Molecule(sim.molfile)
    except IOError as e:
        raise RuntimeError('simFilter: {}. Cannot create filtered.pdb due to problematic pdb: {}'.format(e, sim.molfile))

    if not path.isfile(path.join(outfolder, 'filtered.pdb')):
        if mol.coords.size == 0:  # If we read for example psf or prmtop which have no coords, just add 0s everywhere
            mol.coords = np.zeros((mol.numAtoms, 3, 1), dtype=np.float32)
        mol.write(path.join(outfolder, 'filtered.pdb'), filtsel)


def _autoDetectTrajectories(folder):
    from htmd.molecule.readers import _TRAJECTORY_READERS
    for tt in _TRAJECTORY_READERS:
        trajectories = glob(path.join(folder, '*.{}'.format(tt)))
        if len(trajectories) > 0:
            return natsort.natsorted(trajectories)


def _autoDetectTopology(folder):
    from htmd.molecule.readers import _TOPOLOGY_READERS
    topotypes = ['pdb'] + list(_TOPOLOGY_READERS.keys())  # Prepending PDB so that it's the default
    topo = None
    for tt in topotypes:
        files = glob(path.join(folder, '*.{}'.format(tt)))
        if len(files) > 0:
            if len(files) > 1:
                logger.warning('Multiple "{}" files were found in folder {}. '
                               'Picking {} as the topology'.format(tt, folder, files[0]))
            topo = files[0]
            break
    if topo is None:
        raise RuntimeError('No topology file found in folder {}. '
                           'Supported extensions are {}'.format(folder, list(_TOPOLOGY_READERS.keys())))
    return topo


def _simName(foldername):
    # Uses the name of the last folder as the simulation name
    if os.path.isdir(foldername):
        name = os.path.basename(os.path.normpath(foldername))
    else:
        name = os.path.basename(os.path.dirname(foldername))
    return name
