"""
HTMD can handle a large amount of simulations.
Simulation lists allow to create a simple list containing all relevant information about the simulations to later
perform any type of analysis.
"""

# (c) 2015-2016 Acellera Ltd http://www.acellera.com
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
    molfile : str
        The path to the structural information about the simulation. Usually a PDB file
    """
    # __slots__ = ['simid', 'parent', 'input', 'trajectory', 'molfile']

    def __init__(self, simid, parent, input, trajectory, molfile):
        self.simid = simid
        self.parent = parent
        self.input = input
        self.trajectory = trajectory
        self.molfile = molfile

    def __repr__(self):
        if self.parent is None:
            return 'simid = {}\nparent = {}\ninput = {}\ntrajectory = {}\nmolfile = {}'.format(self.simid, self.parent, self.input, self.trajectory, self.molfile)
        else:
            return 'simid = {}\nparent = {}\ninput = {}\ntrajectory = {}\nmolfile = {}'.format(self.simid, self.parent.simid, self.input, self.trajectory, self.molfile)

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


def simlist(datafolders, molfiles, inputfolders=None):
    """Creates a list of simulations

    Parameters
    ----------
    datafolders : str list
        A list of directories, each containing a single trajectory
    molfiles : str list
        A list of pdb files corresponding to the trajectories in dataFolders. Can also be a single string to a single
        structure which corresponds to all trajectories.
    inputfolders : optional, str list
        A list of directories, each containing the input files used to produce the trajectories in dataFolders

    Return
    ------
    sims : np.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A list of simulations

    Examples
    --------
    >>> simlist(glob('./test/data/*/'), glob('./test/input/*/*.pdb'), glob('./test/input/*/'))
    """

    if not datafolders:
        raise FileNotFoundError('No data folders were given, check your arguments.')
    if not molfiles:
        raise FileNotFoundError('No molecule files were given, check your arguments.')
    if isinstance(molfiles, str):
        molfiles = [molfiles]
    if isinstance(datafolders, str):
        datafolders = [datafolders]
    if inputfolders and isinstance(inputfolders, str):
        inputfolders = [inputfolders]

    #Sim = namedtuple('Sim', ['id', 'parent', 'input', 'trajectory', 'molfile'])

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
    for mol in molfiles:
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
        trajectories = _listTrajectories(datanames[k])

        if not trajectories:
            bar.progress()
            continue

        if len(molfiles) > 1:
            if k not in molnames:
                raise FileNotFoundError('Did not find molfile with folder name ' + k + ' in the given glob')
            molfile = molnames[k]
        else:
            molfile = molfiles[0]

        inputf = []
        if inputfolders:
            if k not in inputnames:
                raise FileNotFoundError('Did not find input with folder name ' + k + ' in the given glob')
            inputf = inputnames[k]

        sims.append(Sim(simid=i, parent=None, input=inputf, trajectory=trajectories, molfile=molfile))
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
        ftrajectory = _listTrajectories(path.join(outFolder, name))
        return Sim(simid=sims[i].simid, parent=sims[i], input=None, trajectory=ftrajectory, molfile=fmolfile)

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

    ftrajectory = _listTrajectories(path.join(outFolder, name))
    #bar.progress()
    return Sim(simid=sims[i].simid, parent=sims[i], input=None, trajectory=ftrajectory, molfile=fmolfile)


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
        raise NameError('simFilter: {}. Cannot create filtered.pdb due to problematic pdb: {}'.format(e, sim.molfile))

    if not path.isfile(path.join(outfolder, 'filtered.pdb')):
        if mol.coords.size == 0:  # If we read for example psf or prmtop which have no coords, just add 0s everywhere
            mol.coords = np.zeros((mol.numAtoms, 3, 1), dtype=np.float32)
        mol.write(path.join(outfolder, 'filtered.pdb'), filtsel)


def _listTrajectories(folder):
    trajtypes = ['xtc', 'trr', 'dcd', 'nc', 'dtr', 'crd', 'gro']
    for tt in trajtypes:
        trajectories = glob(path.join(folder, '*.{}'.format(tt)))
        if len(trajectories) > 0:
            return natsort.natsorted(trajectories)


def _simName(foldername):
    # Uses the name of the last folder as the simulation name
    if os.path.isdir(foldername):
        name = os.path.basename(os.path.normpath(foldername))
    else:
        name = os.path.basename(os.path.dirname(foldername))
    return name
