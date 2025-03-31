"""
HTMD can handle a large amount of simulations.
Simulation lists allow to create a simple list containing all relevant information about the simulations to later
perform any type of analysis.
"""

# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.readers import _TOPOLOGY_READERS
import os
import h5py
import os.path as path
from os import makedirs
from glob import glob as glob
import numpy as np
import unittest
import logging

logger = logging.getLogger(__name__)


class Frame(object):
    """Class used for referencing a specific frame of a :class:`Sim <htmd.simlist.Sim>` object object.

    Parameters
    ----------
    sim : :class:`Sim <htmd.simlist.Sim>` object
        The simulation.
    piece : int
        Which trajectory piece the frame corresponds to.
    frame : int
        The frame of the specific trajectory piece.
    """

    __slots__ = ["sim", "piece", "frame"]

    def __init__(self, sim, piece, frame):
        self.sim = sim
        self.piece = piece
        self.frame = frame

    def __repr__(self):
        return f"sim = {self.sim}\npiece = {self.piece}\nframe = {self.frame}"


class Sim(object):
    """Information on a single simulation trajectory.

    This class is used for storing information related to simulations. This information includes the path to the simulation,
    the path to a structure file (pdb) which corresponds to the simulation, the folder containing the input files used
    to generate the simulation (useful for adaptive), the parent of the simulation (if it was filtered it will point to
    the original simulation) and a unique simulation id.
    """

    def __init__(
        self,
        trajf,
        topology=None,
        inputf=None,
        _simid=None,
        _parent=None,
        _input=None,
        _trajectory=None,
        _molfile=None,
        _numframes=None,
    ):
        import uuid
        from htmd.util import ensurelist

        if trajf is not None:
            if isinstance(trajf, (list, tuple, np.ndarray)):
                _trajectory = [f for f in trajf]
            elif os.path.isdir(trajf):
                _trajectory = _autoDetectTrajectories(trajf)
            elif os.path.isfile(trajf):
                _trajectory = [trajf]
            else:
                raise ValueError(f"Invalid input structure: {trajf}")

            if topology is None:
                if os.path.isfile(trajf):
                    trajf = os.path.dirname(os.path.abspath(trajf))
                _molfile = _autoDetectTopology(trajf)
            else:
                if isinstance(topology, (list, tuple, np.ndarray)):
                    _molfile = [f for f in topology]
                elif os.path.isdir(topology):
                    _molfile = _autoDetectTopology(topology)
                else:
                    _molfile = [topology]

            if inputf is not None:
                _input = inputf
            else:
                _input = os.path.dirname(os.path.abspath(_molfile[0]))

            _numframes = [_readNumFrames(f) for f in _trajectory]
            if _simid is None:
                _simid = str(uuid.uuid4())

        self.simid = str(_simid)
        self.parent = _parent
        self.input = _input
        self.trajectory = _trajectory
        self.molfile = _molfile
        self.numframes = _numframes

    def toHDF5(self, h5group: h5py.Group):
        from moleculekit.util import ensurelist

        h5group.attrs["simid"] = self.simid
        h5group.create_dataset("input", data=ensurelist(self.input))
        h5group.create_dataset("trajectory", data=ensurelist(self.trajectory))
        h5group.create_dataset("molfile", data=ensurelist(self.molfile))
        h5group.create_dataset("numframes", data=ensurelist(self.numframes))
        if self.parent:
            parentgroup = h5group.create_group("parent")
            self.parent.toHDF5(parentgroup)

    @staticmethod
    def fromHDF5(h5group: h5py.Group):
        return Sim(
            trajf=None,
            _simid=str(h5group.attrs["simid"]),
            _input=[x.decode("utf-8") for x in h5group["input"]],
            _trajectory=[x.decode("utf-8") for x in h5group["trajectory"]],
            _molfile=[x.decode("utf-8") for x in h5group["molfile"]],
            _numframes=list(h5group["numframes"]),
            _parent=(
                None if "parent" not in h5group else Sim.fromHDF5(h5group["parent"])
            ),
        )

    def __repr__(self):
        parent = self.parent
        if self.parent is not None:
            parent = self.parent.simid
        return (
            f"\nsimid = {self.simid}\nparent = {parent}\ninput = {self.input}\ntrajectory = {self.trajectory}"
            f"\nmolfile = {self.molfile}\nnumframes = {self.numframes}\n"
        )

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


def simlist(datafolders, topologies, inputfolders=None):
    """Creates a list of simulations

    Parameters
    ----------
    datafolders : str list
        A list of directories, each containing a single trajectory
    topologies : str list
        A list of topology files or folders containing a topology file corresponding to the trajectories in dataFolders.
        Can also be a single path to a single structure which corresponds to all trajectories.
        If the single path is a folder, topology files will be auto-detected in that folder.
    inputfolders : optional, str list
        A list of directories, each containing the input files used to produce the trajectories in dataFolders

    Return
    ------
    sims : np.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A list of simulations

    Examples
    --------
    >>> simlist(glob('./test/data/*/'), './test/data/0/')
    >>> simlist(glob('./test/data/*/'), glob('./test/input/*/'), glob('./test/input/*/'))
    >>> simlist(glob('./test/data/*/'), glob('./test/input/*/*.pdb'), glob('./test/input/*/'))
    """
    from htmd.util import ensurelist
    import natsort

    if not datafolders:
        raise FileNotFoundError("No data folders were given, check your arguments.")
    if not topologies:
        raise FileNotFoundError("No molecule files were given, check your arguments.")
    topologies = ensurelist(topologies)
    datafolders = ensurelist(datafolders)
    for folder in datafolders:
        if not os.path.isdir(folder):
            raise NotADirectoryError(folder)
    if inputfolders:
        inputfolders = ensurelist(inputfolders)
        for folder in inputfolders:
            if not os.path.isdir(folder):
                raise NotADirectoryError(folder)

    # I need to match the simulation names inside the globs given. The
    # reason is that there can be more input folders in the glob than in
    # the data glob due to not having been retrieved. Hence I need to match
    # the folder names.

    # Create a hash map of data folder names
    datanames = dict()
    for folder in datafolders:
        if _simName(folder) in datanames:
            raise RuntimeError(
                "Duplicate simulation name detected. Cannot name-match directories."
            )
        datanames[_simName(folder)] = folder

    molnames = dict()
    for mol in topologies:
        if not os.path.exists(mol):
            raise FileNotFoundError(f"File {mol} does not exist")
        molnames[_simName(mol)] = mol

    if inputfolders:
        inputnames = dict()
        for inputf in inputfolders:
            inputnames[_simName(inputf)] = inputf

    logger.debug("Starting listing of simulations.")
    sims = []

    keys = natsort.natsorted(datanames.keys())
    i = 0
    from tqdm import tqdm

    for k in tqdm(keys, desc="Creating simlist"):
        trajectories = _autoDetectTrajectories(datanames[k])

        if not trajectories:
            continue

        if len(topologies) > 1:
            if k not in molnames:
                raise FileNotFoundError(
                    f"Did not find molfile with folder name {k} in the given glob"
                )
            molfile = molnames[k]
        else:
            molfile = topologies[0]

        if os.path.isdir(molfile):
            molfile = _autoDetectTopology(molfile)

        inputf = []
        if inputfolders:
            if k not in inputnames:
                raise FileNotFoundError(
                    f"Did not find input with folder name {k} in the given glob"
                )
            inputf = inputnames[k]

        numframes = [_readNumFrames(f) for f in trajectories]
        sims.append(
            Sim(
                trajf=None,
                _simid=i,
                _parent=None,
                _input=inputf,
                _trajectory=trajectories,
                _molfile=molfile,
                _numframes=numframes,
            )
        )
        i += 1
    logger.debug("Finished listing of simulations.")
    return sims


def simfilter(sims, outfolder, filtersel, njobs=None):
    """Filters a list of simulations generated by :func:`simlist`

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

    logger.debug("Starting filtering of simulations.")

    from htmd.config import _config
    from htmd.parallelprogress import ParallelExecutor, delayed

    aprun = ParallelExecutor(n_jobs=njobs if njobs is not None else _config["njobs"])
    filtsims = aprun(total=len(sims), desc="Filtering trajectories")(
        delayed(_filtSim)(i, sims, outfolder, filtersel) for i in range(len(sims))
    )

    logger.debug("Finished filtering of simulations")
    return np.array(filtsims)


def simmerge(simlist1, simlist2):
    """Merges two simlists by updating their `simid` fields

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

    logger.debug(f"Processing trajectory {name}")

    fmolfile = [
        path.join(outFolder, "filtered.psf"),
        path.join(outFolder, "filtered.pdb"),
    ]
    (traj, outtraj) = _renameSims(sims[i].trajectory, name, outFolder)
    if not traj:
        ftrajectory = _autoDetectTrajectories(path.join(outFolder, name))
        numframes = _getNumFrames(sims[i], ftrajectory)
        return Sim(
            trajf=None,
            _simid=sims[i].simid,
            _parent=sims[i],
            _input=None,
            _trajectory=ftrajectory,
            _molfile=fmolfile,
            _numframes=numframes,
        )

    try:
        from moleculekit.molecule import Molecule

        mol = Molecule(sims[i].molfile)
    except Exception as e:
        logger.warning(f"Error! Skipping simulation {name} due to error: {e}")
        return

    sel = mol.atomselect(filterSel)

    for j in range(0, len(traj)):
        try:
            mol.read(traj[j])
        except IOError as e:
            logger.warning(f"{e}, skipping trajectory")
            break

        mol.write(outtraj[j], sel)

    ftrajectory = _autoDetectTrajectories(path.join(outFolder, name))
    numframes = _getNumFrames(sims[i], ftrajectory)
    return Sim(
        trajf=None,
        _simid=sims[i].simid,
        _parent=sims[i],
        _input=None,
        _trajectory=ftrajectory,
        _molfile=fmolfile,
        _numframes=numframes,
    )


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
    numframefile = os.path.join(filedir, f".{basename}.numframes")
    numframes = None
    if os.path.exists(numframefile):
        with open(numframefile, "r") as f:
            try:
                numframes = int(f.readline())
            except Exception:
                raise RuntimeError(
                    f"{numframefile} does not contain an integer number of frames. Please delete this file."
                )
    return numframes


def _renameSims(trajectory, simname, outfolder):
    traj = list()
    outtraj = list()

    for t in range(0, len(trajectory)):
        (_, fname) = path.split(trajectory[t])
        (fname, ext) = path.splitext(fname)
        outname = path.join(outfolder, simname, f"{fname}.filtered{ext}")

        if not path.isfile(outname) or (
            path.getmtime(outname) < path.getmtime(trajectory[t])
        ):
            traj.append(trajectory[t])
            outtraj.append(outname)

    return traj, outtraj


def _filterTopology(sim, outfolder, filtsel):
    from htmd.util import ensurelist

    try:
        from moleculekit.molecule import Molecule

        mol = Molecule(sim.molfile)
    except IOError as e:
        raise RuntimeError(f"simFilter: {e}. Cannot read topology file {sim.molfile}")

    if (
        mol.coords.size == 0
    ):  # If we read for example psf or prmtop which have no coords, just add 0s everywhere
        mol.coords = np.zeros((mol.numAtoms, 3, 1), dtype=np.float32)

    extensions = [
        "pdb",
        "psf",
    ]  # Adding pdb and psf to make sure they are always written
    for m in ensurelist(sim.molfile):
        extensions.append(os.path.splitext(m)[1][1:])

    for ext in list(set(extensions)):
        filttopo = path.join(outfolder, f"filtered.{ext}")
        if not path.isfile(filttopo):
            try:
                mol.write(filttopo, filtsel)
            except Exception as e:
                logger.warning(
                    f"Filtering was not able to write {filttopo} due to error: {e}"
                )


def _autoDetectTrajectories(folder):
    from moleculekit.readers import _TRAJECTORY_READERS
    import natsort

    for tt in _TRAJECTORY_READERS:
        if tt in ("xsc",):
            # Some trajectory readers don't actually load trajectories (i.e. xsc files)
            continue
        trajectories = glob(path.join(folder, f"*.{tt}"))
        if len(trajectories) > 0:
            return natsort.natsorted(trajectories)


def _autoDetectTopology(folder):
    topo = {}
    for tt in ("prmtop", "psf", "mae", "pdb"):
        files = glob(path.join(folder, f"*.{tt}"))
        if len(files) > 0:
            if len(files) > 1:
                logger.warning(
                    f'Multiple "{tt}" files were found in folder {folder}. Picking {files[0]} as the topology'
                )
            topo[tt] = files[0]
            break

    if len(topo) == 0:
        raise RuntimeError(
            f"No topology file found in folder {folder}. "
            f"Supported extensions are {list(_TOPOLOGY_READERS.keys())}"
        )

    for tt in ("pdb", "coor"):
        files = glob(path.join(folder, f"*.{tt}"))
        if len(files) > 0:
            if len(files) > 1:
                logger.warning(
                    f'Multiple "{tt}" files were found in folder {folder}. Picking {files[0]} as the topology'
                )
            topo[tt] = files[0]
            break

    return list(topo.values())


def _simName(foldername):
    # Uses the name of the last folder as the simulation name
    if os.path.isdir(foldername):
        name = os.path.basename(os.path.normpath(foldername))
    else:
        name = os.path.basename(os.path.dirname(foldername))
    return name


class _TestSimlist(unittest.TestCase):
    def test_simlist_auto_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            glob(path.join(home(dataDir="adaptive"), "input", "*")),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert len(sims[0].molfile) == 2
        assert _singleMolfile(sims)[0]

    def test_simlist_many_structures(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            glob(path.join(home(dataDir="adaptive"), "input", "*", "structure.pdb")),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]

    def test_simlist_single_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            path.join(home(dataDir="adaptive"), "input", "e1s1_1", "structure.pdb"),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]

    def test_autodetect(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile
        from natsort import natsorted

        datafolders = natsorted(
            glob(path.join(home(dataDir="adaptive"), "data", "*", ""))
        )
        inputfolders = natsorted(
            glob(path.join(home(dataDir="adaptive"), "input", "*"))
        )
        sims = []
        for ff in zip(datafolders, inputfolders):
            sims.append(Sim(ff[0], ff[1], ff[1]))

        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert len(sims[0].molfile) == 2
        assert _singleMolfile(sims)[0]

        sims2 = []
        for ff in zip(datafolders, inputfolders):
            sims2.append(Sim(ff[0], ff[1]))

        for s1, s2 in zip(sims, sims2):
            assert s1.parent == s2.parent
            assert s1.input == s2.input
            assert s1.trajectory == s2.trajectory
            assert s1.molfile == s2.molfile
            assert s1.numframes == s2.numframes

    def test_sim(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            with open(path.join(tmpdir, "test.pdb"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.psf"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.prmtop"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.mae"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.coor"), "w") as f:
                f.write("fake")
            for i in range(11):
                with open(path.join(tmpdir, f"test_{i}.xtc"), "w") as f:
                    f.write("fake")

            # Passing a single directory
            sim = Sim(tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.prmtop"),
                path.join(tmpdir, "test.pdb"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir

            # Passing same directory to all arguments
            sim = Sim(tmpdir, tmpdir, tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.prmtop"),
                path.join(tmpdir, "test.pdb"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir

            # Passing single topology file and arbitrary input folder
            sim = Sim(tmpdir, path.join(tmpdir, "test.pdb"), "/tmp")
            assert sim.molfile == [path.join(tmpdir, "test.pdb")]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == "/tmp"

        # Coor file reading
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(path.join(tmpdir, "test.psf"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.mae"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.coor"), "w") as f:
                f.write("fake")
            for i in range(11):
                with open(path.join(tmpdir, f"test_{i}.xtc"), "w") as f:
                    f.write("fake")

            sim = Sim(tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.psf"),
                path.join(tmpdir, "test.coor"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir


if __name__ == "__main__":
    unittest.main(verbosity=2)
