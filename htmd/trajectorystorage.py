from natsort import natsorted
from htmd.util import ensurelist
from glob import glob
import numpy as np
import os
import h5py
import logging

logger = logging.getLogger(__name__)


class DuplicateTrajectoryError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



class TrajectoryStorage(object):
    def __init__(self, storagelocation):
        self.storagelocation = storagelocation

        if os.path.splitext(storagelocation)[-1] != '.h5' and os.path.splitext(storagelocation)[-1] != '.hdf5':
            raise RuntimeError('Only storage files with .h5 or .hdf5 extensions are allowed')

        with h5py.File(storagelocation, 'w') as f:
            f.create_group('trajectories')

    def __str__(self):
        with h5py.File(self.storagelocation, 'r') as f:
            ntraj = len(list(f['trajectories']))
        if ntraj == 1:
            string = 'TrajectoryStorage with {} trajectory'.format(ntraj)
        else:
            string = 'TrajectoryStorage with {} trajectories'.format(ntraj)
        return string

    def __repr__(self):
        return '<{}.{} object at {}>\n'.format(self.__class__.__module__, self.__class__.__name__, hex(id(self))) \
               + self.__str__()

    @staticmethod
    def _getFileMd5sum(trajectory):
        import hashlib
        return hashlib.md5(open(trajectory, 'rb').read()).hexdigest()

    @staticmethod
    def _getStringMd5sum(string):
        import hashlib
        return hashlib.md5(string.encode('utf-8')).hexdigest()

    def addTrajectory(self, trajectory, topologyFile, inputFolder=None, sortByName=True):
        """ Add a trajectory to the storage

        Parameters
        ----------
        trajectory : file, folder or list of files
            Can be either of:
                - The path to the file of a trajectory
                - A list of paths of trajectory pieces.
                - A folder containing the trajectory
                - A folder containing a trajectory split in multiple pieces
        topologyFile : str
            The path to the topology file corresponding to the trajectory
        inputFolder : str
            The path to the folder which creates all files that were used to create this trajectory
        sortByName : bool
            Set to False to disable automatical sorting of trajectory pieces by their name
        """
        import uuid

        if isinstance(trajectory, list) or isinstance(trajectory, tuple):  # A list of trajectory pieces
            trajpieces = trajectory
        elif os.path.isfile(trajectory):  # A single trajectory piece
            trajpieces = [trajectory]
        elif os.path.isfolder(trajectory):  # A folder containing the trajectories
            trajpieces = _autoDetectTrajectories(trajectory)

        if len(trajpieces) == 0:
            return

        if sortByName:
            trajpieces = natsorted(trajpieces)

        trajpieces = [os.path.abspath(traj) for traj in trajpieces]

        name = TrajectoryStorage._getStringMd5sum(' '.join(trajpieces))
        
        with h5py.File(self.storagelocation, 'a') as f:
            if name in f['trajectories']:
                raise DuplicateTrajectoryError('A trajectory is already in the dataset with files {}'.format(trajpieces))
            trajgroup = f['trajectories'].create_group(name)
            trajgroup.attrs['topology_file'] = topologyFile
            trajgroup.attrs['parent'] = ''
            if inputFolder is not None:
                trajgroup.attrs['input_folder'] = inputFolder

            filesgroup = trajgroup.create_group('trajectory_files')
            for piece in trajpieces:
                piecegroup = filesgroup.create_group(TrajectoryStorage._getStringMd5sum(piece))
                piecegroup.attrs['trajectory_file'] = piece
                piecegroup.attrs['trajectory_md5sum'] = TrajectoryStorage._getFileMd5sum(piece)


    def autoDetectFiles(self, trajectoryFolders, topologyFolders=None, inputFolders=None):
        """ Automatically match trajectories, topologies and input folders by their respective folder names

        This assumes that each trajectory is stored in a separate folder and that it's topology is stored in a folder of the same name
        and that its inputFolder also has the same name. For example any of the following folder structures should work:

        - trajectoryFolders
            - traj1
                - The trajectory file and a topology file for traj1
            - traj2
                - The trajectory file and a topology file for traj2

        Or if trajectories are kept separate from the topologies

        - trajectoryFolders
            - traj1
                - A trajectory file for traj1
            - traj2
                - A trajectory file for traj2
        - topologyFolders
            - traj1
                - A topology file for traj1
            - traj2
                - A topology file for traj1

        Or if the files used to create the trajectories are available as well they can be linked to them with the inputFolders option

        - trajectoryFolders
            - traj1
                - mytraj1.xtc
            - traj2
                - mytraj2.xtc
        - topologyFolders
            - traj1
                - A topology file for traj1
            - traj2
                - A topology file for traj1
        - inputFolders
            - traj1
                - all files used for creating traj1
            - traj2
                - all files used for creating traj2

        Or if the topology file is already stored in the inputFolder of each trajectory

        - trajectoryFolders
            - traj1
                - mytraj1.xtc
            - traj2
                - mytraj2.xtc
        - inputFolders
            - traj1
                - all files used for creating traj1
            - traj2
                - all files used for creating traj2

        Parameters
        ----------
        trajectoryFolders : list
            A list of folders containing trajectories
        topologyFolders : lsit
            A list of folders containing the topologies of the trajectories
        inputFolders : list
            A list of folders containing the input files used to create the trajectories
        """
        from htmd.util import ensurelist
        import natsort

        if not trajectoryFolders:
            raise FileNotFoundError('No trajectory folders were given, check your arguments.')

        topologyFolders = ensurelist(topologyFolders)
        trajectoryFolders = ensurelist(trajectoryFolders)
        for folder in trajectoryFolders:
            if not os.path.isdir(folder):
                raise NotADirectoryError('{}'.format(folder))
        if inputFolders:
            inputFolders = ensurelist(inputFolders)
            for folder in inputFolders:
                if not os.path.isdir(folder):
                    raise NotADirectoryError('{}'.format(folder))

        # I need to match the simulation names inside the globs given. The
        # reason is that there can be more input folders in the glob than in
        # the data glob due to not having been retrieved. Hence I need to match
        # the folder names.

        # Create a hash map of data folder names
        trajectorynames = dict()
        for folder in trajectoryFolders:
            if _simName(folder) in trajectorynames:
                raise RuntimeError('Duplicate simulation name detected. Cannot name-match directories.')
            trajectorynames[_simName(folder)] = folder

        topologynames = dict()
        for mol in topologyFolders:
            if not os.path.exists(mol):
                raise FileNotFoundError('File {} does not exist'.format(mol))
            topologynames[_simName(mol)] = mol

        if inputFolders:
            inputnames = dict()
            for inputf in inputFolders:
                inputnames[_simName(inputf)] = inputf

        keys = natsort.natsorted(trajectorynames.keys())
        notraj = {}  # Folders which don't have trajectories
        from tqdm import tqdm
        for k in tqdm(keys, desc='Detecting trajectories'):
            trajectories = _autoDetectTrajectories(trajectorynames[k])

            if not trajectories:
                notraj[trajectorynames[k]] = True
                continue

            if len(topologyFolders) > 1:
                if k not in topologynames:
                    raise FileNotFoundError('Did not find molfile with folder name ' + k + ' in the given glob')
                topologyfile = topologynames[k]
            else:
                topologyfile = topologyFolders[0]

            if os.path.isdir(topologyfile):
                topologyfile = _autoDetectTopology(topologyfile)

            inputf = None
            if inputFolders:
                if k not in inputnames:
                    raise FileNotFoundError('Did not find input with folder name ' + k + ' in the given glob')
                inputf = inputnames[k]

            try:
                self.addTrajectory(trajectories, topologyfile, inputf)
            except DuplicateTrajectoryError as e:
                logger.info('Skipping trajectory {} since it\'s already in the dataset'.format(trajectorynames[k]))


        logger.info('{} folder(s) were missing trajectories'.format(len(notraj)))


def _autoDetectTrajectories(folder):
    from moleculekit.readers import _TRAJECTORY_READERS
    import natsort
    for tt in _TRAJECTORY_READERS:
        trajectories = glob(os.path.join(folder, '*.{}'.format(tt)))
        if len(trajectories) > 0:
                return trajectories


from moleculekit.readers import _TOPOLOGY_READERS
__readers = list(_TOPOLOGY_READERS.keys())
__defaultReaders = ['pdb', 'prmtop', 'psf']
__otherReaders = list(np.setdiff1d(__readers, __defaultReaders))
__topotypes =  __defaultReaders + __otherReaders  # Prepending PDB, PSF, PRMTOP so that they are the default

def _autoDetectTopology(folder):
    topo = {}
    for tt in __topotypes:
        files = glob(os.path.join(folder, '*.{}'.format(tt)))
        if len(files) > 0:
            if len(files) > 1:
                logger.warning('Multiple "{}" files were found in folder {}. '
                               'Picking {} as the topology'.format(tt, folder, files[0]))
            topo[tt] = files[0]
    if len(topo) == 0:
        raise RuntimeError('No topology file found in folder {}. '
                           'Supported extensions are {}'.format(folder, list(_TOPOLOGY_READERS.keys())))
    return list(topo.values())


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

def _simName(foldername):
    # Uses the name of the last folder as the simulation name
    if os.path.isdir(foldername):
        name = os.path.basename(os.path.normpath(foldername))
    else:
        name = os.path.basename(os.path.dirname(foldername))
    return name