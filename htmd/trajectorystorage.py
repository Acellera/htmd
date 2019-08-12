from natsort import natsorted
from htmd.util import ensurelist
from glob import glob
from tqdm import tqdm
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


class _Trajectory(object):
    def __init__(self, topologyFiles, trajectoryFiles=None, inputFolder=None, parent=''):
        self.topologyFiles = topologyFiles
        self.trajectoryFiles = trajectoryFiles
        self.inputFolder = inputFolder
        self.parent = parent

    @staticmethod
    def _getFileMd5sum(trajectory):
        import hashlib
        return hashlib.md5(open(trajectory, 'rb').read()).hexdigest()

    @staticmethod
    def _getStringMd5sum(string):
        import hashlib
        return hashlib.md5(string.encode('utf-8')).hexdigest()

    def getName(self):
        return _Trajectory._getStringMd5sum(' '.join(self.trajectoryFiles))

    def writeToH5(self, parentgroup):
        trajgroup = parentgroup.create_group(self.getName())
        trajgroup.attrs['parent'] = self.parent
        if self.inputFolder is not None:
            trajgroup.attrs['input_folder'] = self.inputFolder

        filesgroup = trajgroup.create_group('trajectory_files')
        for i, piece in enumerate(self.trajectoryFiles):
            piecegroup = filesgroup.create_group(str(i))
            piecegroup.attrs['trajectory_file'] = piece
            piecegroup.attrs['trajectory_md5sum'] = _Trajectory._getFileMd5sum(piece)

        filesgroup = trajgroup.create_group('topology_files')
        for i, file in enumerate(self.topologyFiles):
            filegroup = filesgroup.create_group(str(i))
            filegroup.attrs['topology_file'] = file
            filegroup.attrs['topology_md5sum'] = _Trajectory._getFileMd5sum(file)

    @staticmethod
    def fromH5(group):
        topologyFiles = None
        inputFolder = None
        trajectoryFiles = None
        parent = ''
        if 'input_folder' in group.attrs:
            inputFolder = group.attrs['input_folder']
        if 'parent' in group.attrs:
            parent = group.attrs['parent']
        if 'trajectory_files' in group:
            trajectoryFiles = []
            for pieceName in natsorted(group['trajectory_files']):
                trajectoryFiles.append(group['trajectory_files'][pieceName].attrs['trajectory_file'])
        if 'topology_files' in group:
            topologyFiles = []
            for pieceName in natsorted(group['topology_files']):
                topologyFiles.append(group['topology_files'][pieceName].attrs['topology_file'])
        return _Trajectory(topologyFiles, trajectoryFiles, inputFolder, parent)

    def toDict(self):
        resdict = self.__dict__.copy()
        resdict['trajName'] = self.getName()
        return resdict



class TrajectoryStorage(object):
    def __init__(self, storagelocation, overwrite=False):
        self.storagelocation = storagelocation

        if os.path.splitext(storagelocation)[-1] != '.h5' and os.path.splitext(storagelocation)[-1] != '.hdf5':
            raise RuntimeError('Only storage files with .h5 or .hdf5 extensions are allowed')

        if os.path.exists(storagelocation) and not overwrite:
            logger.info(f'Loading from {storagelocation} a {self}')

        if not os.path.exists(storagelocation) or overwrite:
            logger.info(f'Creating new trajectory storage at {storagelocation}')
            with h5py.File(storagelocation, 'w') as f:
                f.create_group('trajectories')

    def __str__(self):
        ntraj = self.numTrajectories
        allproj, counts = self.projectionCounts
        nproj = len(allproj)
        if ntraj == 1:
            string = f'TrajectoryStorage with {ntraj} trajectory'
        else:
            string = f'TrajectoryStorage with {ntraj} trajectories'
        if nproj == 1:
            string += f' and 1 projection ({allproj[0]})'
        else:
            string += f' and {nproj} projections'

        return string

    def __repr__(self):
        return '<{}.{} object at {}>\n'.format(self.__class__.__module__, self.__class__.__name__, hex(id(self))) \
               + self.__str__()

    # def __getitem__(self, idx):
    #     with h5py.File(self.storagelocation, 'r') as f:
    #         names = list(f['trajectories'])
    #         return _Trajectory.fromH5(f['trajectories'][names[idx]])

    @property
    def numTrajectories(self):
        with h5py.File(self.storagelocation, 'r') as f:
            return len(list(f['trajectories']))

    @property
    def projectionCounts(self):
        allproj = []
        with h5py.File(self.storagelocation, 'r') as f:
            for trajname in f['trajectories']:
                if 'projections' in f['trajectories'][trajname]:
                    for projname in f['trajectories'][trajname]['projections']:
                        allproj.append(projname)
        vals, counts = np.unique(allproj, return_counts=True)
        return vals, counts
            

    def addTrajectory(self, trajectory, topologyFiles, inputFolder=None, sortByName=True):
        """ Add a trajectory to the storage

        Parameters
        ----------
        trajectory : file, folder or list of files
            Can be either of:
                - The path to the file of a trajectory
                - A list of paths of trajectory pieces.
                - A folder containing the trajectory
                - A folder containing a trajectory split in multiple pieces
        topologyFiles : file, or list of files
            The path to the topology file corresponding to the trajectory
        inputFolder : str
            The path to the folder which creates all files that were used to create this trajectory
        sortByName : bool
            Set to False to disable automatical sorting of trajectory pieces by their name
        """
        import uuid

        if isinstance(trajectory, str) and os.path.isdir(trajectory):  # A folder containing the trajectories
            trajpieces = _autoDetectTrajectories(trajectory)
        else:
            trajpieces = ensurelist(trajectory)

        if len(trajpieces) == 0:
            return

        if sortByName:
            trajpieces = natsorted(trajpieces)

        traj = _Trajectory(topologyFiles=ensurelist(topologyFiles), trajectoryFiles=trajpieces, inputFolder=inputFolder, parent='')

        # trajpieces = [os.path.abspath(traj) for traj in trajpieces]  # I think I prefer relative file paths

        with h5py.File(self.storagelocation, 'a') as f:
            if traj.getName() in f['trajectories']:
                raise DuplicateTrajectoryError('A trajectory is already in the dataset with files {}'.format(trajpieces))
            traj.writeToH5(f['trajectories'])


    def autoDetectFiles(self, trajectoryFolders, topologyFolders=(), inputFolders=()):
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
        trajectoryFolders : list or single file
            A list of folders containing trajectories
        topologyFolders : list or single file
            A list of folders containing the topologies of the trajectories
        inputFolders : list or single file
            A list of folders containing the input files used to create the trajectories
        """
        from htmd.util import ensurelist
        import natsort

        trajectoryFolders = ensurelist(trajectoryFolders)
        topologyFolders = ensurelist(topologyFolders)
        inputFolders = ensurelist(inputFolders)

        for folder in trajectoryFolders:
            if not (os.path.isdir(folder) or (len(trajectoryFolders) == 1 and os.path.isfile(folder))):
                raise NotADirectoryError('{}'.format(folder))
        for folder in topologyFolders:
            if not (os.path.isdir(folder) or (len(topologyFolders) == 1 and os.path.isfile(folder))):
                raise NotADirectoryError('{}'.format(folder))
        for folder in inputFolders:
            if not os.path.isdir(folder):
                raise NotADirectoryError('{}'.format(folder))

        # I need to match the simulation names inside the globs given. The
        # reason is that there can be more input folders in the glob than in
        # the data glob due to not having been retrieved. Hence I need to match
        # the folder names.

        # Create a dictionary of data folder names
        trajectorynames = dict()
        for folder in trajectoryFolders:
            name = _simName(folder)
            if name in trajectorynames:
                raise RuntimeError(f'Duplicate topology name {name} detected. Cannot name-match directories.')
            trajectorynames[name] = folder

        topologynames = dict()
        for folder in topologyFolders:
            name = _simName(folder)
            if name in topologynames:
                raise RuntimeError(f'Duplicate topology name {name} detected. Cannot name-match directories.')
            topologynames[name] = folder

        inputnames = dict()
        for folder in inputFolders:
            name = _simName(folder)
            if name in inputnames:
                raise RuntimeError(f'Duplicate topology name {name} detected. Cannot name-match directories.')
            inputnames[name] = folder

        keys = natsort.natsorted(trajectorynames.keys())
        notraj = {}  # Folders which don't have trajectories
        
        for k in tqdm(keys, desc='Detecting trajectories'):
            trajectories = _autoDetectTrajectories(trajectorynames[k])

            if not trajectories:
                notraj[trajectorynames[k]] = True
                continue

            if len(topologyFolders) == 0:  # If there are no topology folders specified
                # Try to find a topology in the trajectory folder
                topologyfile = _autoDetectTopology(trajectorynames[k])
                if inputFolders and not topologyfile:  # If there are input folders specified search there
                    topologyfile = _autoDetectTopology(inputnames[k])
            elif len(topologyFolders) == 1:
                topologyfile = topologyFolders[0]
                if os.path.isdir(topologyfile):
                    topologyfile = _autoDetectTopology(topologyfile)
            elif len(topologyFolders) > 1:
                if k not in topologynames:
                    raise FileNotFoundError(f'Did not find any topologies with matching folder name {k}.')
                topologyfile = topologynames[k]
                if os.path.isdir(topologyfile):
                    topologyfile = _autoDetectTopology(topologyfile)

            if topologyfile is None:
                raise FileNotFoundError(f'No topology file found in folder {trajectorynames[k]}. Supported extensions are {list(_TOPOLOGY_READERS.keys())}')

            inputfolder = None
            if len(inputFolders) == 1:
                inputfolder = inputFolders[0]
            elif len(inputFolders) > 1:
                if k not in inputnames:
                    raise FileNotFoundError(f'Did not find any input folder with matching name {k}.')
                inputfolder = inputnames[k]

            try:
                self.addTrajectory(trajectories, topologyfile, inputfolder)
            except DuplicateTrajectoryError as e:
                logger.info('Skipping trajectory {} since it\'s already in the dataset'.format(trajectorynames[k]))

        if len(notraj):
            logger.info('{} folder(s) were missing trajectories'.format(len(notraj)))

    def _getAllTrajectoriesAsDicts(self):
        alltraj = []
        with h5py.File(self.storagelocation, 'r') as f:
            for trajname in f['trajectories']:
                traj = _Trajectory.fromH5(f['trajectories'][trajname]).toDict()
                alltraj.append(traj)
        # from IPython.core.debugger import set_trace
        # set_trace()
        return alltraj

    def projectTrajectories(self, projection, projectionName=None, numWorkers=1, skip=1):
        from multiprocessing import Process
        from multiprocessing import JoinableQueue as Queue
        import signal
        from moleculekit.molecule import Molecule

        if projectionName is None:
            projectionName = str(type(projection))
        if len(projectionName) == 0:
            raise RuntimeError('projectionName should not be an empty string')
        logger.info('Calculating projections: {}'.format(projectionName))

        alltraj = self._getAllTrajectoriesAsDicts()

        unique_mol = None
        unique_topo = self._singleTopology()
        if unique_topo is not None:
            unique_mol = Molecule(unique_topo)

        input_queue = Queue()
        output_queue = Queue()

        for item in alltraj:
            input_queue.put(item)

        processes = []

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN) # Set SIGINT to ignore so that child processes inherit it
        processes = [Process(target=_project_worker, args=(input_queue, output_queue, projection, unique_mol, skip)) for _ in range(numWorkers)]
        for p in processes:
            p.start()
        signal.signal(signal.SIGINT, original_sigint_handler)  # Set back the original signal handler

        try:
            with h5py.File(self.storagelocation, 'a') as f:
                for i in tqdm(range(len(alltraj)), desc='Projecting trajectories'):
                    data, ref, trajname = output_queue.get()
                    output_queue.task_done()
                    if data is None:
                        continue
                    if not 'projections' in f['trajectories'][trajname]:
                        f['trajectories'][trajname].create_group('projections')
                    projgroup = f['trajectories'][trajname]['projections']
                    if projectionName in projgroup:
                        del projgroup[projectionName]
                    projgroup.create_group(projectionName)
                    projgroup[projectionName].create_dataset('data', data=data)
                    projgroup[projectionName].create_dataset('ref', data=ref)
        except KeyboardInterrupt:
            for p in processes:
                p.terminate()

        for p in processes:
            input_queue.put(None)

        for p in processes:
            p.join()


    def _singleTopology(self):
        from moleculekit.molecule import mol_equal
        from htmd.util import ensurelist

        with h5py.File(self.storagelocation, 'a') as f:
            alltopo = []
            for trajname in f['trajectories']: 
                currtrajtopos = f['trajectories'][trajname]['topology_files']
                alltopo.append(tuple([currtrajtopos[piece].attrs['topology_file'] for piece in currtrajtopos]))

        uqtopo = list(set(alltopo))

        if len(uqtopo) == 0:
            raise RuntimeError('No topologies found in TrajectoryStorage')
        elif len(uqtopo) == 1:
            return uqtopo[0]
        elif len(uqtopo) > 1:  # If more than one molfile load them and see if they are different Molecules
            ref = Molecule(uqtopo[0], _logger=False)
            for i in range(1, len(uqtopo)):
                mol = Molecule(uqtopo[i], _logger=False)
                if not mol_equal(ref, mol, exceptFields=['coords']):
                    return None
            return uqtopo[0]
        return None


def _calcRef(pieces, fileloc):
    locs = np.array(list([x[0] for x in fileloc]))
    frames = list([x[1] for x in fileloc])
    ref = np.zeros((len(frames), 2), dtype='u4')
    ref[:, 1] = frames
    for i, p in enumerate(pieces):
        ref[locs == p, 0] = i
    return ref

def _project_worker(input_queue, output_queue, projection, unique_mol, skip):
    from moleculekit.projections.projection import Projection
    from moleculekit.molecule import Molecule

    if unique_mol is not None:
        unique_mol = unique_mol.copy()

    while True:
        try:
            item = input_queue.get()
            if item is None:
                break

            mol = unique_mol
            if mol is None:
                mol = Molecule(item['topologyFiles'])

            mol.read(item['trajectoryFiles'], skip=skip)

            if isinstance(projection, Projection):
                res = projection.project(mol)
            elif hasattr(projection, '__call__'): # If it's a function
                res = projection(mol)
            elif isinstance(projection, tuple) and hasattr(projection[0], '__call__'): # If it's a function with extra arguments
                res = projection[0](mol, *projection[1])
            else:
                raise RuntimeError('Invalid projection type {}'.format(type(projection)))

            if isinstance(res, list) or isinstance(res, tuple):
                res = np.array(res)

            _checkProjectionDims(res, mol, 'projection')

            refs = _calcRef(item['trajectoryFiles'], mol.fileloc)

            output_queue.put((res, refs, item['trajName']))
            input_queue.task_done()
        except KeyboardInterrupt:
            break
        except Exception as e:
            logger.error('Failed to project simulation with name {} due to {}'.format(item['trajName'], e))
            input_queue.task_done()
            output_queue.put((None, None, item['trajName']))

def _checkProjectionDims(result, mol, name):
    if (result.ndim == 1 and len(result) != mol.numFrames) or \
       (result.ndim == 2 and result.shape[0] != mol.numFrames):
        raise RuntimeError('The {} produced an incorrect result. It produced a {} shaped array for {} frames. '
                           'The {} should return a numpy array of (nframes, ndim) shape '
                           'where nframes the number of frames in the Molecule it accepts as an '
                           'argument.'.format(name, result.shape, mol.numFrames, name))


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
        return None
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


import unittest
class _TestTrajectoryStorage(unittest.TestCase):
    def test_autoDetectFiles(self):
        import os
        from htmd.home import home
        from htmd.util import tempname

        def compareDicts(dict1, dict2):
            import re
            regex = re.compile('traj(\d)')
            rekeyeddict1 = {int(regex.search(tt['trajectoryFiles'][0])[1]): tt for tt in dict1}
            rekeyeddict2 = {int(regex.search(tt['trajectoryFiles'][0])[1]): tt for tt in dict2}
            for trajname in rekeyeddict1:
                traj1 = rekeyeddict1[trajname]
                traj2 = rekeyeddict2[trajname]
                for key in traj1:
                    for data1, data2 in zip(ensurelist(traj1[key]), ensurelist(traj2[key])):
                        if isinstance(data2, str):
                            assert data2.endswith(data1), f'Mismatch in key: {key}\n\ndict1:\n{rekeyeddict1}\n\ndict2:\n{rekeyeddict2}'
                        else:
                            assert data2 == data1, f'Mismatch in key: {key}\n\ndict1:\n{rekeyeddict1}\n\ndict2:\n{rekeyeddict2}'

        mode1dict = [{'topologyFiles': ['mode1/trajectories/traj1/structure.prmtop'],
                    'trajectoryFiles': ['mode1/trajectories/traj1/output.xtc'],
                    'inputFolder': None,
                    'parent': ''},
                    {'topologyFiles': ['mode1/trajectories/traj2/structure.prmtop'],
                    'trajectoryFiles': ['mode1/trajectories/traj2/output.xtc'],
                    'inputFolder': None,
                    'parent': ''}]

        tmph5 = tempname(suffix='.h5')
        ts = TrajectoryStorage(tmph5)
        ts.autoDetectFiles(natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode1', 'trajectories', '*', ''))))
        alltraj = ts._getAllTrajectoriesAsDicts()
        del ts
        os.remove(tmph5)
        compareDicts(mode1dict, alltraj)


        mode2dict = [{'topologyFiles': ['mode2/topologies/traj1/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj1/output.xtc'],
                    'inputFolder': None,
                    'parent': ''},
                    {'topologyFiles': ['mode2/topologies/traj2/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj2/output.xtc'],
                    'inputFolder': None,
                    'parent': ''}]

        tmph5 = tempname(suffix='.h5')
        ts = TrajectoryStorage(tmph5)
        ts.autoDetectFiles(natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'trajectories', '*', ''))),
                            topologyFolders=natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'topologies', '*', ''))))
        alltraj = ts._getAllTrajectoriesAsDicts()
        del ts
        os.remove(tmph5)
        compareDicts(mode2dict, alltraj)


        mode3dict = [{'topologyFiles': ['mode2/inputs/traj1/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj1/output.xtc'],
                    'inputFolder': 'mode2/inputs/traj1/',
                    'parent': ''},
                    {'topologyFiles': ['mode2/inputs/traj2/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj2/output.xtc'],
                    'inputFolder': 'mode2/inputs/traj2/',
                    'parent': ''}]

        tmph5 = tempname(suffix='.h5')
        ts = TrajectoryStorage(tmph5)
        ts.autoDetectFiles(natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'trajectories', '*', ''))),
                            inputFolders=natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'inputs', '*', ''))))
        alltraj = ts._getAllTrajectoriesAsDicts()
        del ts
        os.remove(tmph5)
        compareDicts(mode3dict, alltraj)

        mode4dict = [{'topologyFiles': ['mode2/topologies/traj1/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj1/output.xtc'],
                    'inputFolder': None,
                    'parent': ''},
                    {'topologyFiles': ['mode2/topologies/traj1/structure.prmtop'],
                    'trajectoryFiles': ['mode2/trajectories/traj2/output.xtc'],
                    'inputFolder': None,
                    'parent': ''}]

        tmph5 = tempname(suffix='.h5')
        ts = TrajectoryStorage(tmph5)
        ts.autoDetectFiles(natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'trajectories', '*', ''))),
                            topologyFolders=os.path.join(home(dataDir='test-trajectory-storage'), 'mode2', 'topologies', 'traj1', 'structure.prmtop'))
        alltraj = ts._getAllTrajectoriesAsDicts()
        del ts
        os.remove(tmph5)
        compareDicts(mode4dict, alltraj)

    def test_projection(self):
        import os
        from htmd.home import home
        from htmd.util import tempname
        from moleculekit.projections.metricdistance import MetricDistance
        
        tmph5 = tempname(suffix='.h5')
        ts = TrajectoryStorage(tmph5)
        ts.autoDetectFiles(natsorted(glob(os.path.join(home(dataDir='test-trajectory-storage'), 'mode1', 'trajectories', '*', ''))))

        proj = MetricDistance('resname BEN and name C4', 'protein and resid 5 and name HE1')
        ts.projectTrajectories(proj, projectionName='C4-HE1 distance')
        projs, counts = ts.projectionCounts
        assert len(projs) == 1
        assert projs[0] == 'C4-HE1 distance'
        assert counts[0] == 2

        proj = MetricDistance('resname BEN and name C4', 'resname BEN and name C7')
        ts.projectTrajectories(proj, projectionName='C4-C7 distance')
        projs, counts = ts.projectionCounts
        assert len(projs) == 2
        assert projs[0] == 'C4-C7 distance'
        assert counts[0] == 2
        assert projs[1] == 'C4-HE1 distance'
        assert counts[1] == 2


if __name__ == '__main__':
    unittest.main(verbosity=2)
