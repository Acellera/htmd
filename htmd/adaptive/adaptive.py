# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import re
import abc
import time
import os.path as path
from glob import glob
from os import makedirs
from shutil import copytree, ignore_patterns
import numpy as np
from natsort import natsorted
from htmd.simlist import _simName
from htmd.molecule.molecule import Molecule
from htmd.userinterface import uisetattr
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS
import logging
logger = logging.getLogger(__name__)


class AdaptiveNew(ProtocolInterface):
    def __init__(self):
        super().__init__()
        from htmd.apps.app import App
        self._cmdObject('app', ':class:`App <htmd.apps.app.App>` object', 'An App class object used to retrieve and submit simulations', None, App)
        self._cmdString('project', 'str', 'The name of the project', 'adaptive')
        self._cmdValue('nmin', 'int', 'Minimum number of running simulations', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('nmax', 'int', 'Maximum number of running simulations', 1, TYPE_INT, RANGE_POS)
        self._cmdValue('nepochs', 'int', 'Maximum number of epochs', 100, TYPE_INT, RANGE_POS)
        self._cmdString('inputpath', 'str', 'The directory used to store input folders', 'input')
        self._cmdString('generatorspath', 'str', 'The directory containing the generators', 'generators')
        self._cmdBoolean('dryrun', 'boolean', 'A dry run means that the adaptive will retrieve and generate a new epoch but not submit the simulations', False)
        self._cmdValue('updateperiod', 'float', 'When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds', 0, TYPE_FLOAT, RANGE_0POS)
        self._running = None

    def run(self):
        """ Runs the adaptive

        Use this command to start the adaptive.

        Example
        -------
        >>> adapt = Adaptive()
        >>> adapt.run()
        """
        while True:
            epoch = self._getEpoch()
            logger.info('Processing epoch ' + str(epoch))

            if epoch == 0 and self.generatorspath:
                logger.info('Epoch 0, generating first batch')
                self._init()
                if not self.dryrun:
                    self.app.submit(natsorted(glob(path.join(self.inputpath, 'e1s*'))))
            else:
                logger.info('Retrieving simulations.')
                self.app.retrieve()

                if epoch >= self.nepochs:
                    logger.info('Reached maximum number of epochs ' + str(self.nepochs))
                    return

                self._running = self.app.inprogress()
                logger.info(str(self._running) + ' simulations in progress')

                if self._running <= self.nmin:
                    self._algorithm()
                    if not self.dryrun:
                        self.app.submit(natsorted(glob(path.join(self.inputpath, 'e' + str(epoch+1) + 's*'))))
                        logger.info('Finished submitting simulations.')

            if self.updateperiod <= 0:
                break
            logger.info('Sleeping for {} seconds.'.format(self.updateperiod))
            time.sleep(self.updateperiod)

    def _init(self):
        folders = natsorted(glob(path.join(self.generatorspath, '*', ''))) # I need the extra ''  to add a finishing /
        if len(folders) == 0:
            logger.info('Generators folder has no subdirectories, using folder itself')
            folders.append(self.generatorspath)

        numF = len(folders)
        numCopies = np.ones(numF, dtype=int) * int(np.floor(self.nmax / numF))
        numExtra = np.mod(self.nmax, numF)
        numCopies = numCopies + np.random.multinomial(numExtra, [1/numF]*numF)  # draw the extra equally from a flat distribution
        if not path.exists(self.inputpath):
            makedirs(self.inputpath)

        # Check if epoch 1 directories already exist in the input folder
        existing = glob(path.join(self.inputpath, 'e1s*'))
        if len(existing) != 0:
            raise NameError('Epoch 1 directories already exist.')

        k = 1
        for i in range(numF):
            for j in range(numCopies[i]):
                name = _simName(folders[i])
                inputdir = path.join(self.inputpath, 'e1s' + str(k) + '_' + name)
                #src = path.join(self.generatorspath, name, '*')
                src = folders[i]
                copytree(src, inputdir, symlinks=True, ignore=ignore_patterns('*.dcd', '*.xtc'))
                k += 1

    def _getEpoch(self):
        """ Compute current epoch of adaptive

        Checks the input folder for the latest epoch inputs and returns the epoch number

        Returns
        -------
        epoch : int
            The current epoch
        """
        folders = glob(path.join(self.inputpath, 'e*', ''))
        epoch = 0
        regex = re.compile('e(\d+)')
        for f in folders:
            res = regex.search(f)
            if res:
                num = res.group(1)
                if int(num) > epoch:
                    epoch = int(num)
        return epoch

    def _writeInputs(self, simsframes, epoch=None):
        if epoch is None:
            epoch = self._getEpoch() + 1

        test = glob(path.join(self.inputpath, 'e' + str(epoch) + '*'))
        if len(test) != 0:
            raise NameError('Input dirs of epoch ' + str(epoch) + ' already exists.')

        if path.exists(path.join(self.inputpath, 'e' + str(epoch) + '_writeinputs.log')):
            raise NameError('Epoch logfile already exists. Cant overwrite it.')

        fid = open(path.join(self.inputpath, 'e' + str(epoch) + '_writeinputs.log'), 'w')

        regex = re.compile('(e\d+s\d+)_')
        for i, f in enumerate(simsframes):
            frameNum = f.frame
            piece = f.piece
            #print(frameNum)
            if f.sim.parent is None:
                currSim = f.sim
            else:
                currSim = f.sim.parent

            traj = currSim.trajectory[piece]
            if currSim.input is None:
                raise NameError('Could not find input folder in simulation lists. Cannot create new simulations.')

            wuName = _simName(traj)
            res = regex.search(wuName)
            if res:  # If we are running on top of adaptive, use the first name part for the next sim name
                wuName = res.group(1)

            # create new job directory
            newName = 'e' + str(epoch) + 's' + str(i+1) + '_' + wuName + 'p' + str(piece) + 'f' + str(frameNum)
            newDir = path.join(self.inputpath, newName, '')
            # copy previous input directory including input files
            copytree(currSim.input, newDir, symlinks=False, ignore=ignore_patterns('*.dcd', '*.xtc', '*.coor'))
            # overwrite input file with new one. frameNum + 1 as catdcd does 1 based indexing
            mol = Molecule()
            mol.read(traj)
            mol.frame = frameNum
            mol.write(path.join(newDir, 'input.coor'))

            # write nextInput file
            fid.write('# {0} \n{1} {2}\n'.format(newName, traj, frameNum))

        fid.close()

    @abc.abstractmethod
    def _algorithm(self):
        return


class Adaptive(object):

    _adaptivecmds = ['app', 'project', 'nmin', 'nmax', 'nepochs', 'inputpath', 'generatorspath',
            'dryrun','updateperiod','running']

    def __init__(self, app=None, project=None, nmin=1, nmax=1, nepochs=1, inputpath='input',
                 generatorspath='generators', dryrun=False, updateperiod=0):
        """ Constructor for the Adaptive class

        Parameters
        ----------
        app : :class:`App <htmd.apps.app.App>` object
            An App class object used to retrieve and submit simulations
        project : str
            The name of the project
        nmin : int
            Minimum number of running simulations
        nmax : int
            Maximum number of running simulations
        nepochs : int
            Maximum number of epochs
        inputpath : str, optional, default='input'
            The directory used to store input folders
        generatorspath : str, optional, default='generators'
            The directory containing the generators
        dryrun : bool, optional, default=False
            A dry run means that the adaptive will retrieve and generate a new epoch but no submit the simulations
        updateperiod : float, optional, default=0
            When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds
        """
        self.app = app
        self.project = project
        self.nmin = nmin
        self.nmax = nmax
        self.nepochs = nepochs
        self.inputpath = inputpath
        self.generatorspath = generatorspath
        self.dryrun = dryrun
        self.updateperiod = updateperiod
        self.running = None

    def __setattr__(self, key, value):
        uisetattr(self, key, value, self._adaptivecmds)

    def run(self):
        """ Runs the adaptive

        Use this command to start the adaptive.

        Example
        -------
        >>> adapt = Adaptive()
        >>> adapt.run()
        """
        while True:
            epoch = self._getEpoch()
            logger.info('Processing epoch ' + str(epoch))

            if epoch == 0 and self.generatorspath:
                logger.info('Epoch 0, generating first batch')
                self._init()
                if not self.dryrun:
                    self.app.submit(natsorted(glob(path.join(self.inputpath, 'e1s*'))))
            else:
                logger.info('Retrieving simulations.')
                self.app.retrieve()

                if epoch >= self.nepochs:
                    logger.info('Reached maximum number of epochs ' + str(self.nepochs))
                    return

                self.running = self.app.inprogress()
                logger.info(str(self.running) + ' simulations in progress')

                if self.running <= self.nmin:
                    self._algorithm()
                    if not self.dryrun:
                        self.app.submit(natsorted(glob(path.join(self.inputpath, 'e' + str(epoch+1) + 's*'))))
                        logger.info('Finished submitting simulations.')

            if self.updateperiod <= 0:
                break
            logger.info('Sleeping for {} seconds.'.format(self.updateperiod))
            time.sleep(self.updateperiod)

    def _init(self):
        folders = natsorted(glob(path.join(self.generatorspath, '*', ''))) # I need the extra ''  to add a finishing /
        if len(folders) == 0:
            logger.info('Generators folder has no subdirectories, using folder itself')
            folders.append(self.generatorspath)

        numF = len(folders)
        numCopies = np.ones(numF, dtype=int) * int(np.floor(self.nmax / numF))
        numExtra = np.mod(self.nmax, numF)
        numCopies = numCopies + np.random.multinomial(numExtra, [1/numF]*numF)  # draw the extra equally from a flat distribution
        if not path.exists(self.inputpath):
            makedirs(self.inputpath)

        # Check if epoch 1 directories already exist in the input folder
        existing = glob(path.join(self.inputpath, 'e1s*'))
        if len(existing) != 0:
            raise NameError('Epoch 1 directories already exist.')

        k = 1
        for i in range(numF):
            for j in range(numCopies[i]):
                name = _simName(folders[i])
                inputdir = path.join(self.inputpath, 'e1s' + str(k) + '_' + name)
                #src = path.join(self.generatorspath, name, '*')
                src = folders[i]
                copytree(src, inputdir, symlinks=True, ignore=ignore_patterns('*.dcd', '*.xtc'))
                k += 1

    def _getEpoch(self):
        """ Compute current epoch of adaptive

        Checks the input folder for the latest epoch inputs and returns the epoch number

        Returns
        -------
        epoch : int
            The current epoch
        """
        folders = glob(path.join(self.inputpath, 'e*', ''))
        epoch = 0
        regex = re.compile('e(\d+)')
        for f in folders:
            res = regex.search(f)
            if res:
                num = res.group(1)
                if int(num) > epoch:
                    epoch = int(num)
        return epoch

    def _writeInputs(self, simsframes, epoch=None):
        if epoch is None:
            epoch = self._getEpoch() + 1

        test = glob(path.join(self.inputpath, 'e' + str(epoch) + '*'))
        if len(test) != 0:
            raise NameError('Input dirs of epoch ' + str(epoch) + ' already exists.')

        if path.exists(path.join(self.inputpath, 'e' + str(epoch) + '_writeinputs.log')):
            raise NameError('Epoch logfile already exists. Cant overwrite it.')

        fid = open(path.join(self.inputpath, 'e' + str(epoch) + '_writeinputs.log'), 'w')

        regex = re.compile('(e\d+s\d+)_')
        for i, f in enumerate(simsframes):
            frameNum = f.frame
            piece = f.piece
            #print(frameNum)
            if f.sim.parent is None:
                currSim = f.sim
            else:
                currSim = f.sim.parent

            traj = currSim.trajectory[piece]
            if currSim.input is None:
                raise NameError('Could not find input folder in simulation lists. Cannot create new simulations.')

            wuName = _simName(traj)
            res = regex.search(wuName)
            if res:  # If we are running on top of adaptive, use the first name part for the next sim name
                wuName = res.group(1)

            # create new job directory
            newName = 'e' + str(epoch) + 's' + str(i+1) + '_' + wuName + 'p' + str(piece) + 'f' + str(frameNum)
            newDir = path.join(self.inputpath, newName, '')
            # copy previous input directory including input files
            copytree(currSim.input, newDir, symlinks=False, ignore=ignore_patterns('*.dcd', '*.xtc', '*.coor'))
            # overwrite input file with new one. frameNum + 1 as catdcd does 1 based indexing
            mol = Molecule()
            mol.read(traj)
            mol.frame = frameNum
            mol.write(path.join(newDir, 'input.coor'))

            # write nextInput file
            fid.write('# {0} \n{1} {2}\n'.format(newName, traj, frameNum))

        fid.close()

    @abc.abstractmethod
    def _algorithm(self):
        return


def reconstructAdaptiveTraj(simlist, trajID):
    """ Reconstructs a long trajectory out of short adaptive runs.

    Parameters
    ----------
    simlist : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects
        A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function
    trajID : int
        The id of the trajectory from which to start going back.

    Returns
    -------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A Molecule object containing the reconstructed trajectory
    chain : np.ndarray
        The simulation IDs of all simulations involved
    pathlist : np.ndarray of str
        The names of all simulations involved.

    Examples
    --------
    >>> mol, chain, pathlist = reconstructAdaptiveTraj(data.simlist, 52)
    """
    from IPython.core.debugger import Tracer

    sim = None
    for s in simlist:
        if s.simid == trajID:
            sim = s
            break
    if sim is None:
        raise NameError('Could not find sim with ID {} in the simlist.'.format(trajID))

    pathlist = []
    pathlist.append(sim.trajectory[0])
    chain = []
    chain.append((sim, -1, -1))

    epo = None
    while epo != 1:
        [sim, piece, frame, epo] = _findprevioustraj(simlist, _simName(sim.trajectory[0]))
        pathlist.append(sim.trajectory[piece])
        chain.append((sim, piece, frame))
    pathlist = pathlist[::-1]
    chain = chain[::-1]

    mol = Molecule(sim.molfile)
    mol.coords = np.zeros((mol.numAtoms, 3, 0), dtype=np.float32)
    mol.fileloc = []
    mol.box = np.zeros((3, 0))
    for i, c in enumerate(chain):
        tmpmol = Molecule(sim.molfile)
        tmpmol.read(c[0].trajectory)
        endpiece = c[1]
        fileloc = np.vstack(tmpmol.fileloc)
        filenames = fileloc[:, 0]
        pieces = np.unique(filenames)
        firstpieceframe = np.where(filenames == pieces[endpiece])[0][0]
        endFrame = firstpieceframe + c[2]
        if endFrame != -1:
            tmpmol.coords = tmpmol.coords[:, :, 0:endFrame + 1]  # Adding the actual respawned frame (+1) since the respawned sim doesn't include it in the xtc
            tmpmol.fileloc = tmpmol.fileloc[0:endFrame + 1]
            tmpmol.box = tmpmol.box[:, 0:endFrame + 1]
        mol.coords = np.concatenate((mol.coords, tmpmol.coords), axis=2)
        mol.box = np.concatenate((mol.box, tmpmol.box), axis=1)
        mol.fileloc += tmpmol.fileloc
        #Tracer()()
    #mol.fileloc[:, 1] = range(np.size(mol.fileloc, 0))

    return mol, chain, pathlist


def _findprevioustraj(simlist, simname):
    regex = re.compile('_(e\d+s\d+)p(\d+)f(\d+)$')
    regex2 = re.compile('_(e\d+s\d+)f(\d+)$')
    m = regex.search(simname)
    if m:
        prevname = m.group(1)
        prevpiece = int(m.group(2)) - 1
        prevframe = int(m.group(3))
    else:
        m2 = regex2.search(simname)
        if m2:
            prevname = m2.group(1)
            prevpiece = 0
            prevframe = int(m2.group(2))
        else:
            raise NameError('Could not match simname: {} with regular expressions.'.format(simname))
    regex = re.compile('e(\d+)s')
    m = regex.match(prevname)
    if not m:
        raise NameError('Could not parse epoch number from name: {}.'.format(prevname))
    epo = int(m.group(1))
    sim = None
    regex = re.compile('{}_'.format(prevname))
    for s in simlist:
        if len(s.trajectory) <= prevpiece:
            continue
        m = regex.search(s.trajectory[prevpiece])
        if m:
            sim = s
            break
    if sim is None:
        raise NameError('Could not find parent of simulation {}.'.format(simname))
    return sim, prevpiece, prevframe, epo

class _AdaptiveTest(Adaptive):
    _cmds = ['datapath']

    def __init__(self, app=None, project=None, nmin=1, nmax=1, nepochs=1, inputpath='input', datapath='data',
                 generatorspath='generators', dryrun=False, updateperiod=0):
        super().__init__(app, project, nmin, nmax, nepochs, inputpath, generatorspath, dryrun, updateperiod)
        self.datapath = datapath

    def __setattr__(self, key, value):
        all = self._adaptivecmds + self._cmds
        uisetattr(self, key, value, all)

    def _algorithm(self):
        """  Select random frames for respawning
        """
        from htmd.projections.metric import Metric
        from htmd.molecule.molecule import Molecule
        from htmd.projections.metriccoordinate import MetricCoordinate
        from htmd.simlist import simlist
        sims = simlist(glob(path.join(self.datapath, '*', '')), glob(path.join(self.inputpath, '*', 'structure.pdb')),
                       glob(path.join(self.inputpath, '*', '')))
        metr = Metric(sims)
        metr.projection(MetricCoordinate(Molecule(sims[0].molfile), 'protein and name CA', 'protein and name CA'))
        data = metr.project()
        simframes = data.abs2sim(np.random.randint(0, data.numFrames, self.nmax-self.running))
        self._writeInputs(simframes)


if __name__ == "__main__":
    """from htmd.apps.acemdlocal import AcemdLocal
    from htmd.adaptive.adaptiverun import AdaptiveRun
    import tempfile
    from htmd.home import home
    import os.path as path
    import shutil
    import logging
    import os
    from glob import glob
    log = logging.getLogger('htmd')
    log.setLevel('DEBUG')

    # Testing epoch 0
    tmpfolder = tempfile.mkdtemp()
    print(tmpfolder)
    md = _AdaptiveTest()
    md.nmin = 0
    md.nmax = 4
    md.nepochs = 2
    md.inputpath = path.join(tmpfolder, 'input')
    md.datapath = path.join(tmpfolder, 'data')
    md.generatorspath = path.join(home(), 'data', 'adaptive', 'generators')
    md.app = AcemdLocal(ngpus=4, datadir=md.datapath, acemd='/shared/acemd/bin/acemd')
    md.dryrun = True
    md.run()
    # Testing epoch 1
    md.inputpath = path.join(home(), 'data', 'adaptive', 'input')
    md.datapath = path.join(home(), 'data', 'adaptive', 'data')
    md.run()
    # Cleaning up
    inputodel = glob(path.join(home(), 'data', 'adaptive', 'input', 'e2*'))
    for i in inputodel:
        shutil.rmtree(i, ignore_errors=True)
    os.remove(path.join(home(), 'data', 'adaptive', 'input', 'e2_writeinputs.log'))""" # No acemd on Travis :(

    '''########## Testing AdaptiveRun. Doesn't work because I don't have enough frames to pass mergeSmall...
    # Testing epoch 0
    tmpfolder = tempfile.mkdtemp()
    print(tmpfolder)
    md = AdaptiveRun()
    md.nmin = 0
    md.nmax = 4
    md.nepochs = 2
    md.metrictype = 'contacts'
    md.datapath = path.join(tmpfolder, 'data')
    md.filteredpath = path.join(tmpfolder, 'filtered')
    md.inputpath = path.join(tmpfolder, 'input')
    md.resultspath = path.join(tmpfolder, 'results')
    md.generatorspath = path.join(home(), 'data', 'adaptive', 'generators')
    md.macronum = 4
    md.ticadim = 0
    md.metricsel1 = 'protein and name CA'
    md.metricsel2 = 'resname BEN and noh'
    md.app = AcemdLocal(ngpus=4, datadir=md.datapath)
    md.dryrun = True
    md.run()
    log.setLevel('WARNING')
    # Testing epoch 1
    md.inputpath = path.join(home(), 'data', 'adaptive', 'input')
    md.datapath = path.join(home(), 'data', 'adaptive', 'data')
    md.run()
    # Cleaning up
    inputodel = glob(path.join(home(), 'data', 'adaptive', 'input', 'e2*'))
    for i in inputodel:
        shutil.rmtree(i, ignore_errors=True, acemd='/shared/acemd/bin/acemd')
    os.remove(path.join(home(), 'data', 'adaptive', 'input', 'e2_writeinputs.log'))'''