# (c) 2015-2018 Acellera Ltd http://www.acellera.com
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
import shutil
from shutil import copytree, ignore_patterns
import numpy as np
from jobqueues.simqueue import RetrieveError, InProgressError, ProjectNotExistError
from htmd.simlist import _simName
from moleculekit.molecule import Molecule
from protocolinterface import ProtocolInterface, val
import logging
logger = logging.getLogger(__name__)

_IGNORE_EXTENSIONS = ('*.dcd', '*.xtc', '*.binpos', '*.trr', '*.nc', '*.h5', '*.lh5', '*.netcdf', '*.vel', '.done', '*.chk', '*restart*')


class AdaptiveBase(abc.ABC, ProtocolInterface):
    def __init__(self):
        super().__init__()
        from jobqueues.simqueue import SimQueue
        self._arg('app', ':class:`SimQueue <jobqueues.simqueue.SimQueue>` object', 'A SimQueue class object used to retrieve and submit simulations', None, val.Object((SimQueue,)))
        self._arg('project', 'str', 'The name of the project', 'adaptive', val.String())
        self._arg('nmin', 'int', 'Minimum number of running simulations', 0, val.Number(int, '0POS'))
        self._arg('nmax', 'int', 'Maximum number of running simulations', 1, val.Number(int, 'POS'))
        self._arg('nepochs', 'int', 'Stop adaptive once we have reached this number of epochs', 1000, val.Number(int, 'POS'))
        self._arg('nframes', 'int', 'Stop adaptive once we have simulated this number of aggregate simulation frames.', 0, val.Number(int, '0POS'))
        self._arg('inputpath', 'str', 'The directory used to store input folders', 'input', val.String())
        self._arg('generatorspath', 'str', 'The directory containing the generators', 'generators', val.String())
        self._arg('dryrun', 'boolean', 'A dry run means that the adaptive will retrieve and generate a new epoch but not submit the simulations', False, val.Boolean())
        self._arg('updateperiod', 'float', 'When set to a value other than 0, the adaptive will run synchronously every `updateperiod` seconds', 0, val.Number(float, '0POS'))
        self._arg('coorname', 'str', 'Name of the file containing the starting coordinates for the new simulations', 'input.coor', val.String())
        self._arg('lock', 'bool', 'Lock the folder while adaptive is ongoing', False, val.Boolean())
        self._running = None

    def run(self):
        """ Runs the adaptive

        Use this command to start the adaptive.

        Example
        -------
        >>> adapt = Adaptive()
        >>> adapt.run()
        """
        from natsort import natsorted
        if self.nmax <= self.nmin:
            raise RuntimeError('nmax option should be larger than nmin.')

        self._setLock()
        while True:
            epoch = self._getEpoch()
            logger.info('Processing epoch ' + str(epoch))

            if epoch == 0 and self.generatorspath:
                logger.info('Epoch 0, generating first batch')
                self._init()
                if not self.dryrun:
                    self.app.submit(natsorted(glob(path.join(self.inputpath, 'e1s*'))))
            else:
                # Retrieving simulations
                logger.info('Retrieving simulations.')
                try:
                    self.app.retrieve()
                except RetrieveError as e:
                    logger.error('Quitting adaptive run due to error in retrieving simulations: {}'.format(e))
                    return
                except ProjectNotExistError:
                    logger.info('Retrieve found no previous simulations for this adaptive. Assuming this is a new adaptive run')

                # Checking how many simulations are in progress (queued/running) on the queue
                try:
                    self._running = self.app.inprogress()
                except InProgressError as e:
                    logger.error('Quitting adaptive run due to error in checking number of simulations in progress: {}'.format(e))
                    return
                except ProjectNotExistError:
                    logger.info('Inprogress found no previous simulations for this adaptive. Assuming this is a new adaptive run')
                    self._running = 0

                logger.info(str(self._running) + ' simulations in progress')

                if epoch >= self.nepochs and self._running == 0:
                    logger.info('Reached maximum number of epochs ' + str(self.nepochs))
                    self._unsetLock()
                    return

                # If currently running simulations are lower than nmin start new ones to reach nmax number of sims
                if self._running <= self.nmin and epoch < self.nepochs:
                    flag = self._algorithm()
                    if flag is False:
                        self._unsetLock()
                        return

                    if not self.dryrun:
                        newsims = glob(path.join(self.inputpath, 'e' + str(epoch+1) + 's*'))
                        try:
                            self.app.submit(natsorted(newsims))
                        except:
                            # If submitting fails delete all simulation inputs to not confuse _getEpoch()
                            for ns in newsims:
                                shutil.rmtree(ns)
                        logger.info('Finished submitting simulations.')

            if self.updateperiod <= 0:
                break
            logger.info('Sleeping for {} seconds.'.format(self.updateperiod))
            time.sleep(self.updateperiod)
        self._unsetLock()

    def _setLock(self):
        import datetime

        if self.lock:
            lockfile = os.path.abspath('./adaptivelock')
            if os.path.exists(lockfile):
                raise FileExistsError('This adaptive folder is locked by a running adaptive application. If this is not'
                                      ' the case, delete the {} file and run adaptive again.'.format(lockfile))

            with open(lockfile, 'w') as f:
                f.write('{}'.format(datetime.datetime.now()))

    def _unsetLock(self):
        if self.lock:
            lockfile = os.path.abspath('./adaptivelock')
            if os.path.exists(lockfile):
                os.remove(lockfile)

    def _init(self):
        from natsort import natsorted
        folders = natsorted(glob(path.join(self.generatorspath, '*', ''))) # I need the extra ''  to add a finishing /
        if len(folders) == 0:
            logger.info('Generators folder has no subdirectories, using folder itself')
            folders.append(self.generatorspath)

        numF = len(folders)
        numCopies = np.ones(numF, dtype=int) * int(np.floor(self.nmax / numF))
        numExtra = np.mod(self.nmax, numF)
        extraChoices = np.random.choice(numF, numExtra, replace=False) # draw the extra
        numCopies[extraChoices] += 1
        # numCopies = numCopies + np.random.multinomial(numExtra, [1/numF]*numF)  # draw the extra equally from a flat distribution
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
                copytree(src, inputdir, symlinks=True, ignore=ignore_patterns(*_IGNORE_EXTENSIONS))
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

        from htmd.parallelprogress import ParallelExecutor
        from htmd.config import _config
        from joblib import delayed

        aprun = ParallelExecutor(n_jobs=_config['njobs'])
        aprun(total=len(simsframes), desc='Writing inputs')(
            delayed(_writeInputsFunction)(i, f, epoch, self.inputpath, self.coorname) for i, f in enumerate(simsframes))

    @abc.abstractmethod
    def _algorithm(self):
        return


def _writeInputsFunction(i, f, epoch, inputpath, coorname):
    regex = re.compile('(e\d+s\d+)_')
    frameNum = f.frame
    piece = f.piece
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
    newName = 'e' + str(epoch) + 's' + str(i + 1) + '_' + wuName + 'p' + str(piece) + 'f' + str(frameNum)
    newDir = path.join(inputpath, newName, '')

    # copy previous input directory including input files
    copytree(currSim.input, newDir, symlinks=False, ignore=ignore_patterns('*.coor', '*.rst', '*.out', *_IGNORE_EXTENSIONS))

    # overwrite input file with new one. frameNum + 1 as catdcd does 1 based indexing

    mol = Molecule(currSim.molfile)  # Always read the mol file, otherwise it does not work if we need to save a PDB as coorname
    mol.read(traj)
    mol.dropFrames(keep=frameNum)  # Making sure only specific frame to write is kept
    mol.write(path.join(newDir, coorname))


def epochSimIndexes(simlist):
    """ Finds the simulation indexes for each epoch.

    Creates a dictionary with the epoch number as key and values the simlist indexes of the simulations corresponding to
    the given epoch.

    Parameters
    ----------
    simlist : list
        A simulation list created using the :func:`simlist <htmd.simlist.simlist>` function

    Returns
    -------

    """
    epochidx = {}
    for i, sl in enumerate(simlist):
        simepoch = getEpochFromName(sl.trajectory[0])
        if simepoch not in epochidx:
            epochidx[simepoch] = []
        epochidx[simepoch].append(i)
    return epochidx


def getEpochFromName(name):
    """ Given a adaptive simulation name, tells you which epoch it belongs to.

    Parameters
    ----------
    name : str
        Simulation name

    Returns
    -------
    epoch : int
        The epoch
    """
    import re
    reg = re.compile('/e(\d+)s\d+_')
    matches = reg.findall(name)
    if len(matches) == 0:
        raise RuntimeError('{} is not an adaptive trajectory'.format(name))
    return int(matches[0])


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
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        A Molecule object containing the reconstructed trajectory
    chain : np.ndarray
        The simulation IDs of all simulations involved
    pathlist : np.ndarray of str
        The names of all simulations involved.

    Examples
    --------
    >>> mol, chain, pathlist = reconstructAdaptiveTraj(data.simlist, 52)
    """

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



if __name__ == "__main__":
    import htmd
    import os
    from htmd.simlist import Frame, simlist
    from htmd.util import tempname

    filedir = htmd.home.home()+'/data/adaptive/'
    sims = simlist(glob(os.path.join(filedir, 'data', '*', '')),
                   glob(os.path.join(filedir, 'input', '*', '')),
                   glob(os.path.join(filedir, 'input', '*', '')))

    outf = tempname()
    os.makedirs(outf)

    f = Frame(sims[0], 0, 5)
    _writeInputsFunction(1, f, 2, outf, 'input.coor')

    mol = Molecule(sims[0])
    mol.read(os.path.join(outf, 'e2s2_e1s1p0f5', 'input.coor'))

    shutil.rmtree(outf)

