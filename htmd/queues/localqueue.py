# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from htmd.queues.simqueue import SimQueue
from protocolinterface import ProtocolInterface, val
import queue
import os
import threading
from subprocess import check_output
from glob import glob as glob
from abc import abstractmethod
import logging
import psutil
from pint import UnitRegistry
logger = logging.getLogger(__name__)

# TODO: Merge CPU and GPU queue into a single one which manages ncpu and ngpu simultaneously


class _LocalQueue(SimQueue, ProtocolInterface):
    """ Support class for local machine queue systems

    Parameters
    ----------
    datadir : str, default=None
        The path in which to store completed trajectories.
    trajext : str, default='xtc'
        Extension of trajectory files. This is needed to copy them to datadir.
    copy : list, default='*.xtc'
        A list of file names or globs for the files to copy to datadir

    """

    def __init__(self):
        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)
        self._arg('datadir', 'str', 'The path in which to store completed trajectories.', None, val.String())
        self._arg('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc',
                  val.String())
        self._arg('copy', 'list', 'A list of file names or globs for the files to copy to datadir', ('*.xtc', ),
                  val.String(), nargs='*')
        self._cmdDeprecated('trajext', 'copy')

        self._states = dict()
        self._queue = None
        self._shutdown = False

    def _setupQueue(self):
        if self._queue is None:
            self._queue = queue.Queue()

            devices = self._getdevices()
            self.memory = self._getmemory()

            self._threads = []
            for d in devices:
                t = threading.Thread(target=self.run_job, args=(d, ))
                t.daemon = True
                t.start()
                self._threads.append(t)

    def run_job(self, gpuid):
        queue = self._queue
        while not self._shutdown:
            path = None
            try:
                path = queue.get(timeout=1)
            except:
                pass

            if path:
                if gpuid is None:
                    logger.info('Running ' + path)
                else:
                    logger.info("Running " + path + " on GPU device " + str(gpuid))
                self._setRunning(path)

                runsh = os.path.join(path, 'run.sh')
                jobsh = os.path.join(path, 'job.sh')
                self._createJobScript(jobsh, path, runsh, gpuid)

                try:
                    ret = check_output(jobsh)
                    logger.debug(ret)
                except Exception as e:
                    logger.info('Error in simulation {}. {}'.format(path, e))
                    self._setCompleted(path)
                    queue.task_done()
                    continue

                logger.info("Completed " + path)
                self._setCompleted(path)
                queue.task_done()

        logger.info("Shutting down worker thread")

    def _createJobScript(self, fname, workdir, runsh, gpudevice=None):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n\n')
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(os.path.normpath(os.path.join(workdir, self._sentinel))))
            f.write('\n')
            if gpudevice is not None:
                f.write('export CUDA_VISIBLE_DEVICES={}\n'.format(gpudevice))
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(os.path.normpath(os.path.join(workdir, self._sentinel))))
            f.write('\n')
            f.write('cd {}\n'.format(os.path.abspath(workdir)))
            f.write('{}'.format(runsh))

            # Move completed trajectories
            if self.datadir is not None:
                datadir = os.path.abspath(self.datadir)
                if not os.path.isdir(datadir):
                    os.mkdir(datadir)
                simname = os.path.basename(os.path.normpath(workdir))
                # create directory for new file
                odir = os.path.join(datadir, simname)
                os.mkdir(odir)
                f.write('\nmv {} {}'.format(' '.join(self.copy), odir))

        os.chmod(fname, 0o700)

    def _setRunning(self, path):
        self._states[path] = 'R'

    def _setCompleted(self, path):
        self._states[path] = 'C'

    def retrieve(self):
        """ Retrieves a list of jobs that have completed since the last call

        Example
        -------
        >>> comp = app.retrieve()
        """
        ret = []
        xx = self._states.copy().keys()
        for i in xx:
            if self._states[i] == 'C':
                del self._states[i]
                ret.append(i)

        return ret

    def submit(self, mydirs):
        """ Queue for execution all of the jobs in the passed list of directories

        Queues all work units in a given directory list with the options given in the constructor opt.

        Parameters
        ----------
        mydirs : list of str
            A list or ndarray of directory paths

        Examples
        --------
        >>> app.submit(glob('input/e2*/'))
        """
        self._setupQueue()

        if isinstance(mydirs, str):
            mydirs = [mydirs]
        self._dirs.extend(mydirs)

        for d in mydirs:
            if not os.path.isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        # if all folders exist, submit
        for d in mydirs:
            dirname = os.path.abspath(d)
            logger.info('Queueing ' + dirname)

            # Clean sentinel files , if existent
            if os.path.exists(os.path.join(d, self._sentinel)):
                try:
                    os.remove(os.path.join(d, self._sentinel))
                except:
                    logger.warning('Could not remove {} sentinel from {}'.format(self._sentinel, d))
                else:
                    logger.info('Removed existing {} sentinel from {}'.format(self._sentinel, d))

            self._states[dirname] = 'Q'
            self._queue.put(dirname)

    def inprogress(self):
        """ Get the number of simulations in progress

        Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Example
        -------
        >>> app.inprogress()
        """
        output_run = sum(x == 'R' for x in self._states.values())
        output_queue = sum(x == 'Q' for x in self._states.values())

        return output_run + output_queue

    def notcompleted(self):
        """Returns the sum of the number of job directories which do not have the sentinel file for completion.

        Returns
        -------
        total : int
            Total number of directories which have not completed
        """
        total = 0
        if len(self._dirs) == 0:
            raise RuntimeError('This method relies on running synchronously.')
        for i in self._dirs:
            if not os.path.exists(os.path.join(i, self._sentinel)):
                total += 1
        return total

    def stop(self):
        self._shutdown = True

    @abstractmethod
    def _getdevices(self):
        return list()

    def _getmemory(self):
        ureg = UnitRegistry()
        total_memory = int(ureg.Quantity(psutil.virtual_memory().total, ureg.byte).to('MiB').magnitude)
        nr_devices = len(self._getdevices())
        if nr_devices != 0:
            return int(total_memory/nr_devices)
        else:
            return None


class LocalGPUQueue(_LocalQueue):
    """ Local machine queue system

    Parameters
    ----------
    datadir : str, default=None
        The path in which to store completed trajectories.
    copy : list, default='*.xtc'
        A list of file names or globs for the files to copy to datadir
    ngpu : int, default=None
        Number of GPU devices that the queue will use. Each simulation will be run on a different GPU. The queue will
        use the first `ngpu` devices of the machine.
    devices : list, default=None
        A list of GPU device indexes on which the queue is allowed to run simulations. Mutually exclusive with `ngpu`
    memory : int, default=None
        The amount of RAM memory available for each job. If None, it will be guessed from total amount of memory and
        the number of devices


    .. rubric:: Methods
    .. autoautosummary:: htmd.queues.localqueue.LocalGPUQueue
       :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.queues.localqueue.LocalGPUQueue
       :attributes:

    """

    def __init__(self):
        super().__init__()
        self._arg('ngpu', 'int', 'Number of GPU devices that the queue will use. Each simulation will be run on '
                                 'a different GPU. The queue will use the first `ngpus` devices of the machine.',
                  None, val.Number(int, '0POS'))
        self._arg('devices', 'list', 'A list of GPU device indexes on which the queue is allowed to run '
                                     'simulations. Mutually exclusive with `ngpus`', None, val.Number(int, '0POS'),
                  nargs='*')
        self._arg('memory', 'int', 'The amount of RAM memory available for each job. If None, it will be guessed from '
                                   'total amount of memory and the number of devices', None,
                  val.Number(int, '0POS'))

    def _getdevices(self):
        ngpu = self.ngpu
        devices = self.devices
        if ngpu is not None and devices is not None:
            raise ValueError('Parameters `ngpu` and `devices` are mutually exclusive.')

        if ngpu is None and devices is None:
            try:
                check_output("nvidia-smi -L", shell=True)
                devices = range(int(check_output("nvidia-smi -L | wc -l", shell=True).decode("ascii")))
            except:
                raise
        elif ngpu is not None:
            devices = range(ngpu)

        if devices is None:
            raise NameError("Could not determine which GPUs to use. "
                            "Specify the GPUs with the `ngpu=` or `devices=` parameters")
        else:
            logger.info("Using GPU devices {}".format(','.join(map(str, devices))))
        return devices

    @property
    def ngpu(self):
        return self.__dict__['ngpu']

    @ngpu.setter
    def ngpu(self, value):
        self.ngpu = value

    @property
    def ncpu(self):
        raise NotImplementedError

    @ncpu.setter
    def ncpu(self, value):
        raise NotImplementedError

    @property
    def memory(self):
        return self.__dict__['memory']

    @memory.setter
    def memory(self, value):
        self.memory = value


class LocalCPUQueue(_LocalQueue):
    """ Local CPU machine queue system

    Parameters
    ----------
    datadir : str, default=None
        The path in which to store completed trajectories.
    copy : list, default='*.xtc'
        A list of file names or globs for the files to copy to datadir
    ncpu : int
        Number of CPU threads that the queue will use. By default, it is the total number of cpus.
    memory : int
        The amount of RAM memory available (in MiB). By default, it will be calculated from total amount
        of memory and the number of devices

    .. rubric:: Methods
    .. autoautosummary:: htmd.queues.localqueue.LocalCPUQueue
       :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.queues.localqueue.LocalCPUQueue
       :attributes:

    """

    def __init__(self):
        super().__init__()
        self._arg('ncpu', 'int', 'Number of CPU threads that the queue will use. If None it will use the `ncpu` '
                                 'configured for HTMD in htmd.configure()', psutil.cpu_count(), val.Number(int, 'POS'))
        self._arg('memory', 'int', 'The amount of RAM memory available for each job.', self._getmemory(),
                  val.Number(int, '0POS'))

    def _getdevices(self):
        return [None] * self.ncpu

    def _getmemory(self):

        memory = psutil.virtual_memory().total/1024**2
        memory *= np.clip(self.ncpu/psutil.cpu_count(), 0, 1)
        memory = int(np.floor(memory))

        return memory

    @property
    def ncpu(self):
        return self.__dict__['ncpu']

    @ncpu.setter
    def ncpu(self, value):
        self.ncpu = value

    @property
    def ngpu(self):
        raise NotImplementedError

    @ngpu.setter
    def ngpu(self, value):
        raise NotImplementedError

    @property
    def memory(self):
        return self.__dict__['memory']

    @memory.setter
    def memory(self, value):
        self.memory = value

if __name__ == "__main__":
    from htmd.home import home

    lo = LocalCPUQueue()

    assert lo.ncpu == psutil.cpu_count()
    assert lo.memory > 1024

    lo.ncpu = 100
    mem1 = lo.memory
    assert lo.ncpu == 100

    lo.ncpu = 1
    mem2 = lo.memory
    assert lo.ncpu == 1

    assert mem1 >= mem2

    folder = os.path.join(home(dataDir='test-localqueue'), 'test_cpu')
    lo.submit([folder] * 2)
    lo.wait(sentinel=False)
    lo.retrieve()

    lo.submit([folder] * 2)
    lo.wait(sentinel=True)
    lo.retrieve()

    lo = LocalGPUQueue()
    try:
        lo._getdevices()
        folders = glob(os.path.join(home(dataDir='test-localqueue'), 'test_gpu*'))
        lo.submit(folders)
        lo.wait()
        for f in folders:
            torem = glob(os.path.join(f, 'output.*')) + glob(os.path.join(f, 'restart.*')) + glob(os.path.join(f, 'log.txt'))
            for r in torem:
                os.remove(r)
    except:
        print('No GPUs detected on this machine')
