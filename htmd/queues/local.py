from htmd.queues.simqueue import SimQueue
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_INT, RANGE_0POS
import queue
import threading
from subprocess import check_output
import os
from glob import glob as glob
import logging
logger = logging.getLogger(__name__)


class LocalGPUQueue(SimQueue, ProtocolInterface):
    """
    Parameters
    ----------
    ngpu : int
        Number of GPU devices that the queue will use. Each simulation will be run on a different GPU. The queue will
        use the first `ngpus` devices of the machine.
    devices : list, default=None
        A list of GPU device indexes on which the queue is allowed to run simulations. Mutually exclusive with `ngpus`
    datadir : str, default=None
        The path in which to store completed trajectories.
    trajext : str, default='xtc'
        Extension of trajectory files. This is needed to copy them to datadir.
    """

    def __init__(self):
        super().__init__()
        self._cmdValue('ngpu', 'int', 'Number of GPU devices that the queue will use. Each simulation will be run on '
                                      'a different GPU. The queue will use the first `ngpus` devices of the machine.',
                       None, TYPE_INT, RANGE_0POS)
        self._cmdList('devices', 'list', 'A list of GPU device indexes on which the queue is allowed to run '
                                         'simulations. Mutually exclusive with `ngpus`', None)
        self._cmdString('datadir', 'str', 'The path in which to store completed trajectories.', None)
        self._cmdString('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc')

        self._states = dict()
        self._queue = None
        self._shutdown = False

    def _getGPUdevices(self):
        ngpu = self.ngpu
        devices = self.devices
        if ngpu is not None and devices is not None:
            raise ValueError('Parameters `ngpu` and `devices` are mutually exclusive.')

        if ngpu is None and devices is None:
            try:
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

    def _setupQueue(self):
        if self._queue is None:
            self._queue = queue.Queue()

            devices = self._getGPUdevices()

            self._threads = []
            for d in devices:
                t = threading.Thread(target=run_job, args=(self, d))
                t.daemon = True
                t.start()
                self._threads.append(t)

    def _createJobScript(self, fname, workdir, runsh, device):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('export CUDA_VISIBLE_DEVICES={}\n'.format(device))
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
                f.write('\nmv *.{} {}'.format(self.trajext, odir))

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

        if isinstance(mydirs, str): mydirs = [mydirs]

        for d in mydirs:
            if not os.path.isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        # if all folders exist, submit
        for d in mydirs:
            dirname = os.path.abspath(d)
            logger.info('Queueing ' + dirname)
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

    def stop(self):
        self._shutdown = True


def run_job(self, gpuid):
    queue = self._queue
    while not self._shutdown:
        path = None
        try:
            path = queue.get(timeout=1)
        except:
            pass

        if path:
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


if __name__ == "__main__":
    from time import sleep

    # TODO: Fix this to work
    """
    a = AcemdLocal(acemd="/shared/acemd/bin/acemd", ngpus=2)
    a.submit("/tmp/job/1")
    a.submit("/tmp/job/2")
    a.submit("/tmp/job/3")
    a.submit("/tmp/job/4")
    a.wait()

    r = a.retrieve()
    if len(r):
        print("Retrieved: ")
        print(r)
    sleep(1)
    a.stop()
    sleep(10)
    print("Done")
    """

