from htmd.apps.app import App
import queue
import threading
from subprocess import check_output
import os
from glob import glob as glob
import logging
logger = logging.getLogger(__name__)


class LocalGPUQueue(App):
    """
    Parameters
    ----------
    ngpus : int
        Number of GPU devices that the queue will use. Each simulation will be run on a different GPU. The queue will
        use the first `ngpus` devices of the machine.
    devices : list
        A list of GPU device indexes on which the queue is allowed to run simulations. Mutually exclusive with `ngpus`
    """

    def __init__(self, jobfunc, jobargs, ngpus=None, devices=None):
        self.states = dict()
        self.queue = queue.Queue()
        self.threads = []
        self.shutdown = False

        if ngpus is not None and devices is not None:
            raise ValueError('Parameters `ngpus` and `devices` are mutually exclusive.')

        if ngpus is None and devices is None:
            try:
                devices = range(int(check_output("nvidia-smi -L | wc -l", shell=True).decode("ascii")))
            except:
                raise
        elif ngpus is not None:
            devices = range(ngpus)

        if devices is None:
            raise NameError("Could not determine which GPUs to use. "
                            "Specify the GPUs with the `ngpus=` or `devices=` parameters")
        else:
            logger.info("Using GPU devices {}".format(','.join(map(str, devices))))

        self.threads = []
        for d in devices:
            t = threading.Thread(target=run_job, args=(self, d, jobfunc, jobargs))
            t.daemon = True
            t.start()
            self.threads.append(t)

    def running(self, path):
        self.states[path] = 'R'

    def completed(self, path):
        self.states[path] = 'C'

    def retrieve(self):
        """ Retrieves a list of jobs that have completed since the last call

        Example
        -------
        >>> comp = app.retrieve()
        """
        ret = []
        xx = self.states.copy().keys()
        for i in xx:
            if self.states[i] == 'C':
                del self.states[i]
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
        >>> app.submit(glob('input/e2*/'));
        """
        if isinstance(mydirs, str): mydirs = [mydirs]

        for d in mydirs:
            if not os.path.isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        # if all folders exist, submit
        for d in mydirs:
            dirname = os.path.abspath(d)
            logger.info('Queueing ' + dirname)
            self.states[dirname] = 'Q'
            self.queue.put(dirname)

    def inprogress(self):
        """ Get the number of simulations in progress

        Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Example
        -------
        >>> app.inprogress()
        """
        output_run = sum(x == 'R' for x in self.states.values())
        output_queue = sum(x == 'Q' for x in self.states.values())

        return output_run + output_queue

    def wait(self):
        """ Blocks script execution until all queued work completes

        Examples
        --------
        >>> app.wait()
        """
        from time import sleep
        import sys
        while self.inprogress() != 0:
            sys.stdout.flush()
            sleep(1)
            # self.queue.join()

    def stop(self):
        self.shutdown = True


def run_job(obj, gpuid, jobfun, jobargs):
    queue = obj.queue
    while not obj.shutdown:
        path = None
        try:
            path = queue.get(timeout=1)
        except:
            pass

        if path:
            try:
                logger.info("Running " + path + " on GPU device " + str(gpuid))
                obj.running(path)

                try:
                    jobfun(*jobargs, path=path, gpuid=gpuid)
                except:
                    obj.completed(path)
                    queue.task_done()
                    continue

                logger.info("Completed " + path)
                obj.completed(path)
                queue.task_done()
            except:
                logger.error("Error running job {}".format(path))
                obj.completed(path)
                queue.task_done()
                continue
    logger.info("Shutting down worker thread")


def _executeMDcommand(cmd, path, datadir, enginename, trajext):
    """ Executes command line MD simulation. Can move finished simulations to datadir.

    Parameters
    ----------
    cmd : str
        The command-line string to execute
    path : str
        The path containing the simulation files. Used only to report errors
    datadir : str
        The path in which to store completed trajectories. Default is None and doesn't move them.
    enginename : str
        Name of the MD engine to report nice error
    trajext : str
        Extension of completed trajectories to be able to move them to datadir
    """
    from subprocess import PIPE, Popen, TimeoutExpired, CalledProcessError
    proc = None
    try:
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        proc.communicate()
    except CalledProcessError:
        logger.error('Error in {} for path: {}. Check the {} file.'.format(enginename, path, os.path.join(path, 'log.txt')))
        if proc:
            proc.kill()
        raise
    except TimeoutExpired:
        if proc:
            proc.kill()
        raise

    import shutil
    # If a datadir is provided, copy finished trajectories there.
    # Only works for nc files.
    if datadir is not None:
        if not os.path.isdir(datadir):
            os.mkdir(datadir)
        simname = os.path.basename(os.path.normpath(path))
        # create directory for new file
        odir = os.path.join(datadir, simname)
        os.mkdir(odir)
        finishedtraj = glob(os.path.join(path, trajext))
        logger.info('Moving simulation {} to {}.'.format(finishedtraj[0], odir))
        shutil.move(finishedtraj[0], odir)


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

