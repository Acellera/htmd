# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
from htmd.apps.app import App
import queue
import threading
from shutil import which, move
from subprocess import call, check_output, CalledProcessError
import os
from glob import glob as glob
import logging
logger = logging.getLogger(__name__)


class AcemdLocal(App):
    def __init__(self, acemd=None, ngpus=None, datadir=None):
        self.states = dict()
        self.queue = queue.Queue()
        self.threads = []
        self.shutdown = False

        if not ngpus:
            try:
                ngpus = int(check_output("nvidia-smi -L | wc -l", shell=True).decode("ascii"))
                logger.info("Found " + str(ngpus) + " GPUs")
            except:
                raise

        if not ngpus:
            raise NameError("Could not determine the number of GPUs. Specify it with the 'ngpus=' named argument")

        if not acemd:
            try:
                # Try to find acemd in the path
                acemd = which("acemd", mode=os.X_OK)
                if acemd:
                    logger.info("Found ACEMD at '" + acemd + "'")
            except:
                pass

        if not acemd:
            raise NameError("Cannot find 'acemd' in the PATH. Set its location with the 'acemd=' named argument")

        if not os.access(acemd, os.X_OK):
            raise NameError("ACEMD file '" + acemd + "' not executable")
        
       # logger.info("Executing ACEMD")

        self.threads = []
        for i in range(ngpus):
            t = threading.Thread(target=run_job, args=(self, i, acemd, datadir))
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
        >>> comp = grid.retrieve()
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
        >>> grid.submit(glob('input/e2*/'));
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
        >>> grid.inprogress()
        """
        output_run = sum(x == 'R' for x in self.states.values())
        output_queue = sum(x == 'Q' for x in self.states.values())

        return output_run + output_queue

    def wait(self):
        """ WAIT - Block until all queued work complete
        """
        from time import sleep
        import sys
        while self.inprogress() != 0:
            sys.stdout.flush()
            sleep(1)
        #self.queue.join()

    def stop(self):
        self.shutdown = True


def run_job(obj, ngpu, acemd, datadir):
    import sys
    queue = obj.queue
    while not obj.shutdown:
        path = None
        try:
            path = queue.get(timeout=1)
        except:
            pass

        if path:
            try:
                logger.info("Running " + path + " on GPU " + str(ngpu))
                obj.running(path)
                cmd = 'cd {}; {} --device {} input > log.txt 2>&1'.format(os.path.normpath(path), acemd, ngpu)
                try:
                    check_output(cmd, shell=True)
                except CalledProcessError:
                    logger.error('Error in ACEMD for path: {}. Check the {} file.'.format(path, os.path.join(path, 'log.txt')))
                    obj.completed(path)
                    queue.task_done()
                    continue

                # If a datadir is provided, copy finished trajectories there. Only works for xtc files.
                if datadir is not None:
                    if not os.path.isdir(datadir):
                        os.mkdir(datadir)
                    simname = os.path.basename(os.path.normpath(path))
                    odir = os.path.join(datadir, simname)
                    os.mkdir(odir)
                    finishedtraj = glob(os.path.join(path, '*.xtc'))
                    logger.info("Moving simulation {} to {}.".format(finishedtraj[0], odir))
                    move(finishedtraj[0], odir)

                logger.info("Completed " + path)
                obj.completed(path)
                queue.task_done()
            except:
                logger.error("Error running job")
                obj.completed(path)
                queue.task_done()
                continue
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
