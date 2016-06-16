from htmd.apps.app import App
from htmd.apps.acemdlocal import AcemdLocal
from subprocess import check_output, CalledProcessError
from glob import glob as glob
from shutil import which, move, copy
import threading
import logging
import os
import queue
import sys
logger = logging.getLogger(__name__)


class PmemdLocal(App):
    """
    Parameters
    ----------
    pmemd_cuda: str
        Path to any executable pmemd binary. If None, it will try to detect
        the pmemd.cuda_SPFP binary
    ngpus : int
        Number of GPU devices that PmemdLocal will use. Each simulation will be
        run on a different GPU. PmemdLocal will use the first `ngpus` devices
        of the machine.
    devices : list
        A list of GPU device indexes on which PmemdLocal is allowed to run
        simulations. Mutually exclusive with `ngpus`
    datadir : str
        A folder to which completed simulations will be moved. If None they
        will be written in the input directory.
    """
    def __init__(self, pmemd_cuda=None, ngpus=None, devices=None, datadir=None):
        self.states = dict()
        self.queue = queue.Queue()
        self.threads = []
        self.shutdown = False
        self.system_name = 'Test'

        if ngpus is not None and devices is not None:
            raise ValueError('Parameters `ngpus` and `devices` are mutually exclusive.')

        if ngpus is None and devices is None:
            try:
                devices = range(int(check_output("nvidia-smi -L | wc -l",
                                shell=True).decode("ascii")))
            except:
                raise
        elif ngpus is not None:
            devices = range(ngpus)

        if devices is None:
            raise NameError("""Could not determine which GPUs to use.
                            Specify the GPUs with the `ngpus=`
                            or `devices=` parameters""")
        else:
            logger.info("Using GPU devices {}".format(','.join(map(str, devices))))

        if not pmemd_cuda:
            try:
                # Try to find pmemd.cuda_SPFP in the path
                pmemd_cuda = which("pmemd.cuda_SPFP", mode=os.X_OK)
                if pmemd_cuda:
                    logger.info("Found pmemd.cuda_SPFP at '" + pmemd_cuda + "'")
            except:
                raise NameError("""Cannot find 'pmemd.cuda_SPFP' in the PATH
                                Set its location with the 'pmemd_cuda=' named
                                argument""")

        if not os.access(pmemd_cuda, os.X_OK):
            raise NameError("pmemd_cuda file '" + pmemd_cuda + "' not executable")

        self.threads = []
        for d in devices:
            t = threading.Thread(target=run_job, args=(self, d, pmemd_cuda,
                                                       datadir, self.system_name))
            t.daemon = True
            t.start()
            self.threads.append(t)

    # Use the methods from the AcemdLocal class
    retrieve = AcemdLocal.retrieve
    inprogress = AcemdLocal.inprogress
    running = AcemdLocal.running

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

            # FIXME copying files to queued directory
            copy('/home/je714/htmd/htmd/protocols/Amber_protocols/Production.in', dirname)

            self.states[dirname] = 'Q'
            self.queue.put(dirname)


def run_job(obj, ngpu, pmemd_cuda, datadir, system_name):
    queue = obj.queue
    while not obj.shutdown:
        path = None
        try:
            path = queue.get(timeout=1)
        except:
            pass
        if path:
            try:
                logger.info('Running {} on GPU device {}'.format(path, str(ngpu)))
                obj.running(path)
                # TODO: Pre-production steps

                # AMBER automatically runs the calculation on the GPU with the
                # most memory even if that GPU is already in use
                # See: http://ambermd.org/gpus/#Running
                # We assume the machine doesn't have the exhaustive mode

                cmd = 'export CUDA_VISIBLE_DEVICES="{}"\n'.format(str(ngpu))
                cmd += """cd {} && {} -O -i Production.in -o {}.out \\
                         -c {} -p {}.prmtop -r {}.rst \\
                         -x {}.nc > log.txt 2>&1""".format(os.path.normpath(path),
                                                           pmemd_cuda,
                                                           system_name,
                                                           system_name,
                                                           system_name,
                                                           system_name,
                                                           system_name)
                try:
                    check_output(cmd, shell=True)
                except CalledProcessError:
                    logger.error('Error in pmemd.cuda for path: {}. Check the {} file.'.format(path, os.path.join(path, 'log.txt')))
                    obj.completed(path)
                    queue.task_done()
                    continue
                # If a datadir is provided, copy finished trajectories there.
                # Only works for nc files.
                if datadir is not None:
                    if not os.path.isdir(datadir):
                        os.mkdir(datadir)
                    simname = os.path.basename(os.path.normpath(path))
                    odir = os.path.join(datadir, simname)
                    os.mkdir(odir)
                    finishedtraj = glob(os.path.join(path, '*.nc'))
                    logger.info('Moving simulation {} to {}.'.format(finishedtraj[0], odir))
                    move(finishedtraj[0], odir)
                logger.info('Completed {}'.format(path))
                obj.completed(path)
                queue.task_done()
            except:
                logger.error('Error running job')
                obj.completed(path)
                queue.task_done()
                continue
    logger.info('Shutting down worker thread')
