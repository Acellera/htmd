# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
import random
import string
import numpy as np
from subprocess import check_output, CalledProcessError
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_ANY, RANGE_0POS, RANGE_POS
from htmd.queues.simqueue import SimQueue
import logging
logger = logging.getLogger(__name__)


class SlurmQueue(SimQueue, ProtocolInterface):
    """ Queue system for SLURM

    Parameters
    ----------
    jobname : str, default=None
        Job name (identifier)
    partition : str, default=None
        The queue (partition) to run on
    priority : str, default='gpu_priority'
        Job priority
    ngpu : int, default=1
        Number of GPUs to use for a single job
    ncpu : int, default=1
        Number of CPUs to use for a single job
    memory : int, default=1000
        Amount of memory per job (MB)
    walltime : int, default=None
        Job timeout (s)
    environment : str, default='ACEMD_HOME,HTMD_LICENSE_FILE'
        Envvars to propagate to the job.
    mailtype : str, default=None
        When to send emails. Separate options with commas like 'END,FAIL'.
    mailuser : str, default=None
        User email address.
    outputstream : str, default='slurm.%N.%j.out'
        Output stream.
    errorstream : str, default='slurm.%N.%j.err'
        Error stream.

    Examples
    --------
    >>> from htmd import *
    >>> s = SlurmQueue()
    >>> s.partition = 'multiscale'
    >>> s.submit('/my/runnable/folder/')  # Folder containing a run.sh bash script
    """
    def __init__(self):
        super().__init__()
        self._cmdString('jobname', 'str', 'Job name (identifier)', None)
        self._cmdString('partition', 'str', 'The queue (partition) to run on', None)
        self._cmdString('priority', 'str', 'Job priority', 'gpu_priority')
        self._cmdValue('ngpu', 'int', 'Number of GPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('ncpu', 'int', 'Number of CPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('memory', 'int', 'Amount of memory per job (MB)', 1000, TYPE_INT, RANGE_0POS)
        self._cmdValue('walltime', 'int', 'Job timeout (s)', None, TYPE_INT, RANGE_POS)
        self._cmdString('environment', 'str', 'Envvars to propagate to the job.', 'ACEMD_HOME,HTMD_LICENSE_FILE')
        self._cmdString('mailtype', 'str', 'When to send emails. Separate options with commas like \'END,FAIL\'.', None)
        self._cmdString('mailuser', 'str', 'User email address.', None)
        self._cmdString('outputstream', 'str', 'Output stream.', 'slurm.%N.%j.out')
        self._cmdString('errorstream', 'str', 'Error stream.', 'slurm.%N.%j.err')  # Maybe change these to job name
        self._cmdString('datadir', 'str', 'The path in which to store completed trajectories.', None)
        self._cmdString('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc')

        # Find executables
        self._qsubmit = SlurmQueue._find_binary('sbatch')
        self._qinfo = SlurmQueue._find_binary('sinfo')
        self._qcancel = SlurmQueue._find_binary('scancel')
        self._qstatus = SlurmQueue._find_binary('squeue')

        self._sentinel = 'htmd.queues.done'
        # For synchronous
        self._dirs = []

    @staticmethod
    def _find_binary(binary):
        ret = shutil.which(binary, mode=os.X_OK)
        if not ret:
            raise FileNotFoundError("Could not find required executable [{}]".format(binary))
        ret = os.path.abspath(ret)
        return ret

    def _createJobScript(self, fname, workdir, runsh):
        workdir = os.path.abspath(workdir)
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#\n')
            f.write('#SBATCH --job-name={}\n'.format(self.jobname))
            f.write('#SBATCH --partition={}\n'.format(self.partition))
            if self.ngpu != 0:
                f.write('#SBATCH --gres=gpu:{}\n'.format(self.ngpu))
            f.write('#SBATCH --cpus-per-task={}\n'.format(self.ncpu))
            f.write('#SBATCH --mem={}\n'.format(self.memory))
            f.write('#SBATCH --priority={}\n'.format(self.priority))
            f.write('#SBATCH --workdir={}\n'.format(workdir))
            f.write('#SBATCH --output={}\n'.format(self.outputstream))
            f.write('#SBATCH --error={}\n'.format(self.errorstream))
            if self.environment is not None:
                f.write('#SBATCH --export={}\n'.format(self.environment))
            if self.walltime is not None:
                f.write('#SBATCH --time={}\n'.format(self.walltime))
            if self.mailtype is not None and self.mailuser is not None:
                f.write('#SBATCH --mail-type={}\n'.format(self.mailtype))
                f.write('#SBATCH --mail-user={}\n'.format(self.mailuser))
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(self._sentinel))
            f.write('\ncd {}\n'.format(workdir))
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

    def retrieve(self):
        # Nothing to do
        pass

    def submit(self, dirs):
        """ Submits all directories

        Parameters
        ----------
        dirs : list
            A list of executable directories.
        """
        if isinstance(dirs, str):
            dirs = [dirs, ]
        self._dirs.extend(dirs)

        # Automatic partition
        if self.partition is None:
            ret = check_output(self._qinfo)
            self.partition = ','.join(np.unique([i.split()[0].strip('*')
                                                 for i in ret.decode('ascii').split('\n')[1:-1]]))

        # if all folders exist, submit
        for d in dirs:
            logger.info('Queueing ' + d)

            # Automatic jobname
            if self.jobname is None:
                self.jobname = os.path.basename(os.path.abspath(d)) + '_' + \
                               ''.join([random.choice(string.digits) for _ in range(5)])

            runscript = os.path.abspath(os.path.join(d, 'run.sh'))

            # Clean sentinel files , if existent
            if os.path.exists(os.path.join(d, self._sentinel)):
                try:
                    os.remove(os.path.join(d, self._sentinel))
                except:
                    logger.warning('Could not remove {} sentinel from {}'.format(self._sentinel, d))
                else:
                    logger.info('Removed existing {} sentinel from {}'.format(self._sentinel, d))

            if not os.path.exists(runscript):
                raise FileExistsError('File {} does not exist.'.format(runscript))
            if not os.access(runscript, os.X_OK):
                raise PermissionError('File {} does not have execution permissions.'.format(runscript))

            jobscript = os.path.abspath(os.path.join(d, 'job.sh'))
            self._createJobScript(jobscript, d, runscript)
            try:
                ret = check_output([self._qsubmit, jobscript])
                logger.debug(ret)
            except:
                raise

    def inprogress(self, debug=False):
        """ Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Returns
        -------
        total : int
            Total running and queued workunits
        """
        import time
        import getpass
        if self.partition is None:
            raise ValueError('The partition needs to be defined.')
        user = getpass.getuser()
        cmd = [self._qstatus, '-n', self.jobname, '-u', user, '-p', self.partition]
        logger.debug(cmd)

        # This command randomly fails so I need to allow it to repeat or it crashes adaptive
        tries = 0
        while tries < 3:
            try:
                ret = check_output(cmd)
            except CalledProcessError:
                if tries == 2:
                    raise
                tries += 1
                time.sleep(3)
                continue
            break

        logger.debug(ret.decode("ascii"))

        # TODO: check lines and handle errors
        l = ret.decode("ascii").split("\n")
        l = len(l) - 2
        if l < 0:
            l = 0  # something odd happened
        return l

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

    def stop(self):
        """ Cancels all currently running and queued jobs
        """
        import getpass
        if self.partition is None:
            raise ValueError('The partition needs to be defined.')
        user = getpass.getuser()
        cmd = [self._qcancel, '-n', self.jobname, '-u', user, '-p', self.partition]
        logger.debug(cmd)
        ret = check_output(cmd)
        logger.debug(ret.decode("ascii"))

    def wait(self, sentinel=False):
        """ Blocks script execution until all queued work completes

        Parameters
        ----------
        sentinel : bool, default=False
            If False, it relies on the queueing system reporting to determine the number of running jobs. If True, it
            relies on the filesystem, in particular on the existence of a sentinel file for job completion.

        Examples
        --------
        >>> SlurmQueue.wait()
        """
        from time import sleep
        import sys

        while (self.inprogress() if not sentinel else self.notcompleted()) != 0:
            sys.stdout.flush()
            sleep(5)


if __name__ == "__main__":
    """
    s=Slurm( name="testy", partition="gpu")
    s.submit("test/dhfr1" )
    ret= s.inprogress( debug=False)
    print(ret)
    print(s)
    pass
    """
