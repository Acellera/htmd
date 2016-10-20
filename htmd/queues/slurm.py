# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
# (c) 2015 Acellera Ltd
# All Rights Reserved
# Distributed under HTMD Software Academic License Agreement v1.0
# No redistribution in whole or part
#
import os
import pwd
import shutil
from os.path import isdir
from subprocess import check_output
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_ANY, RANGE_0POS, RANGE_POS
from htmd import UserInterface
from htmd.queues.queue import Queue
import logging
logger = logging.getLogger(__name__)


class SlurmQueue(Queue, ProtocolInterface):
    def __init__(self):
        super().__init__()
        self._cmdString('jobname', 'str', 'Job name (identifier)', None)
        self._cmdString('partition', 'str', 'The queue (partition) to run on', None)
        self._cmdString('priority', 'str', 'Job priority', 'gpu_priority')
        self._cmdValue('ngpu', 'int', 'Number of GPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('memory', 'int', 'Amount of memory per job (MB)', 4000, TYPE_INT, RANGE_0POS)
        self._cmdValue('walltime', 'int', 'Job timeout (s)', None, TYPE_INT, RANGE_POS)
        self._cmdString('environment', 'str', 'Envvars to propagate to the job.', 'ACEMD_HOME,HTMD_LICENSE_FILE')
        self._cmdString('mailtype', 'str', 'When to send emails. Separate options with commas like \'END,FAIL\'.', None)
        self._cmdString('mailuser', 'str', 'User email address.', None)
        self._cmdString('outputstream', 'str', 'Output stream.', 'slurm.%N.%j.out')
        self._cmdString('errorstream', 'str', 'Error stream.', 'slurm.%N.%j.err')  # Maybe change these to job name

        # Find executables
        self._sbatch = SlurmQueue._find_binary('sbatch')
        self._squeue = SlurmQueue._find_binary('squeue')

    @staticmethod
    def _find_binary(binary):
        ret = shutil.which(binary, mode=os.X_OK)
        if not ret:
            raise FileNotFoundError("Could not find required executable [{}]".format(binary))
        ret = os.path.abspath(ret)
        return ret

    def _createJobScript(self, fname, workdir, runsh):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#\n')
            f.write('#SBATCH --job-name={}\n'.format(self.jobname))
            f.write('#SBATCH --partition={}\n'.format(self.partition))
            f.write('#SBATCH --gres=gpu:{}\n'.format(self.ngpu))
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
            f.write('\ncd {}\n'.format(workdir))
            f.write('{}'.format(runsh))
        os.chmod(fname, 0o700)

    def retrieve(self):
        # Nothing to do
        pass

    def submit(self, dirs):
        """ Submits all directories

        Parameters
        ----------
        dist : list
            A list of executable directories.
        """
        import time

        # if all folders exist, submit
        for d in dirs:
            logger.info('Queueing ' + d)

            runscript = os.path.join(d, 'run.sh')
            if not os.path.exists(runscript):
                raise FileExistsError('File {} does not exist.'.format(runscript))
            if not os.access(runscript, os.X_OK):
                raise PermissionError('File {} does not have execution permissions.'.format(runscript))

            jobscript = os.path.join(d, 'job.sh')
            self._createJobScript(jobscript, d, runscript)
            try:
                ret = check_output([self._sbatch, jobscript])
                logger.debug(ret)
            except:
                raise

    def inprogress(self):
        """ Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Returns
        -------
        total : int
            Total running and queued workunits
        """
        if self.name is None:
            raise ValueError('Project name needs to be defined.')
        if self.queue is None:
            raise ValueError('Project queue needs to be defined.')
        user = pwd.getpwuid(os.getuid()).pw_name
        cmd = [self._squeue, "-n", self.name, "-u", user, "-p", self.queue]
        logger.debug(cmd)
        ret = check_output(cmd)
        logger.debug(ret.decode("ascii"))

        # TODO: check lines and handle errors
        l = ret.decode("ascii").split("\n")
        l = len(l) - 2
        if l < 0:
            l = 0  # something odd happened
        return l


class AcemdSlurm(Slurm):
    def __init__(self):
        super().__init__()
        self._cmdString('acemd', 'str', 'Path to acemd executable', None)
        self._cmdString('inputfile', 'str', 'Name of input files', 'input')
        self._cmdValue('gpuid', 'int', 'GPU id', 0)

    def _findAcemd(self):
        from shutil import which
        if not self.acemd:
            try:
                # Try to find acemd in the path
                acemd = which("acemd", mode=os.X_OK)
                if acemd:
                    logger.info("Found ACEMD at '" + acemd + "'")
                else:
                    acemd = which("acemdhtmd", mode=os.X_OK)
                    if acemd:
                        logger.info("Found ACEMD at '" + acemd + "'")
            except:
                pass
        else:
            acemd = self.acemd

        if not acemd:
            raise NameError("Cannot find 'acemd' in the PATH. Set its location with the 'acemd=' named argument")

        if not os.access(acemd, os.X_OK):
            raise NameError("ACEMD file '" + acemd + "' not executable")

    def retrieve(self):
        # TODO: Make it move the completed trajectories
        pass

    def submit(self, mydirs):
        """ Submits all work units in a given directory list to the engine.

        Parameters
        ----------
        mydirs : list
            A list of job directories
        """
        if isinstance(mydirs, str): mydirs = [mydirs]

        for d in mydirs:
            if not os.path.isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        acemd = SlurmQueue._find_binary(self.acemd)

        for d in mydirs:
            logger.info('Queueing ' + dirname)
            dirname = os.path.abspath(d)
            self._make_jobscript(dirname, self.gpuid)

        super().submit(mydirs)

    def _make_jobscript(self, path, gpuid, acemd, inputfile):
        fn = os.path.join(path, "run.sh")
        with open(fn, 'w') as f:
            f.write('#!/bin/sh')
            f.write('{} --device {} {} > log.txt 2>&1'.format(acemd, gpuid, inputfile))
        os.chmod(fn, 0o700)
        return os.path.abspath(fn)


if __name__ == "__main__":
    """
    s=Slurm( name="testy", partition="gpu")
    s.submit("test/dhfr1" )
    ret= s.inprogress( debug=False)
    print(ret)
    print(s)
    pass
    """
