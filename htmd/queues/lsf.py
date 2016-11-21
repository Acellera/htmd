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
from subprocess import check_output
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_ANY, RANGE_0POS, RANGE_POS
from htmd.queues.simqueue import SimQueue
import logging
logger = logging.getLogger(__name__)


class LsfQueue(SimQueue, ProtocolInterface):
    """ Queue system for LSF

    Parameters
    ----------
    jobname : str, default=None
        Job name (identifier)
    queue : str, default=None
        The queue to run on
    ngpu : int, default=1
        Number of GPUs to use for a single job
    memory : int, default=4000
        Amount of memory per job (MB)
    walltime : int, default=None
        Job timeout (hour:min or min)
    environment : list of strings, default=None
        Things to run before the job (sourcing envs).
    resources : str, default=None
        Resources of the queue
    outputstream : str, default='slurm.%N.%j.out'
        Output stream.
    errorstream : str, default='slurm.%N.%j.err'
        Error stream.

    Examples
    --------
    >>> from htmd import *
    >>> s = LsfQueue()
    >>> s.jobname = 'simulation1'
    >>> s.queue = 'multiscale'
    >>> s.submit('/my/runnable/folder/')  # Folder containing a run.sh bash script
    """
    def __init__(self):
        super().__init__()
        self._cmdString('jobname', 'str', 'Job name (identifier)', None)
        self._cmdString('queue', 'str', 'The queue to run on', None)
        self._cmdValue('ngpu', 'int', 'Number of GPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('memory', 'int', 'Amount of memory per job (MB)', 4000, TYPE_INT, RANGE_0POS)
        self._cmdValue('walltime', 'int', 'Job timeout (hour:min or min)', None, TYPE_INT, RANGE_POS)
        self._cmdString('resources', 'str', 'Resources of the queue', None)
        self._cmdList('environment', 'list', 'Things to run before the job (sourcing envs).', None)
        self._cmdString('outputstream', 'str', 'Output stream.', 'lsf.%J.out')
        self._cmdString('errorstream', 'str', 'Error stream.', 'lsf.%J.err')
        self._cmdString('datadir', 'str', 'The path in which to store completed trajectories.', None)
        self._cmdString('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc')

        # Find executables
        self._qsubmit = LsfQueue._find_binary('bsub')
        self._qstatus = LsfQueue._find_binary('bqueues')
        self._qcancel = LsfQueue._find_binary('bkill')

        self._dirs = []

        # Specific automatic guessing
        ret = check_output(self._qstatus)
        if 'phase6_normal' in ret.decode('ascii'):
            if self.environment is None:
                self.environment = ['source /home/model/MD-SOFTWARE/model_md.bashrc', 'source /home/model/miniconda3/htmd.bashrc']
                logger.info('environment set to {}'.format(self.environment))
            if self.resources is None:
                self.resources = '"{}"'.format('rusage[ngpus_excl_p=1],span[hosts=1]')
                logger.info('resources set to {}'.format(self.resources))
        if 'gpu_priority' in ret.decode('ascii'):
            if self.environment is None:
                self.environment = ['module load acemd', 'module load acellera/test', 'module load gaussian']
                logger.info('environment set to {}'.format(self.environment))
            if self.resources is None:
                self.resources = '"{}"'.format('select[ngpus>0] rusage[ngpus_excl_p=1]')
                logger.info('resources set to {}'.format(self.resources))

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
            f.write('#BSUB -J {}\n'.format(self.jobname))
            f.write('#BSUB -q {}\n'.format(self.queue))
            f.write('#BSUB -n {}\n'.format(self.ngpu))
            f.write('#BSUB -M {}\n'.format(self.memory))
            f.write('#BSUB {}\n'.format(workdir))
            f.write('#BSUB -o {}\n'.format(self.outputstream))
            f.write('#BSUB -e {}\n'.format(self.errorstream))
            if self.walltime is not None:
                f.write('#BSUB -W {}\n'.format(self.walltime))
            if self.resources is not None:
                f.write('#BSUB -R {}\n'.format(self.resources))
            f.write('\n')
            if self.environment is not None:
                for call in self.environment:
                    f.write('{}\n'.format(call))
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

            f.write('\ntouch .done')
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
        if isinstance(dirs, str):
            dirs = [dirs, ]
        self._dirs.extend(dirs)

        # if all folders exist, submit
        for d in dirs:
            logger.info('Queueing ' + d)

            runscript = os.path.abspath(os.path.join(d, 'run.sh'))
            if not os.path.exists(runscript):
                raise FileExistsError('File {} does not exist.'.format(runscript))
            if not os.access(runscript, os.X_OK):
                raise PermissionError('File {} does not have execution permissions.'.format(runscript))

            jobscript = os.path.abspath(os.path.join(d, 'job.sh'))
            self._createJobScript(jobscript, d, runscript)
            try:
                ret = check_output(self._qsubmit + " < " + jobscript, shell=True)
                logger.debug(ret)
            except:
                raise

    def inprogress(self, debug=False):
        inprogress = 0
        for i in self._dirs:
            if not os.path.exists(os.path.join(i, '.done')):
                inprogress += 1
        return inprogress

    def __inprogress(self):
        """ Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Returns
        -------
        total : int
            Total running and queued workunits
        """
        if self.queue is None:
            raise ValueError('The queue needs to be defined.')
        user = pwd.getpwuid(os.getuid()).pw_name
        cmd = [self._qlist, '-J', self.jobname, '-u', user, '-q', self.queue]
        logger.debug(cmd)
        ret = check_output(cmd)
        logger.debug(ret.decode("ascii"))

        # TODO: check lines and handle errors
        l = ret.decode("ascii").split("\n")
        l = len(l) - 2
        if l < 0:
            l = 0  # something odd happened
        return l

    def stop(self):
        """ Cancels all currently running and queued jobs
        """
        if self.queue is None:
            raise ValueError('The queue needs to be defined.')
        user = pwd.getpwuid(os.getuid()).pw_name
        cmd = [self._qcancel, '-J', self.jobname, '-u', user, '-q', self.queue]
        logger.debug(cmd)
        ret = check_output(cmd)
        logger.debug(ret.decode("ascii"))



if __name__ == "__main__":
    """
    s=Slurm( name="testy", partition="gpu")
    s.submit("test/dhfr1" )
    ret= s.inprogress( debug=False)
    print(ret)
    print(s)
    pass
    """
