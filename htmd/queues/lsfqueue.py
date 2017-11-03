# (c) 2015-2017 Acellera Ltd http://www.acellera.com
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
from protocolinterface import ProtocolInterface, val
from htmd.queues.simqueue import SimQueue
from htmd.config import _config
import yaml
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
    app : str, default=None
        The application profile
    ngpu : int, default=1
        Number of GPUs to use for a single job
    ncpu : int, default=1
        Number of CPUs to use for a single job
    memory : int, default=4000
        Amount of memory per job (MiB)
    walltime : int, default=None
        Job timeout (hour:min or min)
    environment : list of strings, default=None
        Things to run before the job (sourcing envs).
    resources : list of strings, default=None
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

    _defaults = {'queue': None, 'app': None, 'gpu_queue': None, 'cpu_queue': None, 'ngpu': 1, 'ncpu': 1,
                 'memory': 4000, 'walltime': None, 'resources': None, 'environment': None}

    def __init__(self, _configapp=None):
        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)
        self._arg('jobname', 'str', 'Job name (identifier)', None, val.String())
        self._arg('queue', 'str', 'The queue to run on', self._defaults['queue'], val.String())
        self._arg('app', 'str', 'The application profile', self._defaults['app'], val.String())
        self._arg('ngpu', 'int', 'Number of GPUs to use for a single job', self._defaults['ngpu'],
                  val.Number(int, '0POS'))
        self._arg('ncpu', 'int', 'Number of CPUs to use for a single job', self._defaults['ncpu'],
                  val.Number(int, '0POS'))
        self._arg('memory', 'int', 'Amount of memory per job (MB)', self._defaults['memory'], val.Number(int, '0POS'))
        self._arg('walltime', 'int', 'Job timeout (hour:min or min)', self._defaults['walltime'], val.Number(int, '0POS'))
        self._arg('resources', 'list', 'Resources of the queue', self._defaults['resources'], val.String(), nargs='*')
        self._arg('environment', 'list', 'Things to run before the job (sourcing envs).', self._defaults['environment'],
                  val.String(), nargs='*')
        self._arg('outputstream', 'str', 'Output stream.', 'lsf.%J.out', val.String())
        self._arg('errorstream', 'str', 'Error stream.', 'lsf.%J.err', val.String())
        self._arg('datadir', 'str', 'The path in which to store completed trajectories.', None, val.String())
        self._arg('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc', val.String())

        # Load LSF configuration profile
        lsfconfig = _config['lsf']
        profile = None
        if _configapp is not None:
            if lsfconfig is not None:
                if os.path.isfile(lsfconfig) and lsfconfig.endswith(('.yml', '.yaml')):
                    try:
                        with open(lsfconfig, 'r') as f:
                            profile = yaml.load(f)
                        logger.info('Loaded LSF configuration YAML file {}'.format(lsfconfig))
                    except:
                        logger.warning('Could not load YAML file {}'.format(lsfconfig))
                else:
                    logger.warning('{} does not exist or it is not a YAML file.'.format(lsfconfig))
                if profile:
                    try:
                        properties = profile[_configapp]
                    except:
                        raise RuntimeError('There is no profile in {} for configuration '
                                           'app {}'.format(lsfconfig, _configapp))
                    for p in properties:
                        self.__dict__[p] = properties[p]
                        logger.info('Setting {} to {}'.format(p, properties[p]))
            else:
                raise RuntimeError('No LSF configuration YAML file defined for the configapp')
        else:
            if lsfconfig is not None:
                logger.warning('LSF configuration YAML file defined without configuration app')

        # Find executables
        self._qsubmit = LsfQueue._find_binary('bsub')
        self._qinfo = LsfQueue._find_binary('bqueues')
        self._qcancel = LsfQueue._find_binary('bkill')
        self._qstatus = LsfQueue._find_binary('bjobs')

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
            f.write('#BSUB -n {}\n'.format(self.ncpu))
            f.write('#BSUB -app {}\n'.format(self.app))
            if self.ngpu != 0:
                f.write('#BSUB -R "select[ngpus>0] rusage[ngpus_excl_p={}]"\n'.format(self.ngpu))
            f.write('#BSUB -M {}\n'.format(self.memory))
            f.write('#BSUB -cwd {}\n'.format(workdir))
            f.write('#BSUB -outdir {}\n'.format(workdir))
            f.write('#BSUB -o {}\n'.format(self.outputstream))
            f.write('#BSUB -e {}\n'.format(self.errorstream))
            if self.walltime is not None:
                f.write('#BSUB -W {}\n'.format(self.walltime))
            if self.resources is not None:
                for resource in self.resources:
                    f.write('#BSUB -R "{}"\n'.format(resource))
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(os.path.normpath(os.path.join(workdir, self._sentinel))))
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

        os.chmod(fname, 0o700)

    def retrieve(self):
        # Nothing to do
        pass

    def _autoJobName(self, path):
        return os.path.basename(os.path.abspath(path)) + '_' + ''.join([random.choice(string.digits) for _ in range(5)])

    def submit(self, dirs):
        """ Submits all directories

        Parameters
        ----------
        dirs : list
            A list of executable directories.
        """
        dirs = self._submitinit(dirs)

        if self.queue is None:
            raise ValueError('The queue needs to be defined.')

        # if all folders exist, submit
        for d in dirs:
            logger.info('Queueing ' + d)

            if self.jobname is None:
                self.jobname = self._autoJobName(d)

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
                ret = check_output(self._qsubmit + " < " + jobscript, shell=True)
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
        import time
        import getpass
        if self.queue is None:
            raise ValueError('The queue needs to be defined.')
        if self.jobname is None:
            raise ValueError('The jobname needs to be defined.')
        user = getpass.getuser()
        cmd = [self._qstatus, '-J', self.jobname, '-u', user, '-q', self.queue]
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

    def stop(self):
        """ Cancels all currently running and queued jobs
        """
        import getpass
        if self.queue is None:
            raise ValueError('The queue needs to be defined.')
        user = getpass.getuser()
        cmd = [self._qcancel, '-J', self.jobname, '-u', user, '-q', self.queue]
        logger.debug(cmd)
        ret = check_output(cmd)
        logger.debug(ret.decode("ascii"))

    @property
    def ncpu(self):
        return self.__dict__['ncpu']

    @ncpu.setter
    def ncpu(self, value):
        self.ncpu = value

    @property
    def ngpu(self):
        return self.__dict__['ngpu']

    @ngpu.setter
    def ngpu(self, value):
        self.ngpu = value

    @property
    def memory(self):
        return self.__dict__['memory']

    @memory.setter
    def memory(self, value):
        self.memory = value


if __name__ == "__main__":
    # TODO: Create fake binaries for instance creation testing
    """
    q = LsfQueue()
    """
