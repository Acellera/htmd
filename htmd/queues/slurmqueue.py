# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
import random
import string
from htmd.config import _config
import yaml
from subprocess import check_output, CalledProcessError
from protocolinterface import ProtocolInterface, val
from htmd.queues.simqueue import SimQueue
import logging
logger = logging.getLogger(__name__)


class SlurmQueue(SimQueue, ProtocolInterface):
    """ Queue system for SLURM

    Parameters
    ----------
    jobname : str, default=None
        Job name (identifier)
    partition : str or list of str, default=None
        The queue (partition) or list of queues to run on. If list, the one offering earliest initiation will be used.
    priority : str, default='gpu_priority'
        Job priority
    ngpu : int, default=1
        Number of GPUs to use for a single job
    ncpu : int, default=1
        Number of CPUs to use for a single job
    memory : int, default=1000
        Amount of memory per job (MiB)
    gpumemory : int, default=None
        Only run on GPUs with at least this much memory. Needs special setup of SLURM. Check how to define gpu_mem on
        SLURM.
    walltime : int, default=None
        Job timeout (s)
    envvars : str, default='ACEMD_HOME,HTMD_LICENSE_FILE'
        Envvars to propagate from submission node to the running node (comma-separated)
    prerun : list of strings, default=None
        Shell commands to execute on the running node before the job (e.g. loading modules).
    mailtype : str, default=None
        When to send emails. Separate options with commas like 'END,FAIL'.
    mailuser : str, default=None
        User email address.
    outputstream : str, default='slurm.%N.%j.out'
        Output stream.
    errorstream : str, default='slurm.%N.%j.err'
        Error stream.
    datadir : str, default=None
        The path in which to store completed trajectories.
    trajext : str, default='xtc'
        Extension of trajectory files. This is needed to copy them to datadir.
    nodelist : list, default=None
        A list of nodes on which to run every job at the *same time*! Careful! The jobs will be duplicated!
    exclude : list
        A list of nodes on which *not* to run the jobs. Use this to select nodes on which to allow the jobs to run on.

    Examples
    --------
    >>> s = SlurmQueue()
    >>> s.partition = 'multiscale'
    >>> s.submit('/my/runnable/folder/')  # Folder containing a run.sh bash script
    """

    _defaults = {'partition': None, 'priority': None, 'ngpu': 1, 'ncpu': 1, 'memory': 1000, 'walltime': None,
                 'envvars': 'ACEMD_HOME,HTMD_LICENSE_FILE', 'prerun': None}

    def __init__(self, _configapp=None):
        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)
        self._arg('jobname', 'str', 'Job name (identifier)', None, val.String())
        self._arg('partition', 'str', 'The queue (partition) or list of queues to run on. If list, the one offering '
                                      'earliest initiation will be used.',
                  self._defaults['partition'], val.String(), nargs='*')
        self._arg('priority', 'str', 'Job priority', self._defaults['priority'], val.String())
        self._arg('ngpu', 'int', 'Number of GPUs to use for a single job', self._defaults['ngpu'],
                  val.Number(int, '0POS'))
        self._arg('ncpu', 'int', 'Number of CPUs to use for a single job', self._defaults['ncpu'],
                  val.Number(int, 'POS'))
        self._arg('memory', 'int', 'Amount of memory per job (MB)', self._defaults['memory'], val.Number(int, 'POS'))
        self._arg('gpumemory', 'int', 'Only run on GPUs with at least this much memory. Needs special setup of SLURM. '
                                      'Check how to define gpu_mem on SLURM.', None, val.Number(int, '0POS'))
        self._arg('walltime', 'int', 'Job timeout (s)', self._defaults['walltime'], val.Number(int, 'POS'))
        self._cmdDeprecated('environment', 'envvars')
        self._arg('mailtype', 'str', 'When to send emails. Separate options with commas like \'END,FAIL\'.', None,
                  val.String())
        self._arg('mailuser', 'str', 'User email address.', None, val.String())
        self._arg('outputstream', 'str', 'Output stream.', 'slurm.%N.%j.out', val.String())
        self._arg('errorstream', 'str', 'Error stream.', 'slurm.%N.%j.err'), val.String()
        self._arg('datadir', 'str', 'The path in which to store completed trajectories.', None, val.String())
        self._arg('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc',
                  val.String())
        self._arg('nodelist', 'list', 'A list of nodes on which to run every job at the *same time*! Careful! The jobs'
                                      ' will be duplicated!', None, val.String(), nargs='*')
        self._arg('exclude', 'list', 'A list of nodes on which *not* to run the jobs. Use this to select nodes on '
                                     'which to allow the jobs to run on.', None, val.String(), nargs='*')
        self._arg('envvars', 'str', 'Envvars to propagate from submission node to the running node (comma-separated)',
                  self._defaults['envvars'], val.String())
        self._arg('prerun', 'list', 'Shell commands to execute on the running node before the job (e.g. '
                                    'loading modules)', self._defaults['prerun'], val.String(), nargs='*')

        # Load Slurm configuration profile
        slurmconfig = _config['slurm']
        profile = None
        if _configapp is not None:
            if slurmconfig is not None:
                if os.path.isfile(slurmconfig) and slurmconfig.endswith(('.yml', '.yaml')):
                    try:
                        with open(slurmconfig, 'r') as f:
                            profile = yaml.load(f)
                        logger.info('Loaded Slurm configuration YAML file {}'.format(slurmconfig))
                    except:
                        logger.warning('Could not load YAML file {}'.format(slurmconfig))
                else:
                    logger.warning('{} does not exist or it is not a YAML file.'.format(slurmconfig))
                if profile:
                    try:
                        properties = profile[_configapp]
                    except:
                        raise RuntimeError('There is no profile in {} for configuration '
                                           'app {}'.format(slurmconfig, _configapp))
                    for p in properties:
                        self.__dict__[p] = properties[p]
                        logger.info('Setting {} to {}'.format(p, properties[p]))
            else:
                raise RuntimeError('No Slurm configuration YAML file defined for the configapp')
        else:
            if slurmconfig is not None:
                logger.warning('Slurm configuration YAML file defined without configuration app')

        # Find executables
        self._qsubmit = SlurmQueue._find_binary('sbatch')
        self._qinfo = SlurmQueue._find_binary('sinfo')
        self._qcancel = SlurmQueue._find_binary('scancel')
        self._qstatus = SlurmQueue._find_binary('squeue')

    @staticmethod
    def _find_binary(binary):
        ret = shutil.which(binary, mode=os.X_OK)
        if not ret:
            raise FileNotFoundError("Could not find required executable [{}]".format(binary))
        ret = os.path.abspath(ret)
        return ret

    def _createJobScript(self, fname, workdir, runsh):
        from htmd.util import ensurelist
        workdir = os.path.abspath(workdir)
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#\n')
            f.write('#SBATCH --job-name={}\n'.format(self.jobname))
            f.write('#SBATCH --partition={}\n'.format(','.join(ensurelist(self.partition))))
            if self.ngpu != 0:
                f.write('#SBATCH --gres=gpu:{}'.format(self.ngpu))
                if self.gpumemory is not None:
                    f.write(',gpu_mem:{}'.format(self.gpumemory))
                f.write('\n')
            f.write('#SBATCH --cpus-per-task={}\n'.format(self.ncpu))
            f.write('#SBATCH --mem={}\n'.format(self.memory))
            f.write('#SBATCH --priority={}\n'.format(self.priority))
            f.write('#SBATCH --workdir={}\n'.format(workdir))
            f.write('#SBATCH --output={}\n'.format(self.outputstream))
            f.write('#SBATCH --error={}\n'.format(self.errorstream))
            if self.envvars is not None:
                f.write('#SBATCH --export={}\n'.format(self.envvars))
            if self.walltime is not None:
                f.write('#SBATCH --time={}\n'.format(self.walltime))
            if self.mailtype is not None and self.mailuser is not None:
                f.write('#SBATCH --mail-type={}\n'.format(self.mailtype))
                f.write('#SBATCH --mail-user={}\n'.format(self.mailuser))
            if self.nodelist is not None:
                f.write('#SBATCH --nodelist={}\n'.format(','.join(ensurelist(self.nodelist))))
            if self.exclude is not None:
                f.write('#SBATCH --exclude={}\n'.format(','.join(ensurelist(self.exclude))))
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(os.path.normpath(os.path.join(workdir, self._sentinel))))
            f.write('\n')
            if self.prerun is not None:
                for call in self.prerun:
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

        if self.partition is None:
            raise ValueError('The partition needs to be defined.')

        # if all folders exist, submit
        for d in dirs:
            logger.info('Queueing ' + d)

            if self.jobname is None:
                self.jobname = self._autoJobName(d)

            runscript = self._getRunScript(d)
            self._cleanSentinel(d)

            jobscript = os.path.abspath(os.path.join(d, 'job.sh'))
            self._createJobScript(jobscript, d, runscript)
            try:
                ret = check_output([self._qsubmit, jobscript])
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
        if self.partition is None:
            raise ValueError('The partition needs to be defined.')
        if self.jobname is None:
            raise ValueError('The jobname needs to be defined.')
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

    def stop(self):
        """ Cancels all currently running and queued jobs
        """
        import getpass
        if self.partition is None:
            raise ValueError('The partition needs to be defined.')
        if self.jobname is None:
            raise ValueError('The jobname needs to be defined.')
        user = getpass.getuser()
        cmd = [self._qcancel, '-n', self.jobname, '-u', user, '-p', self.partition]
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
    q = SlurmQueue()
    """