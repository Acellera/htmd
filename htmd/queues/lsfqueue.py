# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
import random
import string
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
    version : [9, 10], int, default=9
        LSF major version
    jobname : str, default=None
        Job name (identifier)
    queue : list, default=None
        The queue or list of queues to run on. If list, it attempts to submit the job to the first queue listed
    app : str, default=None
        The application profile
    ngpu : int, default=1
        Number of GPUs to use for a single job
    gpu_options : dict, default=None
        Number of GPUs to use for a single job(dict format: {['mode': ['shared', 'exclusive_process'],]['mps': ['yes',
        'no'],]['j_exclusive': ['yes', 'no']]})
    ncpu : int, default=1
        Number of CPUs to use for a single job
    memory : int, default=4000
        Amount of memory per job (MiB)
    walltime : int, default=None
        Job timeout (hour:min or min)
    resources : list, default=None
        Resources of the queue
    outputstream : str, default='lsf.%J.out'
        Output stream.
    errorstream : str, default='lsf.%J.err'
        Error stream.
    datadir : str, default=None
        The path in which to store completed trajectories.
    trajext : str, default='xtc'
        Extension of trajectory files. This is needed to copy them to datadir.
    envvars : str, default='ACEMD_HOME'
        Envvars to propagate from submission node to the running node (comma-separated)
    prerun : list, default=None
        Shell commands to execute on the running node before the job (e.g. loading modules)


    Examples
    --------
    >>> s = LsfQueue()
    >>> s.jobname = 'simulation1'
    >>> s.queue = 'multiscale'
    >>> s.submit('/my/runnable/folder/')  # Folder containing a run.sh bash script
    """

    _defaults = {'version': 9, 'queue': None, 'app': None, 'gpu_queue': None, 'cpu_queue': None,
                 'ngpu': 1, 'gpu_options': None, 'ncpu': 1, 'memory': 4000, 'walltime': None, 'resources': None,
                 'envvars': 'ACEMD_HOME', 'prerun': None}

    def __init__(self, _configapp=None):
        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)
        self._arg('version', 'int', 'LSF major version', self._defaults['version'], valid_values=[9, 10])
        self._arg('jobname', 'str', 'Job name (identifier)', None, val.String())
        self._arg('queue', 'list', 'The queue or list of queues to run on. If list, it attempts to submit the job to '
                                   'the first queue listed', self._defaults['queue'], val.String(), nargs='*')
        self._arg('app', 'str', 'The application profile', self._defaults['app'], val.String())
        self._arg('ngpu', 'int', 'Number of GPUs to use for a single job', self._defaults['ngpu'],
                  val.Number(int, '0POS'))
        self._arg('gpu_options', 'dict', 'Number of GPUs to use for a single job', self._defaults['gpu_options'],
                  val.Dictionary(key_type=str, valid_keys=['mode', 'mps', 'j_exclusive'], val_type=str,
                                 valid_vals={'mode': ['shared', 'exclusive_process'],
                                             'mps': ['yes', 'no'], 'j_exclusive': ['yes', 'no']}
                                 )
                  )
        self._arg('ncpu', 'int', 'Number of CPUs to use for a single job', self._defaults['ncpu'],
                  val.Number(int, '0POS'))
        self._arg('memory', 'int', 'Amount of memory per job (MiB)', self._defaults['memory'], val.Number(int, '0POS'))
        self._arg('walltime', 'int', 'Job timeout (hour:min or min)', self._defaults['walltime'],
                  val.Number(int, '0POS'))
        self._arg('resources', 'list', 'Resources of the queue', self._defaults['resources'], val.String(), nargs='*')
        self._cmdDeprecated('environment', 'prerun')
        self._arg('outputstream', 'str', 'Output stream.', 'lsf.%J.out', val.String())
        self._arg('errorstream', 'str', 'Error stream.', 'lsf.%J.err', val.String())
        self._arg('datadir', 'str', 'The path in which to store completed trajectories.', None, val.String())
        self._arg('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc',
                  val.String())
        self._arg('envvars', 'str', 'Envvars to propagate from submission node to the running node (comma-separated)',
                  self._defaults['envvars'], val.String())
        self._arg('prerun', 'list', 'Shell commands to execute on the running node before the job (e.g. '
                                    'loading modules)', self._defaults['prerun'], val.String(), nargs='*')

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
        from htmd.util import ensurelist
        workdir = os.path.abspath(workdir)
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#\n')
            f.write('#BSUB -J {}\n'.format(self.jobname))
            f.write('#BSUB -q "{}"\n'.format(' '.join(ensurelist(self.queue))))
            f.write('#BSUB -n {}\n'.format(self.ncpu))
            if self.app is not None:
                f.write('#BSUB -app {}\n'.format(self.app))
            if self.ngpu != 0:
                if self.version == 9:
                    if self.gpu_options is not None:
                        logger.warning('gpu_options argument was set while it is not needed for LSF version 9')
                    f.write('#BSUB -R "select[ngpus>0] rusage[ngpus_excl_p={}]"\n'.format(self.ngpu))
                elif self.version == 10:
                    if not self.gpu_options:
                        self.gpu_options = {'mode': 'shared', 'mps': 'no', 'j_exclusive': 'no'}
                    gpu_requirements = list()
                    gpu_requirements.append('num={}'.format(self.ngpu))
                    for i in self.gpu_options:
                        gpu_requirements.append('{}={}'.format(i, self.gpu_options[i]))
                    f.write('#BSUB -gpu "{}"\n'.format(':'.join(gpu_requirements)))
                else:
                    raise AttributeError('Version not supported')
            f.write('#BSUB -M {}\n'.format(self.memory))
            f.write('#BSUB -cwd {}\n'.format(workdir))
            f.write('#BSUB -outdir {}\n'.format(workdir))
            f.write('#BSUB -o {}\n'.format(self.outputstream))
            f.write('#BSUB -e {}\n'.format(self.errorstream))
            if self.envvars is not None:
                f.write('#BSUB --env {}\n'.format(self.envvars))
            if self.walltime is not None:
                f.write('#BSUB -W {}\n'.format(self.walltime))
            if self.resources is not None:
                for resource in self.resources:
                    f.write('#BSUB -R "{}"\n'.format(resource))
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

        if self.queue is None:
            raise ValueError('The queue needs to be defined.')

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
