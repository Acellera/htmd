# (c) 2015-2017 Acellera Ltd http://www.acellera.com
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




class RBoincPsi4Queue(SimQueue, ProtocolInterface):
    """ Queue system for RBoinc systems (e.g. GPUGRID).
    
    Only supports psi4.

    Parameters
    ----------
    jobname : str, default=None
        Job name (RBoinc group)
    priority : str, default='medium'
        Job priority
    memory : int, default=1000
        Amount of memory per job (MiB)
    url : str
        RBoinc contact URL
    app : str, default='QC'
        RBoinc application ID

    Examples
    --------
    To Do
    """

    _defaults = {'priority': 'medium',
                 'url': None,
                 'app': 'QC' }

    # Copy some parameters from acemdboinc.
    _prioritydefaults = {'low': 200, 'normal': 300, 'high': 400}


    def __init__(self, _configapp=None):
        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)

        self._arg('jobname', 'str', 'Job name (RBoinc group)', None, val.String())
        self._arg('priority', 'str', 'Job priority', self._defaults['priority'], val.String())
        self._arg('url', 'str', 'RBoinc contact URL', self._defaults['url'], val.String())
        self._arg('app', 'str', 'RBoinc application ID', self._defaults['app'], val.String())

        self.ncpu = 0
        self.memory = 1024

        # Find executables
        self._submit_pl = _find_binary('rboinc_submit.pl')
        self._retrieve_pl = _find_binary('rboinc_submit.pl')



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


    def inprogress(self):
        """ Calculates the number of simulations running on BOINC

        Returns
        -------
        inpr : int
            The number of simulations of the `project` currently running on the BOINC server

        """
        cmd = [self._retrieve_pl, '-group', self.project, '-url', self._getURL()]
        if self.dryrun:
            print(cmd)
            output = 0
        else:
            try:
                output = check_output(cmd)
            except CalledProcessError:
                output = 0
        # TODO: Capture  error code and nan/inf cases
        return int(output)


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



def _find_binary(binary):
    ret = shutil.which(binary, mode=os.X_OK)
    if not ret:
        raise FileNotFoundError("Could not find required executable [{}]".format(binary))
    ret = os.path.abspath(ret)
    return ret


if __name__ == "__main__":
    # TODO: Create fake binaries for instance creation testing
    """
    q = SlurmQueue()
    """