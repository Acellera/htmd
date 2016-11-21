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
import shutil
import random
import string
import numpy as np
from subprocess import check_output, CalledProcessError
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_ANY, RANGE_0POS, RANGE_POS
from htmd.queues.simqueue import SimQueue
import logging
logger = logging.getLogger(__name__)


class PBSQueue(SimQueue, ProtocolInterface):
    """ Queue system for PBS

    Parameters
    ----------
    jobname : str, default=None
        Job name (identifier)
    queue : str, default=None
        The queue (partition) to run on
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

    Examples
    --------
    >>> from htmd import *
    >>> s = PBSQueue()
    >>> s.queue = 'multiscale'
    >>> s.submit('/my/runnable/folder/')  # Folder containing a run.sh bash script
    """
    def __init__(self):
        super().__init__()
        self._cmdString('jobname', 'str', 'Job name (identifier)', None)
        self._cmdString('queue', 'str', 'The queue to run on', None)
        self._cmdValue('ngpu', 'int', 'Number of GPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('ncpu', 'int', 'Number of CPUs to use for a single job', 1, TYPE_INT, RANGE_0POS)
        self._cmdValue('memory', 'int', 'Amount of memory per job (MB)', 1000, TYPE_INT, RANGE_0POS)
        self._cmdValue('walltime', 'int', 'Job timeout (s)', None, TYPE_INT, RANGE_POS)
        self._cmdString('environment', 'str', 'Envvars to propagate to the job.', 'ACEMD_HOME,HTMD_LICENSE_FILE')
        self._cmdString('mailtype', 'str', 'When to send emails. Separate options with commas like \'END,FAIL\'.', None)
        self._cmdString('mailuser', 'str', 'User email address.', None)
#        self._cmdString('outputstream', 'str', 'Output stream.', 'pbs.%N.%j.out')
#        self._cmdString('errorstream', 'str', 'Error stream.', 'pbs.%N.%j.err')  # Maybe change these to job name
        self._cmdString('datadir', 'str', 'The path in which to store completed trajectories.', None)
        self._cmdString('trajext', 'str', 'Extension of trajectory files. This is needed to copy them to datadir.', 'xtc')

        # Find executables
        self._qsubmit = PBSQueue._find_binary('qsub')
        self._qcancel = PBSQueue._find_binary('qsdel')

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
            f.write('#PBS -N={}\n'.format(self.jobname))
            f.write('#PBS -lselect=%d:ncpus=%d:ngpus=%d:mem=%dMB\n' % ( 1, self.ncpus, self.ngpus, self.mem ) )
            if( self.queue ):
               f.write("#PBS -q  %s\n" % ( self.queue ) )
            f.write('#PBS -lwalltime=%d:0:0\n' % ( self.walltime ) )
            if self.environment is not None:
                f.write('#PBS -v=%s\n' % ( ",".join( self.environment )) )

            f.write('\ncd {}\n'.format(workdir))
            f.write('{}'.format(runsh))
            f.write('\ntouch .done')

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
        dist : list
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
                self.jobname = os.path.basename(os.path.abspath(d)) + '_' + ''.join([random.choice(string.digits) for _ in range(5)])

            runscript = os.path.abspath(os.path.join(d, 'run.sh'))
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
        inprogress = 0
        for i in self._dirs:
            if not os.path.exists(os.path.join(i, '.done')):
                inprogress += 1
        return inprogress


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



