# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import os
from htmd.queues.simqueue import SimQueue
from protocolinterface import ProtocolInterface, val
import logging
import random
import string

logger = logging.getLogger(__name__)
logging.getLogger("botocore").setLevel(logging.CRITICAL)
logging.getLogger("boto3").setLevel(logging.CRITICAL)


class AceCloudQueue(SimQueue, ProtocolInterface):
    """  Class which handles the acecloud queueing.

    Before using this class you should have set up acecloud already by calling `acecloud --setup` from command line.

    Parameters
    ----------
    groupname : str, default=None
        The name of the group of simulations you want to submit. If none is given, a randomly generated string will be used instead.
    datadir : str, default=None
        The directory in which to retrieve your results.
    verbose : bool, default=False
        Turn verbosity mode on or off.

    Examples
    --------
    >>> ac = AceCloudQueue()
    >>> ac.groupname = 'mygroup'
    >>> ac.submit(glob('./mysimfolders/*/'))
    >>> ac.datadir = './mysimfolders/'  # The directory where the output will be stored. Here we use the same as the input
    >>> ac.retrieve()
    >>> ac.inprogress()
    >>> ac.stop() # Kill all simulations of the group

    """
    def __init__(self):

        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)
        self._arg('groupname', 'str', 'The name of the group of simulations you want to submit. If none is given, '
                                      'a randomly generated string will be used instead.', None, val.String())
        self._arg('datadir', 'str', 'The directory in which to retrieve your results.', None, val.String())
        self._arg('instancetype', 'str', 'Instance type', 'g2.2xlarge', val.String(), valid_values=('g2.2xlarge', 'r4.large', 'p2.xlarge'))
        self._arg('hashnames', 'bool', 'If True, each job will have a name created from the hash of its directory '
                                       'instead of using the directory name.', False, val.Boolean())
        self._arg('verbose', 'bool', 'Turn verbosity mode on or off.', False, val.Boolean())
        self._arg('ngpu', 'int', 'Number of GPUs to use for a single job', 0, val.Number(int, '0POS'))
        self._arg('ncpu', 'int', 'Number of CPUs to use for a single job', 1, val.Number(int, '0POS'))
        self._arg('memory', 'int', 'Amount of memory per job (MB)', 8000, val.Number(int, '0POS'))

        self._cloud = None

    def _createCloud(self):
        if self._cloud is not None:
            return
        from acecloud.cloud import Cloud
        self._cloud = Cloud(verbose=self.verbose)

    def submit(self, dirs):
        self._createCloud()

        from acecloud.job import Job
        dirs = self._submitinit(dirs)

        if self.groupname is None:
            self.groupname = ''.join(
                random.SystemRandom().choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in
                range(8))

        for d in dirs:
            if not os.path.isdir(d):
                raise FileExistsError('Submit: directory ' + d + ' does not exist.')
            # runsh = os.path.join(d, 'run.sh')
            # if not os.path.exists(runsh):
            #     raise FileExistsError("File %s does not exist." % runsh)
            # if not os.access(runsh, os.X_OK):
            #     raise FileExistsError("File %s does not have execution permissions." % runsh)

        for d in dirs:
            logger.info("Queueing " + d)
            name = os.path.basename(d)
            if self.hashnames:
                import hashlib
                name = hashlib.sha256(os.path.abspath(d).encode('utf-8')).hexdigest()[:10]

            runscript = self._getRunScript(d)
            self._cleanSentinel(d)

            jobscript = os.path.abspath(os.path.join(d, 'job.sh'))
            self._createJobScript(jobscript, d, runscript)

            from acecloud.requesttype import RequestType
            Job(
                ngpus=1,
                ncpus=1,
                path=d,
                group=self.groupname,
                name=name,
                cloud=self._cloud,
                requesttype=RequestType.SPOT,
                instance_type=self.instancetype,
                verbose=self.verbose
            )

    def inprogress(self):
        self._createCloud()
        from acecloud.status import Status
        jj = self._cloud.getJobs(group=self.groupname)
        count = 0
        for j in jj:
            s = j.status()
            if s == Status.RUNNING or s == Status.PENDING:
                count += 1

        return count

    def retrieve(self):
        self._createCloud()
        from acecloud.status import Status, CloudError
        jj = self._cloud.getJobs(group=self.groupname)
        currdir = os.getcwd()
        for j in jj:
            if self.datadir is not None:
                outf = os.path.join(self.datadir, j.name)
            else:
                outf = j.path
            if not os.path.exists(outf):
                os.makedirs(outf)
            os.chdir(outf)
            try:
                if j.status() == Status.COMPLETED:
                    j.retrieve(directory=outf)  # Duplicate code to avoid race condition with slow retrieve
                    j.delete()
                else:
                    j.retrieve(directory=outf)
            except CloudError as e:
                pass
            except Exception as e:
                logger.warning(e)
                pass
            os.chdir(currdir)

    def stop(self):
        # TODO: This not only stops the job, but also deletes the S3. Not exactly like the stop of other queues
        self._createCloud()
        jj = self._cloud.getJobs(group=self.groupname)
        for j in jj:
            try:
                j.delete()
            except Exception as e:
                logger.warning(e)
                pass

    def _createJobScript(self, fname, workdir, runsh):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('\n')
            # Trap kill signals to create sentinel file
            f.write('\ntrap "touch {}" EXIT SIGTERM\n'.format(self._sentinel))
            f.write('\n')
            f.write('bash run.sh')
        os.chmod(fname, 0o700)

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
