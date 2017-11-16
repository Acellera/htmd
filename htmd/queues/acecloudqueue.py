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
        self._arg('verbose', 'bool', 'Turn verbosity mode on or off.', False, val.Boolean())

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
        if self.datadir is None:
            raise RuntimeError('Please set the datadir parameter before retrieving jobs.')

        self._createCloud()
        from acecloud.status import Status
        jj = self._cloud.getJobs(group=self.groupname)
        currdir = os.getcwd()
        for j in jj:
            if j.status() == Status.COMPLETED:
                try:
                    outf = os.path.join(self.datadir, j.name)
                    if not os.path.exists(outf):
                        os.makedirs(outf)
                    os.chdir(outf)
                    j.retrieve()
                    j.delete()
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

    @property
    def ngpu(self):
        raise NotImplementedError

    @ngpu.setter
    def ngpu(self, value):
        raise NotImplementedError

    @property
    def ncpu(self):
        raise NotImplementedError

    @ncpu.setter
    def ncpu(self, value):
        raise NotImplementedError

    @property
    def memory(self):
        raise NotImplementedError

    @memory.setter
    def memory(self, value):
        raise NotImplementedError
