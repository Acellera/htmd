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
    requesttype : ('spot', 'ondemand'), str, default='spot'
        Choose between "spot" or "ondemand"
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

        super().__init__()
        self._arg('groupname', 'str', 'The name of the group of simulations you want to submit. If none is given, '
                                      'a randomly generated string will be used instead.', None, val.String())
        self._arg('datadir', 'str', 'The directory in which to retrieve your results.', None, val.String())
        self._arg('requesttype', 'str', 'Choose between "spot" or "ondemand"', 'spot', val.String(), valid_values=('spot', 'ondemand'))
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
        if isinstance(dirs, str):
            dirs = [dirs, ]

        if self.groupname is None:
            self.groupname = ''.join(
                random.SystemRandom().choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in
                range(8))

        for d in dirs:
            if not os.path.isdir(d):
                raise FileExistsError('Submit: directory ' + d + ' does not exist.')
            runsh = os.path.join(d, 'run.sh')
            if not os.path.exists(runsh):
                raise FileExistsError("File %s does not exist." % (runsh))
            if not os.access(runsh, os.X_OK):
                raise FileExistsError("File %s does not have execution permissions." % (runsh))

        for d in dirs:
            logger.info("Queueing " + d)
            runsh = os.path.join(d, 'run.sh')
            jobsh = os.path.join(d, 'job.sh')
            self._createJobScript(jobsh, runsh)
            name = os.path.basename(d)

            Job(
                ngpus=1,
                ncpus=1,
                path=d,
                group=self.groupname,
                name=name,
                cloud=self._cloud,
                requesttype=self.requesttype
            )

    def inprogress(self):
        self._createCloud()
        from acecloud.status import Status
        jj = self._cloud.getJobs(group=self.groupname)
        count = 0
        for j in jj:
            s = j.status()
            if s == Status.RUNNING or s == Status.PENDING:
                count = count + 1

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
        self._createCloud()
        jj = self._cloud.getJobs(group=self.groupname)
        for j in jj:
            try:
                j.delete()
            except Exception as e:
                logger.warning(e)
                pass

    def _createJobScript(self, fname, runsh):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('{}'.format(runsh))

        os.chmod(fname, 0o700)
