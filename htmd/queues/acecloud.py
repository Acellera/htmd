# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import os
import pwd
import shutil
from os.path import isdir
from subprocess import check_output
from htmd.protocols.protocolinterface import ProtocolInterface
from htmd.queues.simqueue import SimQueue
from acecloud import Cloud, Job, Status
import logging
import random
import string

logger = logging.getLogger(__name__)
logging.getLogger("botocore").setLevel(logging.CRITICAL)
logging.getLogger("boto3").setLevel(logging.CRITICAL)


class AceCloudQueue(SimQueue, ProtocolInterface):
    def __init__(self, groupname=None):
        super().__init__()
        if not groupname:
            groupname = ''.join(
                random.SystemRandom().choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in
                range(8))

        self._groupname = groupname
        self._cloud = Cloud()
        self._index = 0
        self._jobs = []

    def submit(self, dirs):
        if isinstance(dirs, str):
            dirs = [dirs, ]

        for d in dirs:
            logger.info("Queueing " + d)
            runscript = os.path.abspath(os.path.join(d, "run.sh"))
            if not os.path.exists(runscript):
                raise FileExistsError("File %s does not exist." % (runscript))
            if not os.access(runscript, os.X_OK):
                raise FileExistsError("File %s does not have execution permissions." % (runscript))

            j = Job(
                ngpus=1,
                ncpus=1,
                path=d,
                group=self._groupname,
                name="%6d" % (self._index),
                cloud=self._cloud
            )
            self._index = self._index + 1
            self._jobs.append(j)

    def inprogress(self):
        count = 0
        for j in self._jobs:
            s = j.status()
            if s == Status.RUNNING or s == Status.QUEUED:
                count = count + 1

        return count

    def retrieve(self):
        for j in self._jobs:
            if j.status() == Status.COMPLETED:
                try:
                    j.retrieve()
                    j.delete()
                except:
                    pass

    def stop(self):
        for j in self._jobs:
            try:
                j.stop()
                j.delete()
            except:
                pass
