# (c) 2015-2019 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import json
import logging
import os
import shutil
import sys
import tempfile
import time
import zipfile

from jobqueues.simqueue import SimQueue
from protocolinterface import ProtocolInterface
from protocolinterface.validators import Number, String

logger = logging.getLogger(__name__)

class PlayQueue(SimQueue, ProtocolInterface):

    def __init__(self):

        SimQueue.__init__(self)
        ProtocolInterface.__init__(self)

        self._arg('ngpu', 'int', 'Number of GPUs', default=0, validator=Number(int, '0POS'))
        self._arg('ncpu', 'int', 'Number of CPUs', default=1, validator=Number(int, '0POS'))
        self._arg('memory', 'int', 'Amount of memory (MB)', default=1000, validator=Number(int, 'POS'))
        self._arg('max_jobs', 'int', 'Maximum number of concurent jobs', default=sys.maxsize, validator=Number(int, 'POS'))

        self._arg('token', 'str', 'PM token', required=True, validator=String())
        self._arg('app', 'str', 'App name', required=True, validator=String())

        self._dirs = {}

    def _getSession(self):

        from playmolecule import Session

        return Session(self.token)

    @staticmethod
    def _getCurrentJobID():
        # TODO this could be part of PM API

        configFile = '/data/in/config'

        if not os.path.exists(configFile):
            return None

        with open(configFile) as fd:
            return json.load(fd)['execid']

    def _makeZIP(self, dir_):

        zipFile = os.path.join(dir_, 'input.zip')

        files = os.listdir(dir_)
        with zipfile.ZipFile(zipFile, 'w') as zf:
            for file in files:
                zf.write(os.path.join(dir_, file), arcname=file)

        return zipFile

    def submit(self, dirs):

        logger.info(f'Maximum number of concurrent jobs: {self.max_jobs}')

        dirs = [dirs, ] if isinstance(dirs, str) else dirs
        for dir_ in dirs:

            # Delay submission
            while self.inprogress() >= self.max_jobs:
                time.sleep(5)

            # Submit a job
            job = self._getSession().startApp(self.app)
            job.input = self._makeZIP(dir_) # TODO this is Psi4 specific
            job.submit(_logger=False, childOf=self._getCurrentJobID())
            self._dirs[job._execid] = dir_
            logger.info(f'Submitted job {job._execid}:')
            logger.info(f'    App name: {self.app}')
            logger.info(f'    Input directory: {dir_}')
            logger.info(f'    Resources: {self.ngpu} GPUs, {self.ncpu} CPUs, {self.memory} MB of memory')

    def inprogress(self):

        # TODO updated to use playmolecule.Job.getChildren
        #      see https://github.com/Acellera/htmd/pull/910

        counter = 0
        for jobID in self._dirs:
            job = self._getSession().getJob(id=jobID)
            status = job.getStatus(_logger=False)
            if status in (0, 1, 2, 3, 6, 7): # Queuing, running, etc.
                counter += 1
            elif status in (4, 5): # Completed or errored
                pass
            else:
                raise ValueError('Unknown job status')

        return counter

    def retrieve(self):

        for jobID, dir_ in list(self._dirs.items()):
            job = self._getSession().getJob(id=jobID)
            with tempfile.TemporaryDirectory() as tmpDir:
                status = job.getStatus(_logger=False)
                if status == 4:
                    logger.info(f'Job {jobID} completed (status: {status}):')
                    outDir = job.retrieve(_logger=False, path=tmpDir)
                    for file in os.listdir(outDir):
                        shutil.copy(os.path.join(outDir, file), dir_)
                    self._dirs.pop(jobID)
                    logger.info(f'    Retrieved results to {dir_}')
                elif status == 5:
                    logger.info(f'Job {jobID} failed (status: {status})')
                    self._dirs.pop(jobID)
                elif status in (0, 1, 2, 3, 6, 7):
                    pass
                else:
                    raise ValueError('Unknown job status')

    def stop(self):
        raise NotImplemented()

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
