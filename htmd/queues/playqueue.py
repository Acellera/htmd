# (c) 2015-2019 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import os
import shutil
import tempfile
import zipfile

from htmd.queues.simqueue import SimQueue
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

        self._arg('token', 'str', 'PM token', required=True, validator=String())
        self._arg('app', 'str', 'App name', required=True, validator=String())

        self._jobIDs = {}

    def _getSession(self):

        from playmolecule import Session

        return Session(self.token)

    def _makeZIP(self, directory):

        zipFile = os.path.join(directory, 'input.zip')

        files = os.listdir(directory)
        with zipfile.ZipFile(zipFile, 'w') as zf:
            for file in files:
                zf.write(os.path.join(directory, file), arcname=file)

        return zipFile

    def submit(self, directories):

        self._dirs = self._submitinit(directories)

        for directory in self._dirs:
            job = self._getSession().startApp(self.app)
            job.input = self._makeZIP(directory)
            job.submit()
            self._jobIDs[directory] = job._execid

    def inprogress(self):

        counter = 0
        for directory in self._dirs:
            job = self._getSession().getJob(id=self._jobIDs[directory])
            status = job.getStatus(_logger=False)
            if status in (0, 1, 2, 3, 6, 7): # Queuing, running, etc.
                counter += 1
            elif status in (4, 5): # Completed or errored
                pass
            else:
                raise ValueError('Unknow job status')

        return counter

    def retrieve(self):

        for directory in self._dirs:
            job = self._getSession().getJob(id=self._jobIDs[directory])
            with tempfile.TemporaryDirectory() as tmpDir:
                outDir = job.retrieve(path=tmpDir)
                if outDir:
                    for file in os.listdir(outDir):
                        shutil.copy(os.path.join(outDir, file), directory)

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
