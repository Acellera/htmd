# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

class RetrieveError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SubmitError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InProgressError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ProjectNotExistError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SimQueue(metaclass=ABCMeta):

    def __init__(self):
        super().__init__()
        self._sentinel = 'htmd.queues.done'
        # For synchronous
        self._dirs = None

    @abstractmethod
    def retrieve(self):
        """ Subclasses need to implement this method """
        pass

    @abstractmethod
    def submit(self, dirs):
        """ Subclasses need to implement this method """
        pass

    def _submitinit(self, dirs):
        if isinstance(dirs, str):
            dirs = [dirs, ]
        if self._dirs is None:
            self._dirs = []
        self._dirs.extend(dirs)
        return dirs

    @abstractmethod
    def inprogress(self):
        """ Subclasses need to implement this method """
        pass

    def wait(self, sentinel=False):
        """ Blocks script execution until all queued work completes

        Parameters
        ----------
        sentinel : bool, default=False
            If False, it relies on the queueing system reporting to determine the number of running jobs. If True, it
            relies on the filesystem, in particular on the existence of a sentinel file for job completion.

        Examples
        --------
        >>> self.wait()
        """
        from time import sleep
        import sys

        while (self.inprogress() if not sentinel else self.notcompleted()) != 0:
            self.retrieve()
            sys.stdout.flush()
            sleep(5)

    def notcompleted(self):
        """Returns the sum of the number of job directories which do not have the sentinel file for completion.

        Returns
        -------
        total : int
            Total number of directories which have not completed
        """
        import os
        total = 0
        if self._dirs is None:
            raise RuntimeError('This method relies on running synchronously.')
        for i in self._dirs:
            if not os.path.exists(os.path.join(i, self._sentinel)):
                total += 1
        return total

    def _cleanSentinel(self, d):
        import os
        if os.path.exists(os.path.join(d, self._sentinel)):
            try:
                os.remove(os.path.join(d, self._sentinel))
            except:
                logger.warning('Could not remove {} sentinel from {}'.format(self._sentinel, d))
            else:
                logger.info('Removed existing {} sentinel from {}'.format(self._sentinel, d))

    def _getRunScript(self, d):
        import os
        runscript = os.path.abspath(os.path.join(d, 'run.sh'))
        if not os.path.exists(runscript):
            raise FileExistsError('File {} does not exist.'.format(runscript))
        if not os.access(runscript, os.X_OK):
            raise PermissionError('File {} does not have execution permissions.'.format(runscript))
        return runscript

    @abstractmethod
    def stop(self):
        """ Subclasses need to implement this method """
        pass

    @property
    @abstractmethod
    def ncpu(self):
        """ Subclasses need to have this property """
        pass

    @ncpu.setter
    @abstractmethod
    def ncpu(self, value):
        """ Subclasses need to have this setter """
        pass

    @property
    @abstractmethod
    def ngpu(self):
        """ Subclasses need to have this property """
        pass

    @ngpu.setter
    @abstractmethod
    def ngpu(self, value):
        """ Subclasses need to have this setter """
        pass

    @property
    @abstractmethod
    def memory(self):
        """ Subclasses need to have this property. This property is expected to return a integer in MiB"""
        pass

    @memory.setter
    @abstractmethod
    def memory(self, value):
        """ Subclasses need to have this setter """
        pass
