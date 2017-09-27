# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from abc import ABCMeta, abstractmethod


class SimQueue(metaclass=ABCMeta):

    def __init__(self):
        super().__init__()
        self._sentinel = 'htmd.queues.done'
        # For synchronous
        self._dirs = []

    @abstractmethod
    def retrieve(self):
        """ Subclasses need to implement this method """
        pass

    @abstractmethod
    def submit(self, dirs):
        """ Subclasses need to implement this method """
        pass

    @abstractmethod
    def inprogress(self):
        """ Subclasses need to implement this method """
        pass

    @abstractmethod
    def notcompleted(self):
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
            sys.stdout.flush()
            sleep(5)

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
