# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import abc
from abc import ABCMeta


class SimQueue(metaclass=ABCMeta):

    @abc.abstractmethod
    def retrieve(self):
        """ Subclasses need to implement this method """
        pass

    @abc.abstractmethod
    def submit(self, dirs):
        """ Subclasses need to implement this method """
        pass

    @abc.abstractmethod
    def inprogress(self):
        """ Subclasses need to implement this method """
        pass

    def wait(self):
        """ Blocks script execution until all queued work completes

        Examples
        --------
        >>> queue.wait()
        """
        from time import sleep
        import sys
        while self.inprogress() != 0:
            sys.stdout.flush()
            sleep(1)

    @abc.abstractmethod
    def stop(self):
        """ Subclasses need to implement this method """
        pass