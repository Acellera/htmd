# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import abc


class Projection:
    """
    Parent class for all trajectory projecting classes. Defines abstract functions.
    """

    @abc.abstractmethod
    def project(self, *args, **kwargs):
        """ Subclasses need to implement and overload this method """
        return

    @abc.abstractmethod
    def getMapping(self, mol):
        return

    @abc.abstractmethod
    def _precalculate(self, mol):
        return