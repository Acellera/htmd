# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import os
import inspect


def home(dataDir=None):
    """Return the pathname of the HTMD root directory (or a data subdirectory).

    Parameters
    ----------
    dataDir : str
        If not None, return the path to a specific data directory

    Example
    -------
        >>> htmd.home()                                 # doctest: +ELLIPSIS
        '.../htmd/htmd'
        >>> htmd.home(dataDir="dhfr")                   # doctest: +ELLIPSIS
        '.../htmd/data/dhfr'
        >>> os.path.join(htmd.home(dataDir="dhfr"),"dhfr.pdb")  # doctest: +ELLIPSIS
        '.../htmd/data/dhfr/dhfr.pdb'

    """

    homeDir=os.path.dirname(inspect.getfile(htmd))
    if dataDir:
        return os.path.join(homeDir,"data",dataDir)
    else:
        return homeDir


#Don't know how to do this
# def modulehome(modname):
#    return os.path.dirname(os.path.dirname(inspect.getfile(modname)))

if __name__ == "__main__":
    h = htmd.home()
    print(h)
