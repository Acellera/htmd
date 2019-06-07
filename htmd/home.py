# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import os
import sys
import inspect
import platform


def home(dataDir=None, libDir=False, shareDir=False):
    """Return the pathname of the HTMD root directory (or a data subdirectory).

    Parameters
    ----------
    dataDir : str
        If not None, return the path to a specific data directory
    libDir : bool
        If True, return path to the lib directory

    Returns
    -------
    dir : str
        The directory

    Example
    -------
    >>> from htmd.home import home
    >>> home()                                 # doctest: +ELLIPSIS
    '.../htmd'
    >>> home(dataDir="dhfr")                   # doctest: +ELLIPSIS
    '.../data/dhfr'
    >>> os.path.join(home(dataDir="dhfr"),"dhfr.pdb")  # doctest: +ELLIPSIS
    '.../data/dhfr/dhfr.pdb'
    """

    homeDir=os.path.dirname(inspect.getfile(htmd))
    try:
      if sys._MEIPASS:
         homeDir = sys._MEIPASS
    except:
      pass

    if dataDir:
        return os.path.join(homeDir, 'data', dataDir)
    elif libDir:
        libdir = os.path.join(homeDir, 'lib', platform.system())
        if not os.path.exists(libdir):
            raise FileNotFoundError('Could not find libs.')
        return libdir
    elif shareDir:
        sharedir = os.path.join(homeDir, 'share')
        if not os.path.exists(sharedir):
            raise FileNotFoundError('Could not find HTMD share directory.')
        return sharedir
    else:
        return homeDir


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    h = home()
    print(h)
