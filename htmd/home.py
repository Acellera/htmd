# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import os
import sys
import inspect
import platform


def home(
    dataDir: str | None = None,
    libDir: bool = False,
    shareDir: bool = False,
) -> str:
    """Return the pathname of the HTMD root directory (or a subdirectory).

    Parameters
    ----------
    dataDir : str, optional
        If provided, return the path to the named data subdirectory.
    libDir : bool, optional
        If True, return the path to the platform-specific lib directory.
    shareDir : bool, optional
        If True, return the path to the share directory.

    Returns
    -------
    dir : str
        The requested directory path.

    Raises
    ------
    FileNotFoundError
        If ``libDir`` is True and the lib directory does not exist, or if
        ``shareDir`` is True and the share directory does not exist.

    Examples
    --------
    >>> from htmd.home import home
    >>> home()                                 # doctest: +ELLIPSIS
    '.../htmd'
    >>> home(dataDir="dhfr")                   # doctest: +ELLIPSIS
    '.../data/dhfr'
    >>> os.path.join(home(dataDir="dhfr"),"dhfr.pdb")  # doctest: +ELLIPSIS
    '.../data/dhfr/dhfr.pdb'
    """

    homeDir = os.path.dirname(inspect.getfile(htmd))
    try:
        if sys._MEIPASS:
            homeDir = sys._MEIPASS
    except Exception:
        pass

    if dataDir:
        return os.path.join(homeDir, "data", dataDir)
    elif libDir:
        libdir = os.path.join(homeDir, "lib", platform.system())
        if not os.path.exists(libdir):
            raise FileNotFoundError("Could not find libs.")
        return libdir
    elif shareDir:
        sharedir = os.path.join(homeDir, "share")
        if not os.path.exists(sharedir):
            raise FileNotFoundError("Could not find HTMD share directory.")
        return sharedir
    else:
        return homeDir


if __name__ == "__main__":
    import doctest

    doctest.testmod()

    h = home()
    print(h)
