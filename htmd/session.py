# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import pickle
import logging

logger = logging.getLogger(__name__)


def htmdsave(fname, varnames=None):
    """Saves the workspace to a file

    Parameters
    ----------
    fname : str
        The file in which to save the workspace
    varnames : list
        A list of variables to save. Default: saves all variables.
    """
    from types import ModuleType, FunctionType
    import __main__ as _main_module
    import inspect

    mydict = dict()
    maindict = _main_module.__dict__.copy()
    bannedvars = ["In", "Out", "exit", "quit", "get_ipython"]

    if varnames is None:
        logger.warning(
            "Save will ignore all variables starting with _ and following variable names {}".format(
                bannedvars
            )
        )
    logger.warning(
        "Save will duplicate data which is stored in various pointers or variables. Take care!"
    )

    if varnames is None:
        for k in maindict:
            if (
                k[0] != "_"
                and k not in bannedvars
                and not isinstance(maindict[k], ModuleType)
                and not isinstance(maindict[k], FunctionType)
                and not inspect.isclass(maindict[k])
            ):
                try:
                    pickle.dumps(maindict[k])  # This probably slows down everything
                    mydict[k] = maindict[k]
                except Exception:
                    logger.debug("Skipping variable: {}".format(k))
    else:
        if not isinstance(varnames, list):
            varnames = [varnames]
        for k in varnames:
            mydict[k] = maindict[k]

    logger.debug("Saving variables: {}".format(list(mydict.keys())))
    f = open(fname, "wb")
    pickle.dump(mydict, f)
    f.close()


def htmdload(fname):
    """Loads a previously saved workspace from a file.

    Parameters
    ----------
    fname : str
        The location of the file containing the saved workspace.
    """
    import __main__ as _main_module

    f = open(fname, "rb")
    mydict = pickle.load(f)
    f.close()

    # Inserting all variables into the main dictionary
    _main_module.__dict__.update(mydict)
