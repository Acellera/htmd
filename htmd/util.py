# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import tempfile
import logging
import requests
import os

logger = logging.getLogger(__name__)


def tempname(suffix='', create=False):
    if create:
        file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    else:
        file = tempfile.NamedTemporaryFile(delete=True, suffix=suffix)
    file.close()
    return file.name


def opm(pdb):
    """Download a molecule from the OPM.

    Removes DUM atoms.

    Parameters
    ----------
    pdb: str
        The 4-letter code

    Returns
    -------
    mol: Molecule
        The oriented molecule

    thickness: float
        The bilayer thickness (both layers)

    Examples
    --------
    >>> mol, thickness = opm("1z98")
    >>> mol.numAtoms
    10272
    >>> thickness
    28.2

    """

    from htmd.molecule.support import string_to_tempfile
    from htmd.molecule.molecule import Molecule
    # http://opm.phar.umich.edu/pdb/1z98.pdb
    r = requests.get("http://opm.phar.umich.edu/pdb/{:s}.pdb".format(pdb.lower()))

    if r.status_code != 200:
        raise NameError('PDB code not found in the OPM database')

    tempfile = string_to_tempfile(r.content.decode('ascii'), "pdb")
    mol = Molecule(tempfile)
    mol.filter("not resname DUM")

    # Assuming the half-thickness is the last word in the first line
    # REMARK      1/2 of bilayer thickness:   14.1
    f = open(tempfile)
    h = f.readline()
    f.close()
    os.unlink(tempfile)

    hs = h.split()
    thickness = 2.0*float(hs[-1])

    return mol, thickness


_issuedDeprecationWarnings = {}


def issueDeprecationWarning(msg):
    """ Issue a deprecation warning, only once.

    Parameters
    ----------
    msg : str
        The message.

    """

    import inspect
    caller = inspect.stack()[1][3]
    if not caller in _issuedDeprecationWarnings:
        _issuedDeprecationWarnings[caller] = 1
        logger.warning(
            "Deprecation warning (%s). The function is obsolete and will be removed in a future version! Additional info: %s" % (
            caller, msg))


def _testDeprecation():
    issueDeprecationWarning("Please ignore this first warning")
    issueDeprecationWarning("This second warning should not appear")


if __name__ == "__main__":
    from htmd import *
    import doctest

    doctest.testmod()

    _testDeprecation()
