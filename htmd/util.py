# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import tempfile
import logging
import os
import sys
import numpy as np


logger = logging.getLogger(__name__)


def _getNjobs():
    from htmd.config import _config
    njobs = _config['njobs']
    if njobs < 0:
        import multiprocessing
        njobs = multiprocessing.cpu_count() + njobs + 1
    return njobs


def tempname(suffix='', create=False):
    if create:
        file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    else:
        file = tempfile.NamedTemporaryFile(delete=True, suffix=suffix)
    file.close()
    return file.name


def ensurelist(tocheck, tomod=None):
    """Convert np.ndarray and scalars to lists.

    Lists and tuples are left as is. If a second argument is given,
    the type check is performed on the first argument, and the second argument is converted.
    """
    if tomod is None:
        tomod = tocheck
    if isinstance(tocheck, np.ndarray):
        return list(tomod)
    if isinstance(tocheck, range):
        return list(tocheck)
    if not isinstance(tocheck, list) and not isinstance(tocheck, tuple):
        return [tomod, ]
    return tomod


def diffMolecules(mol1, mol2, sel=None):
    """Check that name, resname, resid, insertion codes match between two molecules.

    Coordinates are not checked.

    Parameters
    ----------
    mol1 : Molecule
        first structure to compare
    mol2 : Molecule
        second structure to compare
    sel: str
        compare only after filtering to this Atom selection string.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__

    Returns
    -------
    diff: list
        a list of differences, as human-readable strings (empty if structures are equal).

    Examples
    --------
    >>> m=Molecule("3PTB")
    >>> m2=m.copy()
    >>> m2.set("resname","HIE","resid 91")
    >>> diffMolecules(m,m2,sel="name CA")
    ['CA   HIS    91     vs   CA   HIE    91  ']
    """

    from numpy import logical_and as l_and

    diff = []

    m1 = mol1.copy()
    m1.filter(sel)
    m2 = mol2.copy()
    m2.filter(sel)

    eq_name = np.equal(m1.name, m2.name)
    eq_resid = np.equal(m1.resid, m2.resid)
    eq_resname = np.equal(m1.resname, m2.resname)
    eq_insertion = np.equal(m1.insertion, m2.insertion)

    eq_all = l_and(eq_name,
                   l_and(eq_resid,
                         l_and(eq_resname, eq_insertion)))

    neq = np.logical_not(eq_all).nonzero()
    for i in neq[0]:
        diff.append("{:4s} {:4s} {:4d} {:1s}   vs   {:4s} {:4s} {:4d} {:1s}".format(
            m1.name[i], m1.resname[i], m1.resid[i], m1.insertion[i],
            m2.name[i], m2.resname[i], m2.resid[i], m2.insertion[i],
        ))
    return diff


def getPdbStrings(mol, sel=None, onlyAtom=True):
    """Return the PDB corresponding to molecule and selection, as a list of strings.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object
    sel : str
        Atom selection string for what to be outputted.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    onlyAtom : bool
        Only return ATOM/HETATM records (default True)

    Examples
    --------
    >>> m = Molecule("3PTB")
    >>> getPdbStrings(m, "resname BEN")         # doctest: +NORMALIZE_WHITESPACE
    ['HETATM    1  C1  BEN A   1      -1.853  14.311  16.658  1.00 19.86      1    C  ',
     'HETATM    2  C2  BEN A   1      -2.107  15.653  16.758  1.00 19.86      1    C  ',
     'HETATM    3  C3  BEN A   1      -1.774  16.341  17.932  1.00 19.86      1    C  ',
     'HETATM    4  C4  BEN A   1      -1.175  15.662  19.005  1.00 19.86      1    C  ',
     'HETATM    5  C5  BEN A   1      -0.914  14.295  18.885  1.00 19.86      1    C  ',
     'HETATM    6  C6  BEN A   1      -1.257  13.634  17.708  1.00 19.86      1    C  ',
     'HETATM    7  C   BEN A   1      -2.193  13.627  15.496  1.00 19.86      1    C  ',
     'HETATM    8  N1  BEN A   1      -2.797  14.235  14.491  1.00 19.86      1    N  ',
     'HETATM    9  N2  BEN A   1      -1.762  12.391  15.309  1.00 19.86      1    N  ']
    """

    f = tempfile.NamedTemporaryFile(suffix=".pdb")
    mol.write(f.name, sel)
    f.seek(0, 0)
    r = f.read()
    f.close()

    rl = r.decode("ascii").split("\n")

    if onlyAtom:
        rl = [x for x in rl if (x.startswith("ATOM") or x.startswith("HETATM"))]

    return rl


def assertSameAsReferenceDir(compareDir, outdir="."):
    """Check if files in refdir are present in the directory given as second argument AND their content matches.

    Raise an exception if not."""

    import filecmp
    import os

    toCompare = os.listdir(compareDir)
    match, mismatch, error = filecmp.cmpfiles(outdir, compareDir, toCompare, shallow=False)
    if len(mismatch) != 0 or len(error) != 0 or len(match) != len(toCompare):
        logger.error("Mismatch while checking directory {} versus reference {}".format(outdir,compareDir))
        logger.error("Files being checked: {}".format(toCompare))
        for f in mismatch:
            logger.error("    diff {} {}".format(os.path.join(outdir, f),
                                                 os.path.join(compareDir, f)   ))
        raise Exception('Mismatch in regression testing.')


def testDHFR():
    import conda
    import shutil
    from jobqueues.localqueue import LocalGPUQueue

    dhfrdir = os.path.abspath(os.path.join(conda.__file__, '../../../../../acemd-examples/dhfr'))

    if not os.path.isdir(dhfrdir):
        raise logger.error('Could not find acemd-examples directory. Do `conda install acemd-examples -c acellera`')

    tmpdir = tempname()
    print(tmpdir)
    shutil.copytree(dhfrdir, tmpdir)

    runsh = os.path.join(tmpdir, 'run.sh')
    with open(runsh, 'w') as f:
        f.write('#!/bin/bash\nacemd >log.txt 2>&1')
    os.chmod(runsh, 0o700)

    logger.info('Starting to run the DHFR test')
    md = LocalGPUQueue()
    try:
        md.submit(tmpdir)
        md.wait()
    except:
        logger.error('Some error occurred. Check {}'.format(tmpdir))
    else:
        shutil.rmtree(tmpdir)
        logger.info('Successfully ran the DHRF test. Temporary dir cleaned.')
    return


if __name__ == "__main__":
    from moleculekit.molecule import Molecule
    import doctest

    sys.exit(doctest.testmod().failed)
