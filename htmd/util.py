# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import tempfile
import logging
import os
import sys
import numpy as np

if TYPE_CHECKING:
    from moleculekit.molecule import Molecule


logger = logging.getLogger(__name__)


def _getNjobs():
    from htmd.config import _config

    njobs = _config["njobs"]
    if njobs < 0:
        import multiprocessing

        njobs = multiprocessing.cpu_count() + njobs + 1
    return njobs


def tempname(suffix: str = "", create: bool = False) -> str:
    """Return a path to a temporary file, optionally creating it on disk.

    Parameters
    ----------
    suffix : str, optional
        File name suffix (including the dot, e.g. ``".pdb"``).
    create : bool, optional
        If True, the file is created on disk. If False, the name is reserved
        but the file is not kept.

    Returns
    -------
    name : str
        Absolute path to the temporary file.
    """
    if create:
        file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    else:
        file = tempfile.NamedTemporaryFile(delete=True, suffix=suffix)
    file.close()
    return file.name


def ensurelist(tocheck, tomod=None) -> list:
    """Convert np.ndarray and scalars to lists.

    Lists and tuples are left as is. If a second argument is given,
    the type check is performed on the first argument, and the second
    argument is converted.

    Parameters
    ----------
    tocheck : object
        The value whose type is inspected to decide whether conversion is needed.
    tomod : object, optional
        The value to convert. If None, ``tocheck`` is converted instead.

    Returns
    -------
    result : list
        ``tomod`` (or ``tocheck``) as a list.
    """
    if tomod is None:
        tomod = tocheck
    if isinstance(tocheck, np.ndarray):
        return list(tomod)
    if isinstance(tocheck, range):
        return list(tocheck)
    if not isinstance(tocheck, list) and not isinstance(tocheck, tuple):
        return [
            tomod,
        ]
    return tomod


def diffMolecules(
    mol1: "Molecule",
    mol2: "Molecule",
    sel: str | np.ndarray | None = None,
) -> list:
    """Check that name, resname, resid, and insertion codes match between two molecules.

    Coordinates are not checked.

    Parameters
    ----------
    mol1 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        First structure to compare.
    mol2 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        Second structure to compare.
    sel : str or np.ndarray, optional
        An atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`). If
        provided, only the selected atoms are compared.

    Returns
    -------
    diff : list
        A list of differences as human-readable strings. Empty if the structures
        are equal.

    Examples
    --------
    >>> m = Molecule("3PTB")
    >>> m2 = m.copy()
    >>> m2.set("resname", "HIE", "resid 91")
    >>> diffMolecules(m, m2, sel="name CA")
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

    eq_all = l_and(eq_name, l_and(eq_resid, l_and(eq_resname, eq_insertion)))

    neq = np.logical_not(eq_all).nonzero()
    for i in neq[0]:
        diff.append(
            "{:4s} {:4s} {:4d} {:1s}   vs   {:4s} {:4s} {:4d} {:1s}".format(
                m1.name[i],
                m1.resname[i],
                m1.resid[i],
                m1.insertion[i],
                m2.name[i],
                m2.resname[i],
                m2.resid[i],
                m2.insertion[i],
            )
        )
    return diff


def getPdbStrings(
    mol: "Molecule",
    sel: str | np.ndarray | None = None,
    onlyAtom: bool = True,
) -> list:
    """Return the PDB corresponding to a molecule and selection as a list of strings.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object.
    sel : str or np.ndarray, optional
        An atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`) for
        what to output.
    onlyAtom : bool, optional
        If True, only ATOM/HETATM records are returned.

    Returns
    -------
    lines : list
        PDB record lines as strings.

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


def assertSameAsReferenceDir(compareDir: str, outdir: str = ".") -> None:
    """Check that files in a reference directory match those in an output directory.

    Subdirectories are compared recursively. Raises an exception if any file is
    missing or its content differs.

    Parameters
    ----------
    compareDir : str
        Path to the reference directory.
    outdir : str, optional
        Path to the output directory to compare against the reference.

    Raises
    ------
    Exception
        If any file is missing, has differing content, or cannot be compared.
    """

    import filecmp
    import os

    def _compare(ref, out, prefix=""):
        entries = os.listdir(ref)
        files = [e for e in entries if os.path.isfile(os.path.join(ref, e))]
        subdirs = [e for e in entries if os.path.isdir(os.path.join(ref, e))]

        match, mismatch, error = filecmp.cmpfiles(out, ref, files, shallow=False)
        bad = [os.path.join(prefix, m) for m in mismatch]
        bad += [os.path.join(prefix, e) for e in error]
        if len(match) != len(files):
            missing = set(files) - set(match) - set(mismatch) - set(error)
            bad += [os.path.join(prefix, m) for m in missing]

        for d in subdirs:
            ref_sub = os.path.join(ref, d)
            out_sub = os.path.join(out, d)
            if not os.path.isdir(out_sub):
                bad.append(os.path.join(prefix, d) + "/")
                continue
            bad += _compare(ref_sub, out_sub, prefix=os.path.join(prefix, d))
        return bad

    bad = _compare(compareDir, outdir)
    if bad:
        raise Exception(
            f"Different results in files {bad}: Compare {outdir} and {compareDir}"
        )


def testDHFR():
    import conda
    import shutil
    from jobqueues.localqueue import LocalGPUQueue

    dhfrdir = os.path.abspath(
        os.path.join(conda.__file__, "../../../../../acemd-examples/dhfr")
    )

    if not os.path.isdir(dhfrdir):
        raise logger.error(
            "Could not find acemd-examples directory. Do `conda install acemd-examples -c acellera`"
        )

    tmpdir = tempname()
    print(tmpdir)
    shutil.copytree(dhfrdir, tmpdir)

    runsh = os.path.join(tmpdir, "run.sh")
    with open(runsh, "w") as f:
        f.write("#!/bin/bash\nacemd >log.txt 2>&1")
    os.chmod(runsh, 0o700)

    logger.info("Starting to run the DHFR test")
    md = LocalGPUQueue()
    try:
        md.submit(tmpdir)
        md.wait()
    except Exception:
        logger.error("Some error occurred. Check {}".format(tmpdir))
    else:
        shutil.rmtree(tmpdir)
        logger.info("Successfully ran the DHRF test. Temporary dir cleaned.")
    return


if __name__ == "__main__":
    from moleculekit.molecule import Molecule
    import doctest

    sys.exit(doctest.testmod().failed)
