# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import tempfile
import logging
import requests
import io

import os
import numpy as np

logger = logging.getLogger(__name__)


def tempname(suffix='', create=False):
    if create:
        file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    else:
        file = tempfile.NamedTemporaryFile(delete=True, suffix=suffix)
    file.close()
    return file.name


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
        compare after filtering with the given atomselection

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
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object
    sel : str
        Atom selection to be output.
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


def view3DHull(m, style="", preamble="", solid=False):
    """Display the convex hull of the given molecule.

    For preamble and color, see http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node128.html

    The function returns an ID. To delete the objects, send the "htmd_delete_graphics <ID>"
    command to VMD (uninmplemented).

    Parameters
    ----------
    m: Molecule
        The object of which to show the hull (only 1 frame)
    style: str
        Style for wireframe lines
    preamble: str
        Commands (material, color) to be prefixed to the output.
        E.g.: "draw color red; graphics top materials on; graphics top material Transparent".
        Note that this affects later preamble-less commands. See
    solid: bool
        Solid or wireframe

    Examples
    --------
    >> m=Molecule("3PTB")
    >> m.view()
    >> m.filter("protein ")
    >> view3DHull(m)

    """
    from scipy.spatial import ConvexHull
    from htmd.vmdviewer import getCurrentViewer

    if m.coords.shape[2] != 1:
        raise Exception("Only one frame is supported")

    r = io.StringIO()
    r.write(preamble + "\n")

    cc = m.coords[:, :, 0]
    hull = ConvexHull(cc)
    for i in range(hull.nsimplex):
        v1, v2, v3 = hull.simplices[i,]
        c1 = cc[v1, :]
        c1s = str(c1).strip('[]')
        c2 = cc[v2, :]
        c2s = str(c2).strip('[]')
        c3 = cc[v3, :]
        c3s = str(c3).strip('[]')
        if solid:
            r.write("draw triangle {{ {:s} }} {{ {:s} }} {{ {:s} }}\n".format(c1s, c2s, c3s))
        else:
            r.write("draw line {{ {:s} }} {{ {:s} }} {:s} \n".format(c1s, c2s, style))
            r.write("draw line {{ {:s} }} {{ {:s} }} {:s} \n".format(c1s, c3s, style))
            r.write("draw line {{ {:s} }} {{ {:s} }} {:s} \n".format(c2s, c3s, style))
            r.write("\n")

    vmdcmd = r.getvalue()

    vmd = getCurrentViewer()
    vmd.send(vmdcmd)


def opm(pdb):
    """Download a molecule from the OPM.

    Removes DUM atoms.

    Parameters
    ----------
    pdb: str
        The 4-letter PDB code

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
    7902
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
    thickness = 2.0 * float(hs[-1])

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

    _testDeprecation()

    import doctest

    doctest.testmod()
