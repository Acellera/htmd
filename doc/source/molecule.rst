Molecule
========

The Molecule class is a central object in HTMD. Most of HTMD functionalities are implemented via Molecules.
A Molecule can actually contain a molecular System (e.g. read from a PDB file), which can of course be composed of
several independent "chemical" molecules (e.g. protein + solvent + ions). This nomenclature is drawn from VMD molecules
and, in order to define "chemical" molecules, one uses segment IDs.

Molecule can read many input formats (PDB, PSF, PRMTOP, etc.) and some trajectory files (XTC, etc). Molecules can be
viewed (VMD or WebGL), aligned, selected, rotated, truncated, appended, and so on.

A very important feature is atom selection. This is identical to `VMD's Atom Selection Language`_, so that it is
possible to verify an atom selection visually and then apply it programmatically.

.. _VMD's Atom Selection Language: http://www.ks.uiuc.edu/Research/vmd/current/ug/node89.html

Contents:

.. toctree::
    :maxdepth: 1

    Molecule module <htmd.molecule.molecule>
