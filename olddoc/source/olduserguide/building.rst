########
Building
########

Atom Selections in HTMD
=======================

Atom Properties
---------------

A :class:`~htmd.molecule.molecule.Molecule` object contains several atoms. Each atom has several properties:

- They most obviously mirror PDB fields
- Some are computed via heuristics (!)
- Constant properties do not depend on the frame (e.g. type, name, mass)
- Variable properties depend on the frame (e.g. coords[serial,axis,frame])

========= ================================================================
Property  Description
========= ================================================================
`serial`   Unique identifier starting on 1
`name`     Atom name from PDB, unique in the residue
`resname`  3-letter residue type, from PDB
`resid`    The usual residue number, **possibly not unique**
`chain`    Chain identifier, also from PDB
`segid`    Use to define *chemical molecule*: important for system building
`element`  Element of the atom, also from PDB
`coords`   These depend on the trajectory frame
========= ================================================================

Atom Selections: Quick Start
----------------------------

Properties of atoms in a :class:`~htmd.molecule.molecule.Molecule` object can be used to build **boolean expressions**.

================ =======
Example          Meaning
================ =======
`chain A`        true `if` the atom's `chain` property equals to string `A`
`resid 40 to 50` true `if` the atom's `resid` property is between 40 and 50
`mass < 15`      true `if` the atom's `mass` property is inferior to 15
`name CA`        true `if` the atom's `name` is `CA`
`water`          VMD expression which is true `if` the atom's `resname` is H2O, HH0, OHH, HOH, OH2, SOL, WAT, TIP, TIP2, TIP3 or TIP4
================ =======

As each expression is evaluated, it can match (or not) a set of atoms. The rest are **boolean operations**.

============================ =========
Example                      Meaning
============================ =========
`chain A and resid 40 to 50` true `if` **both** the atom's `chain` is `A` and its `resid` is between 40 and 50
`basic or acidic`            true `if` **either** the atom is basic (`resname` ARG HIS LYS) or acidic (`resname` ASP GLU)
============================ =========

Therefore, a selection expression yields an atom selection, i.e. a subset of the original set.
These subsets can then be manipulated, e.g. :func:`~htmd.molecule.molecule.Molecule.get` and
:func:`~htmd.molecule.molecule.Molecule.set` methods of :class:`~htmd.molecule.molecule.Molecule`.

Atom Selections: Advanced Expressions and Keywords
--------------------------------------------------

=================== ========
Type                Example
=================== ========
Chemistry           `protein and not hydrogen`
Secondary structure `helix or sheet`
Sequence            `sequence "N..T"`
Maths               `x^2 + y^2 < 40`
... within N of ... `water within 4 of protein`
same ... as ...     `same chain as resid 15`
=================== ========

::

    water     within 4 of protein
    \___/ and \___________\_____/

- selects all `water` atoms within 4 Å of `protein` atoms.

::

    same chain as  resid 15
    \___________/  \______/

- selects all atoms that are in the same chain as any atoms belonging to residue 15

