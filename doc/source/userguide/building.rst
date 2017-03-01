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

General Considerations for System Building
==========================================

PDB format considerations
-------------------------

The PDB format is very old. In an effort to handle its legacy shortcomings, several versions have been made over the
years, they are not all readily interchangeable, and not all software can handle each version perfectly. The most
important things to watch out for are:

- Columns: the PDB format has very rigid rules about what values can go in each space. Keep in mind that it is not a
  space/tab/comma delimited format, but rather has rigid definitions of what should be in each space/column.

- The PDB format as originally designed cannot handle more than 9,999 resids or 99,999 atoms (due to the column format
  issue). Several workarounds have been devised, such as using hexadecimal numbers or other compact number formats. VMD
  has no trouble saving more atoms/residues.

Physical and chemical considerations
------------------------------------

One needs to know well the working system, thus:

- Always review your PDB file: inspect the REMARK sections of the PDB file. You can often find key specific information
  regarding the structure (e.g. disulphide bonds, missing atoms, etc.).

- Disulphide bonds present in the system must be identified. This is automatically done by HTMD in both CHARMM and Amber.

- Metalloproteins: if the metal ion is not an active part of an interaction it may be acceptable to just allow it to act
  as a cation perhaps restraining it with some harmonic constraints if necessary.

- Duplicate atoms in the PDB file: typically simply delete one of the duplicated groups. However, if both conformations
  are potentially important (e.g. such loops involved in molecular recognition) it might be necessary to simulate both
  conformations separately.

Protonation/pH
--------------

The protonation state of the system is critical. Since MD simulations typically don't allow for bond breaking, the
initial protonation of the system must be accurate. Knowing what pH you are trying to reproduce is therefore important
to obtain the correct results. If you suspect changing protonation is important to your system and you still want to use
classical mechanics, consider simulating both states (protonated and not protonated).

Histidine residues can have three different protonations states even at pH 7, therefore, a correct protonation of this
residue is particularly critical. This residue can be protonated at either delta (most common; HSD/HID), epsilon (very
common also; HSE/HIE) or at both nitrogens (special situations and low pH; HSP/HIP).

.. image:: http://docs.htmd.org/img/histidines.png

The best way to determine how histidine should be protonated is to look at the the structure. Typically, a histidine
residue is protonated if it is close enough to an electron donor (e.g. a glutamic acid), thus creating a hydrogen bond.
Since histidines are frequently present at protein active sites, a correct protonation state is particularly important
in ligand binding simulations.

In HTMD, one can use :class:`~htmd.builder.preparation.proteinPrepare` to help with protonation.

List of useful tools
--------------------

====================================================== ==================
:class:`~htmd.molecule.molecule.Molecule` class        Building functions
====================================================== ==================
:func:`~htmd.molecule.molecule.Molecule.append`        :func:`~htmd.molecule.util.maxDistance`
:func:`~htmd.molecule.molecule.Molecule.center`        :func:`~htmd.molecule.util.uniformRandomRotation`
:func:`~htmd.molecule.molecule.Molecule.mutateResidue` :func:`~htmd.molecule.util.boundingBox`
:func:`~htmd.molecule.molecule.Molecule.moveBy`        :func:`~htmd.molecule.util.sequenceID`
:func:`~htmd.molecule.molecule.Molecule.rotateBy`      :func:`~htmd.builder.builder.embed`
---                                                    :func:`~htmd.builder.builder.autoSegment`
====================================================== ==================

List of common CHARMM patches
-----------------------------

- C-terminal patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
CTER -1     standard C-terminus
CT1  0      methylated C-terminus from methyl acetate
CT2  0      amidated C-terminus
CT3  0      N-Methylamide C-terminus
==== ====== ===========

- N-terminal patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
NTER +1     standard N-terminus
ACE  0      acetylated N-terminus (to create dipeptide)
ACP  0      acetylated N-terminus (for proline dipeptide)
PROP +1     Proline N-Terminal
GLYP +1     Glycine N-terminus
==== ====== ===========

- Side-chain patches

==== ====== ===========
Name Charge Description
==== ====== ===========
ASPP 0      patch for protonated aspartic acid, proton on OD2
GLUP 0      patch for protonated glutamic acid, proton on OE2
CYSD -1     patch for deprotonated CYS
DISU +1     patch for disulfides. Patch must be 1-CYS and 2-CYS
HS2  +1     patch for neutral His, move proton from ND1 to NE2
TP1  -1     patch to convert tyrosine to monoanionic phosphotyrosine
TP1A -1     patch to convert tyrosine to monoanionic phenol-phosphate model compound when generating tyr, use first none last none for terminal patches
TP2  -2     patch to convert tyrosine to dianionic phosphotyrosine
TP2A -2     patch to convert tyrosine to dianionic phosphotyrosine when generating tyr, use first none last none for terminal patches this converts a single tyrosine to a phenol phosphate
TMP1 -1     patch to convert tyrosine to monoanionic phosphonate ester O -> methylene (see RESI BMPH)
TMP2 -2     patch to convert tyrosine to dianionic phosphonate ester O -> methylene (see RESI BMPD)
TDF1 -1     patch to convert tyrosine to monoanionic difluoro phosphonate ester O -> methylene (see RESI BDFH)
==== ====== ===========

- Circular protein chain patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
LIG1 0      linkage for cyclic peptide, 1 refers to the C terminus which is a glycine , 2 refers to the N terminus
LIG2 0      linkage for cyclic peptide, 1 refers to the C terminus, 2 refers to the N terminus which is a glycine
LIG3 0      linkage for cyclic peptide, 1 refers to the C terminus which is a glycine, 2 refers to the N terminus which is a glycine
==== ====== ===========
