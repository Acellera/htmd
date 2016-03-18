
Atom selections
===============

**Toni Giorgino**

Institute of Neurosciences (IN-CNR)
Consiglio Nazionale delle Ricerche
Padua, Italy

A *Molecule* has...
-------------------

-  Frames (≥ 0), usually corresponding to a trajectory's time
-  Atoms (a constant number, ≥ 0)
-  Representations
-  Name, top-ness, visibility, etc.
-  (Connectivity & topology)
-  (Volumetric data)

Note
----

The word *Molecule* is somewhat a misnomer. It is more like a *System*:
in general it contains several "chemical" molecules.

E.g. protein + solvent + ions are often **one** "VMD Molecule".

(The closest analog to molecules in the chemical sense is having same
*segment ID*s)

This applies to both VMD and HTMD.

Atoms
-----

A Molecule contains several *atom*s

An *atom* has several **properties**
------------------------------------

-  Most obviously mirror PDB fields
-  Some are computed via heuristics (!)
-  *Constant* properties do not depend on the frame. Eg:
-  ``type``, ``name``, ``mass``, ...
-  *Variable* properties do. Eg:
-  ``x(t)``, ``y(t)``, ``z(t)``, ``user(t)``

\| Property \| Description \| \|-----\|-----\| \| ``serial`` \| Unique
identifier from 1\| \| ``name`` \| Atom name from PDB, unique in the
residue \| \| ``type`` \| Defined by the forcefield, if any \| \|
``resname``\| 3-letter residue type, from PDB \| \| ``resid`` \| The
usual residue number, **possibly not unique** \| \| ``residue``\| Do not
confuse with ``resid``! \| \| ``chain`` \| Also from PDB \| \| ``segid``
\| Must be unique *per chemical molecule*: important for system building
\| \| ``x,y,z`` \| these depend on the trajectory frame \|

What for?
---------

Mostly, we'll build **boolean expressions** based on properties.

    E.g: ``"water"`` or ``"name CA and helix"``

Expressions are evaluated on a set of atoms; each of them either matches
or not.

Therefore, selection expression yield **atom selections**, i.e. subsets
of the original set.

We can then manipulate these subsets, e.g. to *get* and *set*
properties.

\| Example \| Meaning \| \|---------\|---------\| \| ``chain A`` \| true
*iff* the atom's ``chain`` property equals string ``A`` \| \|
``resid 40 to 50`` \| as intended\| \| ``mass < 15`` \| you guessed it\|
\| ``name CA`` \| atom is a Cα \| \| ``water`` \| true iff the atom
belongs to a water molecule \|

Important: all are **boolean** expressions. Each atom in a molecule
either matches, or it doesn't.

Advanced expressions
--------------------

\| Type \| Example \| \|---------\|---------\| \| Chemistry \|
``protein and not hydrogen`` \| \| Secondary structure \|
``helix or sheet`` \| \| Sequence \| ``sequence "N..T"`` \| (ex) Within
\| ``water within 4 of protein`` \| \| Same \|
``same chain as resid 15`` \| \| Maths \| ``x^2 + y^2 < 40`` \|

Common operations
-----------------

-  Set properties, in particular ``resname`` and ``segid``.
-  Copy properties, in some cases.
-  Deleting unwanted atoms.
-  Extraction of ligands for parametrization.
-  Rotation of ligands.
-  In-silico mutation residues: change ``resname`` + delete sidechain
   atoms.

Representations
---------------

A molecule contains several *Representation*s

A *Representation* has...
-------------------------

-  *Selection*: which part of the molecule to show
-  *Style*: ribbon, lines, VdW, ...
-  *Color*: either a number, or a property (e.g: secondary structure)
-  Other style-dependent properties
-  *Material*
-  *Visibility*: hidden or shown (double-click)

Case study
----------

-  Load PDB: **3PTB**
-  Identify non-protein
-  List residues in contact with ligand
-  Are residue numbers in protein contiguous?
-  Are there duplicated residues?

Residues in contact
-------------------

Use an atomselection similar to this

::

    protein     within 4 of resname BEN 
    \_____/ and \___________\_________/

To see whole residues

::

    same residue as  protein within 4 of resname BEN
    \_____________/  ...

To get some context, add a line like

::

    name N CA C O 

Are there breaks in the peptide chain?
--------------------------------------

-  Get ``resid`` of Cαs, compute finite differences; or
-  Use the differece between ``resid`` and ``residue``; or
-  Use ``unique()``, set operators, and the like

Discontinuities in ``resid`` may or may not be breaks in the chain!

Are there duplicate residues?
-----------------------------

-  Get ``resid`` of Cαs, compute finite differences as above

Answers
-------

-  Residues 184, 188, 221 are duplicated
-  Residues in contact (any atom ≤ 4 Å) with BEN :
   ``189 190 191 192 195 213 215 216 219 220 226``
-  Jumps at
   ``{ASN 34} {LEU 67} {THR 125} {SER 130} {LYS 204} {SER 217}``
    These are not actual breaks in the peptide chain
