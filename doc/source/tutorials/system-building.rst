
System building intro
=====================

by Stefan Doerr

You can watch the presentation here:

|image0|

.. |image0| image:: http://pub.htmd.org/73hboiwia98hdj209jq0/opioid_youtube.png
   :target: https://youtu.be/DF9cHKBX19A

Contents
--------

1. Building
2. Tools
3. Workflow
4. Considerations

Building
--------

What is system building:

-  Starting from molecular structures
-  Modify them, position them and combine them to
-  Prepare a biological system for simulation
-  For a given forcefield

Tools
-----

.. figure:: http://pub.htmd.org/73hboiwia98hdj209jq0/molbuilding.png
   :alt: 

Workflow
--------

1. Obtain structures
2. Clean structures
3. Define segments
4. Combine structures
5. Solvate
6. Build and ionize

Before building your system (preliminary considerations)
--------------------------------------------------------

The PDB format is very old. In an effort to handle its legacy
shortcomings, several versions have been made over the years, they are
not all readily interchangeable, and not all software can handle each
version perfectly. The most important things to watch out for are: \*
Columns: the PDB format has very rigid rules about what values can go in
each space. Keep in mind that it is not a space/tab/comma delimited
format, but rather has rigid definitions of what should be in each
space/column. \* The PDB format as originally designed cannot handle
more than 9,999 resids or 99,999 atoms (due to the column format issue).
Several workarounds have been devised, such as using hexadecimal numbers
or other compact number formats. VMD has no trouble saving more
atoms/residues.

In addition, one needs to know well the working system, thus: \* Always
review your PDB file: inspect the REMARK sections of the PDB file. You
can often find keyspecific information regarding the structure (e.g.
disulfide bridges, mising atoms, etc.).

-  Protonation/pH: the protonation state of the system is critical.
   Since molecular dynamics simulations typically don't allow for bond
   breaking, the initial protonation of the system must be accurate.
   Knowing what pH you are trying to reproduce is therefore important to
   obtain the correct results. If you suspect changing protonation is
   important to your system and you still want to use classical
   mechanics, consider simulating both states (protonated and not
   protonated). Histidine residues can have three different protonations
   states even at pH 7, therefore, a correct protonation of this residue
   is particularly critical. This residue can be protonated at either
   delta (most common), epsilon (very common also) or at both nitrogens
   (special situations and low pH).

.. figure:: http://docs.htmd.org/img/histidines.png
   :alt: 

The best way to determine how histidine should be protonated is to look
at the the structure. Typically, a histidine residue is protonated if it
is close enough to an electron donor (e.g. a glutamic acid), thus
creating a hydrogen bond. Certain automated tools predict the
protonation state of histidines based on their surrounding environment
(e.g. Autodock tools). Since histidines are frequently present at
protein active sites, a correct protonation state is particularly
important in ligand binding simulations.

-  Disulfide bonds present in the system must be identified. As shown
   below, this is automatically done by htmd
-  Metalloproteins: if the metal ion is not an active part of an
   interaction it may be acceptable to just allow it to act as a cation
   perhaps restraining it with some harmonic constraints if neccesary.
-  Duplicate atoms in the PDB file: typically simply delete one of the
   duplicated groups. However, if both conformations are potentially
   important (e.g. such loops involved in molecular recognition) it
   might be necessary to simulate both conformations separately.

List of common patches
----------------------

C-terminal patches:

.. raw:: html

   <table class="summarytable">
       <thead>
           <tr>
               <th>

Name

.. raw:: html

   </th>
               <th>

Class

.. raw:: html

   </th>
               <th>

Description

.. raw:: html

   </th>
           </tr>
       </thead>
       <tbody>
           <tr>
               <td>

CTER

.. raw:: html

   </td>
               <td>

-1.00

.. raw:: html

   </td>
               <td>

standard C-terminus

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

CT1

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

methylated C-terminus from methyl acetate

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

CT2

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

amidated C-terminus

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

CT3

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

N-Methylamide C-terminus

.. raw:: html

   </td>
           </tr>
       </tbody>
   </table>

N-terminal patches:

.. raw:: html

   <table class="summarytable">
       <thead>
           <tr>
               <th>

Name

.. raw:: html

   </th>
               <th>

Class

.. raw:: html

   </th>
               <th>

Description

.. raw:: html

   </th>
           </tr>
       </thead>
       <tbody>
           <tr>
               <td>

NTER

.. raw:: html

   </td>
               <td>

1.00

.. raw:: html

   </td>
               <td>

standard N-terminus

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

ACE

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

acetylated N-terminus (to create dipeptide)

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

ACP

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

acetylated N-terminus (for proline dipeptide)

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

PROP

.. raw:: html

   </td>
               <td>

1.00

.. raw:: html

   </td>
               <td>

Proline N-Terminal

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

GLYP

.. raw:: html

   </td>
               <td>

1.00

.. raw:: html

   </td>
               <td>

Glycine N-terminus

.. raw:: html

   </td>
           </tr>
       </tbody>
   </table>

Side chain patches:

.. raw:: html

   <table class="summarytable">
       <thead>
           <tr>
               <th>

Name

.. raw:: html

   </th>
               <th>

Class

.. raw:: html

   </th>
               <th>

Description

.. raw:: html

   </th>
           </tr>
       </thead>
       <tbody>
           <tr>
               <td>

ASPP

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

patch for protonated aspartic acid, proton on od2

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

GLUP

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

patch for protonated glutamic acid, proton on oe2

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

CYSD

.. raw:: html

   </td>
               <td>

-1.0

.. raw:: html

   </td>
               <td>

patch for deprotonated CYS

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

DISU

.. raw:: html

   </td>
               <td>

-0.36

.. raw:: html

   </td>
               <td>

patch for disulfides. Patch must be 1-CYS and 2-CYS

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

HS2

.. raw:: html

   </td>
               <td>

0.00

.. raw:: html

   </td>
               <td>

Patch for neutral His, move proton from ND1 to NE2

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TP1

.. raw:: html

   </td>
               <td>

-1.00

.. raw:: html

   </td>
               <td>

convert tyrosine to monoanionic phosphotyrosine

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TP1A

.. raw:: html

   </td>
               <td>

-1.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to monoanionic phenol-phosphate model compound
when generating tyr, use first none last none for terminal patches

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TP2

.. raw:: html

   </td>
               <td>

-2.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to dianionic phosphotyrosine

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TP2A

.. raw:: html

   </td>
               <td>

-2.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to dianionic phosphotyrosine when generating
tyr, use first none last none for terminal patches this converts a
single tyrosine to a phenol phosphate

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TMP1

.. raw:: html

   </td>
               <td>

-1.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to monoanionic phosphonate ester O ->
methylene (see RESI BMPH)

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TMP2

.. raw:: html

   </td>
               <td>

-2.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to dianionic phosphonate ester O -> methylene
(see RESI BMPD)

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

TDF1

.. raw:: html

   </td>
               <td>

-1.00

.. raw:: html

   </td>
               <td>

patch to convert tyrosine to monoanionic difluoro phosphonate ester O ->
methylene (see RESI BDFH)

.. raw:: html

   </td>
           </tr>
       </tbody>
   </table>

Circular protein chain patches:

.. raw:: html

   <table class="summarytable">
       <thead>
           <tr>
               <th>

Name

.. raw:: html

   </th>
               <th>

Class

.. raw:: html

   </th>
               <th>

Description

.. raw:: html

   </th>
           </tr>
       </thead>
       <tbody>
           <tr>
               <td>

LIG1

.. raw:: html

   </td>
               <td>

0.00000

.. raw:: html

   </td>
               <td>

linkage for cyclic peptide, 1 refers to the C terminus which is a
glycine , 2 refers to the N terminus

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

LIG2

.. raw:: html

   </td>
               <td>

0.00000

.. raw:: html

   </td>
               <td>

linkage for cyclic peptide, 1 refers to the C terminus, 2 refers to the
N terminus which is a glycine

.. raw:: html

   </td>
           </tr>
           <tr>
               <td>

LIG3

.. raw:: html

   </td>
               <td>

0.00000

.. raw:: html

   </td>
               <td>

linkage for cyclic peptide, 1 refers to the C terminus which is a
glycine, 2 refers to the N terminus which is a glycine

.. raw:: html

   </td>
           </tr>
       </tbody>
   </table>
