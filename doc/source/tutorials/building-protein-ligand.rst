
Building a system (globular protein and ligand)
===============================================

by Stefan Doerr

This tutorial shows how to automatically build a system using the HTMD
code. In this example, a system comprised of a globular protein
(trypsin) and its inhibitor (benzamidine) is prepared for ligand binding
simulations using ACEMD (Buch et al. 2011 PNAS 108(25) 10184-89).

Working files
-------------

To generate a system ready for simulation, HTMD needs the topology,
parameters and coordinates (i.e. PDB file or PDB code of the protein) of
all elements in the simulations.

Download all tutorial files from the following
`link <http://docs.htmd.org/download/building_files.zip>`_.

Getting started
---------------

First we import the modules we are going to need for the tutorial

.. code:: python

    from htmd import *
    import numpy as np
    import os
    %pylab inline
    matplotlib.rcParams.update({'font.size': 12})


.. parsed-literal::

    HTMD. 26-27 November 2015 HTMD workshop in Barcelona (fully booked)
    
    You are on the latest HTMD version (unpackaged).
    Populating the interactive namespace from numpy and matplotlib


.. parsed-literal::

    WARNING: pylab import has clobbered these variables: ['random']
    `%matplotlib` prevents importing * from pylab and numpy


Clean and split PDB
-------------------

In this example, the structure is retrieved from the PDB database and
saved in the working directory.

.. code:: python

    prot = Molecule('3PTB')
    prot.write('/tmp/protein.pdb', sel='chain A and protein and noh')

Extract only chain A from the PDB file

.. code:: python

    prot.filter('chain A and (protein or water or resname CA)')

Define segments
---------------

Crystallized water molecules and one calcium ion present in the crystal
structure are also obtained from this PDB.

.. code:: python

    prot.set('segid', 'P', sel='protein and noh')
    prot.set('segid', 'W', sel='water')
    prot.set('segid', 'CA', sel='resname CA')

Center the system to the origin
-------------------------------

.. code:: python

    prot.center()

and calculate the radius of the protein

.. code:: python

    from htmd.molecule.util import maxDistance, uniformRandomRotation
    D = maxDistance(prot, 'all')

Load the ligand and merge molecules
-----------------------------------

Load up the ligand, calculate its geometric center and randomly rotate
it around itself

.. code:: python

    tutfiles = home() + '/data/building-protein-ligand/'
    ligand = Molecule(tutfiles + 'benzamidine.pdb')
    ligand.center()
    ligand.rotateBy(uniformRandomRotation())

Now place the ligand randomly around the protein at the distance defined
above

.. code:: python

    ligand.moveBy([D+10, 0, 0])  # Move the ligand 10 Angstrom away from the furthest protein atom in X dimension
    ligand.rotateBy(uniformRandomRotation())

Set resname and segid of the ligand

.. code:: python

    ligand.set('resname','MOL')
    ligand.set('segid','L')

Join all

.. code:: python

    all = Molecule()
    all.append(prot)
    all.append(ligand)

Solvate the system
------------------

Define the size of the solvation box and solvate the system

.. code:: python

    D = D + 15
    allsol = solvate(all, minmax=[[-D, -D, -D], [D, D, D]])


.. parsed-literal::

    2015-11-23 18:44:22,004 - htmd.builder.solvate - INFO - Using water pdb file at: /shared/sdoerr/Work/pyHTMD/htmd/builder/wat.pdb
    2015-11-23 18:44:22,626 - htmd.builder.solvate - INFO - Replicating 8 water segments, 2 by 2 by 2
    Solvating: 100% (8/8) [############################################] eta 00:01 /


Building the system for CHARMM
------------------------------

Check for the available charmm parameter and topology files

.. code:: python

    charmm.listFiles()


.. parsed-literal::

    ---- Topologies files list: /shared/sdoerr/Work/pyHTMD/htmd/builder/charmmfiles/top/ ----
    top/top_all22star_prot.rtf
    top/top_all36_carb.rtf
    top/top_all36_lipid.rtf
    top/top_all36_prot.rtf
    top/top_water_ions.rtf
    top/top_all36_cgenff.rtf
    top/top_all36_na.rtf
    ---- Parameters files list: /shared/sdoerr/Work/pyHTMD/htmd/builder/charmmfiles/par/ ----
    par/par_all22star_prot.prm
    par/par_all36_carb.prm
    par/par_all36_lipid.prm
    par/par_all36_prot.prm
    par/par_all36_cgenff.prm
    par/par_all36_na.prm
    par/par_water_ions.prm


Indicate the location of the CHARMM topology and parameter files as well
are your own custom parameter and topology files. The CHARMM files can
be included without their full path, using just the name indicated in
the previous list command. Then we build the system for CHARMM.

.. code:: python

    topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf', tutfiles + 'benzamidine.rtf']
    params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', tutfiles + 'benzamidine.prm']
    
    molbuilt = charmm.build(allsol, topo=topos, param=params, outdir='/tmp/build', saltconc=0.15)


.. parsed-literal::

    2015-11-23 18:44:40,576 - htmd.builder.charmm - INFO - Writing out segments.
    Bond between A: [serial 48 resid 22 resname CYS chain A segid P]
                 B: [serial 1007 resid 157 resname CYS chain A segid P]
    
    Bond between A: [serial 185 resid 42 resname CYS chain A segid P]
                 B: [serial 298 resid 58 resname CYS chain A segid P]
    
    Bond between A: [serial 811 resid 128 resname CYS chain A segid P]
                 B: [serial 1521 resid 232 resname CYS chain A segid P]
    
    Bond between A: [serial 853 resid 136 resname CYS chain A segid P]
                 B: [serial 1327 resid 201 resname CYS chain A segid P]
    
    Bond between A: [serial 1084 resid 168 resname CYS chain A segid P]
                 B: [serial 1190 resid 182 resname CYS chain A segid P]
    
    Bond between A: [serial 1265 resid 191 resname CYS chain A segid P]
                 B: [serial 1422 resid 220 resname CYS chain A segid P]
    
    2015-11-23 18:45:12,563 - htmd.builder.builder - INFO - 6 disulfide bonds were added
    2015-11-23 18:45:12,759 - htmd.builder.charmm - INFO - Starting the build.
    2015-11-23 18:45:14,687 - htmd.builder.charmm - INFO - Finished building.
    2015-11-23 18:45:15,869 - htmd.builder.ionize - INFO - Adding 9 anions + 0 cations for neutralizing and 108 ions for the given salt concentration.
    2015-11-23 18:45:16,182 - htmd.builder.ionize - INFO - Min distance of ions from molecule: 5A
    2015-11-23 18:45:16,183 - htmd.builder.ionize - INFO - Min distance between ions: 5A
    2015-11-23 18:45:16,183 - htmd.builder.ionize - INFO - Placing 117 ions.
    2015-11-23 18:46:25,468 - htmd.builder.charmm - INFO - Writing out segments.
    2015-11-23 18:46:37,921 - htmd.builder.charmm - INFO - Starting the build.
    2015-11-23 18:46:39,717 - htmd.builder.charmm - INFO - Finished building.


Note regarding ions: the build command will by default just try to
neutralize the system. To add a specific salt concentration the option
``saltconc`` needs to be used. In the previous command, a 150mM NaCl
salt concentration was used.

Visualize the built system

.. code:: python

    molbuilt.view(sel='water',style='Lines',hold=True)
    molbuilt.view(sel='resname MOL',style='Licorice',hold=True)
    molbuilt.view(sel='ions',style='VDW',hold=True)
    molbuilt.view(sel='protein',style='NewCartoon',color='Secondary Structure')

Building the system for AMBER
-----------------------------

Check for available AMBER forcefield files

.. code:: python

    amber.listFiles()


.. parsed-literal::

    ---- Forcefield files list: /shared/lab/software/AmberTools14/amber14/dat/leap/cmd/ ----
    leaprc.phosaa10
    leaprc.GLYCAM_06j-1
    leaprc.ff14ipq
    leaprc.lipid11
    leaprc.gaff
    leaprc.lipid14
    leaprc.modrna08
    leaprc.GLYCAM_06EPb
    leaprc.ff12SB
    leaprc.ff03.r1
    leaprc.constph
    leaprc.ff14SBonlysc
    leaprc.ff03ua
    leaprc.ffAM1
    leaprc.ffPM3
    leaprc.ff14SB


Indicate the desired forcefield files and build the system for AMBER

.. code:: python

    ffs = ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']  # Missing the parameters for Benzamidine in AMBER
    
    #molbuilt = amber.build(allsol, ff=ffs, outdir='/tmp/build', saltconc=0.15)

Visualize the built system

.. code:: python

    molbuilt.view()

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

