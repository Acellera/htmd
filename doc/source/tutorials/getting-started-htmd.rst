
Getting started with HTMD
=========================

Assuming that you have already download it and installed htmd, this
tutorial guides you through the basic language features.

Any object and function defined by htmd is available in the workspace.
Let's look more carefully at the Molecule class features.

.. code:: python

    from htmd import *


.. parsed-literal::

    HTMD. All material of HTMD2015 will be soon made available
    
    You are on the latest HTMD version (unpackaged).


Molecule objects
----------------

First create an empty molecule object by either \* Fetching it from the
Protein Data Bank by using its PDB code,

.. code:: python

    mol = Molecule('3PTB')

or just use your own PDB file

.. code:: python

    #mol = Molecule('yourprotein.pdb')

Inspect your molecule
---------------------

Molecule.get can be used to check and retreive specific PDB fields, for
example: \* Check the resIds of the triptophan residues present in your
protein,

.. code:: python

    np.unique(mol.get('resid','resname TRP'))




.. parsed-literal::

    array([ 51, 141, 215, 237])



Store the coordinates of a specific atom,

.. code:: python

    mol.get('coords','resname TRP and resid 51 and name CA')




.. parsed-literal::

    array([ 10.06900024,   3.77600002,  34.72100067], dtype=float32)



Display the number of chains or segments present in your PDB file,

.. code:: python

    np.unique(mol.get('chain'))




.. parsed-literal::

    array(['A'], dtype=object)



Note that unless you select a specific atom, unique() is needed to avoid
duplicates in the output.

Duplicate/modify objects and fields
-----------------------------------

Use Molecule.copy to duplicate the molecule into a different object,

.. code:: python

    newmol = mol.copy()

Molecule.writePDB can be used to output a PDB file of your whole
molecule (or just a selection) The following command use the above
copied molecule to write out a PDB file of the ligand atoms present in
the fetched PDB file except for hydrogen

.. code:: python

    newmol.write('/tmp/ligand.pdb','resname BEN and noh')

Alternatively, Molecule.filter can be used to clean/select/remove
specific parts such as chains, segments, etc. For example, clean all
except for protein atoms in chain A

.. code:: python

    mol.filter('chain A and protein')

Molecule.set is instead used to change/name/rename specific fields. For
example, set can create a segid called 'P' out of of the protein atoms,

.. code:: python

    mol.set('segid','P','protein');

or rename all HIS residues to HSN

.. code:: python

    mol.set('resname','HSN','resname HIS')

Joining molecules/segments
--------------------------

Molecule.append append two separated Molecule objects (e.g. ligand,
water or ion segments, etc.) For example, to append the pdb of the
ligand (saved above) to the molecule we are working with, simply do

.. code:: python

    ligand=Molecule('/tmp/ligand.pdb')
    mol.append(ligand)

Playing with coordinates
------------------------

Coordinates can be used to perform geometric tasks on your molecule: \*
Calculate the geometric center of your molecule

.. code:: python

    coo=mol.get('coords')
    c = np.mean(coo,axis=0)

Use Molecule.moveBy to translate and center your molecule to the origin
[0, 0, 0]

.. code:: python

    mol.moveBy(-c)

Check the new center

.. code:: python

    np.mean(mol.get('coords'),axis=0)




.. parsed-literal::

    array([  2.15420940e-07,  -2.36497249e-06,   4.63504330e-06], dtype=float32)



Perform translation and/or rotations (e.g. to build a protein - ligand
system):

Load up your ligand and calculate its geometric center as above

.. code:: python

    ligand = Molecule('/tmp/ligand.pdb')
    ligcenter = np.mean(ligand.coords,axis=0)

Use Molecule.rotateBy to rotate your ligand

.. code:: python

    ligand.rotate([1, 0, 0],math.pi/2)

Note that the uniformRandomRotation() function provide the random
coordinates needed for this rotation.
