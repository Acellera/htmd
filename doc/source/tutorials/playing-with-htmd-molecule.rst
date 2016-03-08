
.. code:: python

    from htmd import *
    from numpy import *
    htmd.config(viewer='ngl')

Playing with Molecule objects in HTMD
=====================================

Assuming that you have already download it and installed htmd, this
tutorial guides you through the basic language features.

Any object and function defined by htmd is available in the workspace.
Let's look more carefully at the Molecule class features.

Molecule objects
----------------

First create an empty molecule object by either

-  Fetching it from the Protein Data Bank by using its PDB code...
-  ...or local files (pdb, mol2, xtc, psf, prmtop)

.. code:: python

    mol = Molecule('3PTB')

Inspect your molecule
---------------------

Printing the instance shows its properties.

.. code:: python

    print(mol)


.. parsed-literal::

    Molecule with 1701 atoms and 1 frames
    PDB field - altloc shape: (1701,)
    PDB field - beta shape: (1701,)
    PDB field - chain shape: (1701,)
    PDB field - charge shape: (1701,)
    PDB field - coords shape: (1701, 3, 1)
    PDB field - element shape: (1701,)
    PDB field - insertion shape: (1701,)
    PDB field - name shape: (1701,)
    PDB field - occupancy shape: (1701,)
    PDB field - record shape: (1701,)
    PDB field - resid shape: (1701,)
    PDB field - resname shape: (1701,)
    PDB field - segid shape: (1701,)
    PDB field - serial shape: (1701,)
    bonds shape: (42, 2)
    box shape: (3, 1)
    fileloc shape: (1, 2)
    frame: 0
    masses shape: (1701,)
    reps: 0
    ssbonds shape: (0,)
    step shape: (0,)
    time shape: (0,)
    viewname: 3PTB


Methods and properties of Molecule instances
--------------------------------------------

Each *molecule* instance has a number of methods (operations that you
can perform on it) and properties (data associated to the molecule).
Properties (among others) correspond to data which is found in PDB
files.

+----------------+--------------+
| Methods        | Properties   |
+================+==============+
| read()         | record       |
+----------------+--------------+
| write()        | serial       |
+----------------+--------------+
| get()          | name         |
+----------------+--------------+
| set()          | resid        |
+----------------+--------------+
| atomselect()   | chain        |
+----------------+--------------+
| filter()       | coords       |
+----------------+--------------+
| remove()       | box          |
+----------------+--------------+
| insert()       | reps         |
+----------------+--------------+
| view()         | ...          |
+----------------+--------------+
| wrap()         |              |
+----------------+--------------+
| align()        |              |
+----------------+--------------+

Properties can be accessed either directly, or via the
``Molecule.get()`` method. Similarly, they can be modified directly, or
via the ``Molecule.set()`` method. This pair of methods is known as
"getter/setter" methods in the object-oriented jargon.

For example:

.. code:: python

    mol.serial




.. parsed-literal::

    array([   1,    2,    3, ..., 1700, 1701, 1702])



.. code:: python

    mol.get("serial")




.. parsed-literal::

    array([   1,    2,    3, ..., 1700, 1701, 1702])



The following sections will show property getters and setters in use in
number of real-world tasks.

Check the resIds of the cystein residues present in your protein
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    mol.get('resid',sel='resname CYS')




.. parsed-literal::

    array([ 22,  22,  22,  22,  22,  22,  42,  42,  42,  42,  42,  42,  58,
            58,  58,  58,  58,  58, 128, 128, 128, 128, 128, 128, 136, 136,
           136, 136, 136, 136, 157, 157, 157, 157, 157, 157, 168, 168, 168,
           168, 168, 168, 182, 182, 182, 182, 182, 182, 191, 191, 191, 191,
           191, 191, 201, 201, 201, 201, 201, 201, 220, 220, 220, 220, 220,
           220, 232, 232, 232, 232, 232, 232])



Note how residue IDs are duplicated. This is due to the fact that one
value is returned per matched atom, and this PDB file has approximately
6 atoms resolved per cystein residue:

.. code:: python

    mol.get('name','resname CYS and resid 58')




.. parsed-literal::

    array(['N', 'CA', 'C', 'O', 'CB', 'SG'], dtype=object)



To obtain one residue ID per residue, we can either further restrict the
selection to carbon α atoms...

.. code:: python

    mol.get('resid',sel='name CA and resname CYS')




.. parsed-literal::

    array([ 22,  42,  58, 128, 136, 157, 168, 182, 191, 201, 220, 232])



...or use Python's ``unique`` function to remove duplicates...

.. code:: python

    unique(mol.get('resid',sel='resname CYS'))




.. parsed-literal::

    array([ 22,  42,  58, 128, 136, 157, 168, 182, 191, 201, 220, 232])



Retrieve the coordinates of a specific atom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is done accessing the "coords" property. It is special, in the
sense that it returns a 3-column vector (for the three coordinates).
Also note how its precision is restricted to the one in the PDB file.

.. code:: python

    mol.get('coords','resname CYS and resid 58 and name CA')




.. parsed-literal::

    array([  4.23999977,  16.49500084,  27.98600006], dtype=float32)



What is returned if more than one atom is selected? A matrix.

.. code:: python

    mol.get('coords','resname CYS and resid 58')




.. parsed-literal::

    array([[  5.12200022,  16.71899986,  26.86300087],
           [  4.23999977,  16.49500084,  27.98600006],
           [  4.87400007,  16.95800018,  29.29999924],
           [  4.23799992,  16.76399994,  30.36199951],
           [  3.94099998,  14.9989996 ,  28.07099915],
           [  2.79200006,  14.45199966,  26.72200012]], dtype=float32)



Display the number of chains or segments present in your PDB file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    unique(mol.get('chain'))




.. parsed-literal::

    array(['A'], dtype=object)



Which means that everything is assigned to the same chain.

List atoms recognized as water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    # Get their indices
    mol.get("serial",sel="water")




.. parsed-literal::

    array([1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649, 1650, 1651,
           1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1661, 1662,
           1663, 1664, 1665, 1666, 1667, 1668, 1669, 1670, 1671, 1672, 1673,
           1674, 1675, 1676, 1677, 1678, 1679, 1680, 1681, 1682, 1683, 1684,
           1685, 1686, 1687, 1688, 1689, 1690, 1691, 1692, 1693, 1694, 1695,
           1696, 1697, 1698, 1699, 1700, 1701, 1702])



.. code:: python

    # Count them
    len(mol.get("serial",sel="water"))




.. parsed-literal::

    62



The atomselect() returns a vector of boolean values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    mol.atomselect("water")




.. parsed-literal::

    array([False, False, False, ...,  True,  True,  True], dtype=bool)



We use the fact that True counts as 1 in sum(), and obtain the same
result.

.. code:: python

    print(mol.atomselect("water"))
    sum(mol.atomselect("water"))


.. parsed-literal::

    [False False False ...,  True  True  True]




.. parsed-literal::

    62



Duplicate/modify objects and fields
-----------------------------------

Use Molecule.copy to duplicate the molecule into a different object,

.. code:: python

    newmol = mol.copy()

Alternatively, Molecule.filter can be used to clean/select/remove
specific parts such as chains, segments, etc. For example, clean all
except for protein atoms in chain A

.. code:: python

    mol.filter('chain A and protein')

Molecule.set is instead used to change/name/rename specific fields. For
example, set can create a segid called 'P' out of of the protein atoms,

.. code:: python

    mol.set('segid','P',sel='protein');

or rename all HIS residues to HSN

.. code:: python

    mol.set('beta',1,sel='resname HIS')

Joining molecules/segments
--------------------------

Molecule.append append two separated Molecule objects (e.g. ligand,
water or ion segments, etc.) For example, to append the pdb of the
ligand (saved above) to the molecule we are working with, simply do

.. code:: python

    ligand=Molecule('3PTB')
    ligand.filter('resname BEN')
    mol.append(ligand)

You can add an atom in a similar fashion.

.. code:: python

    atom = Molecule()
    atom.name = 'CA'
    atom.resid = 0
    atom.chain = 'X'
    atom.coords = [6, 3, 2]
    
    mol.insert(atom, 0)

Writing
-------

The ``writePDB()`` method can be used to output a PDB file of your whole
molecule (or just a selection).

The following command use the above copied molecule to write out a PDB
file of the ligand atoms present in the fetched PDB file except for
hydrogen.

.. code:: python

    ligand.write('/tmp/ligand.pdb','resname BEN and noh')

Playing with coordinates
------------------------

Coordinates can be used to perform geometric tasks on your molecule

**Calculate the geometric center of your molecule**

.. code:: python

    coo=mol.get('coords')
    print(coo)
    c = mean(coo,axis=0)
    print(c)


.. parsed-literal::

    [[  6.           3.           2.        ]
     [ -8.09599972   9.59899998  20.30900002]
     [ -8.08399963   8.70699978  19.11199951]
     ..., 
     [ -2.19300008  13.62699986  15.49600029]
     [ -2.79699993  14.23499966  14.49100018]
     [ -1.76199996  12.39099979  15.30900002]]
    [  2.52720404   7.56192255  23.71602821]


Use ``Molecule.moveBy()`` to translate and center your molecule to the
origin [0, 0, 0]

.. code:: python

    mol.moveBy(-c)

Check the new center

.. code:: python

    mean(mol.get('coords'),axis=0)




.. parsed-literal::

    array([ -1.53262852e-06,  -2.38971347e-06,   3.02510853e-06], dtype=float32)



You can also rotate with ``Molecule.rotateBy``, which requires an axis
and an angle. Note that the ``uniformRandomRotation()`` function provide
the random coordinates needed for this rotation.

.. code:: python

    ligand.rotate([1, 0, 0],math.pi/2)

Visualization
-------------

The Molecule objects can be visualized either in VMD or in the Notebook
integrated javascript NGL viewer.

.. code:: python

    mol = Molecule('3PTB')
    mol.view()

Representations
---------------

It is possible to apply multiple representations to a Molecule as in
VMD. Representations use the same names as in VMD, even when using the
NGL viewer. Important parameters are: **style**, **color**, and **sel**.

There are two ways of applying representations.

The "quick" or "transient" view
-------------------------------

Use the ``view()`` method, specifying the representation as arguments.
Use the ``hold`` parameter so overlay. Representations will be cleared
on every call.

.. code:: python

    mol.view(sel='protein', style='NewCartoon', color='Index', hold=True)
    mol.view(sel='resname BEN', style='Licorice', color=1)

The "explicit" way, for which representations are "sticky"
----------------------------------------------------------

One directly manipulates elements in the ``reps`` property. Views are
stored *in* the molecule object.

.. code:: python

    mol.reps.remove()   # Clear representations
    mol.reps.add(sel='protein', style='NewCartoon', color='Index')
    mol.reps.add(sel='resname BEN', style='Licorice', color=1)
    print(mol.reps)     # Show list of representations
    mol.view()


.. parsed-literal::

    rep 0: sel='protein', style='NewCartoon', color='Index'
    rep 1: sel='resname BEN', style='Licorice', color='1'
    


Atom selection expressions work as in VMD
-----------------------------------------

This removes a slab 6 Å thick (-3 Å :math:`\le x \le` +3 Å).

.. code:: python

    mol.reps.remove()
    mol.view(sel='x*x>9')

Working with trajectories
-------------------------

Molecule provides wrapping and aligning functionallity for working with
MD trajectories and improving the visualization.

.. code:: python

    # molTraj = Molecule('data/filtered.pdb')
    # molTraj.read('data/traj.xtc')
    # molTraj.view()

A realistic case study
----------------------

.. code:: python

    # Load the 'clean' molecule once again
    mol=Molecule('3PTB')

.. code:: python

    # Identify residues in contact with the ligand BEN
    mol.get("resid",sel="name CA and same residue as protein within 4 of resname BEN")




.. parsed-literal::

    array([189, 190, 191, 192, 195, 213, 215, 216, 219, 220, 226])



.. code:: python

    # Identify duplicate residues, based on PDB's insertion attribute
    
    # The quick way
    unique(mol.get('resid',sel='insertion A'))




.. parsed-literal::

    array([184, 188, 221])



.. code:: python

    # Same operation, more explicit steps and pretty-print
    ia=mol.copy()
    ia.filter("insertion A and name CA")
    rid=ia.get('resid')        # ia.resid also works!
    rn=ia.get('resname')
    
    for f, b in zip(rn, rid):
        print(f, b)


.. parsed-literal::

    GLY 184
    GLY 188
    ALA 221


.. code:: python

    # Or, if we don't want to rely on the attribute
    dups=mol.copy()
    dups.filter("name CA and protein")
    
    rid=dups.get('resid')
    rn=dups.get('resname')
    
    nrid, count= unique(rid,return_counts=True)
    nrid[count>1]




.. parsed-literal::

    array([184, 188, 221])



.. code:: python

    count




.. parsed-literal::

    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])



.. code:: python

    # Check whether there are (numeric) holes in the sequence
    ch=mol.copy()
    ch.filter("name CA and protein")
    rid=ch.get('resid')
    rn=ch.get('resname')
    
    # 0 means duplicate residues; >1 means jumps
    deltas=diff(rid)
    print(deltas)
    rid[deltas!=1]


.. parsed-literal::

    [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
     5 1 1 1 1 1 1 1 1 2 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]


.. parsed-literal::

    /home/toni/Apps/miniconda3/lib/python3.4/site-packages/ipykernel/__main__.py:10: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 223 but corresponding boolean dimension is 222




.. parsed-literal::

    array([ 34,  67, 125, 130, 184, 188, 204, 217, 221])



.. code:: python

    # Pretty-print, more explicit
    
    # Iterate over all atoms
    for i in range(size(rid)-1):
        # If at a break...
        if(deltas[i]>1):
            # Remember that deltas[i]=rid[i+1]-rid[i]
            print(rid[i],rn[i],' followed by ',rid[i+1],rn[i+1])


.. parsed-literal::

    34 ASN  followed by  37 SER
    67 LEU  followed by  69 GLY
    125 THR  followed by  127 SER
    130 SER  followed by  132 ALA
    204 LYS  followed by  209 LEU
    217 SER  followed by  219 GLY


