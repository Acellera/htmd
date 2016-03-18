
.. code:: python

    from htmd import *
    htmd.config(viewer='ngl')


.. parsed-literal::

    HTMD. All material of HTMD2015 will be soon made available
    
    You are on the latest HTMD version (unpackaged).


Barnase - Barstar. Building, simulation and adaptive setup.
===========================================================

.. figure:: http://pub.htmd.org/73hboiwia98hdj209jq0/barnasebarstar.png
   :align: center
   :alt: 

Download the two molecules and view them
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can look at their PDB information here and find out their PDB IDs.
Then download them using those IDs.

-  `Barnase <http://www.rcsb.org/pdb/explore.do?structureId=2f4y>`_
-  `Barstar <http://www.rcsb.org/pdb/explore/explore.do?structureId=2hxx>`_

You will need to create Molecule objects. Check the documentation on the
`Molecule <https://www.htmd.org/docs/htmd.molecule.molecule.html>`_
class.

.. code:: python

    # Load the molecules here
    barnase = Molecule('2f4y')
    barstar = Molecule('2hxx')

Filter the structures to keep only one chain of each
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize the structures in VMD. Check the names of the chains and pick
only one chain of each. You can keep the crystalized waters if you want.

.. code:: python

    # Filter here
    barnase.filter('chain A')
    barstar.filter('chain A')

Visualize the filtered structures

Assign a different chain (A, B) and segment to each protein (BRN, STR).
If you kept the waters, assign them to segid W1 and W2 for both
molecules.

Assigning waters of both molecules to the same segid can cause problems
as they have same resids.

.. code:: python

    # Assign the chains for the proteins here
    barnase.set('chain', 'A', 'protein')
    barstar.set('chain', 'B', 'protein')

.. code:: python

    # Assign the segments for the proteins and water here
    barnase.set('segid', 'BRN')
    barstar.set('segid', 'STR')
    barnase.set('segid', 'W1', 'water')
    barstar.set('segid', 'W2', 'water')

`Barstar <http://www.rcsb.org/pdb/explore/explore.do?structureId=2hxx>`_
has a modified residue for which we lack the parametrization (check
under "small molecules"). Mutate the modified Tryptophan in Barstar
(resname 4IN) to a normal Tryptophan (TRP)

.. code:: python

    # Mutate the residue here
    barstar.mutateResidue('resname 4IN', 'TRP')

Combine the proteins and center them
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a new molecule which will contain both other molecules.

.. code:: python

    # Combine here
    mol = Molecule()
    mol.append(barnase)
    mol.append(barstar)

.. code:: python

    # Center here
    mol.center()

Solvate the combined system
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find the maximum distance of the atoms from the center point.

Create a 2D minmax array.

Subtract 5 A from the min coordinates and add 5 A to the max coordinates
to add some space in the box.

.. code:: python

    from htmd.molecule.util import maxDistance
    # Calculate the maximum distance here.
    D = maxDistance(mol)
    print(D)


.. parsed-literal::

    41.7287543262


Solvate (no need to add a salt concentration)

.. code:: python

    # Solvate here
    D += 5
    smol = solvate(mol, minmax=[[-D, -D, -D],[D, D, D]])


.. parsed-literal::

    2016-01-05 09:59:59,991 - htmd.builder.solvate - INFO - Using water pdb file at: /shared/sdoerr/Work/pyHTMD/htmd/builder/wat.pdb
    2016-01-05 10:00:00,289 - htmd.builder.solvate - INFO - Replicating 8 water segments, 2 by 2 by 2
    Solvating: 100% (8/8) [############################################] eta 00:00 /


View the solvated system (can take a minute to load in VMD).

Build the solvated system in CHARMM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    # Build here using charmm
    molbuilt = charmm.build(smol, outdir='./build/')


.. parsed-literal::

    2016-01-05 10:00:10,595 - htmd.builder.charmm - INFO - Writing out segments.
    2016-01-05 10:00:27,406 - htmd.builder.charmm - INFO - Starting the build.
    2016-01-05 10:00:29,051 - htmd.builder.charmm - INFO - Finished building.
    2016-01-05 10:00:30,599 - htmd.builder.ionize - INFO - Adding 0 anions + 4 cations for neutralizing and 0 ions for the given salt concentration.
    2016-01-05 10:00:30,985 - htmd.builder.ionize - INFO - Min distance of ions from molecule: 5A
    2016-01-05 10:00:30,985 - htmd.builder.ionize - INFO - Min distance between ions: 5A
    2016-01-05 10:00:30,986 - htmd.builder.ionize - INFO - Placing 4 ions.
    2016-01-05 10:00:34,019 - htmd.builder.charmm - INFO - Writing out segments.
    2016-01-05 10:00:50,758 - htmd.builder.charmm - INFO - Starting the build.
    2016-01-05 10:00:52,321 - htmd.builder.charmm - INFO - Finished building.


Prepare the equilibration folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have built our system in a folder, we can use the
equilibration protocol to create a new directory containing all files
needed for equilibrating the system.

The number of equilibration steps is set very low here to speed up the
tutorial. In a real case you should use a larger number.

.. code:: python

    from htmd.protocols.equilibration_v1 import Equilibration
    md = Equilibration()
    md.numsteps = 1000
    md.temperature = 300
    md.write('./build', './equil')

Run the equilibration on the local GPU. Takes roughly 5 minutes.

.. code:: python

    mdx = AcemdLocal()
    mdx.submit('./equil')
    mdx.wait()

Prepare the production folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from htmd.protocols.production_v1 import Production
    md = Production()
    md.acemd.show()

.. code:: python

    md.acemd.bincoordinates = 'output.coor'
    md.acemd.extendedsystem  = 'output.xsc'
    md.acemd.binvelocities=None
    md.acemd.binindex=None
    md.acemd.run='50ns'

.. code:: python

    md.temperature = 300

.. code:: python

    md.write('./equil', './generators/s1')

Prepare adaptive
~~~~~~~~~~~~~~~~

.. code:: python

    md = AdaptiveRun()
    md.nmin=2
    md.nmax=4
    md.nepochs = 30
    md.app = AcemdLocal()
    md.metricsel1 = 'name CA and chain A'
    md.metricsel2 = 'name CA and chain B'
    md.metrictype = 'contacts'
    md.ticadim = 3
    md.run()
