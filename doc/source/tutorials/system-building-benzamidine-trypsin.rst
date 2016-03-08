
.. code:: python

    from htmd import *
    htmd.config(viewer='ngl')
    os.chdir('/webdata/73hboiwia98hdj209jq0/')  # Skip this command.

System building Benzamidine Trypsin
===================================

by Stefan Doerr

Download all the required files for the tutorial from this `link <http://pub.htmd.org/73hboiwia98hdj209jq0/building.tar.gz>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can watch the presentation here:

|image0|

.. |image0| image:: http://pub.htmd.org/73hboiwia98hdj209jq0/bentryp_youtube.png
   :target: https://youtu.be/DF9cHKBX19A?t=5m42s

Obtain structures
~~~~~~~~~~~~~~~~~

.. code:: python

    glob('bentryp/*')

.. code:: python

    prot = Molecule('3PTB')
    prot = Molecule('bentryp/trypsin.pdb')
    prot.view()

Clean structures
~~~~~~~~~~~~~~~~

Crystallized water molecules and one calcium ion present in the crystal
structure were also obtained from this PDB. Remove the ligand.

.. code:: python

    prot.filter('chain A and (protein or water or resname CA)')
    prot.view()

Alternatively

.. code:: python

    prot.remove('resname BEN')

Define segments
~~~~~~~~~~~~~~~

To build a system we need to separate the molecules into separate
segments. This prevents the builder from accidentally bonding molecules
and allows us to add caps to them.

.. code:: python

    prot.set('segid', 'P', sel='protein')
    prot.set('segid', 'W', sel='water')
    prot.set('segid', 'CA', sel='resname CA')

Center the protein to the origin

.. code:: python

    prot.center()

Let's add a ligand!
~~~~~~~~~~~~~~~~~~~

.. code:: python

    ligand = Molecule('./bentryp/benzamidine.pdb')
    ligand.center()
    ligand.view(sel='resname MOL', style='Licorice')

But the ligand is located inside the protein... We would like the ligand
to: \* Be at a certain distance from the protein \* Be rotated randomly
to provide different starting conditions

Let's reposition the ligand then:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from htmd.molecule.util import uniformRandomRotation
    ligand.rotateBy(uniformRandomRotation())

This took care of the ligand rotation. But now we still need to position
it far from the protein. We need to find out the radius of the protein:

.. figure:: http://pub.htmd.org/73hboiwia98hdj209jq0/maxdist.png
   :alt: 

.. code:: python

    from htmd.molecule.util import maxDistance
    D = maxDistance(prot, 'all')
    print(D)

.. code:: python

    D += 10
    ligand.moveBy([D, 0, 0])  # Move the ligand 10 Angstrom away from the furthest protein atom in X dimension
    ligand.rotateBy(uniformRandomRotation())

Don't forget the segments
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ligand.set('segid','L')
    ligand.set('resname','MOL')

Mix it all together
~~~~~~~~~~~~~~~~~~~

.. code:: python

    mol = Molecule(name='combo')
    mol.append(prot)
    mol.append(ligand)
    mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
    mol.reps.add(sel='resname MOL', style='Licorice')
    mol.view()

Solvate
~~~~~~~

    Water is the driving force of all nature.

    -Leonardo da Vinci

.. figure:: http://pub.htmd.org/73hboiwia98hdj209jq0/waterbox.png
   :alt: 

.. code:: python

    D = D + 5
    smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
    smol.reps.add(sel='water', style='Lines')
    smol.view()

Build
~~~~~

.. code:: python

    charmm.listFiles()

Build and ionize
~~~~~~~~~~~~~~~~

.. code:: python

    topos  = ['top/top_all22star_prot.rtf', './bentryp/benzamidine.rtf']
    params = ['par/par_all22star_prot.prm', './bentryp/benzamidine.prm']
    
    molbuilt = charmm.build(smol, topo=topos, param=params, outdir='/tmp/build')
    
    molbuilt.view(sel='protein', style='NewCartoon', 
                  color='Secondary Structure', hold=True)
    molbuilt.view(sel='resname MOL', style='Licorice', hold=True)
    molbuilt.view(sel='ions', style='VDW', hold=True)
    molbuilt.view(sel='water', style='Lines')
