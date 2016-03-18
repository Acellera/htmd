
.. code:: python

    from htmd import *
    htmd.config(viewer='ngl')
    os.chdir('/webdata/73hboiwia98hdj209jq0/')  # Skip this command.

System building μ-opioid receptor
=================================

by Stefan Doerr

.. figure:: http://pub.htmd.org/73hboiwia98hdj209jq0/membrane_GPCR.jpg
   :align: center
   :alt: 
Download all the required files for the tutorial from this
`link <http://pub.htmd.org/73hboiwia98hdj209jq0/building.tar.gz>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can watch the presentation here:

`|image0| <https://youtu.be/DF9cHKBX19A?t=22m17s>`_

.. |image0| image:: http://pub.htmd.org/73hboiwia98hdj209jq0/opioid_youtube.png

Building pure proteins
~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from htmd.builder.builder import *
    topos  = ['mor/ff.rtf','top/top_all36_prot.rtf','top/top_all36_lipid.rtf','top/top_water_ions.rtf']
    params = ['mor/ff.prm','par/par_all36_prot.prm','par/par_all36_lipid.prm','par/par_water_ions.prm']
    prot = Molecule('mor/4dkl.pdb')
    prot.filter('protein and noh and chain B or water within 5 of (chain B and protein)')

Detecting segments

.. code:: python

    prot = segmentgaps(prot,'protein','P') 

Building the protein
~~~~~~~~~~~~~~~~~~~~

.. code:: python

    prot = charmm.build(prot, topo=topos, param=params, 
                        outdir='./morbuild/prot/', ionize=False)
    prot.reps.add(sel='segid P1', style='NewCartoon', color=1)
    prot.reps.add(sel='segid P2', style='NewCartoon', color=2)
    prot.view()

Add a sodium in the receptor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    sod = Molecule('mor/sod.pdb')
    sod.set('segid','S1')
    prot.append(sod)
    prot.reps.add(sel='ions', style='VDW')
    prot.view()

Embed the protein into a membrane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    memb = Molecule('mor/membrane80by80C36.pdb')

Center the membrane onto the protein center

.. code:: python

    pcenter = np.mean(prot.get('coords','protein'),axis=0)
    mcenter = np.mean(memb.get('coords'),axis=0)
    memb.moveBy(pcenter-mcenter)

And now embed

.. code:: python

    mol = embed(prot,memb)

Visualize the embedded system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
    mol.reps.add(sel='ions', style='VDW')
    mol.reps.add(sel='lipids', style='Lines')
    mol.view()

Add a ligand
~~~~~~~~~~~~

.. code:: python

    lig = Molecule('mor/QM-min.pdb') 
    lig.set('segid','L');
    lcenter = np.mean(lig.get('coords'),axis=0)
    newlcenter=[random.uniform(-10, 10), random.uniform(-10, 10),  43 ]
    lig.rotateBy(uniformRandomRotation(), lcenter)
    lig.moveBy(newlcenter-lcenter)
    mol.append(lig) 

Put it in a water box
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    coo = mol.get('coords','noh and (lipids or protein)')
    m = np.min(coo, axis=0) + [0, 0, -5]
    M = np.max(coo, axis=0) + [0, 0, 20]
    smol = solvate(mol, minmax=np.vstack((m,M)))
    smol.reps.add(sel='segid L', style='Licorice')
    smol.reps.add(sel='water', style='Lines')
    smol.view()

Build!
~~~~~~

.. code:: python

    molbuilt = charmm.build(smol, topo=topos, param=params, 
                            outdir='./morbuild/build/', saltconc=0.15)

.. code:: python

    molbuilt.view(sel='protein', style='NewCartoon', color='Secondary Structure', hold=True)
    molbuilt.view(sel='segid L', style='Licorice', hold=True)
    molbuilt.view(sel='lipids', style='Lines', hold=True)
    molbuilt.view(sel='ions', style='VDW', hold=True)
    molbuilt.view(sel='water', style='Lines')
