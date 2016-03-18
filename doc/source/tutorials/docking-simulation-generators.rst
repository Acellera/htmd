
.. code:: python

    from htmd import *
    htmd.config(viewer='ngl')
    os.chdir('/webdata/nc983hu3brda/') # Don't use this command

Using docking to generate starting poses for simulations
========================================================

Download the files for this tutorial from this
`link <http://pub.htmd.org/nc983hu3brda/bentryp.tar.gz>`_

Dock the protein with the ligand
--------------------------------

.. code:: python

    prot = Molecule('bentryp/trypsin.pdb')
    prot.center()
    lig = Molecule('bentryp/benzamidine.pdb')
    poses, scores = dock(prot, lig)

Visualize the docked poses
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    mol = Molecule()
    mol.append(prot)
    for i, p in enumerate(poses):
        mol.append(p)
    mol.view(sel='protein', style='NewCartoon', hold=True)
    mol.view(sel='resname MOL', style='Licorice', color=1)

Build systems from docked poses
-------------------------------

.. code:: python

    molbuilt = []
    for i, p in enumerate(poses):
        prot = Molecule('bentryp/trypsin.pdb')
        prot.filter('chain A and (protein or water or resname CA)')
        prot.set('segid', 'P', sel='protein and noh')
        prot.set('segid', 'W', sel='water')
        prot.set('segid', 'CA', sel='resname CA')
        prot.center()
        from htmd.molecule.util import maxDistance
        D = maxDistance(prot, 'all')
        
        ligand = p
        ligand.set('segid','L')
        ligand.set('resname','MOL')
        
        mol = Molecule(name='combo')
        mol.append(prot)
        mol.append(ligand)
        
        D = D + 15
        smol = solvate(mol, minmax=[[-D, -D, -D], [D, D, D]])
        topos  = ['top/top_all22star_prot.rtf', 'top/top_water_ions.rtf', './bentryp/benzamidine.rtf']
        params = ['par/par_all22star_prot.prm', 'par/par_water_ions.prm', './bentryp/benzamidine.prm']
    
        molbuilt.append(charmm.build(smol, topo=topos, param=params, outdir='./docked/build/{}/'.format(i+1), saltconc=0.15))
        if i==1: # For time purposes lets only build the two first
            break

.. code:: python

    from ipywidgets.widgets import Box
    w = []
    for i, m in enumerate(molbuilt):
        m.view(sel='protein', style='NewCartoon', hold=True)
        m.view(sel='water', style='Lines', hold=True)
        h = m.view(sel='resname MOL', style='Licorice', color=0)
        w.append(h)
    b = Box(children=(w[0],w[1]))
    b

Equilibrate the build systems
-----------------------------

.. code:: python

    from htmd.protocols.equilibration_v1 import Equilibration
    md = Equilibration()
    md.numsteps = 1000
    md.temperature = 298
    
    builds = natsort(glob('docked/build/*/'))
    for i, b in enumerate(builds):
        md.write(b, 'docked/equil/{}/'.format(i+1))

.. code:: python

    mdx = AcemdLocal()
    mdx.submit(glob('./docked/equil/*/'))

.. code:: python

    mdx.wait()

Create the production folder
----------------------------

.. code:: python

    from htmd.protocols.production_v1 import Production
    md = Production()
    md.acemd.run = '50ns'
    md.temperature = 298
    
    equils = natsort(glob('docked/equil/*/'))
    for i, b in enumerate(equils):
        md.write(b, 'docked/generators/{}/'.format(i+1))

.. code:: python

    mdx = AcemdLocal()
    mdx.submit(glob('./docked/generators/*/'))

.. code:: python

    mdx.wait()
