
Preparation of the :math:`\mu` opioid receptor with ligand
==========================================================

This is a complex build system as it has several components, the
protein, a sodium ion, the ligand and of course the membrane.

.. code:: python

    from htmd import *
    #get the files
    shutil.copytree(htmd.home()+'/data/mor','/tmp/testmor/pdb')
    os.chdir('/tmp/testmor')
    path='./01_prepare/'

.. code:: python

    %ls /tmp/testmor/pdb


.. parsed-literal::

    4dkl.pdb  ff.prm  ff.rtf  membrane80by80C36.pdb  QM-min.pdb  sod.pdb


Build
-----

.. code:: python

    #Protein 4dkl is taken from opm
    
    topos  = ['top/top_all36_prot.rtf','top/top_all36_lipid.rtf', 'top/top_water_ions.rtf','pdb/ff.rtf']
    params = ['par/par_all36_prot.prm','par/par_all36_lipid.prm', 'par/par_water_ions.prm','pdb/ff.prm']
    prot = Molecule('pdb/4dkl.pdb')
    prot.filter('protein and noh and chain B or water within 5 of (chain B and protein)')
    pcenter = np.mean(prot.get('coords','protein'),axis=0)
    prot = segmentgaps(prot,'protein','P') 
    np.unique(prot.get('segid'))
    prot = charmm.build(prot, topo=topos, param=params, outdir= path+'prot',ionize=False)
    # no need to change protonations
    #prot.view()


.. parsed-literal::

    Found  segment between resid  65  and  263
    Found  segment between resid  270  and  352
    2015-12-16 10:25:28,805 - htmd.builder.charmm - INFO - Writing out segments.
    Bond between A: [serial 3005 resid 140 resname CYS chain B segid P1]
                 B: [serial 3615 resid 217 resname CYS chain B segid P1]
    
    2015-12-16 10:25:29,269 - htmd.builder.builder - INFO - One disulfide bond was added
    2015-12-16 10:25:29,451 - htmd.builder.charmm - INFO - Starting the build.
    2015-12-16 10:25:34,239 - htmd.builder.charmm - INFO - Finished building.


.. code:: python

    #Add sodium in the receptor
    sod = Molecule('pdb/sod.pdb')
    sod.set('segid','S1')
    prot.append(sod)
    
    #Use a POPC membrane created with vmd and C36
    memb = Molecule('pdb/membrane80by80C36.pdb')
    mcenter = np.mean(memb.get('coords'),axis=0)
    memb.moveBy(pcenter-mcenter)
    mol = embed(prot,memb)
    
    #Add ligand, previously parametrized using gaamp
    lig = Molecule('pdb/QM-min.pdb') 
    lig.set('segid','L');
    lcenter=np.mean(lig.get('coords'),axis=0)
    newlcenter=[np.random.uniform(-10, 10), np.random.uniform(-10, 10),  43 ]
    lig.rotateBy(uniformRandomRotation(), lcenter)
    lig.moveBy(newlcenter-lcenter)
    mol.append(lig) 
    
    #Add water
    coo = mol.get('coords','lipids or protein')
    m = np.min(coo,axis=0) + [0,0,-5]
    M = np.max(coo,axis=0) + [0,0,20]
    mol = solvate(mol, minmax=np.vstack((m,M)))
    
    #Build
    mol = charmm.build(mol, topo=topos, param=params, outdir=path+'/build',saltconc=0.15)


.. parsed-literal::

    2015-12-16 10:26:21,835 - htmd.builder.solvate - INFO - Using water pdb file at: /shared/sdoerr/Work/pyHTMD/htmd/builder/wat.pdb
    2015-12-16 10:26:22,270 - htmd.builder.solvate - INFO - Replicating 8 water segments, 2 by 2 by 2
    Solvating: 100% (8/8) [############################################] eta 00:00 /
    2015-12-16 10:26:37,844 - htmd.builder.charmm - INFO - Writing out segments.
    Bond between A: [serial 22800 resid 140 resname CYS chain B segid P1]
                 B: [serial 24036 resid 217 resname CYS chain B segid P1]
    
    2015-12-16 10:27:04,376 - htmd.builder.builder - INFO - One disulfide bond was added
    2015-12-16 10:27:04,518 - htmd.builder.charmm - INFO - Starting the build.
    2015-12-16 10:27:05,943 - htmd.builder.charmm - INFO - Finished building.
    2015-12-16 10:27:07,021 - htmd.builder.ionize - INFO - Adding 14 anions + 0 cations for neutralizing and 70 ions for the given salt concentration.
    2015-12-16 10:27:07,318 - htmd.builder.ionize - INFO - Min distance of ions from molecule: 5A
    2015-12-16 10:27:07,318 - htmd.builder.ionize - INFO - Min distance between ions: 5A
    2015-12-16 10:27:07,319 - htmd.builder.ionize - INFO - Placing 84 ions.
    2015-12-16 10:27:33,622 - htmd.builder.charmm - INFO - Writing out segments.
    2015-12-16 10:27:56,219 - htmd.builder.charmm - INFO - Starting the build.
    2015-12-16 10:27:57,604 - htmd.builder.charmm - INFO - Finished building.


Equilibrate
-----------

.. code:: python

    from htmd.protocols.equilibration_v1 import Equilibration
    md = Equilibration()
    md.numsteps = 10000000
    md.temperature = 300
    md.reference = 'protein and resid 293'
    md.selection = 'segname L and noh'
    md.box = [-25, 25, -25, 25, 43, 45]
    md.k = 5
    md.useconstantratio = True
    md.write(path+'/build',path+'/equil')

.. code:: python

    mdx = AcemdLocal()
    mdx.submit(path+'/equil')
    mdx.wait()

Production
----------

.. code:: python

    from htmd.protocols.production_v1 import Production
    md = Production()
    md.acemd.bincoordinates = 'output.coor'
    md.acemd.extendedsystem  = 'output.xsc'
    md.acemd.binvelocities=None
    md.acemd.binindex=None
    md.acemd.run='50ns'
    md.temperature = 300
    md.reference = 'protein and resid 293'
    md.selection = 'segname L and noh'
    md.k = 5
    md.box = [-25, 25, -25, 25, -10, 45]
    md.write(path +'/equil','gen/s1')
