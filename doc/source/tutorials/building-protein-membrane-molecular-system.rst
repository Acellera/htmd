
Building a protein-membrane molecular system (GPCR)
===================================================

Specifically we are building the :math:`\mu` opiod receptor using the
pdbid 4dkl from the OPM database.

.. code:: python

    %load_ext autoreload 
    %autoreload 2
    from htmd.console import *
    path = home() + '/data/building-protein-membrane/'


.. parsed-literal::

    HTMD Announcement
    16-Oct-2015: New version (0.0.9) of HTMD available please update with 
    	     conda update htmd
    


.. parsed-literal::

    /shared/sdoerr/Software/anaconda3/lib/python3.4/site-packages/pyEMMA-2.0.1-py3.4-linux-x86_64.egg/pyemma/coordinates/util/stat.py:31: DeprecationWarning: Call to deprecated function hist. Called from pyemma.coordinates.util.stat line 31. Please use pyemma.coordinates.histogram()
      def hist(transform, dimensions, nbins):


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


.. code:: python

    prj = 'mornap'
    salt = 0.15
    topos  = ['top/top_all36_prot.rtf','top/top_all36_lipid.rtf', 'top/top_water_ions.rtf']
    params = ['par/par_all36_prot.prm','par/par_all36_lipid.prm', 'par/par_water_ions.prm']

.. code:: python

    prot = Molecule(path + '4dkl.pdb')#from opm
    prot.filter('protein and noh and chain B or water within 5 of (chain B and protein)')
    prot = segmentgaps(prot,'protein','P') 
    pcenter = mean(prot.get('coords','protein'),axis=0)
    unique(prot.get('segid'))


.. parsed-literal::

    Found  segment between resid  65  and  263
    Found  segment between resid  270  and  352




.. parsed-literal::

    array(['A2', 'B2', 'P1', 'P2'], dtype=object)



.. code:: python

    memb = Molecule(path + 'membrane.pdb')
    mcenter = mean(memb.get('coords'),axis=0)
    memb.moveBy(pcenter-mcenter)
    mol = embed(prot,memb)

.. code:: python

    coo = mol.get('coords','lipids or protein')
    m = amin(coo,axis=0) + [0,0,-15]
    M = amax(coo,axis=0) + [0,0,15]
    mol = solvate(mol, minmax=vstack((m,M)))

.. code:: python

    topos  = ['top/top_all36_prot.rtf','top/top_all36_lipid.rtf', 'top/top_water_ions.rtf']
    params = ['par/par_all36_prot.prm','par/par_all36_lipid.prm', 'par/par_water_ions.prm']
    mol = charmm.build(mol, topo=topos, param=params, outdir='/tmp/buil',saltconc=0.15)


.. parsed-literal::

    Bond between A: [serial 17227 resid 140 resname CYS chain B segid P1]
                 B: [serial 17837 resid 217 resname CYS chain B segid P1]
    


.. code:: python

    mol.view()
