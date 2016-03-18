
Conformational analysis of CXCL12
=================================

by Gerard Martinez

Download all the required files for the tutorial from this
`link <http://pub.htmd.org/confana1036hbl2450olw/filtered.tar.gz>`_

(warning: 2.6 Gb)

You can watch the presentation `here <https://youtu.be/I9VISC29Gc4>`_

.. code:: python

    %pylab inline
    from htmd import *
    htmd.config(viewer='ngl')
    os.chdir('/webdata/confana1036hbl2450olw/')  # Skip this command.


.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib
    Videos from the HTMD2015 workshops are available on the Acellera youtube channel: https://www.youtube.com/user/acelleralive
    
    You are on the latest HTMD version (unpackaged).


1. Introduction
---------------

CXCL12 is a chemokine involved in...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  many types of cancer
-  inflammatory diseases
-  early development events

We performed...
~~~~~~~~~~~~~~~

-  ~300 simulations x 200ns each (~60 microseconds in total)
-  we filtered out the water from our trajectories

.. figure:: http://pub.htmd.org/confana1036hbl2450olw/system-protein2.png
   :align: center
   :alt: 

2. Sampling major conformational states
---------------------------------------

In our lab we have tried many different metrics to assess the overall
conformational changes of a protein. From them all, phi and psi angles
of the protein backbone (dihedrals) have been the most successful
descriptors in "blindly" capturing the major protein conformations.

In this section we will project our trajectories on the backbone
dihedrals and, we will reduce the dimensionality by using tICA and then
we will build a Markov Model to asses the major protein conformations in
equilibrium.

.. figure:: http://pub.htmd.org/confana1036hbl2450olw/conformations.png
   :align: center
   :alt: 

Calculate metrics: protein backbone dihedrals
---------------------------------------------

.. code:: python

    fsims = simlist(glob('./filtered/*/'), './filtered/filtered.pdb')


.. parsed-literal::

    2016-03-14 11:31:47,709 - htmd.simlist - INFO - Starting listing of simulations.
    Creating simlist: 100% (289/289) [#################################] eta 00:00 \
    2016-03-14 11:31:48,226 - htmd.simlist - INFO - Finished listing of simulations.


CXCL12 has a very flexible C-terminus loop as well as a transiently
disorderable N-terminal alfa helix. In this study we are not interested
in them but in the core of the chemokine. For this reason, we will
select residues from 10 to 54.

.. code:: python

    metr = Metric(fsims)
    metr.projection(MetricDihedral(protsel='protein and resid 10 to 54', sincos=True))
    data = metr.project()
    data.fstep = 0.1


.. parsed-literal::

    2016-03-14 11:31:48,274 - htmd.projections.metricdiherdal - INFO - Precalculating phi and psi angle atom selections
    2016-03-14 11:31:48,828 - htmd.projections.metricdiherdal - INFO - Finished precalculating phi and psi.
    2016-03-14 11:31:48,833 - htmd.projections.metric - INFO - Metric: Starting projection of trajectories.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done   1 out of 289 | elapsed:    2.1s remaining:  9.9min
    [Parallel(n_jobs=-2)]: Done  21 out of 289 | elapsed:    3.7s remaining:   46.7s
    [Parallel(n_jobs=-2)]: Done  48 out of 289 | elapsed:    6.7s remaining:   33.8s
    [Parallel(n_jobs=-2)]: Done  75 out of 289 | elapsed:   10.3s remaining:   29.3s
    [Parallel(n_jobs=-2)]: Done 102 out of 289 | elapsed:   13.8s remaining:   25.3s
    [Parallel(n_jobs=-2)]: Done 129 out of 289 | elapsed:   18.6s remaining:   23.1s
    [Parallel(n_jobs=-2)]: Done 156 out of 289 | elapsed:   21.7s remaining:   18.5s
    [Parallel(n_jobs=-2)]: Done 183 out of 289 | elapsed:   25.9s remaining:   15.0s
    [Parallel(n_jobs=-2)]: Done 210 out of 289 | elapsed:   29.9s remaining:   11.2s
    [Parallel(n_jobs=-2)]: Done 237 out of 289 | elapsed:   33.3s remaining:    7.3s
    [Parallel(n_jobs=-2)]: Done 264 out of 289 | elapsed:   36.4s remaining:    3.5s


.. parsed-literal::

    2016-03-14 11:32:28,183 - htmd.projections.metric - INFO - Finished projecting the trajectories.
    2016-03-14 11:32:28,184 - htmd.projections.metric - WARNING - Multiple framesteps were read from the simulations. Taking the statistical mode: 0.1ns. If it looks wrong, you can modify it by manually setting the MetricData.fstep property.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done 289 out of 289 | elapsed:   39.1s finished


.. code:: python

    data.plotTrajSizes()



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_8_0.png


Dimensionality reduction
------------------------

.. code:: python

    tica = TICA(data, 20)
    dataTica = tica.project(3)

Clustering
----------

.. code:: python

    dataTica.cluster(MiniBatchKMeans(n_clusters=200), mergesmall=5)


.. parsed-literal::

    2016-03-14 11:32:57,064 - htmd.metricdata - INFO - Mergesmall removed 0 clusters. Original ncluster 200, new ncluster 200.


MSM analysis and visualization
------------------------------

.. code:: python

    model = Model(dataTica)
    model.plotTimescales(lags=list(range(1,1000,50)))



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_14_0.png


.. code:: python

    model.markovModel(600, 8)
    eqDist = model.eqDistribution()
    print(eqDist)


.. parsed-literal::

    2016-03-14 11:33:18,149 - htmd.model - INFO - 99.7% of the data was used
    2016-03-14 11:33:18,190 - htmd.model - INFO - Number of trajectories that visited each macrostate:
    2016-03-14 11:33:18,190 - htmd.model - INFO - [ 31  30  28 289   5  86  11  10]



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_15_1.png


.. parsed-literal::

    [  3.13484342e-04   4.31193525e-03   1.41614615e-02   2.22383619e-02
       2.41110676e-02   1.50397854e-01   3.17462353e-01   4.67003483e-01]


.. code:: python

    #we can now visualize representatives for each of the equilibrium species 
    model.numsamples=1
    model.viewStates(protein=True)


.. parsed-literal::

    [Parallel(n_jobs=1)]: Done   1 jobs       | elapsed:    0.5s
    [Parallel(n_jobs=1)]: Done   2 jobs       | elapsed:    1.2s
    [Parallel(n_jobs=1)]: Done   3 jobs       | elapsed:    1.6s
    [Parallel(n_jobs=1)]: Done   4 jobs       | elapsed:    2.4s
    [Parallel(n_jobs=1)]: Done   5 jobs       | elapsed:    3.0s
    [Parallel(n_jobs=1)]: Done   6 jobs       | elapsed:    3.6s
    [Parallel(n_jobs=1)]: Done   7 jobs       | elapsed:    4.1s
    [Parallel(n_jobs=1)]: Done   8 jobs       | elapsed:    4.7s
    [Parallel(n_jobs=1)]: Done   8 out of   8 | elapsed:    4.7s finished


.. parsed-literal::

    2016-03-14 11:34:10,377 - htmd.molecule.pdbparser - WARNING - Field "serial" of PDB overflows. Your data will be truncated to 5 characters.


.. parsed-literal::

    /shared/sdoerr/Software/anaconda3/lib/python3.5/site-packages/IPython/html.py:14: ShimWarning: The `IPython.html` package has been deprecated. You should import from `notebook` instead. `IPython.html.widgets` has moved to `ipywidgets`.
      "`IPython.html.widgets` has moved to `ipywidgets`.", ShimWarning)


Statistics
----------

What are the major differences between the states X and Y?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    means = getStateStatistic(model, data, range(model.macronum))
    plt.figure()
    plt.bar(range(len(means[0])), means[7] - means[6])
    idx = np.where(np.abs(means[7] - means[6]) > 0.6)[0]
    print(data.map[idx])


.. parsed-literal::

    ['Sine of angle of resid 16 atoms: N CA C resid 17 atoms: N '
     'Cosine of angle of resid 16 atoms: N CA C resid 17 atoms: N '
     'Sine of angle of resid 16 atoms: C resid 17 atoms: N CA C '
     'Cosine of angle of resid 24 atoms: N CA C resid 25 atoms: N '
     'Sine of angle of resid 32 atoms: N CA C resid 33 atoms: N '
     'Cosine of angle of resid 32 atoms: N CA C resid 33 atoms: N '
     'Sine of angle of resid 32 atoms: C resid 33 atoms: N CA C '
     'Sine of angle of resid 44 atoms: N CA C resid 45 atoms: N '
     'Cosine of angle of resid 44 atoms: N CA C resid 45 atoms: N '
     'Sine of angle of resid 44 atoms: C resid 45 atoms: N CA C ']



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_18_1.png


.. code:: python

    # we can visualize which residues are different between states
    filtered = Molecule('./filtered/filtered.pdb')
    filtered.view(sel='protein',style='NewCartoon',hold=True)
    filtered.view(sel='resid 16 17 24 25 32 33 44 45',style='Licorice')

Mapping back
------------

Which trajectory originated the state X?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    np.where(model.macro_ofmicro == 6)




.. parsed-literal::

    (array([  4, 124]),)



.. code:: python

    _,rel = model.sampleStates([10],[10],statetype='micro')
    print(rel)


.. parsed-literal::

    [array([[ 112,  760],
           [ 204, 1704],
           [ 241,  766],
           [ 112,  601],
           [ 204, 1723],
           [ 128,  636],
           [ 241,  809],
           [  48, 1684],
           [ 204, 1819],
           [ 128, 1321]])]


.. code:: python

    print(model.data.simlist[232])


.. parsed-literal::

    id = 232
    parent = None
    input = []
    trajectory = ['./filtered/9x9/9x9-GERARD_VERYLONG_CXCL12_confAna-0-1-RND2283_9.filtered.xtc']
    molfile = ./filtered/filtered.pdb


3. Studying a defined reaction coordinate
-----------------------------------------

Revising the literature related to CXCL12, we find a paper published by
Andrea Bernini et al. (2014) where they describe the opening of a pocket
in CXCL12 located between the 2nd and 3rd beta sheet (see pictures
attached). To try to capture this phenomenon in our simulations, we will
project our trajectories along the 2nd and 3rd beta-sheet distance.

|image0| |image1|

*Figures extracted from "Searching for protein binding sites from
Molecular Dynamics simulations and paramagnetic fragment-based NMR
studies", Andrea bernini et al., 2014 Mar;1844(3):561-6. doi:
10.1016/j.bbapap.2013.12.012. Epub 2013 Dec 27*

.. |image0| image:: http://pub.htmd.org/confana1036hbl2450olw/openclose_struc.jpg
.. |image1| image:: http://pub.htmd.org/confana1036hbl2450olw/openclose_asa.png

.. code:: python

    # The first selection corresponds to beta-sheet 2 carbons alpha, the second one to beta-sheet 3 CA.
    # We specify metric='contacts' to create contact maps instead of proper distances,
    # this means: create an interatom matrix and put 1 if the distance is below cutoff; 0 otherwise. 
    metr = Metric(fsims)
    metr.projection(MetricDistance('resid 38 to 42 and noh', 'resid 22 to 28 and noh', metric='contacts'))
    data3 = metr.project()
    data3.fstep = 0.1


.. parsed-literal::

    2016-03-14 11:35:56,924 - htmd.projections.metric - INFO - Metric: Starting projection of trajectories.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done   1 out of 289 | elapsed:    1.0s remaining:  4.7min
    [Parallel(n_jobs=-2)]: Done  21 out of 289 | elapsed:    1.9s remaining:   23.8s
    [Parallel(n_jobs=-2)]: Done  48 out of 289 | elapsed:    4.0s remaining:   20.2s
    [Parallel(n_jobs=-2)]: Done  75 out of 289 | elapsed:    6.2s remaining:   17.7s
    [Parallel(n_jobs=-2)]: Done 102 out of 289 | elapsed:    8.2s remaining:   15.0s
    [Parallel(n_jobs=-2)]: Done 129 out of 289 | elapsed:   10.3s remaining:   12.7s
    [Parallel(n_jobs=-2)]: Done 156 out of 289 | elapsed:   12.7s remaining:   10.8s
    [Parallel(n_jobs=-2)]: Done 183 out of 289 | elapsed:   15.2s remaining:    8.8s
    [Parallel(n_jobs=-2)]: Done 210 out of 289 | elapsed:   16.9s remaining:    6.3s
    [Parallel(n_jobs=-2)]: Done 237 out of 289 | elapsed:   19.3s remaining:    4.2s
    [Parallel(n_jobs=-2)]: Done 264 out of 289 | elapsed:   21.4s remaining:    2.0s


.. parsed-literal::

    2016-03-14 11:36:20,817 - htmd.projections.metric - INFO - Finished projecting the trajectories.
    2016-03-14 11:36:20,818 - htmd.projections.metric - WARNING - Multiple framesteps were read from the simulations. Taking the statistical mode: 0.1ns. If it looks wrong, you can modify it by manually setting the MetricData.fstep property.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done 289 out of 289 | elapsed:   23.2s finished


.. code:: python

    # tICA projection (dimensionality reduction along the slow process)
    tica3 = TICA(data3, 20)
    dataTica3 = tica3.project(3)

.. code:: python

    # Clustering
    dataTica3.cluster(MiniBatchKMeans(n_clusters=200), mergesmall=5)


.. parsed-literal::

    2016-03-14 11:39:21,305 - htmd.metricdata - INFO - Mergesmall removed 0 clusters. Original ncluster 200, new ncluster 200.


.. code:: python

    # Plot timescales
    model3 = Model(dataTica3)
    model3.plotTimescales(lags=list(range(1,1000,50)))



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_28_0.png


.. code:: python

    # Make Markov Model. we want to pick a lagtime where the timescales are converged (timescale is flat).
    # 600 is the lagtime we want to use (600 frames is equivalent to 60ns). 4 is the number of macrostates.
    model3.markovModel(600, 4)
    eqDist = model3.eqDistribution()
    print(eqDist)


.. parsed-literal::

    2016-03-14 11:39:24,994 - htmd.model - INFO - 100.0% of the data was used
    2016-03-14 11:39:25,022 - htmd.model - INFO - Number of trajectories that visited each macrostate:
    2016-03-14 11:39:25,023 - htmd.model - INFO - [  4   3  14 289]
    2016-03-14 11:39:25,024 - htmd.model - INFO - Take care! Macro 1 has been visited only in 3 trajectories:
    2016-03-14 11:39:25,024 - htmd.model - INFO - id = 59
    parent = None
    input = []
    trajectory = ['./filtered/3x9/3x9-GERARD_VERYLONG_CXCL12_confAna-0-1-RND1251_9.filtered.xtc']
    molfile = ./filtered/filtered.pdb
    2016-03-14 11:39:25,025 - htmd.model - INFO - id = 277
    parent = None
    input = []
    trajectory = ['./filtered/10x23/10x23-GERARD_VERYLONG_CXCL12_confAna-0-1-RND9861_9.filtered.xtc']
    molfile = ./filtered/filtered.pdb
    2016-03-14 11:39:25,025 - htmd.model - INFO - id = 280
    parent = None
    input = []
    trajectory = ['./filtered/10x27/10x27-GERARD_VERYLONG_CXCL12_confAna-0-1-RND0101_9.filtered.xtc']
    molfile = ./filtered/filtered.pdb



.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_29_1.png


.. parsed-literal::

    [ 0.0061427   0.02472505  0.0253045   0.94382775]


.. code:: python

    # Visualize states
    model3.numsamples = 1
    model3.viewStates(protein=True)


.. parsed-literal::

    [Parallel(n_jobs=1)]: Done   1 jobs       | elapsed:    0.5s
    [Parallel(n_jobs=1)]: Done   2 jobs       | elapsed:    0.9s
    [Parallel(n_jobs=1)]: Done   3 jobs       | elapsed:    1.5s
    [Parallel(n_jobs=1)]: Done   4 jobs       | elapsed:    2.3s
    [Parallel(n_jobs=1)]: Done   4 out of   4 | elapsed:    2.3s finished


.. figure:: http://pub.htmd.org/confana1036hbl2450olw/conformation_open.png
   :align: center
   :alt: 

Did you see any macrostate where the pocket is open? what is the
equilibrium population probability? Let's try to find the trajectory
that produced the state...

.. code:: python

    # Map back the trajectory/ies that originated the macro. Substitute 1 for the macro that showed the pocket opening.
    # This function is giving you the microclusters that are inside a given macrocluster
    np.where(model3.macro_ofmicro ==1)




.. parsed-literal::

    (array([ 3, 25]),)



.. code:: python

    # substitute 48 for the micro number from the previous step
    # This function gives you trajectory-frame pairs that visited a given micro
    _,rel = model3.sampleStates([48],[5],statetype='micro')
    print(rel)


.. parsed-literal::

    [array([[ 260, 1881],
           [ 165, 1591],
           [ 193,  152],
           [ 249, 1555],
           [ 174,  703]])]


.. code:: python

    print(model3.data.simlist[277])


.. parsed-literal::

    id = 277
    parent = None
    input = []
    trajectory = ['./filtered/10x23/10x23-GERARD_VERYLONG_CXCL12_confAna-0-1-RND9861_9.filtered.xtc']
    molfile = ./filtered/filtered.pdb


.. code:: python

    # Calculate RMSD of the site of interest for a selected trajectory
    simus = simlist(glob('./filtered/10x23/'), './filtered/filtered.pdb')


.. parsed-literal::

    2016-03-14 11:40:14,656 - htmd.simlist - INFO - Starting listing of simulations.
    Creating simlist: 100% (1/1) [#####################################] eta --:-- /
    2016-03-14 11:40:14,659 - htmd.simlist - INFO - Finished listing of simulations.


.. code:: python

    refmol = Molecule('./filtered/filtered.pdb')
    metr = Metric(simus)
    metr.projection(MetricRmsd(refmol, 'resid 38 to 42 or resid 22 to 28 and noh', trajalnstr='protein'))
    rmsd = metr.project()


.. parsed-literal::

    2016-03-14 11:40:14,856 - htmd.projections.metric - INFO - Metric: Starting projection of trajectories.
    2016-03-14 11:40:17,703 - htmd.projections.metric - INFO - Finished projecting the trajectories.
    2016-03-14 11:40:17,705 - htmd.projections.metric - INFO - Frame step 0.1ns was read from the trajectories. If it looks wrong, redefine it by manually setting the MetricData.fstep property.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done   1 out of   1 | elapsed:    1.9s remaining:    0.0s
    [Parallel(n_jobs=-2)]: Done   1 out of   1 | elapsed:    1.9s finished


.. code:: python

    # Do you see the pocket opening at 50ns?
    plt.plot(rmsd.dat[0])
    plt.xlabel('Simulation length (frames; 0.1ns)', fontsize=10)
    plt.ylabel('RMSD (Angstroms)', fontsize=10)




.. parsed-literal::

    <matplotlib.text.Text at 0x7f284ecfa400>




.. image:: conformational-analysis-CXCL12_files/conformational-analysis-CXCL12_38_1.png


.. code:: python

    # You can also visualize the trajectory from your browser
    refmol.read('./filtered/10x23/10x23-GERARD_VERYLONG_CXCL12_confAna-0-1-RND9861_9.filtered.xtc')
    refmol.align('protein')
    refmol.view()

.. figure:: http://pub.htmd.org/confana1036hbl2450olw/view_trajectory.png
   :align: center
   :alt: 
