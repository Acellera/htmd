
Example MSM for a protein folding process
=========================================

by Stefan Doerr

We demonstrate how to use the HTMD code for analysing a protein folding
process in the case of the protein Villin.

You can download the data and analysis file from the following links:

-  `Datasets <http://pub.htmd.org/kcjlb6b8c7v6djhb/datasets.tar.gz>`_.
   Warning: 3GB filesize.
-  `Analysis script <./protein-folding.py>`_.

Getting started
---------------

First we import the modules we are going to need for the tutorial

.. code:: python

    %pylab inline
    from htmd import *


.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib
    Videos from the HTMD2015 workshops are available on the Acellera youtube channel: https://www.youtube.com/user/acelleralive
    
    You are on the latest HTMD version (unpackaged).


Creating a simulation list
--------------------------

For the purposes of an analysis, full atom information is required by
HTMD. Therefore each simulation trajectory has to be associated with a
corresponding PDB file. Additionally, if you plan to run adaptive, each
trajectory needs to be associated with an input directory which contains
the files used to run the simulations. The simList function allows us to
make these associations

::

    sims = simlist(glob('data/*/'), glob('input/*/structure.pdb'), glob('input/*/'))

Filtering trajectories
----------------------

Typically waters are removed from the trajectories as they are not used
in the analysis and it helps speed up further calculations. These
filtered trajectories are written into a new directory. The advantage of
filtering is that if your input structures only differ in the number of
water molecules, you can produce a single unique pdb file for all
simulations once you remove the waters, which will speed up
calculations. This step is not necessary for the analysis and you can
skip it if you don't mind the slowdown in calculations. In that case
replace in the following commands the fsims variable with sims.

::

    fsims = simfilter(sims, './filtered/', filtersel='not water')

In this tutorial due to space limitations we provide only the filtered
trajectories which are stored in three separate dataset folders.
Therefore we will skip the above two commands and construct the
simulation list directly from the filtered trajectories.

.. code:: python

    os.chdir('/webdata/kcjlb6b8c7v6djhb/')
    sets = glob('datasets/*/')
    sims = []
    for s in sets:
        fsims = simlist(glob(s + '/filtered/*/'), 'datasets/1/filtered/filtered.pdb')
        sims = simmerge(sims, fsims)


.. parsed-literal::

    2016-02-25 15:09:33,806 - htmd.simlist - INFO - Starting listing of simulations.
    Creating simlist: 100% (719/719) [#################################] eta 00:01 /
    2016-02-25 15:09:35,494 - htmd.simlist - INFO - Finished listing of simulations.
    2016-02-25 15:09:35,967 - htmd.simlist - INFO - Starting listing of simulations.
    Creating simlist: 100% (710/710) [#################################] eta 00:01 -
    2016-02-25 15:09:37,775 - htmd.simlist - INFO - Finished listing of simulations.
    2016-02-25 15:09:38,225 - htmd.simlist - INFO - Starting listing of simulations.
    Creating simlist: 100% (708/708) [#################################] eta 00:00 -
    2016-02-25 15:09:39,660 - htmd.simlist - INFO - Finished listing of simulations.


Calculating metrics
-------------------

To build a Markov state model we need to project the atom coordinates
onto a lower dimensional space which can be used for clustering the
conformations into a set of states. For protein-ligand binding systems
we typically use the binary contact map between one or more atoms of the
ligand with the carbon alpha atoms of the protein. Here we have selected
to use the atoms C7 and C4 of Benzamidine and all carbon alpha atoms of
the protein. This will calculate contacts between all atoms of the two
sets.

.. code:: python

    metr = Metric(sims)
    metr.projection(MetricSelfDistance('protein and name CA', metric='contacts'))
    data = metr.project()


.. parsed-literal::

    2016-02-25 15:09:39,717 - htmd.projections.metric - INFO - Metric: Starting projection of trajectories.


.. parsed-literal::

    [Parallel(n_jobs=-2)]: Done   1 out of 2137 | elapsed:    0.4s remaining: 12.8min
    [Parallel(n_jobs=-2)]: Done 189 out of 2137 | elapsed:    8.9s remaining:  1.5min
    [Parallel(n_jobs=-2)]: Done 384 out of 2137 | elapsed:   16.7s remaining:  1.3min
    [Parallel(n_jobs=-2)]: Done 579 out of 2137 | elapsed:   24.0s remaining:  1.1min
    [Parallel(n_jobs=-2)]: Done 774 out of 2137 | elapsed:   32.3s remaining:   56.9s
    [Parallel(n_jobs=-2)]: Done 969 out of 2137 | elapsed:   41.3s remaining:   49.8s
    [Parallel(n_jobs=-2)]: Done 1164 out of 2137 | elapsed:   50.7s remaining:   42.4s
    [Parallel(n_jobs=-2)]: Done 1359 out of 2137 | elapsed:   57.7s remaining:   33.0s
    [Parallel(n_jobs=-2)]: Done 1554 out of 2137 | elapsed:  1.1min remaining:   25.0s
    [Parallel(n_jobs=-2)]: Done 1749 out of 2137 | elapsed:  1.2min remaining:   16.5s
    [Parallel(n_jobs=-2)]: Done 1944 out of 2137 | elapsed:  1.4min remaining:    8.1s
    [Parallel(n_jobs=-2)]: Done 2137 out of 2137 | elapsed:  1.5min finished


.. parsed-literal::

    2016-02-25 15:11:11,335 - htmd.projections.metric - INFO - Finished projecting the trajectories.
    2016-02-25 15:11:11,336 - htmd.projections.metric - WARNING - Multiple framesteps were read from the simulations. Taking the statistical mode: 0.1ns. If it looks wrong, you can modify it by manually setting the MetricData.fstep property.


Here we provide the frame-step in nanoseconds i.e. the time that passes
between two consecutive frames in a trajectory. This is automatically
read from the trajectories, however not all trajectories contain the
correct fstep so it can be useful to manually define it like here.

.. code:: python

    data.fstep = 0.1

Removing trajectories
---------------------

Sometimes the set of trajectories can contain trajectories of incorrect
length. These are typically corrupted trajectories and are removed.

plotTrajSizes plots all trajectory lengths sorted

.. code:: python

    data.plotTrajSizes()



.. image:: protein-folding_files/protein-folding_14_0.png


dropTraj has multiple options for removing simulations from the dataset.
Here we use it to remove all trajectories whose length is not equal to
the mode length.

.. code:: python

    data.dropTraj()


.. parsed-literal::

    2016-02-25 15:11:17,128 - htmd.metricdata - INFO - Dropped 7 trajectories from 2137 resulting in 2130




.. parsed-literal::

    array([  89,  183,  682,  693,  720, 1597, 1901])



TICA
----

TICA is a method that can be used to improve the clustering of the
conformations. This is done by projecting the data onto a
lower-dimensional space which separates well the metastable minima and
places clusters on the transition regions.

.. code:: python

    tica = TICA(data, 20)
    dataTica = tica.project(3)

Bootstrapping
-------------

If we want to bootstrap our calculations we can at this point drop a
random 20% of the trajectories and do the rest of the analysis multiple
times to see if the results are consistent. Alternatively we can keep on
using dataTica in the following commands.

.. code:: python

    dataBoot = dataTica.bootstrap(0.8)

Clustering conformations
------------------------

Once we have cleaned the dataset we proceed to cluster the
conformations.

Here we use the mini-batch kmeans clustering algorithm to procude 1000
clusters. Clusters containing fewer than 5 conformations will get merged
into their next neighbour with more than 5 conformations.

.. code:: python

    dataBoot.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=5)


.. parsed-literal::

    2016-02-25 15:12:26,905 - htmd.metricdata - INFO - Mergesmall removed 0 clusters. Original ncluster 988, new ncluster 988.


.. parsed-literal::

    /shared/sdoerr/Software/anaconda3/lib/python3.4/site-packages/sklearn/cluster/k_means_.py:1273: RuntimeWarning: init_size=300 should be larger than k=1000. Setting it to 3*k
      init_size=init_size)


Building the Markov model
-------------------------

After clustering it is time to build the Markov model.

.. code:: python

    model = Model(dataBoot)

Before constructing the Markov model we need to choose the lag-time at
which it will be built. The lag-time is typically chosen by looking at
the implied timescale (ITS) plot and selecting a lag-time at which the
top timescales start converging. By constructing Markov models at
various lag times HTMD creates a plot which shows the slowest implied
timescales of each Markov model at various lag times. If a model is
Markovian at a specific lag time, the implied timescales should stay
unchanged for any higher lag times. Therefore, given an implied
timescales plot, the user can monitor the convergence and choose the lag
time at which to construct his Markov model, typically the Markov time
which is the shortest lag time at which the timescales converge. Too
large lag times can reduce the temporal resolution of the Markov model
and can create more statistical uncertainty due to fewer transition
counts and thus instability in the implied timescales.

.. code:: python

    model.plotTimescales()


.. parsed-literal::

    25-02-16 15:13:18 pyemma.msm.estimators.implied_timescales.ImpliedTimescales[1] WARNING  Some timescales could not be computed. Timescales array is smaller than expected or contains NaNs


.. parsed-literal::

    /shared/sdoerr/Software/anaconda3/lib/python3.4/site-packages/matplotlib/scale.py:101: RuntimeWarning: invalid value encountered in less_equal
      a[a <= 0.0] = 1e-300



.. image:: protein-folding_files/protein-folding_26_2.png


After seeing the ITS plot we decide on a lag-time of 200 frames (20ns).
Additionally the ITS plot showed us that there is a separation between 4
slow timescales and the rest of the timescales which are fast. Therefore
we choose to lump our microstates together into 5 macrostates.

.. code:: python

    model.markovModel(200, 5)


.. parsed-literal::

    2016-02-25 15:13:29,407 - htmd.model - INFO - 100.0% of the data was used
    2016-02-25 15:13:29,918 - htmd.model - INFO - Number of trajectories that visited each macrostate:
    2016-02-25 15:13:29,919 - htmd.model - INFO - [ 108   39 1659  383   80]


We can also visualize the equilibrium populations of each macrostate
using the following command

.. code:: python

    model.eqDistribution()



.. image:: protein-folding_files/protein-folding_30_0.png




.. parsed-literal::

    array([  9.00269552e-04,   3.22692368e-04,   8.14469876e-01,
             8.43986143e-02,   9.99085483e-02])



Visualizing the states
----------------------

To see what the states look like we use a Matlab integration of VMD. We
load the 3 macrostates and add a ligand representation using the ligand
atomselection.

.. code:: python

    model.viewStates(protein=True)


.. parsed-literal::

    [Parallel(n_jobs=1)]: Done   1 jobs       | elapsed:    5.6s
    [Parallel(n_jobs=1)]: Done   2 jobs       | elapsed:    6.4s
    [Parallel(n_jobs=1)]: Done   3 jobs       | elapsed:    7.8s
    [Parallel(n_jobs=1)]: Done   4 jobs       | elapsed:    8.6s
    [Parallel(n_jobs=1)]: Done   5 jobs       | elapsed:    9.4s
    [Parallel(n_jobs=1)]: Done   5 out of   5 | elapsed:    9.4s finished


Calculating the kinetics
------------------------

One of the major advantages of Markov state models is that they can
provide quantitative results about the kinetics between states.

Provide the Kinetics constructor with the system temperature. It
automatically then calculates the source and sink states.

.. code:: python

    kin = Kinetics(model, temperature=360)
    print(kin.source)
    print(kin.sink)


.. parsed-literal::

    2016-02-25 15:13:46,064 - htmd.kinetics - INFO - Detecting source state...
    2016-02-25 15:13:46,999 - htmd.kinetics - INFO - Guessing the source state as the state with minimum contacts.
    2016-02-25 15:13:46,999 - htmd.kinetics - INFO - Source macro = 2
    2016-02-25 15:13:47,000 - htmd.kinetics - INFO - Detecting sink state...
    2016-02-25 15:13:47,001 - htmd.kinetics - INFO - Sink macro = 4
    2
    4


To see the rates between the source and sink states we use the getRates
method.

.. code:: python

    r = kin.getRates()
    print(r)


.. parsed-literal::

    2016-02-25 15:13:47,006 - htmd.kinetics - INFO - Calculating rates between source: 2 and sink: 4 states.
    mfpton = 2.19E+03 (ns)
    mfptoff = 2.21E+02 (ns)
    kon = 4.56E+05 (1/M 1/s)
    koff = 4.54E+06 (1/s)
    koff/kon = 9.95E+00 (M)
    kdeq = 8.15E+00 (M)
    g0eq = 1.50 (kcal/mol)
    


To plot the free energies and mean first passage times of all state use
the ``plotRates()`` method.

.. code:: python

    kin.plotRates()


.. parsed-literal::

    2016-02-25 15:13:47,116 - htmd.kinetics - INFO - Calculating rates between source: 2 and sink: 0 states.
    2016-02-25 15:13:47,223 - htmd.kinetics - INFO - Calculating rates between source: 2 and sink: 1 states.
    2016-02-25 15:13:47,331 - htmd.kinetics - INFO - Calculating rates between source: 2 and sink: 3 states.
    2016-02-25 15:13:47,444 - htmd.kinetics - INFO - Calculating rates between source: 2 and sink: 4 states.



.. image:: protein-folding_files/protein-folding_38_1.png



.. image:: protein-folding_files/protein-folding_38_2.png



.. image:: protein-folding_files/protein-folding_38_3.png


.. code:: python

    kin.plotFluxPathways()


.. parsed-literal::

    Path flux		%path	%of total	path
    0.010175403883514487	100.0%	100.0%		[2 4]
    3.1813419479578365e-06	0.0%	100.0%		[2 3 4]
    1.8313347228110152e-10	0.0%	100.0%		[2 1 3 4]
    5.834659826955252e-11	0.0%	100.0%		[2 0 3 4]



.. image:: protein-folding_files/protein-folding_39_1.png

