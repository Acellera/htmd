Projections
===========
The simulations trajectories have a very large dimensionality not only in terms of number of frames but also for the intrinsic dimensionality of each frame made by all three-dimensional coordinates. It is therefore necessary to reduce this space by projecting these coordinates into a simpler space. HTMD provides many projection types, e.g. MetricCoordinate to only keep the coordinate of few atoms, or MetricDistance to keep the matrix distance between two sets of atoms.

On top of projection methods it is also highly recommendded to use dimensionality reduction methods to further reduce the space. HTMD provides for instance time independent component analysis (TICA) and Kmeans with triangle inequality. Both perform similarly but TICA is currently the suggested method. KMeansTri is however substantially faster.

Contents:

.. toctree::
   :maxdepth: 2

   Metric used to project simulation lists <htmd.projections.metric>
   MetricData used to store all projected data <htmd.metricdata>
   MetricDistance distance based metrics between different atoms <htmd.projections.metricdistance>
   MetricSelfDistance distance based metrics on the same  atoms <htmd.projections.metricdistance>
   MetricRmsd RMSD based metrics <htmd.projections.metricrmsd>
   MetricCoordinate coordinates of a selection of atoms <htmd.projections.metriccoordinate>
   MetricDihedral internal coordinates of the backbone for proteins <htmd.projections.metricdihedral>
   MetricShell occupancy of a selection of atoms (e.g. water) around another selection <htmd.projections.metricshell>
   MetricPlumed2 access all Plumed2 projections <htmd.projections.metricplumed2>
   MetricSecondaryStructure secondary structure based metrics <htmd.projections.metricsecondarystructure>
   TICA Time independent component analysis dimensionality reduction <htmd.projections.tica>
   KMeansTri kmeans triangle inequality dimensionality reduction <htmd.projections.kmeanstri>
   GWPCA principal component analysis dimensionality reduction <htmd.projections.gwpca>
