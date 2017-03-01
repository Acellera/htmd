Projections
===========
The simulations trajectories have a very large dimensionality not only in terms of number of frames but also for the
intrinsic dimensionality of each frame made by all three-dimensional coordinates. It is therefore necessary to reduce
this space by projecting these coordinates into a simpler space. HTMD provides many projection types,
e.g. MetricCoordinate to only keep the coordinate of few atoms, or MetricDistance to keep the matrix distance between
two sets of atoms.

Contents:

.. toctree::
    :maxdepth: 1

    MetricData - Storage for projected data <htmd.metricdata>
    Metric - Helper class for combining Metrics for projection <htmd.projections.metric>
    MetricCoordinate - coordinates of an atom selection <htmd.projections.metriccoordinate>
    MetricDistance - (Self-)distance-based metrics between atoms selections <htmd.projections.metricdistance>
    MetricDihedral - Dihedral-based metrics <htmd.projections.metricdihedral>
    MetricRmsd - RMSD-based metric <htmd.projections.metricrmsd>
    MetricShell - occupancy of an atom selection (e.g. water) around another selection <htmd.projections.metricshell>
    MetricSecondaryStructure - secondary structure-based metric <htmd.projections.metricsecondarystructure>
    MetricPlumed2 - access all Plumed2 metrics (CVs) <htmd.projections.metricplumed2>
    MetricSasa - Solvent accessible surface area <htmd.projections.metricsasa>
    MetricTMscore - TMscore-based metric <htmd.projections.metrictmscore>
    MetricFluctuation - RMSF-based metric <htmd.projections.metricfluctuation>