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
    MetricCoordinate - coordinates of an atom selection <moleculekit.projections.metriccoordinate>
    MetricDistance - (Self-)distance-based metrics between atoms selections <moleculekit.projections.metricdistance>
    MetricDihedral - Dihedral-based metrics <moleculekit.projections.metricdihedral>
    MetricRmsd - RMSD-based metric <moleculekit.projections.metricrmsd>
    MetricShell - occupancy of an atom selection (e.g. water) around another selection <moleculekit.projections.metricshell>
    MetricSecondaryStructure - secondary structure-based metric <moleculekit.projections.metricsecondarystructure>
    MetricPlumed2 - access all Plumed2 metrics (CVs) <moleculekit.projections.metricplumed2>
    MetricSasa - Solvent accessible surface area <moleculekit.projections.metricsasa>
    MetricTMscore - TMscore-based metric <moleculekit.projections.metrictmscore>
    MetricFluctuation - RMSF-based metric <moleculekit.projections.metricfluctuation>