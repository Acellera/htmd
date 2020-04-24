Projections
===========
Simulation trajectories have very large dimensionality not only in terms of number of frames but also for the
intrinsic dimensionality of each frame made by all three-dimensional coordinates. It is therefore necessary to reduce
this space by projecting these coordinates into a simpler space. HTMD provides many projection types,
e.g. MetricCoordinate to only keep the coordinate of few atoms, or MetricDistance to keep the matrix distance between
two sets of atoms.

Contents:

.. toctree::
    :maxdepth: 1

    MetricData - Storage for projected data <htmd.metricdata>
    Metric - Helper class for combining Metrics for projection <htmd.projections.metric>
    MetricCoordinate - coordinates of an atom selection <../moleculekit/moleculekit.projections.metriccoordinate.rst>
    MetricDistance - (Self-)distance-based metrics between atoms selections <../moleculekit/moleculekit.projections.metricdistance.rst>
    MetricDihedral - Dihedral-based metrics <../moleculekit/moleculekit.projections.metricdihedral.rst>
    MetricRmsd - RMSD-based metric <../moleculekit/moleculekit.projections.metricrmsd.rst>
    MetricShell - occupancy of an atom selection (e.g. water) around another selection <../moleculekit/moleculekit.projections.metricshell.rst>
    MetricSecondaryStructure - secondary structure-based metric <../moleculekit/moleculekit.projections.metricsecondarystructure.rst>
    MetricPlumed2 - access all Plumed2 metrics (CVs) <../moleculekit/moleculekit.projections.metricplumed2.rst>
    MetricSasa - Solvent accessible surface area <../moleculekit/moleculekit.projections.metricsasa.rst>
    MetricTMscore - TMscore-based metric <../moleculekit/moleculekit.projections.metrictmscore.rst>
    MetricFluctuation - RMSF-based metric <../moleculekit/moleculekit.projections.metricfluctuation.rst>