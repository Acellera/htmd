#########
Changelog
#########

1.8.0 (stable) - 2017/07/17
===========================

Concomitant release of 1.9.0 (latest).

**Main new features (relative to 1.6.0):**

- Adds new metric for spherical coordinates (``MetricSphericalCoordinate``).
- ``Model.markovModel``: automatic reduction of macrostates number
- Improves stability of builders, membrane tiling and equilibration protocol.
- Faster workflow from build to ``proteinprepare``.
- Psi4 version updated to 1.1
- Improves treatment of atom alternative locations
- Adds ``MetricData.sampleRegion``, which can return conformations from a specific region of data-space
- Improves file reading (CHARMM CRD CARD file format , ENT PDB files and gunzipped files).

1.6.0 (stable) - 2017/03/01
===========================

Concomitant release of 1.7.0 (latest).

**Main new features (relative to 1.4.0):**

- AdaptiveGoal implemented
- Parameterize tool released
- Apps changed for Queues
- New Metrics implemented: MetricSasa, MetricTMscore, MetricFluctuation
- Equilibration and Production protocols updated
- Documentation improved

