#########
Changelog
#########

1.12.0 (stable) - 2018/03/22
============================

**Main new features (relative to 1.10.0):**

- Improved User Guide: tutorials updated, new refactored projections tutorial available
- Compability with NGLview 1.0 visualization
- Support for multiple queues in Slurm and LSF
- Improvements in MOL2 reading and writing
- Introduced Membrane Builder
- Parameterization tooling gradual improvements

1.10.0 (stable) - 2017/10/25
============================

**Main new features (relative to 1.8.0):**

- Parameterize reloaded
   - More robust Psi4 settings
   - New dihedral angle fitting procedure
   - Improvement on ESP partial charges calculation
   - Refactoring of the back-end
   - Parallelized local runs and improved configurations to run on clusters
- New cofactors and non-canonical residues added for Amber
- Improvements on MOL2 file reading and writing
- Update on AceCloudQueue to launch ACEMD jobs on AceCloud from HTMD

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

