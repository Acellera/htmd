# HTMD and moleculekit

HTMD and [moleculekit](https://software.acellera.com/moleculekit/) are two separate Python packages that work together. Knowing which package owns which feature makes both easier to navigate.

## moleculekit: structures and frames

[moleculekit](https://software.acellera.com/moleculekit/) is the lower-level library. It handles everything that lives inside a single molecular system:

- the {py:class}`~moleculekit.molecule.Molecule` class - per-atom arrays, the bond graph, periodic box info, and per-frame trajectory data;
- file I/O for PDB, mmCIF/BCIF, PSF, PRMTOP, XTC, DCD, TRR, NetCDF, and others;
- atom selection (string language and boolean masks);
- geometry: alignment, distances, dihedrals, RMSD/RMSF, projections, interactions;
- system preparation: protonation, non-standard residues, mutation, missing-sidechain and partial backbone repair, segmentation;
- visualization via the built-in MolStar viewer and VMD/PyMOL bridges.

If you can do it without launching a simulation, it's probably in moleculekit.

## HTMD: campaigns of simulations

HTMD builds on top of moleculekit and adds the machinery for **running and analysing many simulations** at once:

- **System building** ({py:mod}`htmd.builder` and {py:mod}`htmd.membranebuilder`): solvating, ionising, and parameterising systems through CHARMM, AMBER, or OpenMM; placing membranes; closing loops.
- **Queue-aware simulation drivers** (re-exported by {py:mod}`htmd.ui` from the `jobqueues` package - `LocalGPUQueue`, `LocalCPUQueue`, `SlurmQueue`, `LsfQueue`, `PBSQueue`): launching simulations across local GPUs or HPC queueing systems with a uniform interface.
- **Adaptive sampling** ({py:mod}`htmd.adaptive`): driving multi-epoch campaigns whose next starting frames are chosen by an MSM built from earlier epochs - see [Adaptive sampling](adaptive-sampling.md).
- **MSM analysis** ({py:mod}`htmd.simlist`, {py:mod}`htmd.metricdata`, {py:mod}`htmd.model`, {py:mod}`htmd.kinetics`, {py:mod}`htmd.clustering`): the simulation-list → projection → clustering → model → kinetics pipeline; see [MSM workflow](msm-workflow.md).
- **HTMD-specific projections** ({py:mod}`htmd.projections.metric`, {py:mod}`htmd.projections.tica`, {py:mod}`htmd.projections.gwpca`): TICA / GWPCA dimensionality reduction over the per-frame features moleculekit produces.

If it involves more than one trajectory, or a job queue, it's in HTMD.

## How they cooperate

A typical HTMD workflow imports moleculekit objects indirectly through {py:mod}`htmd.ui`, which re-exports {py:class}`~moleculekit.molecule.Molecule`, the moleculekit projection metrics, viewer helpers, and friends:

```python
from htmd.ui import *

mol = Molecule("3PTB")          # moleculekit Molecule
mol.view(viewer="molstar")      # moleculekit viewer

# moleculekit projection used as input to an HTMD AdaptiveMD run
metric = MetricDistance("name CA", "resname BEN", periodic="selections")
```

You can import directly from either package; `htmd.ui` exists for convenience so a single namespace covers everything you need in an analysis script.

## Where to look

- moleculekit topics (reading/writing structures, atom selection, projections, system preparation) - [**moleculekit docs**](https://software.acellera.com/moleculekit/).
- HTMD topics ([building](system-building.md), queues, [adaptive sampling](adaptive-sampling.md), [MSM construction](msm-workflow.md), kinetics) - [**this site**](../index.md).
