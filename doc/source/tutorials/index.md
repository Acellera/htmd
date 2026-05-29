# Tutorials

Step-by-step lessons that take you from zero to a working result. For working with the {py:class}`~moleculekit.molecule.Molecule` class itself (loading structures, atom selection, trajectories, visualization) see the [moleculekit tutorials](https://software.acellera.com/moleculekit/tutorials/index.html).

## System building

The headline capability. Each tutorial walks through building a real, simulation-ready system end-to-end - protein preparation, segmentation, ligand placement, membrane embedding, and parameter assignment under CHARMM / AMBER. The {doc}`system-building overview <../explanation/system-building>` lays out the full feature set, including non-canonical amino acids, stapled peptides, isopeptides, and disulfide handling.

```{toctree}
:maxdepth: 1

system-prep/index
```

## Running simulations and MSM analysis

Three end-to-end MSM analyses on real systems - benchmark trypsin/benzamidine, villin folding, and CXCL12. Each tutorial walks the full simulation-list → projection → clustering → model → kinetics pipeline. The {doc}`MSM workflow concept page <../explanation/msm-workflow>` sketches that pipeline before you start.

```{toctree}
:maxdepth: 1

analysis/index
```

## Adaptive sampling

Multi-epoch campaigns whose next starting frames are picked from an on-the-fly MSM. The {doc}`adaptive sampling concept page <../explanation/adaptive-sampling>` explains the rationale.

```{toctree}
:maxdepth: 1

adaptive/index
```
