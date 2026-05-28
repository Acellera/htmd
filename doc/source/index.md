# HTMD

HTMD is a Python environment for **preparing, parameterising, simulating, and analysing biomolecular systems at scale**. It builds on [moleculekit](https://software.acellera.com/moleculekit/); see [HTMD and moleculekit](explanation/htmd-vs-moleculekit.md) for which package owns which feature.

::::{grid} 2
:gutter: 3

:::{grid-item-card} 🎓 Tutorials
:link: tutorials/index
:link-type: doc

Step-by-step lessons. Start here if you're new.
:::

:::{grid-item-card} 🛠 How-to guides
:link: how-to/index
:link-type: doc

Task-focused recipes. "How do I X?"
:::

:::{grid-item-card} 📖 Reference
:link: reference/index
:link-type: doc

Full API documentation, generated from source.
:::

:::{grid-item-card} 💡 Explanation
:link: explanation/index
:link-type: doc

Concepts and mental models.
:::

::::

## Headline capabilities

```{list-table}
:header-rows: 1
:widths: 22 78

* - Capability
  - What it covers
* - **System building**
  - Solvating, ionising, and parameterising proteins, ligands, membranes, non-canonical amino acids, stapled peptides, isopeptides, cyclic peptides, and disulfide-crosslinked assemblies through AMBER, OpenMM, or CHARMM. See {doc}`explanation/system-building`.
* - **MD simulation & MSM analysis**
  - Queue-aware simulation drivers plus the {py:class}`~htmd.simlist.Sim` → {py:class}`~htmd.projections.metric.Metric` → {py:class}`~htmd.projections.tica.TICA` → clustering → {py:class}`~htmd.model.Model` → {py:class}`~htmd.kinetics.Kinetics` pipeline. See {doc}`explanation/msm-workflow`.
* - **Adaptive sampling**
  - Multi-epoch simulation campaigns whose next starting frames are picked from an on-the-fly MSM. See {doc}`explanation/adaptive-sampling`.
```

## Installation

```bash
pip install acellera-htmd
```

See [Installation](installation.md) for the conda install, `uv`, and licence registration.

## Quick start: build, simulate and analyse

A complete run end to end - prepare a protein, build it under AMBER, equilibrate, produce, then project, cluster, build an MSM, and extract kinetics.

```python
from htmd.ui import *
from moleculekit.tools.preparation import systemPrepare
from acemd.protocols import setup_equilibration, setup_production
from acemd import acemd
from sklearn.cluster import MiniBatchKMeans

# 1. Build Trp-cage for AMBER.
mol = Molecule("1L2Y")
prepared, specs = systemPrepare(mol, pH=7.4)
amber.build(prepared, outdir="./build")

# 2. Equilibrate, then run production.
setup_equilibration("./build", "./equil", run="10ns")
acemd("./equil")

setup_production("./equil", "./prod", run="100ns", temperature=300)
acemd("./prod")

# 3. Project Cα-Cα distances, cluster, build the MSM, and extract kinetics.
sims = simlist(datafolders="./prod", topologies="./prod")
metr = Metric(sims)
metr.set(MetricSelfDistance("protein and name CA"))
data = metr.project()
data.cluster(MiniBatchKMeans(n_clusters=100))

model = Model(data)
model.markovModel(lag=20, macronum=4)

kin = Kinetics(model, temperature=300)
print(kin.getRates())
```

For protein-ligand and non-canonical-residue builds, see the [system-building tutorials](tutorials/system-prep/index.md); for the analysis pipeline in depth, see the [MSM analysis tutorials](tutorials/analysis/index.md).

## Citing

If you use HTMD in published work, please cite:

S. Doerr, M. J. Harvey, F. Noé, and G. De Fabritiis.
*HTMD: High-throughput molecular dynamics for molecular discovery.*
J. Chem. Theory Comput. **2016**, 12 (4), 1845-1852.
[doi:10.1021/acs.jctc.6b00049](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049)

```{toctree}
:maxdepth: 1
:hidden:

installation
tutorials/index
how-to/index
explanation/index
reference/index
```
