---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Ligand binding MSM: trypsin + benzamidine

**You will learn:** how to take unbiased binding trajectories of a small drug-like ligand and a protein receptor and build a Markov-state model that gives you the bound macrostate, the binding pathways, the binding affinity (K<sub>d</sub>), and on / off rates.

**Prerequisites:**
- HTMD installed.
- `wget` available on `$PATH`.
- ~5 GB of disk for the trajectory download.

## The system

Trypsin is a serine protease; benzamidine is a small competitive inhibitor that docks into the S1 specificity pocket. The dataset contains short unbiased trajectories started from many random ligand poses around the protein - some find the binding pocket, most don't. Aggregate sampling is ≈ 17 µs across 852 trajectories.

The point of this analysis is to reconstruct the binding free-energy surface and rates *without ever steering the ligand* - just by counting transitions between metastable states discovered by clustering the trajectories.

## The flow

Same MSM pipeline as the {doc}`villin folding tutorial <villin-folding>` - the only thing that changes is **the projection**: for a binding problem the slow coordinate is the protein-ligand distance pattern, not the protein's own contact map.

1. Download + simlist.
2. Project with {py:class}`~moleculekit.projections.metricdistance.MetricDistance` between protein Cα atoms and ligand heavy atoms.
3. TICA to find the slow binding modes.
4. Cluster, build MSM, lump into macrostates.
5. Read out the FES, MFPTs, k<sub>on</sub> / k<sub>off</sub>, and K<sub>d</sub>.

## Setup

```{code-cell} python
import os
from glob import glob
from pathlib import Path
from htmd.ui import (
    simlist, simmerge,
    Metric, MetricDistance,
    TICA, Model, Kinetics,
)
from sklearn.cluster import MiniBatchKMeans
```

## Step 1 - Build the simlist

The trajectory bundle ships on Figshare ([HTMD tutorial data, DOI 10.6084/m9.figshare.32541291](https://doi.org/10.6084/m9.figshare.32541291)) as [ligand_binding_datasets.zip](https://ndownloader.figshare.com/files/65180823) (~3 GB).

The dataset is split into several "epochs" (adaptive-sampling rounds). {py:func}`~htmd.simlist.simlist` builds one per-epoch list (it dedups by folder name within a single call), and {py:func}`~htmd.simlist.simmerge` stitches them into a combined list while preserving traj IDs:

```{code-cell} python
DATASETS = Path(os.environ["HTMD_TUTORIAL_DATASETS"]) / "ligand_binding_datasets"
topology = str(DATASETS / "1" / "filtered")  # any epoch works - all share the same topology
sims = []
for epoch in sorted(DATASETS.glob("*/")):
    trajs = glob(os.path.join(epoch, "filtered", "*", ""))
    sims = simmerge(sims, simlist(trajs, topology))
len(sims)
```

## Step 2 - Project: protein-ligand contacts

For binding, the relevant coordinate is "which protein residue is the ligand currently in contact with". {py:class}`~moleculekit.projections.metricdistance.MetricDistance` computes the distance matrix between two atom selections; with `metric="contacts"` you get a binary contact map per frame between every protein Cα and every ligand heavy atom (1 if the distance is below the `threshold`, 0 otherwise - default `threshold=8` Å).

```{code-cell} python
metr = Metric(sims)
metr.set(MetricDistance(
    "protein and name CA",
    "resname MOL and noh",
    periodic="selections",
    metric="contacts",
))
data = metr.project()
data.fstep = 0.1
```

```{code-cell} python
data.plotTrajSizes()
data.dropTraj()
```

`periodic="selections"` handles the case where the ligand has wrapped through PBC mid-trajectory; the wrapping is undone per frame so contact distances are correct across the periodic boundary.

## Step 3 - TICA

```{code-cell} python
tica = TICA(data, 2, units="ns")
dataTica = tica.project(3)
```

## Step 4 - Cluster + MSM

```{code-cell} python
dataBoot = dataTica.bootstrap(0.8)
dataBoot.cluster(MiniBatchKMeans(n_clusters=1000))

model = Model(dataBoot)
model.plotTimescales(maxlag=15, units="ns")
```

The slow timescale for binding is much shorter than folding (ligand diffusion happens on the ns-100ns scale, not µs). Read the lag off the ITS plot - usually around 5 ns for this system.

```{code-cell} python
model.markovModel(5, 5, units="ns")
```

Five macrostates: bound + a few "encountered but not docked" intermediates + the bulk solution.

## Step 5 - Free-energy surface

```{code-cell} python
model.plotFES(0, 1, temperature=298)
```

```{code-cell} python
model.plotFES(0, 1, temperature=298, states=True)
```

The bound state usually appears as a deep basin in one corner of TIC1 / TIC2; the bulk-solution state spreads across most of the projected area; the intermediates are shallow basins between them.

To overlay representative protein + ligand snapshots from each macrostate, run `model.viewStates(ligand="resname MOL and noh")` from an interactive session - it launches VMD. (Omitted here because it needs a display.)

## Step 6 - Kinetics + K<sub>d</sub>

```{code-cell} python
kin = Kinetics(model, temperature=298, concentration=0.0037)
```

`concentration=0.0037` mol/L is the effective bulk ligand concentration in this simulation box (one ligand in the periodic box of this volume). {py:class}`~htmd.kinetics.Kinetics` uses it to convert the simulated rates into experimental units - **k<sub>on</sub>** (M⁻¹ s⁻¹), **k<sub>off</sub>** (s⁻¹), and **K<sub>d</sub> = k<sub>off</sub> / k<sub>on</sub>**.

```{code-cell} python
r = kin.getRates()
print(r)
```

Each row is `source → sink`: MFPT, rate, and (for the bulk → bound entry) k<sub>on</sub> / k<sub>off</sub> / K<sub>d</sub>. Match those against the experimental K<sub>d</sub> to validate the model.

```{code-cell} python
kin.plotRates()
```

```{code-cell} python
kin.plotFluxPathways()
```

The flux pathways show which intermediates the ligand visits on the way from solution to the bound pocket - useful for understanding the binding mechanism (e.g. is there a kinetic trap on the surface? Does the ligand approach from a specific direction?).

## Parameters that matter

| Knob | Effect |
| --- | --- |
| `MetricDistance(sel1, sel2, ... metric="contacts")` | Coarse but the right default for binding - thresholding collapses the huge unbound bulk region into a single "no contacts" state. `metric="distances"` would spread that bulk across thousands of useless microstates and overwhelm the bound-state resolution. |
| `threshold` on `MetricDistance` | Default 8 Å. Lower (e.g. 5 Å) tightens what counts as "in contact" and emphasises tighter poses; higher dilates the bound basin. |
| `concentration` on `Kinetics` | Critical for **k<sub>on</sub>** and **K<sub>d</sub>**. Compute it as (n<sub>ligands</sub> / n<sub>waters</sub>) · 55.4 mol/L - the water count tracks the real bulk volume more accurately than the box volume, which over-counts because it includes the protein's excluded volume. |
| `periodic="selections"` on the projection | **Essential** when the ligand wraps through the box during the trajectory. Skipping it produces nonsense contacts at PBC crossings. |
| `model.markovModel(lag, n_macro)` n_macro | More macrostates → more pathway resolution but harder to interpret. 4-6 is typical for binding. |

## Gotchas

- **The bulk state is huge and unstructured.** With `metric="contacts"` the unbound bulk collapses to essentially a single microstate (all-zero contact vector), which is exactly what you want for binding analysis. Don't over-interpret intermediate basins that have very small populations - they may be undersampled.
- **Bound-state validation.** Before trusting K<sub>d</sub>, open `viewStates(ligand=...)` and confirm the bound macrostate actually puts the ligand in the experimental pocket. If it's binding to the wrong site, your model is sampling a metastable mis-pose.
- **Symmetric ligands.** Benzamidine is roughly C₂v-symmetric; for ligands without that symmetry, atom-pair distances can flip when the ligand rotates 180°. Either symmetrise the contact features manually or accept that the model will distinguish the two flipped orientations as separate states.

## See also

- {doc}`MSM workflow explanation <../../explanation/msm-workflow>` - what's happening under the hood.
- {doc}`Protein folding MSM <villin-folding>` - same pipeline, different projection.
- {doc}`Adaptive sampling <../adaptive/index>` - how the binding trajectory set was generated (random starting poses + adaptive sampling targeting under-explored regions).
