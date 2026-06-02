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

# Protein folding MSM: villin headpiece

**You will learn:** how to take a set of unbiased MD trajectories of a small protein and build a Markov-state model that gives you the folded / unfolded macrostates, the folding timescale, and the folding rate.

**Prerequisites:**
- HTMD installed.
- `wget` available on `$PATH` (or substitute your favourite downloader).
- Patience for a one-time ~1 GB dataset download.

## The system

The villin headpiece is a 35-residue protein with a single hydrophobic core - one of the smallest naturally-occurring autonomously folding domains. We use a set of unbiased explicit-solvent trajectories run at 360 K (above the melting point, so both folded and unfolded states are populated). Total aggregate ≈ 107 µs across 2137 short trajectories.

## The flow

The MSM pipeline you'll follow:

1. **Download + index** the filtered trajectories with {py:func}`~htmd.simlist.simlist`.
2. **Project** the geometry onto a low-dimensional distance space with {py:class}`~moleculekit.projections.metricdistance.MetricSelfDistance`.
3. **Reduce dimensions** with TICA - find the slow collective coordinates.
4. **Cluster** the TICA space into ~1000 microstates.
5. **Build** a Markov model: choose a lag time from the implied-timescales plot, then lump microstates into macrostates with PCCA.
6. **Read kinetics**: equilibrium populations, free-energy surface, mean first-passage times.

## Setup

```{code-cell} python
import os
from glob import glob
from pathlib import Path
from htmd.ui import (
    simlist, simmerge,
    Metric, MetricSelfDistance,
    TICA, Model, Kinetics,
)
from sklearn.cluster import MiniBatchKMeans
```

## Step 1 - Build a simulation list

The filtered trajectories (already stripped of waters and ions, aligned for analysis) ship as a zip on Figshare ([HTMD tutorial data, DOI 10.6084/m9.figshare.32541291](https://doi.org/10.6084/m9.figshare.32541291)) - download [protein_folding_datasets.zip](https://ndownloader.figshare.com/files/65180772) and extract it once. ~2.5 GB. The tutorial reads the absolute path from the `HTMD_TUTORIAL_DATASETS` environment variable so the same notebook works regardless of where you extracted to.

The dataset is organised as several "epochs" (independent adaptive-sampling rounds). {py:func}`~htmd.simlist.simlist` builds one list per epoch (it dedups by folder name within a call), and {py:func}`~htmd.simlist.simmerge` stitches them into a single combined list while preserving traj IDs:

```{code-cell} python
DATASETS = Path(os.environ["HTMD_TUTORIAL_DATASETS"]) / "protein_folding_datasets"
topology = str(DATASETS / "1" / "filtered")  # any epoch works - all share the same topology
sims = []
for epoch in sorted(DATASETS.glob("*/")):
    trajs = glob(os.path.join(epoch, "filtered", "*", ""))
    sims = simmerge(sims, simlist(trajs, topology))
len(sims)
```

A {py:class}`~htmd.simlist.Sim` is the unit of analysis: trajectory frames + the matching topology PDB.

## Step 2 - Project geometry to Cα distances

Folding is encoded in the network of Cα-Cα contacts: in the unfolded state the helices are apart, in the folded state they pack against each other. {py:class}`~moleculekit.projections.metricdistance.MetricSelfDistance` computes the upper-triangular distance matrix between every Cα pair, frame by frame. `metric="contacts"` thresholds those distances into 0/1 contacts (1 if below the `threshold`, 0 otherwise - default `threshold=8` Å). Contact features are cheaper and more robust than raw distances for MSMs.

```{code-cell} python
metr = Metric(sims)
metr.set(MetricSelfDistance("protein and name CA", metric="contacts", periodic="chains"))
data = metr.project()
data.fstep = 0.1  # ns per frame - tells downstream code the time axis
```

`data` is a {py:class}`~htmd.metricdata.MetricData` object: per-trajectory projected coordinates plus the metric metadata.

```{code-cell} python
data.plotTrajSizes()
```

Drop any anomalously short trajectories (sub-sampled or crashed) - they confuse the lag-time fit later:

```{code-cell} python
data.dropTraj()
```

## Step 3 - TICA: find the slow coordinates

TICA (Time-lagged Independent Component Analysis) re-projects the high-dimensional contact space onto the directions of slowest decorrelation. For an MSM, those slow directions are exactly what you want to discretise.

```{code-cell} python
tica = TICA(data, 2, units="ns")        # lag time for the TICA fit
dataTica = tica.project(3)              # keep the 3 slowest TIC coordinates
```

The lag time you pick here (2 ns) should be a few times the fastest non-trivial relaxation; for villin folding 2 ns is plenty.

## Step 4 - Cluster the TICA space

Discretise the continuous TICA space into ~1000 microstates with mini-batch k-means. Bootstrap the data first (sample 80% of trajs) - re-running with different bootstrap draws gives you deviations / error bars on the timescales and is the standard check that the model is sample-stable:

```{code-cell} python
dataBoot = dataTica.bootstrap(0.8)
dataBoot.cluster(MiniBatchKMeans(n_clusters=1000))
```

## Step 5 - Build the MSM

```{code-cell} python
model = Model(dataBoot)
model.plotTimescales(maxlag=40, units="ns")
```

The implied-timescales (ITS) plot tells you when the model becomes Markovian: pick the **shortest** lag time at which the top timescales have plateaued. For villin folding that's ~25 ns - the slow folding timescale is clearly separated from the bulk continuum.

```{code-cell} python
model.markovModel(25, 4, units="ns")
```

`25` is the MSM lag time in ns; `4` is the number of macrostates PCCA lumps the microstates into. We pick 4 because the ITS plot shows the top 3 timescales clearly separated from the rest of the spectrum - 3 slow processes partition the microstate space into 3 + 1 = 4 metastable basins. The result captures the folded basin + a couple of unfolded sub-states + the transition region.

## Step 6 - Equilibrium populations and FES

```{code-cell} python
model.eqDistribution()
```

```{code-cell} python
model.plotFES(0, 1, temperature=360)
```

The FES on TIC1 / TIC2 shows the macrostate basins. Overlay the state assignments:

```{code-cell} python
model.plotFES(0, 1, temperature=360, states=True)
```

To eyeball the macrostates as structures, run `model.viewStates(protein=True)` from an interactive session - it launches VMD with one representative frame per macrostate. (Omitted here because it needs a display.)

## Step 7 - Kinetics

```{code-cell} python
kin = Kinetics(model, temperature=360)
r = kin.getRates()
print(r)
```

`Kinetics` reads the macrostate model and computes mean-first-passage times + folding/unfolding rates between every pair of macrostates. The output table has `source → sink` MFPTs and rates in ns and 1/s.

```{code-cell} python
kin.plotRates()
```

```{code-cell} python
kin.plotFluxPathways()
```

Flux pathways visualise the dominant transition routes in TPT (transition-path theory) form - which intermediates the protein actually visits on its way from unfolded to folded.

## Parameters that matter

| Knob | Effect |
| --- | --- |
| `MetricSelfDistance(... metric="contacts")` | Contact features are coarser but more robust than raw distances - the binary thresholding makes clusters insensitive to small Cα-Cα fluctuations within a folded basin. |
| `TICA(data, 2, units="ns")` lag | Set a few times the fastest non-trivial relaxation. Too small → noise; too large → loses fast slow modes. |
| `tica.project(3)` | Number of TIC coordinates kept. 3-5 is typical; check the eigenvalue spectrum to decide. |
| `MiniBatchKMeans(n_clusters=1000)` | More clusters = finer microstate resolution but slower. 500-2000 is the usual range. |
| `model.markovModel(25, 4, units="ns")` lag | Read off the ITS plot - shortest lag at which slow timescales plateau. |
| `model.markovModel(... n_macro=4)` | PCCA macrostates. Look at the timescale gaps in the ITS plot to pick the number. |

## Gotchas

- **ITS plateau, not minimum.** Pick the lag where slow timescales *first level off*, not where they're highest. Too-long lags throw away kinetic resolution for no gain.
- **`bootstrap(0.8)`** lets you get deviations / error bars and see whether the model is sample-stable. Re-run with several bootstrap draws and compare; if different bootstraps give wildly different ITS / states, your dataset is undersampled.
- **Macrostate count is a judgement call.** Look at the timescale gap pattern - if you see a clear gap after k slow timescales, take k+1 macrostates.

## See also

- {doc}`MSM workflow explanation <../../explanation/msm-workflow>` - the conceptual pipeline behind these calls.
- {doc}`Ligand binding MSM <trypsin-benzamidine-binding>` - the same pipeline for a binding system (uses `MetricDistance` between protein and ligand instead of `MetricSelfDistance` on the protein).
- {doc}`Adaptive sampling <../adaptive/index>` - how to *generate* the trajectory set that feeds this analysis.
