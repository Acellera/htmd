# The MSM analysis workflow

HTMD's analysis side is built around a fixed five-step pipeline. Each step has a dedicated class, and the output of one step is the input of the next. Knowing the shape of the pipeline makes the API obvious in retrospect.

```{mermaid}
flowchart LR
    A[SimList] --> B[Projection<br/>(MetricDistance, MetricRmsd, ...)]
    B --> C[Dimensionality<br/>reduction<br/>(TICA / GWPCA)]
    C --> D[Clustering<br/>(MiniBatchKMeans,<br/>KCenters, RegularGrid)]
    D --> E[Markov state model<br/>(Model)]
    E --> F[Kinetics<br/>(rates, MFPTs)]
```

## 1. SimList - enumerate trajectories

{py:func}`~htmd.simlist.simlist` walks a `data/` (or `filtered/`) directory tree and returns a list of {py:class}`~htmd.simlist.Sim` objects, each pointing at one trajectory together with its topology. This is the moral equivalent of `glob` plus topology matching, and it's the only object you carry through the rest of the pipeline.

```python
from htmd.ui import *
sims = simlist(glob("data/*/"), glob("input/*/structure.pdb"))
```

## 2. Projection - per-frame features

A **projection** maps each frame of each simulation to a feature vector. Projections come from moleculekit's `metric*` family - `MetricDistance`, `MetricRmsd`, `MetricDihedral`, `MetricSecondaryStructure`, ... - configured once and applied across the whole `SimList`:

```python
metr = Metric(sims)
metr.set(MetricDistance("name CA", "resname BEN", periodic="selections"))
data = metr.project()
```

The result is a {py:class}`~htmd.metricdata.MetricData` object: a `(n_simulations, n_frames_per_sim, n_features)` tensor plus the mapping back to individual frames.

## 3. Dimensionality reduction - keep what matters

For all but the simplest systems, the raw projection has too many features for clustering to behave well. HTMD provides:

- {py:class}`~htmd.projections.tica.TICA` - time-lagged independent component analysis, the default choice for MD.
- {py:class}`~htmd.projections.gwpca.GWPCA` - generalised weighted PCA.

Both classes take a {py:class}`~htmd.metricdata.MetricData` and return a new {py:class}`~htmd.metricdata.MetricData` in the reduced space.

## 4. Clustering - discretise the state space

Clustering assigns each frame to one of `K` microstates. HTMD wraps several clustering schemes:

- `MiniBatchKMeans` from scikit-learn (the usual default).
- {py:class}`~htmd.clustering.kcenters.KCenters` - good for low-population states.
- {py:class}`~htmd.clustering.regular.RegularGrid` - deterministic grid in feature space.

The result is a per-frame integer label stored on the same `MetricData` object.

## 5. Markov state model and kinetics

{py:class}`~htmd.model.Model` builds the Markov state model from the clustered data: counts transitions at a chosen lag time, computes the transition matrix, identifies macrostates by PCCA+, and exposes equilibrium populations and timescales.

{py:class}`~htmd.kinetics.Kinetics` consumes the `Model` and computes mean first-passage times, rate constants, and free energies between user-specified source and sink states.

```python
model = Model(data_clust)
model.markovModel(lag=20, macronum=4)
kin = Kinetics(model, temperature=300, concentration=1e-3)
kin.getRates()
```

## What this means in practice

Each tutorial in [Analysing trajectories](../tutorials/analysis/index.md) is one trip through this pipeline on a different system. If a tutorial mentions a class you don't recognise, find the box it sits in above and you'll know what its inputs and outputs look like.
