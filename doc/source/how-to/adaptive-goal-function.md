# How to write a custom goal function for AdaptiveGoal

## Goal

Drive {py:class}`~htmd.adaptive.adaptivegoal.AdaptiveGoal` toward a specific objective (a target secondary structure, a known binding pocket geometry, an RMSD to a reference, ...) by writing a goal function that scores each frame against that objective. Adaptive then biases new simulations to spawn from high-scoring states (directed sampling) blended with the usual under-explored-microstate exploration.

## Minimal example

```python
import numpy as np
from htmd.adaptive.adaptivegoal import AdaptiveGoal
from moleculekit.molecule import Molecule
from moleculekit.projections.metricdistance import MetricSelfDistance
from moleculekit.projections.metricsecondarystructure import MetricSecondaryStructure
from jobqueues.localqueue import LocalGPUQueue

# Reference: per-residue secondary structure of the crystal target
crystal = Molecule("crystal.pdb")
crystal_ss = MetricSecondaryStructure().project(crystal)[0]   # shape (n_residues,)


def goal_function(mol, crystal_ss=crystal_ss):
    """Score = fraction of residues whose SS matches the crystal.
    Higher is better. Returns one score per frame.
    """
    ss = MetricSecondaryStructure().project(mol)               # (n_frames, n_residues)
    return (ss == crystal_ss).mean(axis=1)


ad = AdaptiveGoal()
ad.app = LocalGPUQueue()
ad.nmin, ad.nmax, ad.nepochs = 5, 10, 30
ad.projection = MetricSelfDistance("protein and name CA", metric="contacts")
ad.goalfunction = goal_function
ad.ucscale = 0.5                                               # 50% undirected / 50% directed
ad.run()
```

`goal_function(mol)` takes a {py:class}`~moleculekit.molecule.Molecule` (one or many frames) and returns a 1-D NumPy array with one score per frame. Higher scores are "better" - adaptive uses them to pick which sampled microstates to re-spawn from.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `goalfunction` | A callable `mol -> np.ndarray` of shape `(n_frames,)`. Closures / `functools.partial` are fine. |
| `ucscale` | Mixing ratio between **undirected** (exploration of under-sampled microstates) and **directed** (goal-driven) components. `0` = pure goal, `1` = pure exploration, `0.5` = balanced. |
| `statetype` | `"micro"` (default) or `"macro"` - which state granularity scores get aggregated to. (`"cluster"` is inherited from `AdaptiveMD` but **not supported** here - the directed-component code only handles `micro` / `macro`.) |
| `autoscale` | When `True`, adjusts `ucscale` automatically based on how stuck the run is on its goal score. Companion knobs: `autoscalediff` (default 10, epochs window), `autoscalemult` (default 1, ucscale step), `autoscaletol` (default 0.2, goal-improvement tolerance). |
| `savegoal` | Optional path to a `.pkl` file. If set, AdaptiveGoal pickles the projected goal values per trajectory each epoch - useful when iterating on the goal function. |

## Common variations

### RMSD-to-reference goal

```python
from moleculekit.projections.metricrmsd import MetricRmsd

ref = Molecule("target.pdb")

def goal_rmsd(mol, ref=ref):
    # Higher score for lower RMSD: 1 / (1 + RMSD)
    rmsd = MetricRmsd(ref, "protein and name CA").project(mol)
    return 1.0 / (1.0 + rmsd)
```

### Multi-criteria goal

```python
def goal_combined(mol):
    rmsd = MetricRmsd(ref, "protein and name CA").project(mol)
    sasa = MetricSasa(sel="resname LIG").project(mol).sum(axis=1)
    rmsd_score = 1.0 / (1.0 + rmsd)
    sasa_score = 1.0 - np.clip(sasa / 1000.0, 0, 1)     # bury the ligand
    return 0.6 * rmsd_score + 0.4 * sasa_score
```

Combine projections by hand inside the goal; adaptive just sees the final scalar score.

### Use cached reference data

```python
from functools import partial

ref_features = MetricSecondaryStructure().project(Molecule("ref.pdb"))[0]

def goal(mol, ref):
    ss = MetricSecondaryStructure().project(mol)
    return (ss == ref).mean(axis=1)

ad.goalfunction = partial(goal, ref=ref_features)
```

`partial` lets you pre-bind reference state so the goal doesn't recompute it for every adaptive epoch.

## Gotchas

- The goal function must return **one score per frame**. Returning `(n_frames, k)` raises a confusing error inside the spawn logic.
- Adaptive calls the goal once per epoch on **every accumulated frame**, so an O(n²) goal becomes prohibitive after a few epochs. Vectorise / cache.
- Don't make the goal too sharp - if only 0.01% of frames score above zero, the directed component degenerates to the highest-scoring single frame and you lose diversity. Use a smooth scoring function (`1 / (1 + x)`, tanh-shaped, etc.).
- `ucscale=0` (pure goal) tends to over-exploit one basin. `ucscale=0.5` is a good default; `autoscale=True` adapts it dynamically.

## See also

- {doc}`How to configure adaptive sampling <adaptive-configure>` - all the `AdaptiveMD` knobs that `AdaptiveGoal` inherits.
- {doc}`Adaptive sampling explanation <../explanation/adaptive-sampling>` - how directed + undirected components combine in FAST.
- {py:class}`htmd.adaptive.adaptivegoal.AdaptiveGoal` - API reference.
