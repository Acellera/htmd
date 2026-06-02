# How to drop bad trajectories from an MSM dataset

## Goal

Remove crashed, truncated, or otherwise anomalous trajectories from a {py:class}`~htmd.metricdata.MetricData` (or simlist) before clustering, so they don't pollute the implied-timescales fit.

## Minimal example - drop everything that isn't the modal length

```python
from htmd.ui import Metric, MetricSelfDistance, simlist
from glob import glob

sims = simlist(glob("./data/*/"), "./structure.pdb")

metr = Metric(sims)
metr.set(MetricSelfDistance("protein and name CA", metric="contacts"))
data = metr.project()
data.fstep = 0.1

# Drop any traj whose length isn't the most common (modal) length
data.dropTraj()
```

By default {py:meth}`~htmd.metricdata.MetricData.dropTraj` keeps only trajectories of the **statistical mode length**, i.e. the most common frame count. Anything shorter (crashed) or longer (re-runs) is dropped.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `limits` | `[low, high]` frame-count range. Trajectories outside this range are dropped. |
| `multiple` | List of integers - **keeps** trajectories whose length is a multiple of at least one of the given values; drops the rest. Useful when sims write at irregular checkpoints. |
| `idx` | Explicit list of trajectory indices to drop. Use after manual inspection. |
| `keepsims` | List of `Sim` objects to preserve. Drops everything else. Useful for bootstrap-then-recover-bad workflows. |

## Common variations

### Length-bounded drop

```python
data.dropTraj(limits=[500, 5000])    # keep only trajs with 500-5000 frames
```

### Keep only trajectories that are exact multiples of a checkpoint length

```python
data.dropTraj(multiple=[100])        # drop anything that isn't divisible by 100 frames
```

Useful when each sim is supposed to checkpoint every 100 frames - non-multiples indicate a partial / truncated run.

### Manual outlier removal after looking at trajLengths

```python
import numpy as np

# plotTrajSizes scales the x-axis by data.fstep when set - put fstep in ns
# beforehand for an axis in ns; otherwise the axis is in frames.
data.plotTrajSizes()                # eyeball the length distribution
outliers = np.where(data.trajLengths < 100)[0].tolist()
data.dropTraj(idx=outliers)
```

## Gotchas

- `dropTraj` modifies the `MetricData` **in place** - if you want to keep the original, `.copy()` first.
- After `dropTraj`, the cluster assignments (if any) are invalidated - re-cluster explicitly by calling `data.cluster(clusterobj)` again before building a model. Clustering is **not** idempotent; each call refits from scratch with the clusterobj you pass.
- When `dropTraj` runs on the **post-TICA** `MetricData` (the object returned by `tica.project(...)`), it also drops the corresponding trajectories from the linked un-reduced parent. Plain `Metric.project()` output has no `.parent` link, so this auto-sync doesn't apply there.

## See also

- {doc}`How to build a simlist from a non-standard layout <msm-simlist-custom-layout>` - the producer side.
- {doc}`How to filter trajectories with simfilter <msm-filter-trajectories>` - upstream water/ion stripping for faster projection.
- {doc}`How to bootstrap a model for error bars <msm-bootstrap-errors>` - the related sub-sampling workflow.
- {py:meth}`htmd.metricdata.MetricData.dropTraj` - API reference.
