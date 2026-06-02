# How to bootstrap an MSM for error bars on timescales

## Goal

Estimate the uncertainty on implied timescales and macrostate populations by re-fitting the MSM on multiple random sub-samples of the trajectory set. Bootstrap is the standard sample-stability check: if different subsets give wildly different timescales, your dataset is undersampled.

## Minimal example

```python
import numpy as np
from htmd.model import Model
from sklearn.cluster import MiniBatchKMeans

n_boots = 10
timescales = []
for i in range(n_boots):
    boot = dataTica.bootstrap(0.8)                       # 80% random subset
    boot.cluster(MiniBatchKMeans(n_clusters=1000))
    model = Model(boot)
    model.markovModel(lag=25, macronum=4, units="ns")
    # Timescales live on the underlying deeptime MSM at model.msm
    timescales.append(model.msm.timescales() * model.data.fstep)

# Each row has length (n_active_states - 1) - active = the largest strongly-connected
# microstate set kept by counts.submodel_largest() at this lag. The count can differ
# across bootstraps if connectivity changes, so the array may be ragged - wrap in a
# list and pad if you need a rectangular array.
timescales = np.array(timescales)                        # shape (n_boots, n_active_states - 1)
slowest = timescales[:, :3]                              # top three implied timescales
print("mean top-3 timescales (ns):", slowest.mean(axis=0))
print("std  top-3 timescales (ns):", slowest.std(axis=0))
```

{py:meth}`~htmd.metricdata.MetricData.bootstrap(ratio)` returns a new `MetricData` containing a random `ratio` fraction of the trajectories (no replacement by default). The kept count is `int(floor(numtraj * ratio))` - on small N this rounds down (e.g. 7 trajs × 0.8 = 5 kept, not 5.6), so bootstrapping ratios are coarser than they look. Each bootstrap is independent - re-cluster + re-fit from scratch.

`model.msm.timescales()` (deeptime API) returns the implied timescales in **the same units as `model.lag`** - and `markovModel(..., units="ns")` stores the lag in frames after the unit conversion, so the timescales come back in frames. Multiply by `model.data.fstep` to convert to whatever time unit you set `fstep` to (ns by convention).

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `ratio` | Fraction of trajectories kept in the subset (e.g. `0.8`). The remaining trajs are dropped for that bootstrap. |
| `replacement` | `True` to allow drawing the same trajectory multiple times (classic bootstrap). `False` (default) gives a sub-sample. |
| `n_boots` | Number of bootstrap rounds to run. 10-20 is enough for ballpark error bars; 50+ for tight CIs. |

## Common variations

### Bootstrap full distributions, not just means

```python
boot_results = []
for _ in range(n_boots):
    boot = dataTica.bootstrap(0.8)
    boot.cluster(MiniBatchKMeans(n_clusters=1000))
    model = Model(boot)
    model.markovModel(lag=25, macronum=4, units="ns")
    boot_results.append({
        "ts":   model.msm.timescales() * model.data.fstep,   # in ns
        # plot=False suppresses the matplotlib window that eqDistribution opens by default
        "eq":   model.eqDistribution(plot=False),
    })

# Percentile-based 95% CI on the slowest timescale:
slow = np.array([b["ts"][0] for b in boot_results])
print(f"slowest timescale: {np.percentile(slow, 50):.1f} ns "
      f"(95% CI [{np.percentile(slow, 2.5):.1f}, {np.percentile(slow, 97.5):.1f}])")
```

### ITS plot with bootstrap error bars

```python
its_per_boot = []
for _ in range(n_boots):
    boot = dataTica.bootstrap(0.8)
    boot.cluster(MiniBatchKMeans(n_clusters=1000))
    model = Model(boot)
    # plotTimescales returns the ITS table when plot=False, results=True
    its, lags = model.plotTimescales(
        lags=[5, 10, 25, 50, 100], units="frames",
        plot=False, results=True,
    )
    its_per_boot.append(its)                              # (n_lags, n_its)

its = np.array(its_per_boot)                              # (n_boots, n_lags, n_its)
mean_its = its.mean(axis=0)
std_its  = its.std(axis=0)
```

`plot=False` suppresses the matplotlib figure (otherwise `plotTimescales` opens one per bootstrap), and `results=True` makes it return the timescales array + lag list instead of `None`.

Plot `mean_its ± std_its` per slow process vs lag time - that's the **bootstrap-version of the ITS plot**, much more informative than a single deterministic line.

## Gotchas

- `bootstrap(0.8)` samples **at the trajectory level**, not the frame level. Each bootstrap drops or keeps whole trajectories - if you have 10 trajectories, you can't get cleaner-than-10% error bars.
- `bootstrap()` only copies the trajectory subset - independence across bootstraps **only** holds if you also call `boot.cluster(...)` on every draw, so each bootstrap gets its own cluster centroids and transition counts. Reusing cluster assignments from the full data biases the result toward the full model.
- 10 bootstraps with 1000 clusters each on a multi-million-frame dataset is the slow step in this workflow. Cache `dataTica` between bootstraps; only the cluster + Model fit re-runs.
- If different bootstraps give qualitatively different macrostate **count** (PCCA picks different `n_macro`), the underlying dataset doesn't support the macrostate decomposition - reduce `n_macro` until bootstraps converge.

## See also

- {doc}`How to read off and interpret an ITS plot <msm-interpret-its-plot>` - what bootstrapped timescales look like in practice.
- {doc}`How to drop bad trajectories <msm-drop-bad-trajectories>` - relevant pre-bootstrap cleanup.
- {py:meth}`htmd.metricdata.MetricData.bootstrap` - API reference.
