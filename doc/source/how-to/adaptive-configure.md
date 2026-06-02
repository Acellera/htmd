# How to configure adaptive sampling parameters

## Goal

Set up an {py:class}`~htmd.adaptive.adaptiverun.AdaptiveMD` run: choose how many simulations run per epoch, when to stop, and where input / data / filtered trajectories live on disk.

## Minimal example

```python
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.projections.tica import TICA
from moleculekit.projections.metricdistance import MetricSelfDistance
from jobqueues.localqueue import LocalGPUQueue
from htmd.clustering.kcenters import KCenter

ad = AdaptiveMD()
ad.app = LocalGPUQueue()                      # or SlurmQueue() / LsfQueue() / PBSQueue()
ad.nmin = 5                                   # minimum concurrent sims to keep running
ad.nmax = 10                                  # maximum concurrent sims at any time
ad.nepochs = 30                               # stop after this many epochs (waits for all running sims to finish)
ad.generatorspath = "./generators"            # input templates for new sims
ad.inputpath = "./input"                      # where adaptive writes per-epoch input dirs
ad.datapath = "./data"                        # where completed sims land
ad.filteredpath = "./filtered"                # auto-filtered trajectories (water stripped)
ad.projection = MetricSelfDistance("protein and name CA", metric="contacts")
ad.clustmethod = KCenter
ad.macronum = 4
ad.ticadim = 3
ad.ticalag = 20                               # in frames
ad.run()
```

`ad.run()` blocks until either `nepochs` adaptive rounds have completed or you `KeyboardInterrupt`. Each epoch: collect any newly-finished sims, project + cluster + build the model, choose new starting frames from under-explored microstates, write new input directories under `inputpath/`, submit them through `app`.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `app` | The {py:class}`~jobqueues.simqueue.SimQueue` that runs the actual ACEMD jobs. Local GPU, Slurm, LSF, PBS - same interface. |
| `nmin` | Low-water mark on the running fleet. A new epoch only fires when the running count drops to `<= nmin`. Adaptive does **not** actively keep the fleet above `nmin` between epochs. |
| `nmax` | Refill ceiling. When an epoch fires it spawns `nmax - running` new sims. |
| `nepochs`, `nframes` | Independent OR-stop conditions: the run ends when *either* the epoch count reaches `nepochs` *or* the aggregate frame count reaches `nframes`. `nframes=0` (default) disables the frame check. |
| `generatorspath` | Directory containing template input directories. Only seeds **epoch 1**; later epochs copy from parent sim input dirs instead. |
| `inputpath`, `datapath`, `filteredpath` | Per-epoch input dirs / completed sim outputs / water-stripped filtered trajs. |
| `coorname` | Name of the starting-coordinate file each new sim gets (the parent frame). Default `"input.coor"`. |
| `boxname` | Name of the starting-box file (PBC). Default `"input.xsc"`; set to `"none"` to skip writing a box file. |
| `projection` | A moleculekit `Metric*` (or list) - the feature space adaptive clusters on. |
| `ticadim`, `ticalag` | TICA dimensions to keep / lag in frames. Set `ticadim=0` to disable TICA. `ticalag` is silently clamped at runtime to `min(min(trajLengths)/2, ticalag)` (floor 2). |
| `clustmethod`, `macronum` | Clustering algorithm class + number of macrostates. |
| `skip` | Sub-sample frames before clustering - useful when sims write at high frequency. |
| `filter`, `filtersel` | Auto-strip waters from completed sims before clustering. Defaults: `True`, `"not water"`. |
| `dryrun` | Set to `True` to test the spawning logic without submitting jobs. |

## Common variations

### Conservative exploration (few sims, many epochs)

```python
ad.nmin = 2
ad.nmax = 4
ad.nepochs = 100
```

Useful when each simulation is long (e.g. 100 ns) and you want a slow, careful walk through state space.

### Aggressive exploration (large fleet, fewer epochs)

```python
ad.nmin = 50
ad.nmax = 100
ad.nepochs = 10
```

Useful when you have a large GPU cluster and per-sim cost is low - throws many short trajectories at the problem.

### Stop on aggregate sampling rather than epoch count

```python
ad.nframes = 1_000_000   # stop when we've simulated ~100 µs total at 0.1 ns/frame
ad.nepochs = 1000        # safety upper bound
```

## Gotchas

- `generatorspath` must contain **complete, runnable** ACEMD input directories - each subdir is what `acemd` would consume directly. Generators only seed **epoch 1**: adaptive picks one generator per new sim and overwrites `coorname`/`boxname` with the starting frame. From epoch 2 onwards adaptive copies the **parent sim's input directory** instead and overwrites the same two files.
- `nmin` is a **low-water mark**, not a floor. A new epoch fires only when the running fleet drops to `<= nmin`; the epoch then spawns `nmax - running` new sims to refill back to `nmax`. With `nmin=5, nmax=10` and 4 sims running, the next epoch fires (4 ≤ 5) and spawns 6 new sims (10 − 4). Adaptive does **not** actively keep the fleet above `nmin` between epochs.
- The `projection` is recomputed from scratch every epoch on **all** completed sims, so very large feature spaces (e.g. all heavy-atom self-distances on a 500-residue protein) slow each epoch. Project to Cα-only or use TICA aggressively.
- `KeyboardInterrupt` during `ad.run()` is safe - the next `ad.run()` call resumes from the last completed epoch by globbing `inputpath/e*/` to count epochs already on disk.

## See also

- {doc}`How to write a custom goal function <adaptive-goal-function>` - directed adaptive via {py:class}`~htmd.adaptive.adaptivegoal.AdaptiveGoal`.
- {doc}`How to inspect a finished adaptive run <adaptive-inspect-run>` - reconstruct the parent / child tree from `datapath/`.
- {doc}`Adaptive sampling explanation <../explanation/adaptive-sampling>` - why under-explored microstates get re-spawned.
