# How to build a simlist from a non-standard directory layout

## Goal

Construct an {py:class}`htmd.simlist.Sim` list for trajectories whose directory structure doesn't fit the default `data/<simname>/{topology, traj.xtc}` pattern - per-trajectory topology, multiple trajectories per simulation, manually-curated subsets, or trajectories spread across several roots.

## Minimal example

```python
from glob import glob
from htmd.simlist import simlist

sims = simlist(
    glob("./data/*/"),                  # one subdir per simulation
    "./structure.pdb",                  # one shared topology for all sims
)
```

The simplest case: every sim folder uses the same topology PDB. `simlist` discovers a moleculekit-supported trajectory inside each folder (xtc, dcd, nc/netcdf, trr, binpos, h5, ...) and emits one `Sim` per folder.

```{tip}
**When all sims share the same topology, pass a single path.** A single shared topology lets every `Sim` reference the same file, and the downstream projection / Metric pipeline only loads it once for the whole simlist. The per-sim list form is needed only when the topologies actually differ.
```

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `datafolders` | List of directories, each containing one or more trajectory files for that sim. |
| `topologies` | Either a single path (shared topology) or a list with one path per `datafolder` (per-sim topology). Folder names must match for the list form. |
| `inputfolders` | Optional - list of input directories matching `datafolders`. Required for adaptive-sampling traceback. |

## Common variations

### Per-simulation topology

```python
sims = simlist(
    glob("./data/*/"),
    glob("./input/*/structure.pdb"),    # one PDB per sim - matched by folder name
)
```

`simlist` matches by **directory basename** - `data/sim_42/` is paired with the topology under `input/sim_42/structure.pdb`. Missing topologies raise `FileNotFoundError`; **duplicate folder basenames** within a single call raise `RuntimeError`.

### Trajectories spread across multiple roots

Since duplicate folder basenames raise inside a single `simlist` call, multi-root datasets must be split per root and stitched together with {py:func}`~htmd.simlist.simmerge`:

```python
from htmd.simlist import simmerge

sims = []
for root in ["./run1/data", "./run2/data", "./run3/data"]:
    fsims = simlist(glob(f"{root}/*/"), f"{root}/../input/0/structure.pdb")
    sims = simmerge(sims, fsims)
```

This is the pattern the {doc}`villin folding <../tutorials/analysis/villin-folding>` and {doc}`trypsin-benzamidine <../tutorials/analysis/trypsin-benzamidine-binding>` tutorials use to merge multiple adaptive epochs.

### Custom file extension or naming

`simlist` walks moleculekit's trajectory-reader list (xtc, dcd, nc/netcdf, trr, binpos, h5, lh5) and stops at the **first** extension that matches inside each folder - if a folder happens to contain both `.xtc` and `.dcd`, only one set is picked up. If your trajectories live one level deeper (e.g. `data/<sim>/output/traj.xtc`):

```python
sims = simlist(
    glob("./data/*/output/"),           # point at the inner folder
    "./structure.pdb",
)
```

### A hand-curated subset

```python
# Only keep sims whose name starts with "binding_"
keep = [d for d in glob("./data/*/") if "binding_" in d]
sims = simlist(keep, "./structure.pdb")
```

Or pass {py:func}`~htmd.simlist.simfilter` post-hoc to subset by content (see {doc}`How to drop bad trajectories <msm-drop-bad-trajectories>`).

## Gotchas

- Trajectory folder names must be **unique within a single `simlist` call** - the dedup logic raises on duplicates. If you have epochs that reuse names, split per-epoch and `simmerge`.
- The topology must contain the **same number of atoms** as the trajectory's first frame. Stripping waters in one but not the other breaks the projection step downstream.
- `inputfolders` only matters for adaptive-sampling traceback ({doc}`see how-to <adaptive-inspect-run>`); skip it for static analysis.
- `simlist` is cheap: it doesn't open the trajectory files for coordinate data. Frame counts on each `Sim` may show up as `None` until projection actually reads each trajectory.

## See also

- {doc}`How to filter trajectories with simfilter <msm-filter-trajectories>` - strip waters / ions before projection.
- {doc}`How to drop bad trajectories <msm-drop-bad-trajectories>` - prune the resulting simlist.
- {doc}`Villin folding MSM <../tutorials/analysis/villin-folding>` - real example with multi-epoch merge.
- {py:mod}`htmd.simlist` - API reference.
