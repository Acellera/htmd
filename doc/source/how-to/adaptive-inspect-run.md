# How to inspect a finished adaptive run

## Goal

After an {py:class}`~htmd.adaptive.adaptiverun.AdaptiveMD` run finishes (or while one is paused), reconstruct which simulation spawned which child, how aggregate sampling grew across epochs, and how the goal score (if any) progressed.

## Minimal example

```python
from glob import glob
from htmd.simlist import simlist

sims = simlist(glob("./data/*/"),                  # one subdir per completed sim
               glob("./input/*/structure.pdb"),    # matching topology per sim
               glob("./input/*/"))                 # input dirs for adaptive-traceback

# Per-sim metadata
for s in sims[:5]:
    print(s.simid,       # numeric index assigned in load order
          s.input,       # input folder path - encodes epoch/spawn/parent in its name
          s.molfile)     # topology file(s) for this sim
```

Each {py:class}`~htmd.simlist.Sim` records the path to its `input/` directory; the directory **name** is the parent-tracking key (e.g. `e3s7_e1s2p0f120`). The leading `e<n>s<m>` is the epoch / spawn, and everything after the `_` records which parent frame the sim spawned from. See {doc}`Adaptive sampling explanation <../explanation/adaptive-sampling>` ("Naming scheme") for the full breakdown of the pattern.

## Common variations

### Group sims by epoch

```python
from collections import defaultdict
import re

by_epoch = defaultdict(list)
for s in sims:
    m = re.search(r"e(\d+)s(\d+)", s.input)
    epoch = int(m.group(1)) if m else 0
    by_epoch[epoch].append(s)

for epoch in sorted(by_epoch):
    print(f"epoch {epoch:3d}: {len(by_epoch[epoch])} sims")
```

### Reconstruct the parent → child tree

```python
from htmd.adaptive.adaptive import reconstructAdaptiveTraj

# Walk back from a current sim/frame through its adaptive ancestry.
# Signature: reconstructAdaptiveTraj(simlist, trajID) -> (mol, chain, pathlist)
mol, chain, pathlist = reconstructAdaptiveTraj(sims, sims[-1].simid)
mol.view()       # the concatenated parent-chain trajectory leading to this frame
```

{py:func}`~htmd.adaptive.adaptive.reconstructAdaptiveTraj` chains together the parent trajectory + the spawn point + the new sim's trajectory back to the original seed, so you can watch the full multi-epoch path that produced a frame of interest. It returns a 3-tuple: the concatenated {py:class}`~moleculekit.molecule.Molecule`, a list of parent-sim names in order, and the list of trajectory file paths it pieced together.

### Aggregate sampling over epochs

```python
import numpy as np

# s.numframes is a list (one entry per trajectory file); entries can be None
# if simlist didn't manage to read the frame count.
trajlens = np.array([sum(n or 0 for n in s.numframes) for s in sims])  # frames per sim
epochs   = np.array([int(re.search(r"e(\d+)s", s.input).group(1)) for s in sims])

agg_frames = np.zeros(epochs.max() + 1, dtype=int)
for e in range(epochs.max() + 1):
    agg_frames[e] = trajlens[epochs <= e].sum()

print("aggregate frames after each epoch:", agg_frames)
```

Useful for plotting "sampling vs epoch" curves to see whether your adaptive run is still exploring or has flat-lined.

### Replay the goal score across epochs

```python
import numpy as np
from moleculekit.molecule import Molecule

goal_per_epoch = []
for epoch in sorted(by_epoch):
    epoch_mols = []
    for s in by_epoch[epoch][:50]:
        m = Molecule(s.molfile)        # topology
        m.read(s.trajectory)           # all trajectory pieces for this sim
        epoch_mols.append(m)
    scores = np.concatenate([goal_function(m) for m in epoch_mols])
    goal_per_epoch.append(scores.mean())
print(goal_per_epoch)
```

For an {py:class}`~htmd.adaptive.adaptivegoal.AdaptiveGoal` run, plotting the per-epoch mean goal score tells you whether the directed component is improving the objective over time.

## Gotchas

- The epoch / spawn numbers live in the **directory name** (e.g. `e3s7_e1s2p0f120`), not in any database. If you rename or reorganise the `data/` tree post-hoc you lose the lineage.
- `sim.input` and `sim.molfile` use absolute paths from when adaptive ran - if you move the project directory, regenerate the simlist with the new paths.
- `reconstructAdaptiveTraj` walks the **simlist you pass it**: every ancestor sim must be present in that simlist (typically the simlist over `data/` or `filtered/`). If you've dropped trajectories before reconstruction it'll fail to find their parents.
- If adaptive is still running, `data/*/` will gain new subdirs while you inspect - re-run `simlist(...)` to refresh.

## See also

- {doc}`How to configure adaptive sampling <adaptive-configure>` - the producer side of the data this how-to consumes.
- {doc}`Adaptive sampling explanation <../explanation/adaptive-sampling>` - what the epoch / spawn naming actually means.
- {py:func}`htmd.adaptive.adaptive.reconstructAdaptiveTraj` - API reference for ancestor-chain reconstruction.
