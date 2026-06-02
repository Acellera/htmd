# How to filter trajectories with simfilter

## Goal

Strip waters (or any other atom selection) out of a simlist's trajectories *before* they hit the Metric / TICA / clustering pipeline. Filtered trajectories are written to disk as a new XTC set, dramatically reducing per-frame I/O and memory use for downstream analysis - and reproducing the same step the adaptive-sampling pipeline performs automatically on completed sims.

## Minimal example

```python
from glob import glob
from htmd.simlist import simlist, simfilter

sims  = simlist(glob("./data/*/"), "./structure.pdb")
fsims = simfilter(sims, "./filtered/", filtersel="not water")
# fsims now points at water-stripped XTCs under ./filtered/
```

{py:func}`~htmd.simlist.simfilter` writes a stripped XTC per input `Sim` plus a single shared topology pair (`filtered.psf` + `filtered.pdb`) in `outfolder`, then returns a new simlist where each `Sim` points at its filtered trajectory with `molfile=[filtered.psf, filtered.pdb]`. The original trajectories are untouched.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `sims` | Input simlist from {py:func}`~htmd.simlist.simlist`. |
| `outfolder` | Directory where filtered XTCs and the shared topology PDB are written. Created if missing. |
| `filtersel` | Atom-selection string passed to moleculekit's selector. Atoms matching the selection are **kept**. `"not water"` strips waters; `"protein or resname LIG"` keeps only protein + ligand. |
| `njobs` | Number of parallel workers. Defaults to `htmd.config['njobs']`, which is **1** out of the box. Pass `njobs=N` explicitly (or set `htmd.config['njobs']`) to enable parallelism. |

## Common variations

### Strip waters and ions

```python
fsims = simfilter(sims, "./filtered/", filtersel="not (water or ion)")
```

The right default for MSM analysis - water and counter-ions are noise for any conformational / binding projection, and removing them speeds projection 5-20× depending on the box.

### Keep just the protein backbone

```python
fsims = simfilter(sims, "./filtered/", filtersel="protein and backbone")
```

Useful when projection only needs backbone Cα / N / C / O atoms (e.g. RMSD-to-reference, secondary structure). Smaller filtered XTCs, faster every step downstream.

### Reuse a pre-existing filtered directory

```python
# If ./filtered/ already has the filtered trajectories from a prior run
fsims = simlist(glob("./filtered/*/"), "./filtered/filtered.pdb")
```

`simfilter` doesn't keep state - it just writes new files to disk. If you've already filtered once and want to re-analyse, build the simlist directly from `./filtered/` and skip the `simfilter` call.

## Gotchas

- `simfilter` writes filtered XTCs to disk and is the slow path (linear in total trajectory size). Run it **once** and reuse the filtered directory across analyses; don't call it inside an inner loop.
- The filtered topology pair (`filtered.psf` + `filtered.pdb`) is written from **`sims[0]`'s topology only** - `simfilter` does not check that the other sims have matching topologies, so atom-count mismatches downstream can surface silently if you mix incompatible sims. Pre-validate your simlist before filtering.
- `filtersel` applies the same selection to **every** sim. There's no per-sim filter customisation - if some sims need a different selection, split into separate simlists and `simfilter` each.
- Adaptive-sampling auto-filters completed sims as part of its `_algorithm` loop (controlled by `AdaptiveMD.filter` default `True`, `filtersel` default `"not water"`, `filteredpath` default `"filtered"`). You usually don't need to re-run `simfilter` on an adaptive dataset - just build the simlist directly from `filteredpath/`.

## See also

- {doc}`How to build a simlist from a non-standard layout <msm-simlist-custom-layout>` - the producer of the input simlist.
- {doc}`How to drop bad trajectories <msm-drop-bad-trajectories>` - the downstream cleanup once you have filtered XTCs.
- {py:func}`htmd.simlist.simfilter` - API reference.
