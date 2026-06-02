# How to extract a representative structure per macrostate

## Goal

Get one (or more) actual coordinate frames for each macrostate of a {py:class}`~htmd.model.Model` - so you can render them in any viewer, write them to PDB for figures, or use them as starting structures for downstream simulations. {py:meth}`~htmd.model.Model.viewStates` opens VMD interactively; this how-to uses {py:meth}`~htmd.model.Model.getStates` to get the same frames as plain {py:class}`~moleculekit.molecule.Molecule` objects without needing a display.

## Minimal example

```python
# One representative frame per macrostate, returned as a list of Molecule objects
mols = model.getStates(numsamples=1)

for state_idx, m in enumerate(mols):
    m.write(f"./state_{state_idx}.pdb")
```

{py:meth}`~htmd.model.Model.getStates` samples the requested states, loads the corresponding frames out of the trajectories, **wraps and aligns** them, and returns one {py:class}`~moleculekit.molecule.Molecule` per state. Each molecule contains `numsamples` frames. Writing it directly with `m.write("file.pdb")` emits only the current frame (`mol.frame`) - see the ensemble variation below for multi-frame output.

## Parameters that matter

| Parameter on `getStates` | What it does |
| --- | --- |
| `states` | List of state indices to sample. Default: all states of `statetype`. |
| `numsamples` | Number of frames per state. Default 50. Set to 1 for a single representative; larger for an ensemble figure. |
| `statetype` | `"macro"` (default), `"micro"`, or `"cluster"`. |
| `samplemode` | For `statetype="macro"` one of `"weighted"` (default; proportional to microstate eq populations - best for "representative" frames), `"even"` (uniform across all microstates of the macrostate), or `"random"`. For `statetype="micro"` / `"cluster"` only `"random"` is meaningful. |
| `wrapsel` | Atom selection used to recentre the box. Default `"protein"`. |
| `alignsel` | Atom selection used to align all returned frames. Default `"name CA"`. Set to `None` to skip alignment. |
| `alignmol` | Reference {py:class}`~moleculekit.molecule.Molecule` to align to. Defaults to the topology of the first sim. |
| `simlist` | Override the simlist used to look up trajectories (e.g. a filtered simlist for water-stripped frames). |

## Common variations

### Many frames per state for an ensemble figure

```python
mols = model.getStates(numsamples=20)

for state_idx, m in enumerate(mols):
    m.write(f"./state_{state_idx}.pdb")           # topology (first frame)
    m.write(f"./state_{state_idx}.xtc")           # all 20 aligned snapshots
```

By default the PDB writer emits only the current frame (`mol.frame`). To get all 20 in one PDB file you can pass `m.write("state_N.pdb", frames=range(m.numFrames))` (multi-MODEL output), but most viewers handle XTC trajectory files better - write the topology to a single-frame `.pdb` and the 20 frames to a `.xtc` alongside, then load both together (`Molecule("state_0.pdb").read("state_0.xtc")`).

### Sample microstates instead of macrostates

```python
mols = model.getStates(states=[42], statetype="micro", numsamples=10)
mols[0].write("./micro_42.pdb")
```

### Highest-population microstate within each macrostate

```python
import numpy as np

micro_to_macro = model.macro_ofmicro
micro_eq       = model.msm.stationary_distribution

best_micros = []
for macro in range(model.macronum):
    micros_in_macro = np.where(micro_to_macro == macro)[0]
    best_micros.append(micros_in_macro[np.argmax(micro_eq[micros_in_macro])])

mols = model.getStates(states=best_micros, statetype="micro", numsamples=1)
```

Useful when you want the **single most-populated microstate** of each macrostate (a tighter "representative" than weighted-random sampling).

## Gotchas

- `getStates` requires every Sim in the simlist to share a **structurally equivalent topology** - the `_singleMolfile` check compares Molecules via `mol_equal(..., exceptFields=["coords"])`, so the file paths can differ as long as the atoms / bonds / charges match. If your simlist mixes genuinely different topologies, pre-filter to a common one (`simfilter`) and pass the filtered simlist via `simlist=`.
- `samplemode="weighted"` and `"even"` only make sense for `statetype="macro"` - the function falls back to random for micro / cluster.
- The macrostate assignment is fixed at the moment you call `model.markovModel(lag, macronum)`. Re-fitting with a different `macronum` re-numbers macrostates - state 0 in one fit isn't state 0 in another.
- Macrostate **0 is not always the bound / folded basin** - inspect the structures (or check `model.eqDistribution()`) before assuming.

## See also

- {doc}`How to read off and interpret an ITS plot <msm-interpret-its-plot>` - upstream: how you chose `macronum` in the first place.
- {doc}`Villin folding MSM <../tutorials/analysis/villin-folding>` - the canonical macrostate-decomposition flow.
- {py:meth}`htmd.model.Model.getStates`, {py:meth}`htmd.model.Model.sampleStates` - API references.
