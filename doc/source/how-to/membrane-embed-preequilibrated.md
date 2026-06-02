# How to embed a protein in a pre-equilibrated membrane

## Goal

Drop a protein into a membrane bilayer that has already been built and equilibrated separately (e.g. from CHARMM-GUI, a previous {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` run, or a downloaded membrane file). This is the alternative to {py:func}`buildMembrane(solute=...) <htmd.membranebuilder.build_membrane.buildMembrane>`, which builds the membrane *around* the protein from scratch.

## Minimal example

```python
from moleculekit.molecule import Molecule
from htmd.builder.builder import embed

prot = Molecule("./protein-opm-aligned.pdb")    # protein, bilayer centre at z=0
memb = Molecule("./equilibrated-membrane.pdb")  # pre-built bilayer + waters

system = embed(prot, memb, gap=1.3)             # protein into membrane
system.write("./protein-in-membrane.pdb")
```

{py:func}`~htmd.builder.builder.embed` removes any residues of the **second** molecule whose atoms come within `gap` Å of any atom of the first, then appends the first molecule onto the trimmed second. For protein-in-membrane: pass the protein as `mol1` and the membrane as `mol2` so the **lipids** clashing with the protein get removed (carving the hole the protein sits in).

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `mol1`, `mol2` | The two molecules to merge. The function removes residues of `mol2` clashing with `mol1`, then appends `mol1` onto the trimmed `mol2`. For protein-in-membrane: pass the **protein as `mol1`** and the **membrane as `mol2`** so lipids get carved away around the protein. |
| `gap` | Minimum allowed atom-atom distance in Å. Default 1.3 Å - increase to ~2.0 Å for tighter packing if your equilibrated membrane is already dense. |

## Common variations

### Lipid-only removal (keep waters)

If your pre-equilibrated membrane file has waters and you only want to remove clashing lipids (not waters that happen to sit near the protein):

```python
from htmd.builder.builder import removeLipidsInProtein

# Drop lipids that overlap the protein footprint, leave waters in place.
# Returns a (trimmed_mol2, n_removed_fragments) tuple.
trimmed_memb, n_removed = removeLipidsInProtein(prot, memb, lipidsel="lipids")
system = Molecule()
system.append(trimmed_memb)
system.append(prot)
```

{py:func}`~htmd.builder.builder.removeLipidsInProtein` is more surgical than `embed` - it only touches atoms matching `lipidsel`, leaving waters and ions alone.

### Re-centering the membrane on the protein

```python
import numpy as np

pcenter = np.mean(prot.coords[:, :2, 0], axis=0)   # XY centroid of the protein
mcenter = np.mean(memb.coords[:, :2, 0], axis=0)   # XY centroid of the membrane
memb.moveBy([pcenter[0] - mcenter[0],
             pcenter[1] - mcenter[1], 0])

system = embed(prot, memb)
```

`buildMembrane(solute=...)` does this re-centering automatically. With a pre-built membrane you usually need to do it by hand before `embed`.

### Then run the standard build

After `embed` produces a clash-free protein + membrane + (existing) waters, solvate any uncovered extramembrane regions and build under AMBER as usual:

```python
from htmd.builder.solvate import solvate
from htmd.builder import amber

# Drop the membrane's water layer first - we re-solvate to cover the protein's
# extramembrane domains too
system.remove("water", _logger=False)

lipid_mask = system.atomselect("lipid")
xy_min = system.coords[lipid_mask, :2, 0].min(axis=0)
xy_max = system.coords[lipid_mask, :2, 0].max(axis=0)
z_min  = system.coords[:, 2, 0].min() - 15
z_max  = system.coords[:, 2, 0].max() + 15
system = solvate(system, minmax=[[xy_min[0], xy_min[1], z_min],
                                  [xy_max[0], xy_max[1], z_max]])
amber.build(system, outdir="./build", ionize=True, saltconc=0.15)
```

## Gotchas

- The protein and the membrane must be in the **same coordinate frame** before `embed`. If the membrane sits at `z ∈ [50, 90]` but the protein is OPM-aligned (`z=0` is the bilayer centre), `embed` will see no overlap and the protein ends up floating in the water layer.
- `embed` removes **entire residues** of `mol2` whose any-atom comes within `gap` of `mol1`. For a membrane this means whole lipids get dropped, not partial ones — so the carved hole is the union of every clashing lipid residue (no convex-hull expansion; that's a separate function, {py:func}`~htmd.builder.builder.removeLipidsInProtein`).
- After `embed`, the membrane's water layer probably doesn't cover the protein's intracellular / extracellular domains. Always re-solvate with a membrane-aware box (XY from the lipid extent, Z padded above and below the tallest atom) before the build.
- When in doubt, prefer {py:func}`buildMembrane(solute=...) <htmd.membranebuilder.build_membrane.buildMembrane>` - it handles re-centering, carve-out, and waters in one call. Reach for `embed` only when you have a pre-equilibrated membrane you want to preserve.

## See also

- {doc}`Build a membrane-embedded protein <../tutorials/system-prep/07-membrane>` - the canonical `buildMembrane(solute=...)` path.
- {py:func}`htmd.builder.builder.embed` and {py:func}`htmd.builder.builder.removeLipidsInProtein` - API references.
