# How to use a custom force field with `amber.build`

## Goal

Pass non-default force field, topology, or parameter files to {py:func}`htmd.builder.amber.build` so tLeap parameterizes residues that aren't covered by the built-in AMBER defaults (custom ligands, post-translational modifications, modified bases, etc.).

## Minimal example

```python
from htmd.builder import amber

amber.build(
    mol,
    outdir="./build",
    ff=["leaprc.protein.ff14SB", "leaprc.gaff2", "leaprc.water.tip3p"],
    topo=["./ligand.cif"],          # one .prepi / .prep / .in / .cif / .mol2 per residue
    param=["./ligand.frcmod"],      # one .frcmod per residue
)
```

The three arguments override the {py:func}`~htmd.builder.amber.defaultFf`, {py:func}`~htmd.builder.amber.defaultTopo`, and {py:func}`~htmd.builder.amber.defaultParam` lists respectively. The pre-populated defaults (AMBER ff14SB for protein, GAFF2 for organics, TIP3P water, ions, lipid17, etc.) are skipped entirely when you pass any of these arguments - re-add the defaults explicitly if you only want to *append*.

```{tip}
tLeap pairs parameters to atoms by the **atom-type names inside the files** - the topo and param filenames don't have to match each other for the build to succeed. However, a useful convention is to **name each pair after the residue's resname** (`LIG.cif` + `LIG.frcmod` for a residue `LIG`): when HTMD's built-in cofactor / NCAA / PTM library contains an entry with the same name, the auto-load logic detects the user-supplied `LIG.frcmod` and skips re-adding the built-in `LIG` parameters on top - avoiding duplicate-atomtype errors when your custom parameters reuse atom types that collide with the built-in set.
```

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `ff` | List of `leaprc.*` file paths or names. Overrides the entire default protein/water/ion stack. |
| `topo` | List of `.prepi` / `.prep` / `.in` / `.cif` / `.mol2` files - one per non-canonical residue. CIF/MOL2 are auto-converted to prepi via `prepgen` at build time. |
| `param` | List of `.frcmod` files containing the bond / angle / dihedral / VdW parameters for the residues in `topo`. |
| `atomtypes` | Extra GAFF atom types in `(name, mass, comment)` tuples - rarely needed. |
| `offlibraries` | Additional `.off` / `.lib` libraries to load (legacy AMBER format). |

## Common variations

### Combining custom topologies with the defaults

```python
from htmd.builder import amber

amber.build(
    mol,
    outdir="./build",
    ff=amber.defaultFf(),
    topo=amber.defaultTopo() + ["./ligand.cif"],
    param=amber.defaultParam() + ["./ligand.frcmod"],
)
```

This is what {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs`'s `out.topo_paths` / `out.frcmod_paths` give you - and how {py:func}`amber.build` consumes them in the {doc}`protein-with-ligand tutorial <../tutorials/system-prep/02-protein-ligand>`.

### Using a different protein force field

```python
amber.build(mol, outdir="./build", ff=[
    "leaprc.protein.ff19SB",
    "leaprc.gaff2",
    "leaprc.water.opc",
])
```

ff19SB + OPC is the AMBER-recommended modern stack (replacing ff14SB + TIP3P). Make sure `leaprc.water.opc` and `leaprc.protein.ff19SB` are in your tLeap search path - {py:func}`htmd.builder.amber.listFiles` lists what's available.

## Gotchas

- The default `topo` / `param` lists already contain entries for water, ions, GAFF, lipid17, RNA, DNA, etc. If you replace them, re-add the defaults you still need - tLeap will fail on residues whose parameters aren't loaded.
- `.cif` files coming from the `parameterize` tool are auto-converted to `.prepi` via `prepgen`. If `prepgen` isn't on `$PATH`, `amber.build` raises a clear error - install AmberTools or the `antechamber-unofficial` PyPI package.
- The order of `leaprc.*` files in `ff` matters: tLeap applies them in sequence, so later entries can shadow earlier ones. Put `leaprc.water.<model>` last so the water model isn't overwritten by an earlier protein/lipid leaprc.

## See also

- {doc}`Build a protein with a ligand <../tutorials/system-prep/02-protein-ligand>` - the full NCAA pipeline that emits the `topo` / `param` lists you'd pass here.
- {doc}`How to mix AMBER and OpenFF parameters <system-build-amber-openff>` - same problem but using OpenFF small-molecule parameters.
- {py:mod}`htmd.builder.amber` - full API reference.
