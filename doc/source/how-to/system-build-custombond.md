# How to add a covalent crosslink by hand

## Goal

Tell tLeap to form a covalent bond between two specific atoms during the build - useful for disulfide bridges with non-standard cysteines, head-to-tail cyclic peptides, isopeptides, drug-protein covalent adducts, and any other crosslink that {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` can't infer on its own.

```{note}
You usually don't need to write `custombonds` by hand. The canonical NCAA flow ({py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` + {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs`) detects inter-residue crosslinks from the input PDB's connectivity and emits the corresponding `custombonds` list automatically as `out.custombonds`. Reach for this recipe only when detect misses the bond, when the crosslink isn't in the PDB at all (e.g. you're modelling a covalent adduct that wasn't crystallised), or when you want to override what detect found.
```

## Minimal example

```python
from htmd.builder import amber

custombonds = [
    # (atom selection for endpoint 1, atom selection for endpoint 2)
    ('segid "P0" and resid 1 and name "CA"',
     'segid "P0" and resid 50 and name "CA"'),
]

amber.build(
    mol,
    outdir="./build",
    custombonds=custombonds,
)
```

Each entry in `custombonds` is a tuple of two atom-selection strings (each must resolve to exactly one atom). At build time, `amber.build` translates them into tLeap `bond` directives.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `custombonds` | List of `(sel1, sel2)` atom-selection-string pairs. Each pair adds one `bond mol.<pos1>.<name1> mol.<pos2>.<name2>` to the generated `tleap.in`. |
| `disulfide` | Separate list of disulfide bonds with the same shape. `amber.build` auto-detects S-S within 2.5 Å, so you only need this for non-standard cases. |

## Common variations

### Head-to-tail cyclic peptide

```python
custombonds = [
    ('segid "PEP" and resid 1 and name "N"',
     'segid "PEP" and resid 10 and name "C"'),
]
amber.build(mol, outdir="./build",
            custombonds=custombonds,
            caps={"PEP": ("none", "none")})   # suppress ACE/NME caps on the cyclised segment
```

### Isopeptide crosslink (Gln side-chain to Lys ε-NH)

```python
custombonds = [
    ('segid "P0" and resid 18 and name "CD"',   # Gln18 side-chain γ-carbon
     'segid "P0" and resid 42 and name "NZ"'),  # Lys42 ε-nitrogen
]
```

You'll usually also need to template both endpoints under different resnames so each side carries the right partial charges - see {doc}`Build a stapled peptide <../tutorials/system-prep/04-stapled-peptide>` for the full flow.

### Multi-bond entry from `parameterizeFromSpecs`

```python
out = parameterizeFromSpecs(specs, prepared, outdir="./params")
amber.build(
    prepared,
    outdir="./build",
    custombonds=out.custombonds,           # list of tuples emitted by the spec inference
    topo=out.topo_paths,
    param=out.frcmod_paths,
)
```

When you use the canonical NCAA flow, `out.custombonds` already contains the right selections for every inter-residue bond detect found. You only hand-write `custombonds` for cases detect doesn't cover.

## Gotchas

- Each selection in a `custombonds` tuple must resolve to **exactly one atom**. If a selection matches zero or several atoms, `amber.build` raises a clear error. Use `mol.atomselect(sel).sum()` to debug.
- The atoms named in a bond must already exist on the two residues - tLeap won't add atoms, only bonds. If your endpoint is a side-chain atom that's normally stripped (e.g. an `OXT` of a mid-chain residue), template the residue first so the right atoms are present.
- For disulfides between residues `amber.build` already auto-detects (CYS pairs within 2.5 Å, both with free SG), you don't need a `custombonds` entry - it'll just add a duplicate `bond` directive. Use `disulfide=...` if you want to override the detection.

## See also

- {doc}`Build a stapled peptide <../tutorials/system-prep/04-stapled-peptide>` - the standard pattern for a side-chain crosslink, including the upstream NCAA templating.
- {doc}`Build a cyclic peptide <../tutorials/system-prep/05-cyclic-peptide>` - head-to-tail variant with cap suppression.
- {py:func}`htmd.builder.nonstandard.parameterizeFromSpecs` - the function that auto-emits `custombonds` for detected NCAA crosslinks.
