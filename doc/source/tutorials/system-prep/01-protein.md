---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Build a protein

**You will learn:** how to take a PDB file with a canonical protein, prepare it, solvate it, and produce an AMBER topology and coordinates pair ready for MD.

**Prerequisites:**
- HTMD installed.

## The flow

For a protein containing only canonical residues, building is four steps after loading:

1. {py:func}`~moleculekit.tools.autosegment.autoSegment` - split the molecule into one segment per chain.
2. {py:func}`~moleculekit.tools.preparation.systemPrepare` - protonate at the chosen pH, fix protonation states, fill in missing sidechains and some missing backbone heavy atoms.
3. {py:func}`~htmd.builder.solvate.solvate` - wrap a water box around the prepared structure.
4. {py:func}`htmd.builder.amber.build` - run tLeap to produce a `prmtop` + `pdb` pair (and ionise to the requested salt concentration).

No NCAA detection or parameterisation step is needed when every residue is in tLeap's built-in ff14SB library.

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.solvate import solvate
```

## Step 1 - Load the structure

We use Trp-cage (PDB `1L2Y`): a 20-residue mini-protein with no ligands and no non-standard residues. {py:class}`~moleculekit.molecule.Molecule` accepts either a local file path (PDB, mmCIF, PSF, PRMTOP, ...) or a four-character PDB ID it downloads from RCSB on the fly:

```{code-cell} python
mol = Molecule("1L2Y")
```

## Step 2 - Segment the chains

```{code-cell} python
mol = autoSegment(mol, fields=("segid", "chain"))
```

`autoSegment` walks the chain graph and assigns each contiguous protein chain - and each non-protein residue - to its own segment. This is what stops {py:func}`~htmd.builder.amber.build` from extending a protein chain through a HETATM ligand later, and from auto-capping the wrong terminus when an NCAA is at the chain end.

For protein-only systems the practical effect is one segment per chain (here, just `P0`).

## Step 3 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4)
```

{py:func}`~moleculekit.tools.preparation.systemPrepare` runs PDB2PQR under the hood: predicts pKa values, picks protonation states at the requested pH, adds missing hydrogens, and returns the prepared molecule plus a list of non-canonical specs. When no `detect_specs` argument is supplied the function auto-detects and returns that list; when one is supplied it is returned unchanged. For a canonical-only protein the spec list is empty and you can ignore it.

## Step 4 - Solvate

```{code-cell} python
solvated = solvate(prepared, pad=10)
```

{py:func}`~htmd.builder.solvate.solvate` wraps a TIP3P water box around the prepared molecule. `pad=10` adds 10 Å of water in every direction beyond the molecule's bounding box - large enough to prevent the solute from seeing its own periodic image at typical box sizes.

## Step 5 - Build under AMBER

```{code-cell} python
amber.build(solvated, outdir="./build", ionize=True, saltconc=0.15)
```

That single call:

- Writes a tLeap input script consuming ff14SB for protein, TIP3P for water, and standard ion parameters.
- Detects disulfide bridges in `solvated` and feeds them to tLeap as `bond` directives - the default is auto-detect; override with `disulfide=[(sel1, sel2), ...]` if needed.
- Adds Na⁺ / Cl⁻ counter-ions to neutralise the system and reach 0.15 M NaCl.
- Returns a built {py:class}`~moleculekit.molecule.Molecule` and writes the AMBER files into `./build/`:

```text
build/
├── structure.prmtop      # AMBER topology
├── structure.pdb         # built coordinates as PDB
├── structure.crd         # built coordinates as CRD
├── tleap.in              # tLeap input we generated
├── leap.log              # tLeap's log
└── ff*_leaprc.*          # force-field paths sourced by tLeap
```

The `structure.prmtop` + `structure.pdb` pair is what {py:func}`acemd.protocols.setup_equilibration` (or any other MD driver) consumes.

## Parameters that matter

| Argument | Effect |
| --- | --- |
| `outdir` | Where the build files land. The directory is created. |
| `ff` | Force-field overrides, e.g. `ff=["leaprc.protein.ff14SB", "leaprc.water.tip3p"]`. Defaults to {py:func}`amber.defaultFf() <htmd.builder.amber.defaultFf>` (protein + water + DNA + RNA + lipid21 + GAFF2). |
| `ionize` | Add counter-ions and salt. `True` by default. |
| `saltconc` | NaCl concentration in mol/L when `ionize=True`. Defaults to `0` (counter-ions only). |
| `disulfide` | `None` for auto-detect, or a list of `(sel1, sel2)` atom-selection pairs. |
| `caps` | Per-segment caps as `{"P0": ("ace", "nme")}`. Auto by default; pass `("none", "none")` for a free terminus. |

## Gotchas

- Always run `autoSegment` *before* `systemPrepare`. The segmentation reflects covalent connectivity, and `systemPrepare`'s C-terminal handling relies on it.
- Hydrogens in the input PDB are dropped and re-added by `systemPrepare`. Don't try to preserve input Hs - the prep is what guarantees the protonation state is consistent with the chosen pH.
- Solvate *before* `amber.build`, not after. `amber.build` does not wrap a water box on its own.
- `saltconc` defaults to `0` (just neutralising counter-ions). Pass `saltconc=0.15` (or your target) for physiological ionic strength.
- For systems with bound metals like Zn or Ca, you'll usually want to keep them by including them in the input - tLeap has parameters for the common ones.

## See also

- {doc}`Build a protein with a ligand <02-protein-ligand>` - the next step up: non-canonical residues, SMILES templating, parameterisation.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
