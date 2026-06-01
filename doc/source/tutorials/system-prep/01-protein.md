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

1. {py:func}`~moleculekit.tools.autosegment.autoSegment` - automatically split the structure into independent segments by walking the residue connectivity and starting a new segment at every real chain break (residue-number gap validated by backbone distance for proteins).
2. {py:func}`~moleculekit.tools.preparation.systemPrepare` - protonate at the chosen pH, fix protonation states, fill in missing sidechains and some missing backbone heavy atoms.
3. {py:func}`~htmd.builder.solvate.solvate` - wrap a water box around the prepared structure.
4. {py:func}`htmd.builder.amber.build` - run tLeap to produce a `prmtop` + `pdb` pair (and ionise to the requested salt concentration).

No NCAA detection or parameterization step is needed when every residue is in tLeap's built-in ff14SB library.

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.solvate import solvate
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load the structure

We use Trp-cage (PDB `1L2Y`): a 20-residue mini-protein with no ligands and no non-standard residues. {py:class}`~moleculekit.molecule.Molecule` accepts either a local file path (PDB, mmCIF, PSF, PRMTOP, ...) or a four-character PDB ID it downloads from RCSB on the fly:

```{code-cell} python
mol = Molecule("1L2Y")
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

## Step 2 - Segment the chains

```{code-cell} python
mol = autoSegment(mol, fields=("segid", "chain"))
```

{py:func}`~moleculekit.tools.autosegment.autoSegment` walks the structure residue-by-residue and starts a new segment whenever it sees a real chain break: a jump in residue numbers (e.g. resids 22 → 50 with no intermediate residues) that the spatial check confirms is a true gap in the backbone, *or* a transition from one chain identifier to another. Each contiguous run of bonded residues gets its own segid. This is what stops {py:func}`~htmd.builder.amber.build` from extending a protein chain through a HETATM ligand later, and from auto-capping the wrong terminus when a non-canonical residue sits at the chain end.

For Trp-cage there's just one continuous chain, so autoSegment produces a single segment named `P0`.

## Step 3 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4)
```

```{code-cell} python
:tags: [remove-input]
show3d(prepared, ball_and_stick="all")
```

{py:func}`~moleculekit.tools.preparation.systemPrepare` runs PDB2PQR under the hood: predicts pKa values, picks protonation states at the requested pH, adds missing hydrogens, and returns the prepared molecule plus a list of non-canonical specs. When no `detect_specs` argument is supplied the function auto-detects and returns that list; when one is supplied it is returned unchanged. For a canonical-only protein the spec list is empty and you can ignore it.

## Step 4 - Solvate

```{code-cell} python
solvated = solvate(prepared, pad=10)
```

```{code-cell} python
:tags: [remove-input]
show3d(solvated)
```

{py:func}`~htmd.builder.solvate.solvate` wraps a TIP3P water box around the prepared molecule. `pad=10` adds 10 Å of water in every direction beyond the molecule's bounding box. We keep the box small to keep the tutorial fast; for a production run you'd typically use a larger pad (15-20 Å) so that any local unfolding or large-scale motion can't reach across the periodic boundary and interact with the protein's own image.

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

## Force-field defaults and overrides

`amber.build` consumes three categories of force-field input, each with its own argument:

| Argument | Default | Format |
| --- | --- | --- |
| `ff` | {py:func}`amber.defaultFf() <htmd.builder.amber.defaultFf>` | tLeap `leaprc.*` files (the master force-field selectors). |
| `topo` | {py:func}`amber.defaultTopo() <htmd.builder.amber.defaultTopo>` (empty) | Residue topology templates: `.prepi`, `.prep`, `.in`, `.cif`, or `.mol2` (the last two are converted to prepi via `prepgen` automatically). |
| `param` | {py:func}`amber.defaultParam() <htmd.builder.amber.defaultParam>` (empty) | Parameter overlays: `.frcmod`. |

Inspect the active defaults:

```{code-cell} python
print(amber.defaultFf())
```

To **replace** the force-field set, pass your own list to `ff=`:

```python
amber.build(
    solvated,
    outdir="./build",
    ff=["leaprc.protein.ff19SB", "leaprc.water.opc"],   # ff19SB protein + OPC water
)
```

To **augment** the build with extra residue templates or parameters (e.g. an NCAA you parameterized with antechamber, or a custom cofactor), pass them via `topo=` and `param=`:

```python
amber.build(
    solvated,
    outdir="./build",
    topo=["./params/MY_RES.prepi"],
    param=["./params/MY_RES.frcmod"],
)
```

This is exactly what {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` produces for non-canonical residues in {doc}`tutorial 02 <02-protein-ligand>` and onwards - the function returns `topo_paths` and `frcmod_paths` ready to feed straight in.

To browse what's bundled with HTMD's tLeap install (leaprc files, prepi templates, frcmod overlays):

```python
amber.listFiles()
```

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
- You don't have to strip input hydrogens before `systemPrepare`. If `titration=True` (the default) PDB2PQR re-protonates all titratable residues at the chosen pH and effectively overrides whatever Hs the input carried; passing `titration=False` keeps your input protonation. Either way, the prep is what reconciles the hydrogens with the rest of the build.
- Solvate *before* `amber.build`, not after. `amber.build` does not wrap a water box on its own.
- `saltconc` defaults to `0` (just neutralising counter-ions). Pass `saltconc=0.15` (or your target) for physiological ionic strength.
- For systems with bound metals like Zn or Ca, you'll usually want to keep them by including them in the input - tLeap has parameters for the common ones.

## See also

- {doc}`Build a protein with a ligand <02-protein-ligand>` - the next step up: non-canonical residues, SMILES templating, parameterization.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
