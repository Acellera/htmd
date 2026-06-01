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

# Build a bicyclic peptide

**You will learn:** how to build a peptide whose backbone is threaded onto a small-molecule **scaffold** through three side-chain crosslinks, exemplified by PDB `8QFZ` chain B - a 12-mer cyclic peptide whose three cysteines are thioether-linked to an LFI scaffold.

**Prerequisites:**
- HTMD installed.
- You've worked through {doc}`Build a stapled peptide <04-stapled-peptide>` - this tutorial extends the same scaffolded-NCAA pattern.

```{note}
The workflow below is **identical** to {doc}`Build a protein with a ligand <02-protein-ligand>` - the only change is the single SMILES you pass to `templateResidueFromSmiles` (for the LFI scaffold). {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` reads the three `SG-Cn` thioether bonds and the three CYS chain-positions from the input structure's connectivity on its own, and {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` carries them through to the build without any extra wiring.
```

## What the bicycle is

PDB `8QFZ` chain B is a 12-mer peptide containing three cysteine residues - **CYS11 (N-terminal), CYS17 (mid-chain), CYS22 (C-terminal)** - whose sulphur atoms (`SG`) are each thioether-bonded to one of the three anchor carbons (`C10`, `C11`, `C12`) of a small heterocyclic scaffold residue called **`LFI`** (1,3,5-tris(3-bromopropanoyl)hexahydro-1,3,5-triazine). Each `SG-Cn` closure replaces an LFI bromine, so the templating SMILES carries the unbound form with the three `Br` leaving groups.

The combined topology makes the peptide **bicyclic**: the backbone plus the three SG-Cn crosslinks define two non-trivial loops. In the build flow this means:

1. {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` emits one {py:class}`~moleculekit.tools.nonstandard_residues.ScaffoldSpec` (for `LFI`) plus three {py:class}`~moleculekit.tools.nonstandard_residues.ChainResidueSpec` entries (for the three `CYS`), each with a *different* `new_resname` because the three CYS sit in three distinct chain-position buckets (N-terminal, mid-chain, C-terminal).
2. {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` produces one `LFI.cif` for the scaffold plus one `.prepi` per CYS bucket, and three entries in `custombonds` - one per `SG-Cn` thioether closure.

```{note}
This tutorial **skips solvation and ionisation** so the build runs in seconds and the focus stays on the scaffolding. For a production run, either solvate first with {py:func}`~htmd.builder.solvate.solvate` (and keep `ionize=True` on the build) or build with implicit solvent by passing `gbsa=True` to {py:func}`~htmd.builder.amber.build`.
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.nonstandard import parameterizeFromSpecs
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load and segment

```{code-cell} python
mol = Molecule("8QFZ")
mol.filter("chain B")
mol = autoSegment(mol, fields=("segid", "chain"))
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

## Step 2 - Detect

```{code-cell} python
specs = detectNonStandardResidues(mol)
for spec in specs:
    print(spec)
```

You should see one `ScaffoldSpec(resname='LFI', ...)` and three `ChainResidueSpec(resname='CYS', ..., new_resname='XX*')` entries. The three CYS residues each get a **different** `new_resname` so the downstream parameterisation can write a distinct topology per residue. Two reasons for the rename:

1. **They're not canonical anymore.** Each one carries an extra `SG-Cn` thioether bond to the scaffold, so the standard ff14SB `CYS` template no longer matches - the sulphur valence is different, and ff14SB wouldn't know what to do at the new bond. Renaming to a non-canonical bucket lifts these residues out of the built-in `CYS` library and lets `parameterizeFromSpecs` write a per-residue topology that combines ff14SB backbone types with GAFF2 sulphur-side typing.
2. **Each sits in a different chain position.** `CYS11` is N-terminal, `CYS17` is mid-chain, `CYS22` is C-terminal - so the backbone hydrogen counts and formal charges differ between them. `parameterizeFromSpecs` buckets by `(resname, is_n_term, is_c_term)`, so giving the three CYS the same `new_resname` would collapse them into one bucket and lose the chain-position distinction.

Detect also drops the `HG` hydrogen on each CYS sulphur, since those positions are now occupied by the thioether bond.

## Step 3 - Template the LFI scaffold from SMILES

```{code-cell} python
LFI_SMILES = "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr"
mol.templateResidueFromSmiles("resname LFI", LFI_SMILES, addHs=True)
```

The SMILES describes the **unbound LFI** with all three bromopropanoyl arms intact. `templateResidueFromSmiles` matches the central triazinane skeleton against the residue in the structure, sets bond orders and formal charges on the matched atoms, then **strips the three unmatched terminal `Br` atoms** automatically - the same mechanism that strips a terminal `-OH` for mid-chain amino acids. What remains is the bound LFI form with three carbons (`C10`, `C11`, `C12`) primed for the thioether closures.

The three CYS residues are canonical and don't need a SMILES template.

## Step 4 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)
```

PDB2PQR protonates the canonical residues normally. The three CYS specs from step 2 carry the rename + `HG`-drop instructions, which `systemPrepare` applies at the end of the prep. The bond-capture mechanism preserves the three pre-existing `SG-Cn` bonds across the PDB2PQR pass.

## Step 5 - Parameterise

```{code-cell} python
out = parameterizeFromSpecs(
    specs,
    prepared,
    outdir="./params",
    charge_method="gasteiger",
)
print(out)
```

`parameterizeFromSpecs` runs antechamber once per cluster and emits:

- One `LFI.cif` for the scaffold (full GAFF2 typing).
- Three `.prepi` files, one per CYS bucket (each carries ff14SB types on the standard backbone + canonical sidechain atoms, with the sulphur side staying GAFF2 to match the scaffold).
- Three `custombonds` entries - one per `SG-Cn` thioether closure.

The `.frcmod` files include **cross-force-field bond terms** that bridge ff14SB backbone atoms (`N`, `CT`, `C`, `O`, ...) to GAFF2 sidechain types (lowercase atom types) so tLeap can resolve the canonical-to-scaffold junctions.

## Step 6 - Build

```{code-cell} python
amber.build(
    prepared,
    outdir="./build",
    custombonds=out.custombonds,
    topo=out.topo_paths,
    param=out.frcmod_paths,
    caps={"P0": ("none", "none")},
    ionize=False,
)
```

`caps={"P0": ("none", "none")}` disables auto-capping on the peptide segment: the N-terminal CYS is bonded to LFI (no free amine for an ACE cap to attach to) and the C-terminal CYS already carries `OXT` from its template.

## Step 7 - Verify

```{code-cell} python
import numpy as np
built = Molecule("./build/structure.prmtop")
built.read("./build/structure.pdb")

sg_idxs = np.where(built.name == "SG")[0]
lfi_c_idxs = np.where(
    (built.resname == "LFI") & np.isin(built.name, ["C10", "C11", "C12"])
)[0]
n_thioether = sum(
    1 for sg in sg_idxs for c in lfi_c_idxs if built.hasBond(int(sg), int(c))[0]
)
print("SG-Cn thioether bonds:", n_thioether)        # should be 3
```

## Gotchas

- Pass the **unbound** scaffold SMILES (with the leaving groups). `templateResidueFromSmiles` recognises the leaving-group atoms as terminal heavy atoms unmatched by the structure and strips them automatically; templating the *bound* form (no leaving groups) doesn't have enough atoms to disambiguate the MCS match.
- `parameterizeFromSpecs` buckets canonical anchors by `(resname, is_n_term, is_c_term)`, so an N-terminal CYS, a mid-chain CYS, and a C-terminal CYS *each* get their own `.prepi`. Don't be surprised by three CYS-derived topology files - the cluster pipeline needs them distinct because the backbone protonation differs at each position.
- The same flow handles other "constrained peptide" scaffolds (TATA, TBMT, ...) - swap `LFI_SMILES` for the scaffold's SMILES and let detect emit the appropriate `ScaffoldSpec` + canonical-anchor renames.

## See also

- {doc}`Build a stapled peptide <04-stapled-peptide>` - the simpler case of a direct NCAA-NCAA crosslink with no scaffold.
- {doc}`Build a cyclic peptide <03-cyclic-peptide>` - head-to-tail cyclisation as another macrocycle topology.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
