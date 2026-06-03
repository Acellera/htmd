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

# Build a stapled peptide

**You will learn:** how to build a peptide whose two NCAAs are joined by a chemical crosslink (a hydrocarbon staple), exemplified by an NF-Y-derived stapled peptide from PDB `8QU4`.

**Prerequisites:**
- HTMD installed.
- You've worked through {doc}`Build a protein with a ligand <02-protein-ligand>`.

```{note}
The workflow below is **identical** to {doc}`Build a protein with a ligand <02-protein-ligand>` - the only change is the two extra SMILES strings you pass to {py:meth}`~moleculekit.molecule.Molecule.templateResidueFromSmiles` (one per NCAA). {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` reads the staple bond from the input structure's connectivity on its own, and {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` carries it through to the build without any extra wiring.
```

## What the staple is

PDB `8QU4` chain A is a 13-mer designed peptide containing two NCAAs:

- `NLE` (norleucine) at one end of the staple.
- `MK8` (an α-methyl-norleucine variant) at the other end.

The staple is a single CE-CE bond between the two NCAAs - the closure product of a ring-closing metathesis between two olefinic side chains. There is **no canonical anchor between them**: the crosslink joins two NCAAs directly.

For the build flow, this means {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` returns two {py:class}`~moleculekit.tools.nonstandard_residues.ChainResidueSpec` entries (one per NCAA) with `anchor_atom` set to the staple atom (`CE`), and {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` emits exactly one entry in `custombonds` for the staple closure.

```{note}
This tutorial **skips solvation and ionisation** so the build runs in seconds and the focus stays on the staple. For a production run, either solvate first with {py:func}`~htmd.builder.solvate.solvate` (and keep `ionize=True` on the build) or build with implicit solvent by passing `gbsa=True` to {py:func}`~htmd.builder.amber.build`.
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
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
mol = Molecule("8QU4")
mol.filter("chain A")
mol.segid[:] = "P"
```

```{code-cell} python
:tags: [remove-input]
show3d(
    mol,
    highlight_bonds=[("resname NLE and name CE", "resname MK8 and name CE")],
)
```

The orange tube marks the **NLE.CE - MK8.CE staple bond** - the lone covalent crosslink that turns the linear 13-mer into a stapled macrocycle.

## Step 2 - Detect

```{code-cell} python
specs = detectNonStandardResidues(mol)
for spec in specs:
    print(spec)
```

Each NLE / MK8 spec should carry `anchor_atom="CE"`, marking the side-chain carbon that participates in the staple.

## Step 3 - Template both NCAAs

```{code-cell} python
NLE_SMILES = "CCCC[C@@H](C(=O)O)N"
MK8_SMILES = "CCCC[C@](C)(C(=O)O)N"

mol.templateResidueFromSmiles("resname NLE", NLE_SMILES, addHs=True)
mol.templateResidueFromSmiles("resname MK8", MK8_SMILES, addHs=True)
```

These are the free amino-acid SMILES for the two residues, with the neutral amine + free acid form (NLE and MK8 are not ionisable beyond the backbone). For ionisable residues, encode the protonation state at your target pH explicitly. `templateResidueFromSmiles` strips the terminal `-OH` automatically when the residue sits inside a peptide chain.

## Step 4 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, detect_specs=specs)
```

PDB2PQR protonates the canonical residues of the 13-mer peptide as normal; the bond-capture mechanism is what preserves the inter-NCAA staple bond across the prep. Passing `detect_specs=specs` lets you reuse the list we already computed (avoiding a duplicate detect call), edit it before prep if needed, and thread the same list into `parameterizeFromSpecs`. The spec list is returned unchanged; rebinding it keeps the data flow visually obvious.

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

You'll see one `NLE.prepi` and one `MK8.prepi` in `topo_paths` (plus matching `.frcmod` files in `frcmod_paths`) and exactly one entry in `custombonds` - the `NLE.CE - MK8.CE` staple.

## Step 6 - Build

```{code-cell} python
amber.build(
    prepared,
    outdir="./build",
    custombonds=out.custombonds,
    topo=out.topo_paths,
    param=out.frcmod_paths,
    caps={"P": ("none", "none")},
    ionize=False,
)
```

Two non-default arguments here:

- `caps={"P": ("none", "none")}` disables automatic ACE / NME capping on this segment - the staple peptide has its own designed termini and adding caps on top would create overvalent atoms.
- `ionize=False` skips ion placement because we haven't solvated.

After the build, verify the staple is in the topology:

```{code-cell} python
import numpy as np
built = Molecule("./build/structure.prmtop")
built.read("./build/structure.pdb")

nle_ce = np.where((built.resname == "NLE") & (built.name == "CE"))[0][0]
mk8_ce = np.where((built.resname == "MK8") & (built.name == "CE"))[0][0]
print("staple bond present:", built.hasBond(nle_ce, mk8_ce)[0])
```

## Gotchas

- Disable capping (`caps={"<segid>": ("none", "none")}`) for designed peptides whose termini are explicit in the input. Default capping assumes a free protein N/C terminus and adds ACE / NME caps - which clash with hand-crafted termini.
- Other staple chemistries (lactam, thioether, click triazole) use the same flow - detect sees the inter-NCAA bond from the input file's connectivity, regardless of which atoms anchor the staple.

## See also

- {doc}`Build a cyclic peptide <05-cyclic-peptide>` - head-to-tail cyclisation as a special case of inter-NCAA bonding.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
