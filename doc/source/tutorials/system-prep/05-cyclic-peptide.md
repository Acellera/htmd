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

# Build a cyclic peptide

**You will learn:** how to build a head-to-tail cyclic peptide whose residues are non-canonical, exemplified by cyclosporin (PDB `4TOT` chain E).

**Prerequisites:**
- HTMD installed.
- You've worked through {doc}`Build a protein with a ligand <02-protein-ligand>` - this tutorial builds on the same five-step flow.

```{note}
The workflow below is **identical** to {doc}`Build a protein with a ligand <02-protein-ligand>` - the only change is the SMILES dictionary you pass to `templateResidueFromSmiles`. {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` reads the input structure's connectivity to find the non-canonical residues and the ring-closing peptide bond on its own, and {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` handles the cluster parameterization without any extra wiring.
```

## What makes cyclic peptides interesting

Cyclosporin A is a head-to-tail cyclic 11-residue peptide. Almost every residue is N-methylated or otherwise modified, and there are no canonical anchors - every residue is a non-canonical amino acid (NCAA), and the first and last residues are covalently joined to close the ring.

For the build flow, the practical implication: {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` will return one {py:class}`~moleculekit.tools.nonstandard_residues.ChainResidueSpec` per NCAA, plus the inter-residue peptide bonds are added as `custombonds` so tLeap closes the ring correctly. **You don't have to wire the cyclisation by hand** - detect sees the existing peptide bond between the last and first residues and emits the corresponding custombond.

```{note}
This tutorial **skips solvation and ionisation** so the build runs in seconds and the focus stays on the cyclisation. For a production run, either solvate first with {py:func}`~htmd.builder.solvate.solvate` (and keep `ionize=True` on the build) or build with implicit solvent by passing `gbsa=True` to {py:func}`~htmd.builder.amber.build`.
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
mol = Molecule("4TOT")
mol.filter("chain E")          # one of the cyclosporin copies in the crystal
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

You'll see one {py:class}`~moleculekit.tools.nonstandard_residues.ChainResidueSpec` per NCAA. Because every residue in chain E is peptide-bonded on both sides (no free terminus), they all sit in the same chain-position bucket: mid-chain.

## Step 3 - Template every NCAA from SMILES

```{code-cell} python
smiles = {
    "33X": "CC(C=O)NC",
    "34E": "CN[C@@H]([C@H](C)CN1CCN(CCOC)CC1)C=O",
    "ABA": "CC[C@H](C=O)N",
    "BMT": "C/C=C/C[C@@H](C)[C@H]([C@@H](C=O)NC)O",
    "DAL": "C[C@H](C=O)N",
    "MLE": "CC(C)C[C@@H](C=O)NC",
    "MVA": "CC(C)[C@@H](C=O)NC",
}
for resname, smi in smiles.items():
    mol.templateResidueFromSmiles(f'resname "{resname}"', smi, addHs=True)
```

Templating is per *unique* resname, not per occurrence - cyclosporin's three `MLE` residues all share one SMILES. RCSB chemical-component SMILES are a starting point: encode the **protonation state at your target pH** (none of these cyclosporin residues are ionisable, so the neutral form is correct). `templateResidueFromSmiles` strips the terminal `-OH` / `-OXT` automatically when the residue is peptide-bonded on one or both sides, so the same SMILES works in both contexts.

## Step 4 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)
```

For a cyclic peptide with no canonical residues, `systemPrepare`'s PDB2PQR pass has nothing to do on the protonation side - the bond-capture / restore mechanism is what preserves the inter-NCAA peptide bonds (including the ring-closing one) through the prep. Passing `detect_specs=specs` rather than relying on the default auto-detect lets you reuse the list we already computed, edit it before prep (drop entries you don't want spec-handled, tweak a `new_resname`, ...), and thread the same list into `parameterizeFromSpecs`. The spec list is returned unchanged; rebinding it back into `specs` keeps the data flow visually obvious.

## Step 5 - Parameterize

```{code-cell} python
out = parameterizeFromSpecs(
    specs,
    prepared,
    outdir="./params",
    charge_method="gasteiger",
)
print(out)
```

`parameterizeFromSpecs` dedupes by `(resname, is_n_term, is_c_term)`. Three `MLE` residues at mid-chain produce **one** `MLE.prepi`. The `out.custombonds` list carries every NCAA-NCAA peptide bond - including the closing bond between residue N and residue 1, which is what makes this a cycle.

## Step 6 - Build

```{code-cell} python
amber.build(
    prepared,
    outdir="./build",
    custombonds=out.custombonds,
    topo=out.topo_paths,
    param=out.frcmod_paths,
    ionize=False,
)
```

`ionize=False` skips ion placement because we haven't solvated; tLeap still adds `bond` directives for every entry in `custombonds`, including the head-to-tail closure, so the resulting `prmtop` carries a closed ring.

## Gotchas

- For peptides whose ring closes through a *side chain* (lactam bridges, isopeptide cycles, thioether cycles), detect emits a {py:class}`~moleculekit.tools.nonstandard_residues.ChainResidueSpec` with the appropriate anchor atom and the cycle closure is wired automatically - the same pattern that the {doc}`stapled-peptide tutorial <04-stapled-peptide>` shows.

## See also

- {doc}`Build a stapled peptide <04-stapled-peptide>` - chemical crosslink between two NCAAs at non-backbone atoms.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
