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

# Build a protein with a ligand

**You will learn:** how to detect non-standard residues in a structure, template them from SMILES, parameterise them with antechamber, and feed the result to {py:func}`~htmd.builder.amber.build`.

**Prerequisites:**
- HTMD installed.

## When you need this flow

If your structure contains any residue your force field doesn't know - a small-molecule ligand, a non-canonical amino acid, a covalently-bound drug, a phosphorylated residue - the canonical {doc}`protein build <01-protein>` won't work as-is. The builder will raise a *"could not find residue X"* error.

The full flow adds two steps between `systemPrepare` and `amber.build`:

1. {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` - inspect the molecule, return one spec per non-canonical residue.
2. {py:meth}`~moleculekit.molecule.Molecule.templateResidueFromSmiles` - fix bond orders and hydrogens on each non-canonical residue from its reference SMILES.
3. {py:func}`~moleculekit.tools.preparation.systemPrepare` - protonate canonicals, preserve templated non-canonicals.
4. {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` - run antechamber on the cluster of each non-canonical, emit topology / parameter files and a `custombonds` list.
5. {py:func}`~htmd.builder.amber.build` - tLeap consumes those outputs.

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.nonstandard import parameterizeFromSpecs
from htmd.builder.solvate import solvate
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load and segment

We use trypsin in complex with benzamidine (PDB `3PTB`). `BEN` is the non-canonical residue:

```{code-cell} python
mol = Molecule("3PTB")
mol = autoSegment(mol, fields=("segid", "chain"))
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

## Step 2 - Detect non-standard residues

```{code-cell} python
specs = detectNonStandardResidues(mol)
for spec in specs:
    print(spec)
```

You should see one spec for `BEN` - a {py:class}`~moleculekit.tools.nonstandard_residues.LigandSpec`. The spec records what kind of non-canonical we're looking at (free ligand, covalent ligand, chain-resident NCAA, scaffold, ...) and any anchor information needed to handle crosslinks.

## Step 3 - Template the non-canonical residues from SMILES

`templateResidueFromSmiles` matches the SMILES against the atoms in the selection, sets correct bond orders and formal charges, and adds the missing hydrogens:

```{code-cell} python
BEN_SMILES = "[NH2+]=C(N)c1ccccc1"
mol.templateResidueFromSmiles('resname "BEN"', BEN_SMILES, addHs=True)
```

The SMILES carries the **protonated** benzamidinium form (one of the amidine nitrogens is `[NH2+]`), which is the physiologically relevant state at pH 7.4. The RCSB chemical component for `BEN` is stored as the neutral amidine `N=C(N)c1ccccc1` - that's a starting point, but the SMILES you pass to `templateResidueFromSmiles` must encode the **protonation state at your target pH** with explicit formal charges. Templating the neutral form locks the wrong charges into the parameterisation.

You don't have to hand-edit the SMILES for the mid-chain case: when a residue is peptide-bonded on one or both sides, the function automatically strips the terminal `-OH` / `-OXT` that's absent in the bonded copy and retries the match. Full heavy-atom coverage and explicit hydrogens still work best.

## Step 4 - Prepare with the spec list

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)
```

If you don't pass `detect_specs`, `systemPrepare` calls `detectNonStandardResidues` for you. Passing `detect_specs=specs` explicitly is useful when you want to **reuse** the list we already computed in step 2 (avoiding the duplicate detect call), **edit it** before prep (drop specs for residues you want to leave alone, tweak a `new_resname`, etc.), and **thread** the same list into `parameterizeFromSpecs` in step 5. `systemPrepare` returns the spec list unchanged as its second value; the rebind keeps the data flow visually obvious. To opt out of non-standard residue handling entirely, pass `detect_specs=[]`.

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

For each unique `(resname, terminal-position)` bucket, `parameterizeFromSpecs` runs antechamber to assign GAFF2 atom types and per-atom partial charges, then writes:

- `out.topo_paths` - one topology file per unique non-canonical bucket. Free ligands like `BEN` get a `.cif`; chain-resident NCAAs get a `.prepi`.
- `out.frcmod_paths` - the matching `BEN.frcmod` with bond / angle / dihedral parameters.
- `out.custombonds` - atom-selection pairs naming the inter-residue bonds tLeap should add (empty here because `BEN` is a free ligand with no covalent connection to the protein).
- `out.xml_paths` - OpenMM force-field XML(s) for the same residues, in case you later want to switch the build backend.

The default charge method is **AM1-BCC** - antechamber runs an AM1-BCC calculation per residue, producing standard GAFF-quality charges. We pass `charge_method="gasteiger"` here because it's much faster and good enough for a tutorial; for production builds drop the argument (or set it explicitly to `"am1-bcc"`) for higher-quality electrostatics.

## Step 6 - Solvate

```{code-cell} python
solvated = solvate(prepared, pad=10)
```

`amber.build` does not wrap a water box on its own, so we solvate the prepared structure first - same as in the {doc}`canonical protein tutorial <01-protein>`.

## Step 7 - Build

```{code-cell} python
built = amber.build(
    solvated,
    outdir="./build",
    custombonds=out.custombonds,
    topo=out.topo_paths,
    param=out.frcmod_paths,
    ionize=True,
    saltconc=0.15,
)
```

The three new arguments (`custombonds`, `topo`, `param`) wire `parameterizeFromSpecs`'s outputs into tLeap. Ionisation, disulfide detection, and capping happen the same way as in a canonical-only build.

```{code-cell} python
:tags: [remove-input]
show3d(built)
```

## Gotchas

- For ionisable ligands, template with the SMILES of the **protonation state at your pH** (e.g. `[NH2+]=C(N)c1ccccc1` for benzamidinium at 7.4). Templating the neutral form locks the wrong protonation state into the parameterisation.
- If a non-canonical residue is at the C-terminus and its template already carries `OXT`, pass `caps={"<segid>": ("ace", "none")}` to {py:func}`~htmd.builder.amber.build` so tLeap doesn't try to add an NME cap on top.
- `parameterizeFromSpecs` dedupes by `(resname, is_n_term, is_c_term)`. Two MLE residues in the same chain-position bucket share one `MLE.prepi`; an N-terminal MLE plus a mid-chain MLE produce two distinct prepis.

## See also

- {doc}`Build a cyclic peptide <03-cyclic-peptide>` - same flow with a head-to-tail cyclisation.
- {doc}`Build a stapled peptide <04-stapled-peptide>` - chemical crosslink between two NCAAs.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
