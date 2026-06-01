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

# Build a protein with a ligand using OpenFF

**You will learn:** how to parameterize the benzamidine ligand with the **OpenFF Sage** small-molecule force field (instead of the default GAFF2) and build the resulting trypsin + benzamidine system through {py:func}`htmd.builder.openmm.build`.

**Prerequisites:**
- HTMD installed with the OpenFF extra: `pip install acellera-htmd[openff,nagl]`.
- You've worked through {doc}`Build a protein with a ligand <02-protein-ligand>` - this tutorial reuses the same five prep steps and only changes the small-molecule force field.

## Why a different small-molecule force field

{doc}`Tutorial 02 <02-protein-ligand>` parameterizes `BEN` with **GAFF2** through antechamber. GAFF2 is a solid generic default and is what the AMBER ecosystem has long used for small molecules. The [OpenFF / SMIRNOFF family](https://openforcefield.org/force-fields/) ("Sage") is a more recent alternative:

- **Direct SMARTS-based parameter assignment** via the [SMIRNOFF format](https://openforcefield.github.io/standards/standards/smirnoff/): a single ~300-line force-field XML replaces GAFF2's many-line atom-typing rules, with chemistry-aware specificity (sulfonamides, phosphates, and bridgehead nitrogens are typed correctly from Sage 2.1 onwards).
- **vdW parameters retrained** against experimental mass density and enthalpy-of-mixing data from NIST.
- **Valence parameters** (bonds, angles, torsions) fit against quantum-chemical optimised geometries and torsion profiles, with torsions split per central-bond multiplicity from Sage 2.3.
- **Charges via NAGL**, an AM1-BCC graph-neural-network surrogate introduced with Sage 2.3 - much faster than running AM1-BCC per-residue, and what the Sage release was validated against.
- **Native OpenMM XML output**, no tLeap round-trip. Training inputs, scripts, and results are public on the [OpenFF GitHub org](https://github.com/openforcefield).

For any ligand where you want the latest, externally-validated force-field assignment, swap GAFF2 for an OpenFF release. The {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` interface is unchanged - you flip two arguments:

- `forcefield="openff_unconstrained-2.3.0.offxml"` (Sage 2.3, the current production line).
- `charge_method="nagl"` (the GNN surrogate for AM1-BCC that Sage was fit against).

The OpenMM builder ({py:func}`htmd.builder.openmm.build`) consumes the emitted force-field XML directly - no separate prepi/frcmod step.

```{warning}
Parameterizing **any** residue under OpenFF binds the whole system to the OpenMM builder. The `.frcmod` format that tLeap consumes can't represent the full SMIRNOFF feature set, so {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` emits only the OpenMM force-field XML on the SMIRNOFF branch and {py:func}`htmd.builder.amber.build` has nothing to consume. If you need an AMBER `prmtop` produced by tLeap, stay on the GAFF2 path in {doc}`tutorial 02 <02-protein-ligand>`. The OpenMM builder still writes a `prmtop` + `pdb` pair, so the *downstream* MD engine (ACEMD, OpenMM, GROMACS, ...) doesn't care.
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import openmm
from htmd.builder.nonstandard import parameterizeFromSpecs
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Steps 1-4 - Same as the GAFF2 tutorial

The prep pipeline is unchanged from {doc}`tutorial 02 <02-protein-ligand>` - see that tutorial for the per-step prose. Code-only recap:

```{code-cell} python
mol = Molecule("3PTB")
mol = autoSegment(mol, fields=("segid", "chain"))
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

```{code-cell} python
specs = detectNonStandardResidues(mol)

BEN_SMILES = "[NH2+]=C(N)c1ccccc1"
mol.templateResidueFromSmiles('resname "BEN"', BEN_SMILES, addHs=True)

prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)
```

## Step 5 - Parameterize with OpenFF Sage

```{code-cell} python
out = parameterizeFromSpecs(
    specs,
    prepared,
    outdir="./params",
    forcefield="openff_unconstrained-2.3.0.offxml",
    charge_method="nagl",
)
print(out)
```

`forcefield` strings that don't start with `gaff` are treated as SMIRNOFF offxml filenames and dispatched through OpenFF Interchange. The output `ClusterOutputs` will carry `.xml` files in `xml_paths` (per-cluster OpenMM force-field XML); `topo_paths` and `frcmod_paths` may be empty for the SMIRNOFF branch, since the XML is self-contained.

`charge_method="nagl"` runs the OpenFF NAGL GNN surrogate for AM1-BCC, which is the charge model Sage's vdW and torsion parameters were fit against. Mixing Sage with Gasteiger or any other model emits a mismatched-charges warning - it works, but you lose some of the consistency the Sage release was validated for.

## Step 6 - Build under OpenMM

```{code-cell} python
# 3PTB ships with crystallographic waters that systemPrepare keeps; strip
# them so openmm.build's internal solvation step runs (otherwise it sees
# the mol as already-solvated and skips it, leaving only the close-in
# crystal waters that the ioniser can't safely swap for ions).
prepared.remove("water", _logger=False)

built, system = openmm.build(
    prepared,
    outdir="./build",
    extra_xml=list(out.xml_paths),
    custombonds=out.custombonds,
    solvate=True,
    padding=15.0,
    ionize=True,
    saltconc=0.15,
)
```

Argument breakdown:

| Argument | Effect |
| --- | --- |
| `extra_xml` | Additional XML files describing the non-canonical residues - here, the per-cluster Sage XML emitted by `parameterizeFromSpecs`. |
| `custombonds` | Inter-residue bonds (same shape as in `amber.build`). |
| `solvate=True` | The OpenMM builder solvates internally; no separate {py:func}`~htmd.builder.solvate.solvate` call is needed. If you'd rather pre-solvate (with your own padding, water model, or `centersel`), call {py:func}`~htmd.builder.solvate.solvate` on `prepared` first, then pass `solvate=False` here and the builder will skip the internal step. |
| `padding` | Solvation padding in Å. |
| `ionize` + `saltconc` | Same meaning as in `amber.build`. |

`openmm.build` returns a `(built_mol, openmm_system)` tuple. The `built_mol` is a {py:class}`~moleculekit.molecule.Molecule` of the final system; `openmm_system` is the OpenMM `System` object you can hand straight to an OpenMM `Simulation`. The output directory still contains `structure.prmtop` + `structure.pdb` so the result is consumable by any downstream MD driver, including ACEMD via {py:func}`acemd.protocols.setup_equilibration`.

## Gotchas

- `frcmod` (tLeap's input format) can't express the full SMIRNOFF feature set, so the OpenFF path emits OpenMM XML only - `amber.build` cannot consume it. Picking OpenFF for one residue commits the whole build to {py:func}`htmd.builder.openmm.build`.
- Pair `forcefield="openff_*.offxml"` with `charge_method="nagl"`. Sage's parameters were fit alongside AM1-BCC charges, and NAGL is the recommended fast surrogate; other charge models (Gasteiger, antechamber's AM1-BCC) work but emit a mismatched-charges warning and lose some of the consistency the Sage release was validated for.
- `openmm.build` does its own solvation when `solvate=True`. If you want to pre-solvate (e.g. to use a non-default padding or `centersel`), call {py:func}`~htmd.builder.solvate.solvate` first and pass `solvate=False` to the builder. Don't pre-solvate and leave `solvate=True` - you'd end up with a doubly-solvated box.
- For pure-canonical protein systems (no ligand, no NCAA) there's nothing for OpenFF to do at the small-molecule level; you don't need Sage and the GAFF2 path in tutorial 02 is enough. This tutorial pays off when the ligand is the focus.

## See also

- {doc}`Build a protein with a ligand <02-protein-ligand>` - the same system parameterized under GAFF2 via tLeap.
- {doc}`System-building overview <../../explanation/system-building>` - the conceptual map of the whole stack.
