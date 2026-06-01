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

# Build a membrane-embedded protein

**You will learn:** how to fetch a bilayer-aligned protein from the [OPM database](https://opm.phar.umich.edu/), wrap a custom-composition lipid bilayer around it, and produce an AMBER `prmtop` for the full membrane system - all in one script.

**Prerequisites:**
- HTMD installed.
- You've worked through {doc}`Build a protein with a ligand <02-protein-ligand>` - this tutorial reuses the same NCAA detection / templating / parameterization pipeline.

## The system

PDB `5VBL` is the **apelin receptor** (a GPCR, chain B) bound to a **peptide agonist** (chain A) containing five chain-resident non-canonical amino acids - `200`, `ALC`, `HRG`, `NLE`, `OIC`. The structure also carries a co-crystallised free lipid (`OLC`). The receptor is membrane-embedded - OPM ships it with a fitted bilayer thickness of 33.4 Å - so we'll wrap a mixed **POPC / cholesterol** bilayer around it.

```{note}
This tutorial sets `minimize=0` and `equilibrate=0` on {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` so the tutorial executes in seconds. The resulting membrane carries only the **initial lipid placement** - lipid tails will have clashes and the bilayer is not relaxed. For a production system, set `minimize` to a few hundred steps and `equilibrate` to a few nanoseconds (and run on `platform="CUDA"` if you have a GPU).
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.opm import get_opm_pdb
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.nonstandard import parameterizeFromSpecs
from htmd.builder.solvate import solvate
from htmd.membranebuilder.build_membrane import buildMembrane
import numpy as np
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load the RCSB structure and align it to OPM

```{code-cell} python
mol = Molecule("5VBL")
mol.remove("water", _logger=False)

ref, thickness = get_opm_pdb("5VBL", validateElements=False)
mol.align("protein and name CA", refmol=ref, mode="structure")

print(f"bilayer thickness: {thickness} Å, atoms after removing water: {mol.numAtoms}")
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

The workflow is: load the structure you actually want to simulate (from RCSB or a local file), then rotate it onto OPM's membrane-aligned orientation. {py:func}`~moleculekit.opm.get_opm_pdb` fetches OPM's mirror of the same PDB - rotated and translated so that the bilayer center sits at `z=0` and the membrane normal points along `+z`. {py:meth}`~moleculekit.molecule.Molecule.align` with `mode="structure"` runs TM-align between the two, so the two structures don't need to have identical atom counts (OPM often trims residues that don't sit in the bilayer plane). The returned `thickness` is the OPM-fit double-leaflet thickness (33.4 Å for 5VBL). This is the alignment {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` expects from its `solute=` argument.

## Step 2 - Detect the non-standard residues

```{code-cell} python
mol = autoSegment(mol, fields=("segid", "chain"))

specs = detectNonStandardResidues(mol)
for spec in specs:
    print(spec)
```

Look at what detect found:

- Five `ChainResidueSpec` entries for the chain-resident NCAAs - `HRG`, `ALC`, `OIC`, `NLE`, and `200` (the C-terminal one). These are exactly the NCAAs we expect from the peptide's chemistry.
- **Two extra `ChainResidueSpec` entries for canonical residues** - `GLU` (resid 10) renamed to `XX1` with `anchor_atom='CD'`, and `LYS` (resid 13) renamed to `XX2` with `anchor_atom='NZ'`. Those are the canonical residues at the **two ends of an isopeptide bond**: the apelin agonist is a side-chain macrocycle closed by a `GLU.CD - LYS.NZ` γ-glutamyl / ε-lysyl crosslink. Detect saw the inter-residue bond on the side-chain atoms (not the backbone) and emitted the rename + anchor automatically, so the downstream parameterization gives each end of the isopeptide its own per-residue topology. The `custombonds` list emitted by {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` will close the ring at build time.
- One `LigandSpec` for the free `OLC` lipid co-crystallised with the protein.

So we need to template every chain-resident NCAA *and* the free lipid - but **not** `GLU` / `LYS` (they're canonical residues whose only modification is the inter-side-chain bond, which detect handles entirely from connectivity):

```{code-cell} python
smiles = {
    "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
    "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
    "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
    "NLE": "CCCC[C@@H](C=O)N",
    "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    "OLC": "CCCCCCCC(=O)OC[C@H](O)CO",
}
for resname, smi in smiles.items():
    mol.templateResidueFromSmiles(f'resname "{resname}"', smi, addHs=True)

prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)

out = parameterizeFromSpecs(
    specs, prepared, outdir="./params", charge_method="gasteiger"
)
print(out)
```

## Step 3 - Build the membrane around the protein

```{code-cell} python
memb = buildMembrane(
    [80, 80],
    ratioupper={"popc": 0.7, "chl1": 0.3},
    ratiolower={"popc": 0.7, "chl1": 0.3},
    minimize=0,
    equilibrate=0,
    platform="CPU",
    outdir="./memb",
    solute=prepared,
)
print(f"membrane mol: {memb.numAtoms} atoms")
```

What the arguments do:

| Argument | Effect |
| --- | --- |
| `[80, 80]` | XY footprint in Å. Make it ~10 Å wider than the protein in each direction so the bilayer fully surrounds it. |
| `ratioupper` / `ratiolower` | Per-leaflet lipid mole fractions. Different upper/lower compositions are allowed - useful for asymmetric bilayers. Run {py:func}`~htmd.membranebuilder.build_membrane.listLipids` to see what's available. |
| `solute=prepared` | The OPM-aligned protein. {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` carves the protein footprint out of the lipid placement so lipids and protein don't overlap from the start. |
| `minimize=0`, `equilibrate=0` | Skip the OpenMM relaxation step - only the initial lipid placement runs. Set both to positive values for production. |
| `platform="CPU"` | OpenMM platform for the lipid placement step. Switch to `"CUDA"` if you have a GPU. |

The returned `memb` is the **lipid bilayer + a water buffer above and below** - the protein is *not* included in the return; we'll merge them in the next step.

## Step 4 - Merge protein + membrane (lipids only)

```{code-cell} python
memb.remove("water", _logger=False)

system = Molecule()
system.append(prepared)
system.append(memb)
print(f"merged system: {system.numAtoms} atoms")
```

```{code-cell} python
:tags: [remove-input]
show3d(system)
```

Two things going on here:

- `buildMembrane` already carved the protein footprint out of the lipid placement (that's what `solute=prepared` does in the previous step), so the membrane and protein occupy disjoint XY space and no additional clash-removal is needed - a direct append is enough.
- We drop `buildMembrane`'s own water layer (a short z-buffer above and below the bilayer that the builder adds for the LJ packing's sake) and re-solvate the full system in the next step so the water box covers the protein's extramembrane regions too.

## Step 5 - Solvate with a membrane-aware box

```{code-cell} python
lipid_mask = system.atomselect("lipid")
xy_min = system.coords[lipid_mask, :2, 0].min(axis=0)
xy_max = system.coords[lipid_mask, :2, 0].max(axis=0)
z_min = system.coords[:, 2, 0].min() - 15
z_max = system.coords[:, 2, 0].max() + 15

system = solvate(
    system,
    minmax=[[xy_min[0], xy_min[1], z_min],
            [xy_max[0], xy_max[1], z_max]],
)
```

The XY box matches the lipid extent (so water doesn't poke out past the bilayer edge), and Z extends 15 Å above and below the tallest / lowest atom in the system - enough to fully solvate the protein's intracellular and extracellular domains. `solvate` will only place waters where there's space - no waters appear inside the bilayer.

```{code-cell} python
:tags: [remove-input]
show3d(system)
```

## Step 6 - Build under AMBER

```{code-cell} python
amber.build(
    system,
    outdir="./build",
    custombonds=out.custombonds,
    topo=out.topo_paths,
    param=out.frcmod_paths,
    caps={"P0": ("none", "none")},
    ionize=True,
    saltconc=0.15,
)
```

A couple of notes on the arguments:

- `caps={"P0": ("none", "none")}` - the peptide inhibitor's C-terminal NCAA (`200`) carries its own `OXT`, so we suppress the auto-ACE/NME on segment `P0` to avoid a clash. If your autoSegment assigns a different segid to the inhibitor, swap `"P0"` for the correct value (inspect `set(prepared.segid)`).
- `ionize=True, saltconc=0.15` - tLeap places counter-ions and salt in the waters added by `solvate` in the previous step.

The output `./build/structure.prmtop` + `structure.pdb` is now ready for an equilibration protocol (e.g. {py:func}`acemd.protocols.setup_equilibration`).

## Gotchas

- {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` with `solute=` assumes the protein already has the bilayer centered at `z=0` and the membrane normal pointing along `+z`. If you're not sure your structure is in that frame, align it onto an OPM reference with {py:func}`~moleculekit.opm.get_opm_pdb` (for PDBs present in OPM) or {py:func}`~moleculekit.opm.align_to_opm` (for new structures) first.
- {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` returns the bilayer **without** the protein; merge them explicitly. Same goes for any cofactors or co-crystallised ligands.
- `minimize=0, equilibrate=0` produces a raw, unrelaxed membrane straight out of the initial lipid placement. Re-run {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` with `minimize=300, equilibrate=1` (1 ns) before any production simulation; if your machine has a GPU also switch `platform="CUDA"`.
- The default lipid library ({py:func}`~htmd.membranebuilder.build_membrane.listLipids`) covers POPC, POPE, POPG, POPS, POPA, POPI, DOPC, DOPE, ..., plus cholesterol (`chl1`) and a few sphingolipids. Mixed-composition leaflets work; if you need a lipid that isn't listed, add a custom lipid PDB to the lipid directory (or open an issue upstream).

## See also

- {doc}`Build a protein with a ligand <02-protein-ligand>` - the NCAA workflow this tutorial reuses.
- {doc}`System-building overview <../../explanation/system-building>` - where the membrane builder fits in the wider stack.
