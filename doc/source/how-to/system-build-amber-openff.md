# How to mix AMBER protein with OpenFF small-molecule parameters

## Goal

Build a system where the protein is parameterized under AMBER (ff14SB / ff19SB) but the small-molecule ligand uses **OpenFF / SMIRNOFF** parameters instead of GAFF2. The build is driven through OpenMM (which natively supports the SMIRNOFF templating) and emits an output bundle in two interchangeable forms:

- **AMBER form**: `structure.prmtop` + the topology PDB / CIF. The CIF (`structure.cif`) is **always** written; the PDB (`structure.pdb`) is additionally written only when the system fits within PDB's 99999-atom CONECT serial limit. The prmtop runs in ACEMD / OpenMM / GROMACS directly.
- **ForceField-XML form**: `structure.cif` + `structure.pdb` (when it fits) + one `structure.smallmol_<i>.xml` per OpenFF/SMIRNOFF-templated ligand + any `extra_xml` files copied verbatim + `system.yaml` (a manifest listing the structure, the ordered parameter set, and the box). This is the form ACEMD reads natively when you pass `system.yaml` to it - no prmtop round-trip through ParmEd's atom-type bucketing, parameters stay byte-exact as OpenMM built them.

Pick the AMBER form for engine interoperability; pick the ForceField-XML form (via `system.yaml`) when you want OpenMM-exact reproducibility under ACEMD.

## Minimal example

```python
from htmd.builder import openmm as builder
from openff.toolkit import Molecule as OFFMolecule

ligand = OFFMolecule.from_smiles("[NH2+]=C(N)c1ccccc1")  # benzamidinium

builder.build(
    mol,                                    # the moleculekit Molecule with protein + ligand
    outdir="./build",
    ff=["amber14-all.xml", "amber14/tip3p.xml"],
    small_molecule_ff="openff-2.3.0",       # SMIRNOFF force field for the ligand
    molecules=[ligand],                     # OpenFF Molecule object(s) for the ligand(s)
    ionize=True, saltconc=0.15,
)
```

The protein is parameterized with the standard AMBER XML force field; the ligand is matched against `mol`'s atoms and parameterized on-the-fly via the OpenFF SMIRNOFF template generator. The output `./build/structure.prmtop` is fully compatible with `acemd --input` or OpenMM's `app.AmberPrmtopFile`.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `ff` | OpenMM XML force field files - the protein / lipid / nucleic / water stack. Default is {py:func}`~htmd.builder.openmm.defaultFf` = `["amber14/protein.ff14SB.xml", "amber14/lipid17.xml", "amber14/DNA.bsc1.xml", "amber14/RNA.OL3.xml", "amber14/tip3p.xml"]`. Replace, don't append - passing `ff=` overrides the full default. |
| `small_molecule_ff` | Name of the small-molecule force field, e.g. `"openff-2.3.0"` (SMIRNOFF) or `"gaff-2.2.20"` (GAFF). |
| `molecules` | List of `openff.toolkit.Molecule` objects (or paths to SDF files) describing every small molecule that needs `small_molecule_ff` templating. |
| `extra_xml` | Additional OpenMM XML files for non-standard residues that aren't templated. |
| `solvate=True`, `padding=10.0`, `water_model="tip3p"` | Wrap a water box around the system at build time. Set `solvate=False` if you've already pre-solvated via {py:func}`~htmd.builder.solvate.solvate`. |

## Common variations

### Multiple ligands

```python
lig1 = OFFMolecule.from_smiles("[NH2+]=C(N)c1ccccc1")  # benzamidine+
lig2 = OFFMolecule.from_smiles("c1ccc(O)cc1")          # phenol

builder.build(
    mol, outdir="./build",
    small_molecule_ff="openff-2.3.0",
    molecules=[lig1, lig2],
)
```

The template generator matches each `OFFMolecule` against the atoms in `mol` by graph isomorphism, so the SMILES has to be the **exact protonation state** of the residue in `mol` (formal charges and explicit hydrogens included).

### Using GAFF-2.2 with OpenMM instead of OpenFF

```python
builder.build(mol, outdir="./build",
              small_molecule_ff="gaff-2.2.20",
              molecules=[ligand])
```

Same flow but the ligand gets GAFF-2 atom types and charges via `openmmforcefields`. Useful when you want to keep small-molecule parameters consistent with an AMBER-only build elsewhere in the project.

### Loading ligands from SDF

```python
builder.build(mol, outdir="./build",
              small_molecule_ff="openff-2.3.0",
              molecules=["./ligand.sdf"])
```

The SDF must encode the same **chemical graph** (elements + connectivity + bond orders + formal charges) as the corresponding residue in `mol` - atom order doesn't have to match, the SMIRNOFF template generator uses RDKit substructure matching. Useful when the ligand has stereochemistry that's hard to express in a one-liner SMILES.

## Gotchas

- `small_molecule_ff` and `molecules` must be passed **together**. The template generator is only registered when both are non-empty (see `_setup_forcefield`); passing one without the other is a silent no-op and the build fails later with an unhelpful "no template found" error from OpenMM.
- `small_molecule_ff` + `molecules` is the **explicit** templating path - the OpenFF/GAFF parameters are computed once and locked into the prmtop. If you instead pass the ligand as part of `extra_xml` (manually generated XML), that's the *pre-parameterized* path and the small-molecule ff isn't used.
- OpenFF requires correct bond orders and formal charges on the input molecule. Run {py:meth}`~moleculekit.molecule.Molecule.templateResidueFromSmiles` first to fix bond orders if you loaded the ligand from a plain PDB.
- Sage **2.3** (`openff-2.3.0`) assigns partial charges through a pretrained graph-neural-network surrogate (NAGL) — no conformer / AM1-BCC step, so the first build on a new ligand is fast (sub-second per molecule). Older versions (`openff-2.0` / `2.1` / `2.2`) and GAFF (`gaff-2.2.20`) still run real AM1-BCC and take a few seconds per unique ligand. Caching is in-memory per process only - separate Python runs always re-evaluate charges.
- AMBER+OpenFF mixing is the only mixed-stack we recommend. **Don't** mix CHARMM and AMBER protein force fields - their non-bonded parameter conventions are incompatible.

## See also

- {doc}`Build a protein with an OpenFF-parameterized ligand <../tutorials/system-prep/03-openff-ligand>` - the full executable tutorial this how-to compresses.
- {doc}`How to use a custom force field with amber.build <system-build-custom-forcefield>` - the GAFF / `amber.build` equivalent.
- [OpenFF toolkit docs](https://docs.openforcefield.org/) and [openmmforcefields](https://github.com/openmm/openmmforcefields).
