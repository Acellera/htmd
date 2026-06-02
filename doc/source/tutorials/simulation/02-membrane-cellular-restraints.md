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

# Membrane simulations with cellular restraints

**You will learn:** how to keep a ligand on the correct side of the membrane during a long production run, using {py:func}`acemd.protocols.get_cellular_restraints` to derive a set of flat-bottomed restraints that follow the periodic image.

**Prerequisites:**
- HTMD installed.
- ACEMD installed.
- You've worked through {doc}`Build a membrane-embedded protein <../system-prep/07-membrane>` and {doc}`Run an MD simulation with ACEMD <01-acemd-md>`.

## Why this matters

In a membrane-protein simulation the box wraps periodically in all three directions. The bilayer extends across the XY plane, with the protein crossing it once - but in the **+z / -z direction**, a ligand or peptide on the extracellular side can diffuse out of the top of the box and reappear at the bottom, suddenly sitting on the intracellular side. The same atom is one image away; the simulation doesn't know the difference. For binding / unbinding studies this completely scrambles the kinetics.

{py:func}`acemd.protocols.get_cellular_restraints` derives a set of **flat-bottomed positional restraints** that confine a chosen molecule to its native compartment, optionally allowing it to *enter* but not *cross* the bilayer, or to cross only one way. The restraints are anchored to the lipid centre of mass so they track the membrane through translation but assume the box height stays fixed - hence "only use these in NVT" (the docstring spells this out).

## The flow

1. Start from a pre-built membrane-embedded system (see {doc}`tutorial 07 <../system-prep/07-membrane>` for the full build).
2. Compute the restraint list with {py:func}`~acemd.protocols.get_cellular_restraints`.
3. {py:func}`~acemd.protocols.setup_equilibration` with **conventional positional restraints** on the ligand (because equilibration runs as NPT and the box height changes).
4. Run the equilibration via the `acemd` CLI - this produces the `output.coor` (final equilibrated coordinates) and `output.xsc` (final box) that production setup reads.
5. {py:func}`~acemd.protocols.setup_production` with the **cellular restraints** for the NVT production run.
6. Run production via the `acemd` CLI.

```{note}
Steps 1-3 execute in the docs build. Step 5 (`setup_production`) is shown but not executed because it needs the equilibration's `output.coor` + `output.xsc`, and steps 4 / 6 (the `acemd` CLI) run on a GPU - copy the commands into a script when you're ready.
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from acemd.protocols import (
    setup_equilibration,
    setup_production,
    get_cellular_restraints,
)
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load a built apelin receptor + membrane system

The {doc}`membrane tutorial <../system-prep/07-membrane>` walks through the OPM-aligned 5VBL build end to end. We reuse its output: `structure.prmtop` + `structure.pdb` for the apelin receptor (a GPCR) embedded in a POPC + cholesterol bilayer with the agonist peptide bound on the extracellular side. The pre-built pair lives under `tutorials/simulation/_systems/`; the cell below copies it into a local `./build/` and loads it into a Molecule. The box dimensions come from the PDB's `CRYST1` record - `get_cellular_restraints` needs them.

```{code-cell} python
:tags: [remove-cell]
import os, shutil
from pathlib import Path
systems = Path(os.environ["HTMD_TUTORIAL_SYSTEMS"])
Path("./build").mkdir(exist_ok=True)
shutil.copy(systems / "apelin-membrane.prmtop", "./build/structure.prmtop")
shutil.copy(systems / "apelin-membrane.pdb", "./build/structure.pdb")
```

```{code-cell} python
built = Molecule("./build/structure.prmtop")
built.read("./build/structure.pdb")
```

```{code-cell} python
:tags: [remove-input]
show3d(built, focus='segid "P0"', ball_and_stick='segid "P0"')
```

## Step 2 - Compute the cellular restraints

The apelin peptide sits in segment `P0` (chain A of the original PDB). On the OPM-aligned structure, `+z` points toward the extracellular side of the membrane, so the peptide lives in the **extracellular** compartment. We want a restraint that keeps it there throughout production.

`membrane_rel_z` is the fraction of the box height at which the bilayer is centred. The default `0.5` only holds when the bilayer sits exactly mid-box; in practice the build's z-padding above and below the membrane is rarely symmetric, so compute it from the actual lipid coordinates:

```{code-cell} python
lipid_z = built.coords[built.atomselect("lipid"), 2, 0]
membrane_rel_z = float(lipid_z.mean() / built.box[2, 0])
print(f"membrane centred at z = {lipid_z.mean():.1f} / {built.box[2, 0]:.1f} "
      f"-> membrane_rel_z = {membrane_rel_z:.3f}")
```

```{code-cell} python
extracellular_sel = 'segid "P0"'   # apelin peptide
intracellular_sel = None           # no intracellular ligand in this system

restraints = get_cellular_restraints(
    built,
    extracellular_sel=extracellular_sel,
    intracellular_sel=intracellular_sel,
    membrane_rel_z=membrane_rel_z,
    extracellular_crossing=None,   # don't even enter the membrane
)
for r in restraints:
    print(r)
```

Each entry is a dict in ACEMD's `extforces` format - the same shape `setup_equilibration` and `setup_production` accept. The restraint is a one-dimensional flat-bottomed potential along z, centred on the lipid centre of mass with an offset that places the "free zone" entirely in the extracellular compartment. Inside the free zone the ligand feels nothing; cross the boundary and a harmonic force kicks in to push it back.

The `extracellular_crossing` knob controls what's allowed:

| Value | Behaviour |
| --- | --- |
| `"cross"` | Ligand may travel through the membrane and reach the intracellular side (still no PBC wrapping). |
| `"enter"` | Ligand may enter the membrane but not cross it - useful for partition-coefficient studies. |
| `None` | Ligand stays in the extracellular compartment - the default for "keep it on one side". |

`intracellular_crossing` is the mirror knob for any intracellular ligand. Pass `None` for either selection when you don't need a restraint on that side.

## Step 3 - Equilibration (NPT) with positional restraints

Equilibration runs NPT - the box height changes as the lipid bilayer relaxes - so **cellular restraints aren't safe here**. Use simple positional restraints to keep the peptide near its starting coordinates while the rest of the system settles.

```{code-cell} python
equil_restraints = [
    {
        "type": "positionalrestraint",
        "sel": f"{extracellular_sel} and noh",
        "setpoints": ["5@0"],
    }
]

setup_equilibration(
    "./build", "./equil",
    run="4ns",
    extforces=equil_restraints,
)
print(open("./equil/input.yaml").read())
```

The peptide's heavy atoms (`... and noh`) are held with a force constant of 5 kcal/mol/Å² for the whole 4 ns - long enough for the bilayer and water box to settle without the peptide drifting. Hydrogens are left unrestrained so the peptide can still relax internally. The minimisation + protein-backbone defaults from {py:func}`~acemd.protocols.setup_equilibration` still fire on top of this.

## Step 4 - Run the equilibration

Hand the equilibration directory to the `acemd` CLI from your shell. This is the slow step - it needs a GPU and runs for hours - so don't try to execute it inside a notebook.

```bash
acemd --input ./equil
```

When the command exits, `./equil/` contains:

- `output.coor` - the final equilibrated coordinates.
- `output.xsc` - the final equilibrated box.
- `output.vel` - the final velocities.
- `output.xtc` + `log.txt` + restart files - the trajectory and run state.

`setup_production` reads `output.coor` (as the starting coordinates for production) and `output.xsc` (as the production box) from this directory, so you must run equilibration *before* setting up production. The topology (`structure.prmtop`) and force-field parameter files are reused from the equilibration directory as-is.

## Step 5 - Production setup (NVT) with cellular restraints

Once equilibration has fixed the box dimensions, switch to NVT for production and replace the positional restraints with the cellular ones from step 2.

```python
setup_production(
    "./equil", "./prod",
    run="100ns",
    temperature=300,
    barostat=False,         # NVT - lipid restraints assume constant box
    extforces=restraints,   # the cellular restraints we computed above
)
```

`barostat=False` is the NVT switch. The `extforces` list flows straight from `get_cellular_restraints` into the production `input.yaml`.

## Step 6 - Run production

```bash
acemd --input ./prod
```

## Gotchas

- `get_cellular_restraints` reads the box from `mol.box`, the lipid positions from `mol.coords`, and the membrane location from `membrane_rel_z`. Always compute `membrane_rel_z` from `mol.coords[mol.atomselect("lipid"), 2].mean() / mol.box[2]` rather than hardcoding `0.5` - the asymmetric water padding above and below the bilayer shifts the relative centre away from mid-box in practice.
- **NVT only.** The restraints are anchored to the lipid centre and to a fraction of the box height; an NPT run where the box rescales would silently misalign the free-zone with the lipids. Equilibrate the box dimensions under NPT with conventional positional restraints first, then switch to NVT for production.
- `lipidsel` defaults to the standard AMBER lipid resnames. If you're using a custom lipid (or built with CHARMM-GUI conventions), pass `lipidsel="<your selection>"` explicitly.

## See also

- {doc}`Run an MD simulation with ACEMD <01-acemd-md>` - the canonical equilibration → production flow this tutorial extends.
- {doc}`Build a membrane-embedded protein <../system-prep/07-membrane>` - the long-form membrane build whose output we load here.
