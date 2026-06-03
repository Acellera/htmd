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

# Run an MD simulation with ACEMD

**You will learn:** how to take a built, solvated, ionised system and drive it through a complete MD protocol - equilibration then production - using [ACEMD](https://software.acellera.com/acemd/).

**Prerequisites:**
- HTMD installed.
- ACEMD installed (`pip install acemd` matched to your CUDA build).
- The {doc}`Build a protein <../system-prep/01-protein>` tutorial walks through the build steps reused below.

## The flow

ACEMD's Python API exposes two convenience setups - one per protocol stage - plus the runner itself:

1. {py:func}`acemd.protocols.setup_equilibration` - writes the equilibration `input.yaml` and copies the structure / parameter / coordinate files into a fresh equilibration directory.
2. The `acemd` CLI - run from the shell against the equilibration directory, reads the input file and produces trajectory + log + restart files alongside.
3. {py:func}`acemd.protocols.setup_production` - reads the equilibrated state from step 2 and writes a production `input.yaml` into a new production directory.
4. The `acemd` CLI again - run against the production directory.

Equilibration is short (a few ns of NPT with restraints relaxing); production is the long unrestrained NVT/NPT run that generates the trajectory you'll analyse.

```{note}
This tutorial executes the `setup_equilibration` step so you can see the output it generates. The two `acemd` invocations and the production setup are shown but **not** executed - a real run takes hours to days on a GPU, well outside a tutorial. Copy the code blocks into a script and launch it on a workstation or queue when you're ready to simulate for real.
```

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from acemd.protocols import setup_equilibration, setup_production
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Get a built, solvated, ionised system

For this tutorial we reuse a Trp-cage (1L2Y) system that the {doc}`canonical protein build <../system-prep/01-protein>` already produced - load → segment → prepare → solvate → AMBER build with ionisation. ACEMD doesn't care which builder produced the files - an AMBER `structure.prmtop` + `structure.pdb` pair, a CHARMM PSF + coords, or an OpenMM XML force field all work as inputs. We solvate and ionise here because for explicit-solvent runs you'll usually want both - implicit-solvent runs (`gbsa=True` on the build) are also supported if you want to skip the water box.

```python
mol = Molecule("1L2Y")
mol = autoSegment(mol, fields=("segid", "chain"))
prepared, _ = systemPrepare(mol, pH=7.4)
solvated = solvate(prepared, pad=12)
amber.build(solvated, outdir="./build", ionize=True, saltconc=0.15)
```

To keep the docs build fast we skip the actual build here and pick up from the resulting `./build/` directory. Load it back into a Molecule to verify:

```{code-cell} python
:tags: [remove-cell]
import os, shutil
from pathlib import Path
systems = Path(os.environ["HTMD_TUTORIAL_SYSTEMS"])
Path("./build").mkdir(exist_ok=True)
shutil.copy(systems / "trpcage.prmtop", "./build/structure.prmtop")
shutil.copy(systems / "trpcage.pdb", "./build/structure.pdb")
```

```{code-cell} python
mol = Molecule("./build/structure.prmtop")
mol.read("./build/structure.pdb")
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

## Step 2 - Set up the equilibration

```{code-cell} python
setup_equilibration("./build", "./equil", run="4ns")
```

{py:func}`~acemd.protocols.setup_equilibration` reads the structure + parameters from `./build/`, copies them into `./equil/`, and writes an `input.yaml` that configures:

- A short steepest-descent minimisation (`minimize: 500` steps).
- An NPT ensemble: `thermostat: true` at `thermostattemperature: 300` K, `barostat: true` at the default 1 atm.
- `velocities: 300` K - Maxwell-Boltzmann velocities seeded at this temperature when **no checkpoint exists**; on a resumed run the velocity field is loaded from the checkpoint instead.
- `restart: true` - on re-invocation of `acemd --input ./equil`, the run resumes from `./equil/restart.chk` rather than starting over.
- Default positional restraints in `extforces`: four entries by default - protein Cα with `1@0`, protein heavy atoms (non-Cα) with `0.1@0`, nucleic backbone with `1@0`, nucleic non-backbone heavy with `0.1@0`. Each decays linearly to zero by step `500000`, gradually releasing the system into the equilibrated water box. Pass `defaultrestraints=False` to drop them, or add your own via `extforces=[...]`.
- The integration length you passed in (`run: 4ns`), and the periodic box size taken from the build (`boxsize: [...]`).

Inspect the produced files:

```{code-cell} python
sorted(os.listdir("./equil"))
```

`input.yaml` is the file ACEMD reads; the other files are the structure / parameters / starting coordinates it references. The YAML itself is human-readable - inspect it to see (and override) any of the integration parameters:

```{code-cell} python
print(open("./equil/input.yaml").read())
```

## Step 3 - Run the equilibration

Hand the equilibration directory to the `acemd` CLI from your shell. This is the slow step - it needs a GPU and runs for hours - so don't try to execute it inside a notebook.

```bash
acemd --input ./equil
```

ACEMD reads `./equil/input.yaml`, runs minimisation + the requested ns of NPT, and writes the trajectory (`output.xtc`), the final state (`output.coor`, `output.vel`, `output.xsc`), and the checkpoint (`restart.chk`) back into `./equil/`. The bundled `./equil/run.sh` wrapper redirects ACEMD's stdout/stderr into `log.txt`; if you call `acemd --input ./equil` directly the log goes to your terminal instead. When the command exits, the system is ready for production.

## Step 4 - Set up production

Production reuses the last frame of the equilibration as its starting state.

```python
setup_production("./equil", "./prod", run="100ns", temperature=300)
```

{py:func}`~acemd.protocols.setup_production` reads the final coordinates and box (`output.coor` + `output.xsc`) from `./equil/`, copies them into `./prod/`, and writes a production `input.yaml` configured for an unrestrained **NVT** run at the given temperature for the requested length. Velocities are **regenerated** at Maxwell-Boltzmann at `temperature` rather than copied from equilibration (set `barostat=True` to switch to NPT if you need pressure coupling in production).

## Step 5 - Run production

```bash
acemd --input ./prod
```

The trajectory lives in `./prod/output.xtc` and is what you'd feed into the {doc}`MSM analysis tutorials <../analysis/index>` (or any other downstream analysis).

## Running on a queue

For longer campaigns - many starting structures, replicas, or adaptive epochs - drive ACEMD through HTMD's job queues instead of running the `acemd` CLI by hand:

```python
from htmd.ui import SlurmQueue

queue = SlurmQueue()         # or LocalGPUQueue() / LsfQueue() / PBSQueue()
queue.submit("./prod")       # ./prod is the directory setup_production produced
queue.wait()                 # blocks until all submitted runs finish
```

The queues run ACEMD on each directory you submit and stream the resulting trajectories back. This is the pattern the [adaptive-sampling tutorials](../adaptive/index.md) build on.

## Parameters that matter

| Argument | Effect |
| --- | --- |
| `run="4ns"` (equilibration) | Equilibration length. |
| `run="100ns"` (production) | Production length. A single MD trajectory; pick to match the timescale you're studying or use adaptive sampling for orders-of-magnitude reach. |
| `temperature=300` | Production temperature in K. Match the equilibration. |
| `setup_equilibration(..., extforces=[...])` | List of external-force dicts (e.g. positional restraints on a ligand). See {py:func}`~acemd.protocols.setup_equilibration` for the schema. |
| `setup_equilibration(..., minimize=300)` | Override the minimisation step count. Useful when starting from a freshly-relaxed structure. |

## Gotchas

- The two `acemd` invocations dominate wall time; everything else (the build, the two `setup_*` calls) is seconds. Profile your run on a short equilibration before queuing up many production trajectories.
- The output trajectory format follows the extension of `trajectoryfile` (default `output.xtc`); pass e.g. `setup_production(..., trajectoryfile="output.dcd")` to switch. MSM analysis through {py:class}`~htmd.simlist.Sim` and {py:class}`~htmd.projections.metric.Metric` reads both XTC and DCD natively.

## See also

- {doc}`Build a protein <../system-prep/01-protein>` - the canonical build whose output this tutorial consumes.
- {doc}`MSM analysis tutorials <../analysis/index>` - what to do with the trajectory once production finishes.
- {doc}`Adaptive sampling explanation <../../explanation/adaptive-sampling>` - the alternative to single long trajectories for hard-to-sample processes.
