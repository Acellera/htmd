# Installation

HTMD is published on PyPI as `acellera-htmd` and also available as a conda package. Pip is the default; use conda if you need the bundled external tools (CHARMM build helpers, etc.) that conda packages but pip does not.

## With pip (default)

```bash
python -m venv .venv
source .venv/bin/activate
pip install acellera-htmd
```

## With conda

We recommend Miniforge: <https://github.com/conda-forge/miniforge>. Once it's installed:

```bash
conda create -n htmd
conda activate htmd
conda install python=3.10 htmd -c acellera -c conda-forge
```

Some HTMD features (e.g. CHARMM building, [ACEMD](https://software.acellera.com/acemd/) integration) depend on tools that conda can install but pip cannot. If you hit a `FileNotFoundError` for an external binary, switch to the conda install.

## Licence registration

HTMD prints a short copyright reminder on import in interactive sessions. To silence it and let Acellera know who is using HTMD, you can register your install once:

```bash
htmd_register
```

Registration is optional - the package works fully without it.

## Verify

```python
from htmd.ui import *

mol = Molecule("3PTB")
print(mol.numAtoms)
```

Expected output: an integer count of atoms in PDB entry `3PTB`.

## For contributors: `uv`

If you're developing HTMD (running the test suite, building docs, working from a checkout), use [`uv`](https://docs.astral.sh/uv/):

```bash
uv sync --group dev      # install dev + test deps
uv sync --group docs     # add the doc-build deps
uv run pytest            # run the test suite
cd doc && uv run make html
```

See the `pyproject.toml` `[dependency-groups]` table for the full list of optional dependency groups.
