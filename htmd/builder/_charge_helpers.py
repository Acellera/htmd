# (c) 2015-2025 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
"""Charge-assignment helpers shared across the GAFF / SMIRNOFF backends.

These helpers compute partial charges on a :class:`moleculekit.molecule.Molecule`
and write them to ``mol.charge`` in place. Each helper lazy-imports its
own dependency stack so importing this module is cheap and never needs
the heavy ML / QM packages.

Lives in its own module to break what would otherwise be a circular
import between ``_ambertools.py`` (which now uses them in the GAFF path)
and ``openmm.py`` (which uses them in the SMIRNOFF path).
"""

import numpy as np


def _assign_rdkit_gasteiger_charges(mol):
    """Compute Gasteiger (PEOE) partial charges with RDKit and write them
    onto ``mol.charge`` in place.

    antechamber's own ``-c gas`` ignores the net charge and produces a
    charge set summing to zero even for an ion. RDKit's Gasteiger seeds
    from the per-atom formal charges, so the charges sum to the actual
    net charge - it handles charged species correctly. RDKit is already
    a moleculekit dependency and is available under Pyodide.
    """
    from moleculekit.rdkittools import molecule_to_rdkitmol
    from rdkit.Chem import AllChem

    rdmol = molecule_to_rdkitmol(mol, sanitize=True, _logger=False)
    AllChem.ComputeGasteigerCharges(rdmol)
    charges = np.array(
        [float(a.GetDoubleProp("_GasteigerCharge")) for a in rdmol.GetAtoms()],
        dtype=np.float32,
    )
    if not np.all(np.isfinite(charges)):
        raise RuntimeError(
            "RDKit Gasteiger produced non-finite charges; the molecule may "
            "have an atom in an environment its parameters do not cover."
        )
    mol.charge = charges


def _assign_nagl_charges(mol, model_name="openff-gnn-am1bcc-1.0.0.pt"):
    """Compute AM1-BCC-equivalent partial charges with NAGL (a GNN
    surrogate for AM1-BCC) and write them onto ``mol.charge`` in place.

    NAGL is a SMIRNOFF-ecosystem tool that reproduces AM1-BCC charges
    via a graph neural network, avoiding the SQM step. It is several
    orders of magnitude faster than antechamber's AM1-BCC on medium-to-
    large molecules and honours the formal charge. Works with either
    the GAFF or SMIRNOFF backend - the charges are just floats, the
    typing engine is independent.

    NAGL is an opt-in dependency (it pulls PyTorch + pytorch-lightning).
    Install with ``pip install acellera-htmd[nagl]`` (or ``uv sync
    --extra nagl``).
    """
    try:
        from openff.nagl import GNNModel
        from openff.nagl_models import get_model
    except ImportError as e:
        raise ImportError(
            "charge_method='nagl' requires the openff-nagl stack (NAGL + "
            "PyTorch). Install with 'pip install acellera-htmd[nagl]' or "
            "'uv sync --extra nagl'."
        ) from e

    off_mol = mol.toOpenFFMolecule(sanitize=True, assignStereo=True)
    model = GNNModel.load(get_model(model_name))
    charges = model.compute_property(off_mol)
    arr = np.asarray(getattr(charges, "magnitude", charges), dtype=np.float32)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(
            "NAGL produced non-finite charges; the molecule may have an "
            "atom in an environment the GNN was not trained on."
        )
    mol.charge = arr


def _assign_resp_charges(mol, multi_conf=False):
    """Compute RESP charges (Restrained ElectroStatic Potential fit to a
    Psi4-computed QM ESP) and write them onto ``mol.charge`` in place.

    Dispatches to ``parameterize.charge.psi4resp.get_resp_psi4_charges``,
    which runs HF/6-31G* (or def2-SV(P) when heavy elements are present)
    via Psi4, builds the Merz-Singh-Kollman grid via openff-recharge, and
    runs the iterative two-stage RESP fit. With ``multi_conf=True`` the
    fit is performed jointly over up to 10 RDKit-generated conformers.
    Works with either the GAFF or SMIRNOFF backend; the typing engine
    is independent.

    Requires the Acellera ``parameterize`` package (private; not on
    public PyPI) and Psi4 + openff-recharge.
    """
    try:
        from parameterize.charge.psi4resp import get_resp_psi4_charges
    except ImportError as e:
        raise ImportError(
            "charge_method='resp'/'resp-multiconf' requires the Acellera "
            "'parameterize' package (private, not on public PyPI) plus "
            "Psi4 and openff-recharge. Contact the Acellera team for "
            "access, or pick a different charge method."
        ) from e

    charges = get_resp_psi4_charges(mol, multi_conf=multi_conf)
    arr = np.asarray(charges, dtype=np.float32)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(
            "RESP fit produced non-finite charges; the Psi4 ESP "
            "calculation or the iterative solver likely diverged."
        )
    mol.charge = arr
