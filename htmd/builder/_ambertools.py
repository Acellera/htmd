"""AmberTools dispatch helpers shared by the nonstandard-residue builder.

Wraps the ``antechamber`` / ``parmchk2`` / ``prepgen`` invocations the
parameterization pipeline relies on, with a single switch to dispatch via
``antechamber_pyodide`` instead of subprocess when running under Pyodide.

Forcefield (atom typing) is restricted to ``GAFF`` / ``GAFF2`` since those
are the two parameter sets ``parmchk2`` understands. Charge methods are
restricted to the ones antechamber can compute directly from the input
structure with no extra QM input: ``AM1-BCC``, ``Gasteiger``, ``ABCG2``,
or ``None`` to skip charge assignment.

Under Pyodide only the ``Gasteiger`` charge method is available (the
WASM build of antechamber bundled in ``antechamber_pyodide`` ships
without the SQM backend that AM1-BCC / ABCG2 require).
"""

import logging
import os
from copy import deepcopy
from itertools import product

import numpy as np

from moleculekit.molecule import Molecule

logger = logging.getLogger(__name__)


_CHARGE_METHOD_MAP = {
    "am1-bcc": "bcc",
    "gasteiger": "gas",
    "abcg2": "abcg2",
    None: None,
}

_FORCEFIELD_MAP = {
    "gaff": ("gaff", "1"),
    "gaff2": ("gaff2", "2"),
}

_PYODIDE_CHARGE_METHODS = {"gasteiger", None}


def _write_pyodide_output(path, data):
    """Write bytes or str data returned by antechamber_pyodide to a file."""
    with open(path, "wb") as f:
        f.write(data if isinstance(data, bytes) else data.encode())


def _run_ambertools(program, cmd, cwd, use_pyodide=False, input_files=None):
    """Run an AmberTools command via native subprocess or antechamber_pyodide.

    Parameters
    ----------
    program : str
        The program name (e.g. "antechamber", "parmchk2", "prepgen").
    cmd : list of str
        Full command-line arguments (program name first).
    cwd : str
        Working directory. Native runs with cwd; pyodide ignores it.
    use_pyodide : bool
        If True, dispatch via antechamber_pyodide.run.
    input_files : dict or None
        For pyodide: mapping of filename -> bytes content.

    Returns
    -------
    dict or None
        Pyodide returns a dict of output filename -> bytes.
        Native returns None (outputs are written to disk).
    """
    if use_pyodide:
        from antechamber_pyodide import run as _ac_run

        return _ac_run(program, cmd, input_files=input_files or {})

    import subprocess

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(
            f"{program} failed (exit {result.returncode}):\n{result.stderr}"
        )
    return None


def _assign_rdkit_gasteiger_charges(mol):
    """Compute Gasteiger (PEOE) partial charges with RDKit and write them
    onto ``mol.charge`` in place.

    antechamber's own ``-c gas`` ignores the net charge and produces a
    charge set summing to zero even for an ion. RDKit's Gasteiger seeds
    from the per-atom formal charges, so the charges sum to the actual
    net charge - it handles charged species correctly. RDKit is already a
    moleculekit dependency and is available under Pyodide.
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


def _fftype_antechamber(
    mol,
    tmpdir,
    forcefield="GAFF2",
    netcharge=0,
    charge_method=None,
    use_pyodide=False,
    output_mol2="typed.mol2",
    output_frcmod="mol.frcmod",
):
    """Assign GAFF / GAFF2 atom types and generate frcmod using antechamber/parmchk2.

    Parameters
    ----------
    mol : Molecule
        Input molecule.
    tmpdir : str
        Working directory for intermediate files.
    forcefield : str
        Atom-typing target passed to antechamber's ``-at`` flag.
        One of ``"GAFF"`` or ``"GAFF2"`` (case-insensitive).
    netcharge : int
        Net formal charge of the molecule.
    charge_method : str or None
        Charge model. ``"Gasteiger"`` is computed with RDKit (it honours
        the net charge, unlike antechamber's own Gasteiger); ``"AM1-BCC"``
        and ``"ABCG2"`` are passed to antechamber's ``-c`` flag. ``None``
        skips charge assignment. Case-insensitive. Under Pyodide only
        ``"Gasteiger"`` (or ``None``) is supported.
    use_pyodide : bool
        Dispatch via antechamber_pyodide instead of subprocess.
    output_mol2 : str
        Basename for the typed mol2 file written into *tmpdir*.
    output_frcmod : str
        Basename for the frcmod file written into *tmpdir*.

    Returns
    -------
    typed_mol : Molecule
        Molecule with the requested atom types assigned.
    typed_path : str
        Absolute path to the typed mol2 file written into ``tmpdir``.
    frcmod_path : str
        Absolute path to the generated frcmod file.
    """
    ff_key = forcefield.lower() if forcefield else None
    if ff_key not in _FORCEFIELD_MAP:
        raise ValueError(
            f"Unsupported forcefield {forcefield!r}. "
            f"Supported: {sorted(_FORCEFIELD_MAP)}"
        )
    at_flag, parmchk_s = _FORCEFIELD_MAP[ff_key]

    charge_key = charge_method.lower() if charge_method else None
    if charge_key not in _CHARGE_METHOD_MAP:
        raise ValueError(
            f"Unsupported charge_method {charge_method!r}. "
            f"Supported: {sorted(k for k in _CHARGE_METHOD_MAP if k)}"
        )
    if use_pyodide and charge_key not in _PYODIDE_CHARGE_METHODS:
        raise ValueError(
            f"charge_method {charge_method!r} is not supported under Pyodide. "
            f"Use 'Gasteiger' or None."
        )
    ac_charge = _CHARGE_METHOD_MAP[charge_key]

    if charge_key == "gasteiger":
        # antechamber's '-c gas' ignores the net charge (its charges sum
        # to zero even for an ion). Compute Gasteiger with RDKit instead,
        # which seeds from formal charges, and let antechamber only assign
        # atom types, carrying these charges through unchanged.
        mol = mol.copy()
        _assign_rdkit_gasteiger_charges(mol)
        ac_charge = None

    input_mol2 = "input.mol2"
    mol2path = os.path.join(tmpdir, input_mol2)
    mol.write(mol2path)

    cmd = [
        "antechamber",
        "-i", input_mol2,
        "-fi", "mol2",
        "-o", output_mol2,
        "-fo", "mol2",
        "-at", at_flag,
        "-nc", str(netcharge),
        "-dr", "n",
        "-j", "1",
    ]
    if ac_charge is not None:
        cmd += ["-c", ac_charge]

    input_files = None
    if use_pyodide:
        with open(mol2path, "rb") as f:
            input_files = {input_mol2: f.read()}

    result = _run_ambertools(
        "antechamber",
        cmd,
        cwd=tmpdir,
        use_pyodide=use_pyodide,
        input_files=input_files,
    )

    typed_path = os.path.join(tmpdir, output_mol2)
    frcmod_path = os.path.join(tmpdir, output_frcmod)

    if use_pyodide:
        typed_mol2_data = result[output_mol2]
        _write_pyodide_output(typed_path, typed_mol2_data)
        input_files = {output_mol2: typed_mol2_data}
    else:
        with open(typed_path, "rb") as f:
            input_files = {output_mol2: f.read()}

    cmd = [
        "parmchk2",
        "-i", output_mol2,
        "-f", "mol2",
        "-o", output_frcmod,
        "-s", parmchk_s,
        "-a", "Y",
    ]

    result = _run_ambertools(
        "parmchk2",
        cmd,
        cwd=tmpdir,
        use_pyodide=use_pyodide,
        input_files=input_files,
    )

    if use_pyodide:
        _write_pyodide_output(frcmod_path, result[output_frcmod])

    typed_mol = Molecule(typed_path)

    # Verify the assigned charges sum to the requested net charge. This
    # catches a charge method that silently ignores it - antechamber's
    # own Gasteiger does, summing to zero even for an ion. The tolerance
    # absorbs antechamber's per-atom charge rounding.
    if charge_key is not None:
        total = float(np.sum(typed_mol.charge))
        if abs(total - netcharge) > 0.02:
            raise RuntimeError(
                f"{charge_method} charges sum to {total:.4f}, but the net "
                f"charge of the molecule is {netcharge} - the charge method "
                f"did not honour the net charge."
            )

    return typed_mol, typed_path, frcmod_path


def _fix_prepi_atomname_capitalization(mol, prepi):
    """Fix antechamber's case-flips on atom names in a generated prepi.

    Antechamber occasionally rewrites two-letter element atom names with
    a different case (``CL`` becomes ``Cl``, etc.); the prepi has to use
    the original capitalization so tLeap matches it against the mol2 atom
    names. Rewrites the prepi in place, using ``mol.name`` as the
    reference set.
    """
    uqnames = {x.upper(): x for x in np.unique(mol.name)}

    with open(prepi, "r") as f:
        lines = f.readlines()

    section = None
    for i in range(len(lines)):
        if lines[i].strip() == "":
            section = None

        if section == "atoms":
            if len(lines[i].strip().split()) < 4:
                continue
            old_name = lines[i][6:9].strip()
            if old_name.upper() in uqnames and uqnames[old_name.upper()] != old_name:
                logger.info(
                    f"Fixed residue {mol.resname[0]} atom name {old_name} -> "
                    f"{uqnames[old_name.upper()]} to match the input structure."
                )
                lines[i] = (
                    f"{lines[i][:6]}{uqnames[old_name.upper()]:4s}{lines[i][10:]}"
                )
            continue
        if section in ("loop", "improper"):
            n_pieces = len(lines[i].strip().split())
            for j in range(n_pieces):
                piece = lines[i][j * 5 : (j + 1) * 5].strip()
                if piece.upper() in uqnames and uqnames[piece.upper()] != piece:
                    logger.info(
                        f"Fixed residue {mol.resname[0]} atom name {piece} -> "
                        f"{uqnames[piece.upper()]} to match the input structure."
                    )
                    lines[i] = (
                        f"{lines[i][: j * 5]}{uqnames[piece.upper()]: >5}"
                        f"{lines[i][(j + 1) * 5 :]}"
                    )
            continue

        stripped = lines[i].strip()
        if stripped.startswith("CORR"):
            section = "atoms"
        elif stripped.startswith("LOOP"):
            section = "loop"
        elif stripped.startswith("IMPROPER"):
            section = "improper"

    with open(prepi, "w") as f:
        f.writelines(lines)


def _duplicate_parameters(prm, original):
    """Duplicate bond/angle/dihedral/improper entries under all permutations
    of the renamed atom types in ``original`` (``{old_type: [new_type, ...]}``).

    Used after a chain-resident NCAA's backbone GAFF2 types are renamed to
    ff14SB: any frcmod entry that mentions an original GAFF2 type also gets
    written under the corresponding ff14SB-typed permutation so torsions
    spanning the backbone-sidechain boundary still resolve at build time.
    """
    def _gen_permutations(typ, original):
        possibles = []
        for tt in typ:
            if tt in original:
                possibles.append(list(original[tt]) + [tt])
            else:
                possibles.append([tt])
        return list(product(*possibles))

    for prm_typ in ("bond_types", "angle_types", "dihedral_types", "improper_types"):
        new_bt = {}
        for typ in getattr(prm, prm_typ):
            perms = _gen_permutations(typ, original)
            if len(perms) == 1:
                continue
            for perm in perms:
                new_bt[perm] = getattr(prm, prm_typ)[typ]

        for bt, val in new_bt.items():
            prm.__dict__[prm_typ][bt] = deepcopy(val)
