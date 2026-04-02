# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.tools.graphalignment import makeMolGraph, compareGraphs
from moleculekit.molecule import Molecule
from moleculekit.util import ensurelist
from htmd.home import home

import networkx as nx

import shutil
import parmed
import numpy as np
import tempfile
import os
import logging

logger = logging.getLogger(__name__)


triala = Molecule(os.path.join(home(shareDir=True), "builder", "triala_capped.cif"))
triala_g = triala.toGraph()
triala_bb = nx.shortest_path(triala_g, source=1, target=38)
triala_g.remove_edge(14, 16)
triala_g.remove_edge(24, 26)
triala_comp = sorted(
    [list(x) for x in nx.connected_components(triala_g)],
    key=lambda c: min(c),
)
start_map = np.zeros(triala.numAtoms, dtype=bool)
start_map[triala_comp[0]] = True
end_map = np.zeros(triala.numAtoms, dtype=bool)
end_map[triala_comp[2]] = True


def _extend_residue(mol, nterm=True, cterm=True):
    # Adds ALA on either side of the residue for full parameterization
    mol = mol.copy()
    mol.resid[:] = 3

    def _ix(mol, name, resid=3):
        return np.where((mol.name == name) & (mol.resid == resid))[0][0]

    backbone_idx = nx.shortest_path(
        mol.toGraph(), source=_ix(mol, "N"), target=_ix(mol, "C")
    )
    if nterm:
        # Remove backbone formal charge since we'll add ALA
        mol.formalcharge[_ix(mol, "N")] = 0
        try:
            hn2_idx = _ix(mol, "HN2")
        except Exception:
            logger.info(
                "Could not find HN2 atom. Will create peptide bond using H atom of backbone N."
            )
            hn2_idx = _ix(mol, "H")

        startcap = triala.copy()
        startcap.align([14, 16, 18], refmol=mol, refsel=[hn2_idx] + backbone_idx[:2])
        startcap.filter(start_map, _logger=False)

        mol.remove(f"index {hn2_idx} or name HN2", _logger=False)
        mol.insert(startcap, index=0)
        mol.bonds = np.vstack((mol.bonds, [14, _ix(mol, "N")]))
        mol.bondtype = np.hstack((mol.bondtype, ["1"]))

    backbone_idx = nx.shortest_path(
        mol.toGraph(), source=_ix(mol, "N"), target=_ix(mol, "C")
    )
    if cterm:
        # Remove backbone formal charge since we'll add ALA
        mol.formalcharge[_ix(mol, "C")] = 0
        oxt_idx = _ix(mol, "OXT")

        endcap = triala.copy()
        endcap.align([18, 24, 26], refmol=mol, refsel=backbone_idx[-2:] + [oxt_idx])
        endcap.filter(end_map, _logger=False)

        mol.remove("name OXT HXT", _logger=False)
        mol.append(endcap)
        mol.bonds = np.vstack((mol.bonds, [_ix(mol, "C", 3), _ix(mol, "N", 4)]))
        mol.bondtype = np.hstack((mol.bondtype, ["1"]))

    # TODO: Need to add debumping code
    return mol


_CHARGE_METHOD_MAP = {
    "am1-bcc": "bcc",
    "gasteiger": "gas",
    None: None,
    "none": None,
}


def _write_pyodide_output(path, data):
    """Write bytes or str data returned by antechamber_pyodide to a file."""
    with open(path, "wb") as f:
        f.write(data if isinstance(data, bytes) else data.encode())


def _find_antechamber():
    """Return path to native antechamber binary, or None."""
    return shutil.which("antechamber", mode=os.X_OK)


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
    else:
        import subprocess

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
        if result.returncode != 0:
            raise RuntimeError(
                f"{program} failed (exit {result.returncode}):\n{result.stderr}"
            )
        return None


def _fftype_antechamber(
    mol,
    tmpdir,
    method="gaff2",
    netcharge=0,
    charge_method=None,
    use_pyodide=False,
    output_mol2="typed.mol2",
    output_frcmod="mol.frcmod",
):
    """Assign GAFF2 atom types and generate frcmod using antechamber/parmchk2.

    Parameters
    ----------
    mol : Molecule
        Input molecule.
    tmpdir : str
        Working directory for intermediate files.
    method : str
        Force-field atom-type method ("gaff" or "gaff2").
    netcharge : int
        Net formal charge of the molecule.
    charge_method : str or None
        Charge model passed to antechamber's ``-c`` flag.
    use_pyodide : bool
        Dispatch via antechamber_pyodide instead of subprocess.
    output_mol2 : str
        Basename for the typed mol2 file written into *tmpdir*.
    output_frcmod : str
        Basename for the frcmod file written into *tmpdir*.

    Returns
    -------
    typed_mol : Molecule
        Molecule with GAFF2 atom types assigned.
    frcmod_path : str
        Absolute path to the generated frcmod file.
    """
    input_mol2 = "input.mol2"
    mol2path = os.path.join(tmpdir, input_mol2)
    mol.write(mol2path)

    ac_charge = _CHARGE_METHOD_MAP.get(charge_method.lower() if charge_method else None)

    cmd = [
        "antechamber",
        "-i",
        input_mol2,
        "-fi",
        "mol2",
        "-o",
        output_mol2,
        "-fo",
        "mol2",
        "-at",
        method.lower(),
        "-nc",
        str(netcharge),
        "-dr",
        "n",
        "-j",
        "1",
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
        "-i",
        output_mol2,
        "-f",
        "mol2",
        "-o",
        output_frcmod,
        "-s",
        {"gaff": "1", "gaff2": "2"}[method.lower()],
        "-a",
        "Y",
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
    return typed_mol, frcmod_path


def parameterizeNonCanonicalResidues(
    cifs,
    outdir,
    forcefield="GAFF2",
    charge_model="AM1-BCC",
    is_nterm=False,
    is_cterm=False,
    backend=None,
):
    """Parameterize non-canonical residues.

    Parameters
    ----------
    cifs : str or list of str
        Path to the CIF file(s) containing the non-canonical residues.
    outdir : str
        Path to the output directory.
    forcefield : str, optional
        Force field to use for parameterization. Default is "GAFF2".
    charge_model : str, optional, choices: ["AM1-BCC", "Gasteiger", "Espaloma", "ESP", "RESP-PSI4", "None"]
        Charge model to use for parameterization
    is_nterm : bool, optional
        Whether the residue is a terminal N-term
    is_cterm : bool, optional
        Whether the residue is a terminal C-term
    backend : str or None, optional
        Force a specific backend: "parameterize", "antechamber_pyodide", or
        "antechamber_native".  When None (default) the first available backend
        is selected automatically.

    """
    cifs = ensurelist(cifs)
    if forcefield.lower() not in ("sage", "gaff2"):
        raise AttributeError(
            "Parameterization can only be performed with SAGE or GAFF2 forcefields"
        )

    for cif in cifs:
        mol = Molecule(cif)
        _parameterize_non_canonical_residue(
            mol,
            outdir,
            forcefield,
            is_nterm=is_nterm,
            is_cterm=is_cterm,
            charge_model=charge_model,
            backend=backend,
        )


def _parameterize_non_canonical_residue(
    mol,
    outdir,
    forcefield,
    charge_model,
    is_nterm=False,
    is_cterm=False,
    backend=None,
):
    if backend is not None:
        _backend = backend
    else:
        try:
            from parameterize.main import parameterize_molecule  # noqa: F401

            _backend = "parameterize"
        except ImportError:
            try:
                import antechamber_pyodide  # noqa: F401

                _backend = "antechamber_pyodide"
            except ImportError:
                if _find_antechamber() is not None:
                    _backend = "antechamber_native"
                else:
                    raise ImportError(
                        "No parameterization backend available. Install one of: "
                        "parameterize (private, contact info@acellera.com), "
                        "antechamber_pyodide (for Pyodide environments), "
                        "or AmberTools (conda install ambertools -c conda-forge) "
                        "which provides the antechamber command."
                    )

    os.makedirs(outdir, exist_ok=True)

    if len(np.intersect1d(mol.resname, triala.resname)):
        raise RuntimeError(
            f"Your residue cannot be called any of the following: {' '.join(np.unique(triala.resname))}"
        )

    mol = mol.copy()
    resn = mol.resname[0]

    with tempfile.TemporaryDirectory() as tmpdir:
        xmol = _extend_residue(mol, nterm=not is_nterm, cterm=not is_cterm)

        cmol = xmol.copy()
        cmol.resname[:] = resn
        cmol.resid[:] = 1

        if _backend == "parameterize":
            from parameterize.main import parameterize_molecule

            exclude_atoms = list(np.where(xmol.resname != resn)[0])
            cmol.write(os.path.join(tmpdir, f"{resn}.cif"))

            parameterize_molecule(
                outdir=tmpdir,
                molfile=os.path.join(tmpdir, f"{resn}.cif"),
                forcefield=forcefield,
                minimizer="ffeval",
                charge_type=charge_model,
                keep_atomnames=True,
                exclude_atoms=exclude_atoms,
            )

            typed_mol = Molecule(os.path.join(tmpdir, "parameters", f"{resn}-orig.cif"))
            frcmod = os.path.join(tmpdir, "parameters", f"{resn}.frcmod")
        else:
            netcharge = int(np.sum(cmol.formalcharge))
            typed_mol, frcmod = _fftype_antechamber(
                cmol,
                tmpdir,
                method=forcefield,
                netcharge=netcharge,
                charge_method=charge_model,
                use_pyodide=_backend == "antechamber_pyodide",
            )

        shutil.copy(frcmod, outdir)
        use_pyodide = _backend == "antechamber_pyodide"
        _post_process_parameterize(
            outdir,
            typed_mol,
            frcmod,
            resn,
            refmol=xmol,
            use_pyodide=use_pyodide,
        )


def _post_process_parameterize(
    outdir, mol, frcmod, resn, refmol=None, use_pyodide=False
):
    from collections import defaultdict

    padding_atoms = np.zeros(mol.numAtoms, dtype=bool)
    if refmol is not None:
        # TODO: Move this to parameterize (?)
        fields = ("element",)
        g1 = makeMolGraph(refmol, "all", fields)
        g2 = makeMolGraph(mol, "all", fields)
        _, _, matching = compareGraphs(
            g1, g2, fields=fields, tolerance=0.5, returnmatching=True
        )
        for pp in matching:  # Rename atoms in reference molecule
            mol.name[pp[0]] = refmol.name[pp[1]]
            mol.formalcharge[pp[0]] = refmol.formalcharge[pp[1]]
            mol.resname[pp[0]] = refmol.resname[pp[1]]

        # Antechamber can change the case of atom names (e.g. "CL" -> "Cl").
        # The graph matching above may not cover every atom, so fix any
        # remaining case mismatches using the reference molecule's names.
        ref_name_case = {n.upper(): n for n in np.unique(refmol.name)}
        for i in range(mol.numAtoms):
            upper = mol.name[i].upper()
            if upper in ref_name_case and mol.name[i] != ref_name_case[upper]:
                mol.name[i] = ref_name_case[upper]

        padding_atoms[mol.resname != resn] = True

    # Rename backbone atom types
    original = defaultdict(list)
    backbone_at = {"N": "N", "H": "H", "CA": "CT", "HA": "H1", "C": "C", "O": "O"}
    for key, val in backbone_at.items():
        sel = mol.name == key
        for at in mol.atomtype[sel]:
            original[at].append(val)
        mol.atomtype[sel] = val

    prm = parmed.amber.AmberParameterSet(frcmod)

    # Duplicate parameters for renamed atoms which don't exist in the backbone
    _duplicate_parameters(prm, original)

    # Remove unused parameters
    _clean_prm(prm, mol, list(backbone_at.values()), padding_atoms)
    prm.write(os.path.join(outdir, f"{resn}.frcmod"))

    # Remove the caps
    mol.remove(padding_atoms, _logger=False)

    with tempfile.TemporaryDirectory() as tmpdir:
        mol2file = os.path.join(tmpdir, f"{resn}.mol2")
        mol.write(mol2file)

        acfile = os.path.join(tmpdir, f"{resn}_mod.ac")
        ac_cmd = [
            "antechamber",
            "-i",
            f"{resn}.mol2",
            "-fi",
            "mol2",
            "-o",
            f"{resn}_mod.ac",
            "-fo",
            "ac",
            "-nc",
            str(mol.formalcharge.sum()),
            "-pf",
            "y",
            "-dr",
            "n",
            "-j",
            "0",
            "-an",
            "n",
        ]

        input_files = None
        if use_pyodide:
            with open(mol2file, "rb") as f:
                input_files = {f"{resn}.mol2": f.read()}

        result = _run_ambertools(
            "antechamber",
            ac_cmd,
            cwd=tmpdir,
            use_pyodide=use_pyodide,
            input_files=input_files,
        )
        if use_pyodide:
            _write_pyodide_output(acfile, result[f"{resn}_mod.ac"])

        mol_g = mol.toGraph()
        n_idx = np.where(mol.name == "N")[0][0]
        c_idx = np.where(mol.name == "C")[0][0]
        backbone = nx.shortest_path(mol_g, source=n_idx, target=c_idx)

        mainchain = os.path.join(tmpdir, f"mainchain.{resn.lower()}")
        with open(mainchain, "w") as f:
            f.write("HEAD_NAME N\n")
            f.write("TAIL_NAME C\n")
            for bb_idx in backbone[1:-1]:
                f.write(f"MAIN_CHAIN {mol.name[bb_idx]}\n")
            f.write("PRE_HEAD_TYPE C\n")
            f.write("POST_TAIL_TYPE N\n")
            f.write(f"CHARGE {mol.formalcharge.sum():.1f}\n")

        prepi = os.path.join(outdir, f"{resn}.prepi")
        pg_cmd = [
            "prepgen",
            "-i",
            f"{resn}_mod.ac",
            "-o",
            f"{resn}.prepi",
            "-f",
            "prepi",
            "-m",
            f"mainchain.{resn.lower()}",
            "-rn",
            resn,
        ]

        input_files = None
        if use_pyodide:
            with open(acfile, "rb") as f:
                ac_bytes = f.read()
            with open(mainchain, "rb") as f:
                mc_bytes = f.read()
            input_files = {
                f"{resn}_mod.ac": ac_bytes,
                f"mainchain.{resn.lower()}": mc_bytes,
            }

        result = _run_ambertools(
            "prepgen",
            pg_cmd,
            cwd=tmpdir,
            use_pyodide=use_pyodide,
            input_files=input_files,
        )
        if use_pyodide:
            _write_pyodide_output(prepi, result[f"{resn}.prepi"])
        else:
            shutil.copy(os.path.join(tmpdir, f"{resn}.prepi"), prepi)

        _fix_prepi_atomname_capitalization(mol, prepi)


def _fix_prepi_atomname_capitalization(mol, prepi):
    # Fix wrong prepi atom name capitalization
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
            # Fix wrong prepi atom name capitalization
            if old_name.upper() in uqnames and uqnames[old_name.upper()] != old_name:
                logger.info(
                    f"Fixed residue {mol.resname[0]} atom name {old_name} -> {uqnames[old_name.upper()]}"
                    " to match the input structure."
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
                        f"Fixed residue {mol.resname[0]} atom name {piece} -> {uqnames[piece.upper()]}"
                        " to match the input structure."
                    )
                    lines[i] = (
                        f"{lines[i][:j*5]}{uqnames[piece.upper()] : >5}{lines[i][(j+1)*5:]}"
                    )
            continue

        if lines[i].strip().startswith("CORR"):
            section = "atoms"
        if lines[i].strip().startswith("LOOP"):
            section = "loop"
        if lines[i].strip().startswith("IMPROPER"):
            section = "improper"

    with open(prepi, "w") as f:
        for line in lines:
            f.write(line)


def _clean_prm(prm, mol, backbone_at, padding_atoms):
    from moleculekit.util import calculateAnglesAndDihedrals

    all_at = np.unique(mol.atomtype[~padding_atoms]).tolist() + backbone_at
    angles, dihedrals = calculateAnglesAndDihedrals(mol.bonds)

    at_dict = {
        "bond_types": mol.atomtype[mol.bonds].tolist(),
        "angle_types": mol.atomtype[angles].tolist(),
        "dihedral_types": mol.atomtype[dihedrals].tolist(),
    }

    to_delete = []
    for at in prm.atom_types:
        if at not in all_at:
            to_delete.append(at)
    for at in to_delete:
        del prm.atom_types[at]

    for param_t in ("bond_types", "angle_types", "dihedral_types"):
        to_delete = []
        seen = []
        for bt in getattr(prm, param_t):
            if (
                not all(np.isin(bt, all_at))
                or all(np.isin(bt, backbone_at))
                or (
                    list(bt) not in at_dict[param_t]
                    and list(bt)[::-1] not in at_dict[param_t]
                )
                or list(bt)[::-1] in seen
            ):
                to_delete.append(bt)
            seen.append(list(bt))

        for bt in to_delete:
            del prm.__dict__[param_t][bt]

    for param_t in ("improper_types", "improper_periodic_types"):
        to_delete = []
        for bt in getattr(prm, param_t):
            if not all(np.isin(bt, all_at)) or all(np.isin(bt, backbone_at)):
                to_delete.append(bt)
        for bt in to_delete:
            del prm.__dict__[param_t][bt]


def _duplicate_parameters(prm, original):
    from copy import deepcopy

    def _gen_permutations(typ, original):
        import itertools

        possibles = []
        for tt in typ:
            if tt in original:
                possibles.append(original[tt] + [tt])
            else:
                possibles.append([tt])
        return list(itertools.product(*possibles))

    # Duplicate parameters for renamed atoms which don't exist in the backbone
    for prm_typ in ("bond_types", "angle_types", "dihedral_types", "improper_types"):
        new_bt = {}
        for typ in getattr(prm, prm_typ):
            perms = _gen_permutations(typ, original)
            if len(perms) == 1:
                continue  # No replacements
            for perm in perms:
                new_bt[perm] = getattr(prm, prm_typ)[typ]

        for bt in new_bt:
            prm.__dict__[prm_typ][bt] = deepcopy(new_bt[bt])
