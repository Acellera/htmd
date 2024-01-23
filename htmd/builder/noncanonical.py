# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.tools.graphalignment import makeMolGraph, compareGraphs
from moleculekit.molecule import Molecule
from moleculekit.util import ensurelist
from htmd.home import home
from glob import glob
import networkx as nx
import unittest
import shutil
import parmed
import numpy as np
import tempfile
import os
import logging

logger = logging.getLogger(__name__)


triala = Molecule(os.path.join(home(shareDir=True), "builder", "triala_capped.mmtf"))
triala_g = triala.toGraph()
triala_bb = nx.shortest_path(triala_g, source=1, target=38)
triala_g.remove_edge(14, 16)
triala_g.remove_edge(24, 26)
triala_comp = [list(x) for x in list(nx.connected_components(triala_g))]
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


def parameterizeNonCanonicalResidues(
    cifs,
    outdir,
    forcefield="GAFF2",
    calculator="AIMNet2",
    is_nterm=False,
    is_cterm=False,
):
    cifs = ensurelist(cifs)
    if forcefield.lower() not in ("sage", "gaff2"):
        raise AttributeError(
            "Parameterization can only be performed with SAGE or GAFF2 forcefields"
        )

    for cif in cifs:
        mol = Molecule(cif)
        _parameterize_non_canonical_residue(
            mol, outdir, forcefield, calculator, is_nterm=is_nterm, is_cterm=is_cterm
        )


def _parameterize_non_canonical_residue(
    mol, outdir, forcefield, calculator, is_nterm=False, is_cterm=False
):
    try:
        from parameterize.cli import main_parameterize, list_dihedrals
    except ImportError:
        raise ImportError(
            "You are missing the parameterize library. It's only available for private installations. Please contact info@acellera.com"
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
        exclude_atoms = list(np.where(xmol.resname != resn)[0])

        cmol = xmol.copy()
        # So that CIF is written as smallmol instead of macromol
        cmol.resname[:] = resn

        dih = list_dihedrals(cmol, exclude_atoms=exclude_atoms, _logger=False)

        # Delete pure backbone dihedrals to not reparameterize them and thus waste atomtype names and computation time
        for dd in ("CH3-C-N-CA", "CA-C-N-CH3", "C-N-CA-C", "N-CA-C-N", "CA-C-N-CA"):
            if dd in dih:
                del dih[dd]

        main_parameterize(
            cmol,
            forcefield=forcefield,
            charge_type="AM1-BCC",
            min_type="mm",
            dihed_fit_type="iterative",
            dihed_opt_type="mm",
            fit_dihedral=calculator is not None,
            dihedrals=list(dih.values()),
            calculator=calculator,
            outdir=tmpdir,
            exclude_atoms=exclude_atoms,
        )
        shutil.copy(
            glob(os.path.join(tmpdir, "parameters", "*", f"{resn}-orig.cif"))[0],
            os.path.join(tmpdir, f"{resn}.cif"),
        )
        shutil.copy(
            glob(os.path.join(tmpdir, "parameters", "*", f"{resn}.frcmod"))[0],
            os.path.join(outdir, f"{resn}.frcmod"),
        )
        _post_process_parameterize(xmol, tmpdir, outdir, resn)


def _post_process_parameterize(xmol, tmpdir, outdir, resn):
    from subprocess import call
    from collections import defaultdict

    # TODO: Move this to parameterize (?)
    mol = Molecule(os.path.join(tmpdir, f"{resn}.cif"))
    fields = ("element",)
    g1 = makeMolGraph(xmol, "all", fields)
    g2 = makeMolGraph(mol, "all", fields)
    _, _, matching = compareGraphs(
        g1, g2, fields=fields, tolerance=0.5, returnmatching=True
    )
    for pp in matching:  # Rename atoms in reference molecule
        mol.name[pp[0]] = xmol.name[pp[1]]
        mol.formalcharge[pp[0]] = xmol.formalcharge[pp[1]]
        mol.resname[pp[0]] = xmol.resname[pp[1]]

    # Rename backbone atom types
    backbone_at = {"N": "N", "H": "H", "CA": "CT", "HA": "H1", "C": "C", "O": "O"}
    padding_atoms = np.zeros(mol.numAtoms, dtype=bool)
    padding_atoms[xmol.resname != resn] = True
    original = defaultdict(list)

    for key, val in backbone_at.items():
        sel = mol.name == key
        for at in mol.atomtype[sel]:
            original[at].append(val)
        mol.atomtype[sel] = val

    frcmod = os.path.join(outdir, f"{resn}.frcmod")
    prm = parmed.amber.AmberParameterSet(frcmod)

    # Duplicate parameters for renamed atoms which don't exist in the backbone
    _duplicate_parameters(prm, original)

    # Remove unused parameters
    _clean_prm(prm, mol, list(backbone_at.values()), padding_atoms)
    prm.write(frcmod)

    # Remove the caps
    mol.remove(padding_atoms, _logger=False)

    mol2file = os.path.join(tmpdir, f"{resn}.mol2")
    mol.write(mol2file)

    acfile = os.path.join(tmpdir, f"{resn}_mod.ac")
    call(
        [
            "antechamber",
            "-fi",
            "mol2",
            "-fo",
            "ac",
            "-i",
            mol2file,
            "-o",
            acfile,
            "-dr",
            "n",
            "-j",
            "0",
            "-nc",
            str(mol.formalcharge.sum()),
            "-an",  # adjust atom names
            "n",  # no
            "-pf",
            "y",
        ]
    )

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
    call(
        [
            "prepgen",
            "-i",
            acfile,
            "-o",
            prepi,
            "-f",
            "prepi",
            "-m",
            mainchain,
            "-rn",
            resn,
        ]
    )
    print("")  # Add a newline after the output of prepgen

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
                lines[
                    i
                ] = f"{lines[i][:6]}{uqnames[old_name.upper()]:4s}{lines[i][10:]}"
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
                    lines[
                        i
                    ] = f"{lines[i][:j*5]}{uqnames[piece.upper()] : >5}{lines[i][(j+1)*5:]}"
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
    from moleculekit.util import guessAnglesAndDihedrals

    all_at = np.unique(mol.atomtype[~padding_atoms]).tolist() + backbone_at
    angles, dihedrals = guessAnglesAndDihedrals(mol.bonds)

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


# The below is used for testing only
try:
    from parameterize.cli import main_parameterize
except ImportError:
    _PARAMETERIZE_INSTALLED = False
else:
    _PARAMETERIZE_INSTALLED = True


def _compare_frcmod(
    refFile,
    resFile,
    dihedralForceConstAbsTol=1e-6,
    dihedralPhaseAbsTol=1e-6,
):
    import difflib

    # Default tolerances
    defaultAbsTol = 1e-6
    defaultRelTol = 1e-6

    with open(refFile) as ref, open(resFile) as res:
        refLines, resLines = ref.readlines(), res.readlines()

    # Removes the first line with the HTMD version
    refLines = refLines[1:]
    resLines = resLines[1:]

    if len(refLines) != len(resLines):
        raise AssertionError(
            f"Reference file {refFile} has a different number of lines ({len(refLines)}) from file {resFile} ({len(resLines)})"
        )

    # Iterate over lines in the file
    for refLine, resLine in zip(refLines, resLines):
        refFields = refLine.split()
        resFields = resLine.split()

        # Iterate over fields in the line
        for iField, (refField, resField) in enumerate(zip(refFields, resFields)):
            # Set tolerance
            if (
                len(refFields) == 7 and iField == 2
            ):  # Detect dihedral force constant column
                absTol = dihedralForceConstAbsTol
                relTol = 0
            elif len(refFields) == 7 and iField == 3:  # Detect dihedral phase column
                absTol = dihedralPhaseAbsTol
                relTol = 0
            else:
                absTol = defaultAbsTol
                relTol = defaultRelTol

            try:
                refField = float(refField)
                resField = float(resField)
            except ValueError:
                # The fields cannot be converted to floats, so compare them directly
                if refField == resField:
                    continue
            else:
                # The fields can be converted to floats, so compare them with the tolerances
                if np.isclose(refField, resField, atol=absTol, rtol=relTol):
                    continue

            # Print in case of failure
            print(f"Failed: {refField} == {resField}")
            print(f"Absolute tolerance: {absTol}")
            print(f"Relative tolernace: {relTol}")

            diff = difflib.unified_diff(
                refLines, resLines, fromfile=refFile, tofile=resFile, n=1
            )
            if len(list(diff)):
                raise RuntimeError("".join(diff))


def _compare_prepis(refFile, resFile):
    import difflib

    with open(refFile) as f:
        reflines = f.readlines()
    with open(resFile) as f:
        newlines = f.readlines()

    diff = difflib.unified_diff(
        reflines, newlines, fromfile=refFile, tofile=resFile, n=1
    )
    if len(list(diff)):
        raise RuntimeError("".join(diff))


class _TestNCAAResParam(unittest.TestCase):
    @unittest.skipUnless(
        _PARAMETERIZE_INSTALLED, "Can only run with parameterize installed"
    )
    def test_ncaa_residue_parameterization(self):
        # import shutil

        refdir = home(dataDir=os.path.join("test-custom-residue-param"))
        refresdir = os.path.join(refdir, "gaff2-params")
        cifs = glob(os.path.join(refdir, "*.cif"))
        for cif in cifs:
            with tempfile.TemporaryDirectory() as tmpdir, self.subTest(cif):
                parameterizeNonCanonicalResidues(
                    cif, tmpdir, forcefield="GAFF2", calculator=None
                )

                frcmod = glob(os.path.join(tmpdir, "*.frcmod"))[0]
                name = os.path.basename(frcmod)
                _compare_frcmod(os.path.join(refresdir, name), frcmod)

                prepi = glob(os.path.join(tmpdir, "*.prepi"))[0]
                name = os.path.basename(prepi)
                _compare_prepis(os.path.join(refresdir, name), prepi)

        refresdir = os.path.join(refdir, "gaff2-params-terminals")
        cifs = glob(os.path.join(refdir, "*.cif"))
        for cif in cifs:
            with tempfile.TemporaryDirectory() as tmpdir, self.subTest(cif):
                parameterizeNonCanonicalResidues(
                    cif,
                    tmpdir,
                    forcefield="GAFF2",
                    is_cterm=True,
                    is_nterm=True,
                    calculator=None,
                )

                frcmod = glob(os.path.join(tmpdir, "*.frcmod"))[0]
                name = os.path.basename(frcmod)
                # shutil.copy(frcmod, os.path.join(refresdir, name))
                _compare_frcmod(os.path.join(refresdir, name), frcmod)

                prepi = glob(os.path.join(tmpdir, "*.prepi"))[0]
                name = os.path.basename(prepi)
                # shutil.copy(prepi, os.path.join(refresdir, name))
                _compare_prepis(os.path.join(refresdir, name), prepi)

    @unittest.skipUnless(
        _PARAMETERIZE_INSTALLED, "Can only run with parameterize installed"
    )
    def test_ncaa_residue_parameterization_xTB(self):
        refdir = home(dataDir=os.path.join("test-custom-residue-param"))
        refresdir = os.path.join(refdir, "xtb-params")
        cif = os.path.join(refdir, "33X.cif")
        with tempfile.TemporaryDirectory() as tmpdir, self.subTest(cif):
            parameterizeNonCanonicalResidues(
                cif, tmpdir, forcefield="GAFF2", calculator="xTB"
            )
            frcmod = glob(os.path.join(tmpdir, "*.frcmod"))[0]
            name = os.path.basename(frcmod)
            _compare_frcmod(os.path.join(refresdir, name), frcmod)

            prepi = glob(os.path.join(tmpdir, "*.prepi"))[0]
            name = os.path.basename(prepi)
            _compare_prepis(os.path.join(refresdir, name), prepi)


if __name__ == "__main__":
    unittest.main(verbosity=2)
