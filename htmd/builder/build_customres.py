# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.tools.graphalignment import makeMolGraph, compareGraphs
from moleculekit.molecule import Molecule
from htmd.home import home
import unittest
import shutil
import parmed
import numpy as np
import tempfile
import os

from moleculekit.util import ensurelist


ace = Molecule(os.path.join(home(shareDir=True), "builder", "caps", "ACE.cif"))
nme = Molecule(os.path.join(home(shareDir=True), "builder", "caps", "NME.cif"))


def _cap_residue(mol):
    mol = mol.copy()
    resn = mol.resname[0]

    acec = ace.copy()
    acec.align("name N H CA", refmol=mol)
    acec.remove("name N H CA", _logger=False)

    nmec = nme.copy()
    nmec.align("name CA C O", refmol=mol)
    nmec.remove("name CA C O", _logger=False)

    mol.insert(acec, index=0)
    mol.append(nmec)
    # Bond ACE/NME to the residue
    mol.bonds = np.vstack(
        (
            mol.bonds,
            [
                np.where((mol.resname == "NME") & (mol.name == "N"))[0][0],
                np.where((mol.resname == resn) & (mol.name == "C"))[0][0],
            ],
            [
                np.where((mol.resname == "ACE") & (mol.name == "C"))[0][0],
                np.where((mol.resname == resn) & (mol.name == "N"))[0][0],
            ],
        )
    )
    mol.bondtype = np.hstack((mol.bondtype, ["1"], ["1"]))
    return mol


def parameterizeCustomResidues(cifs, outdir, method="gaff2", nnp=None):
    cifs = ensurelist(cifs)
    method = method.lower()
    if method == "ani-2x" and nnp is None:
        raise RuntimeError(
            "The user must provide an NNP calculator to user ANI-2x parameterization"
        )
    if method not in ("gaff2", "ani-2x"):
        raise AttributeError(
            "Parameterization can only be performed with GAFF2 or ANI-2x methods"
        )

    for cif in cifs:
        mol = Molecule(cif)
        _parameterize_custom_residue(mol, outdir, method, nnp=nnp)


def _parameterize_custom_residue(mol, outdir, method, nnp=None):
    try:
        from parameterize.parameterization.cli import main_parameterize
    except ImportError:
        raise ImportError(
            "You are missing the parameterize library. Please install it with conda install parameterize -c acellera -c conda-forge"
        )

    os.makedirs(outdir, exist_ok=True)

    mol = mol.copy()
    # Remove backbone formal charge from templated molecule
    mol.formalcharge[mol.name == "N"] = 0
    resn = mol.resname[0]
    with tempfile.TemporaryDirectory() as tmpdir:
        cmol = _cap_residue(mol)

        sdffile = os.path.join(tmpdir, "mol.sdf")
        cmol.write(sdffile)

        # TODO: Improve it to not parameterize the cap dihedrals by excluding those dihedrals
        main_parameterize(
            cmol,
            user_charge=int(cmol.formalcharge.sum()),
            forcefield="GAFF2",
            charge_type="AM1-BCC",
            min_type="mm",
            dihed_fit_type="iterative",
            dihed_opt_type="mm",
            fit_dihedral=method != "gaff2",
            nnp=nnp,
            outdir=tmpdir,
        )
        shutil.copy(
            os.path.join(tmpdir, "parameters", "GAFF2", "mol-orig.mol2"),
            os.path.join(tmpdir, f"{resn}.mol2"),
        )
        shutil.copy(
            os.path.join(tmpdir, "parameters", "GAFF2", "mol.frcmod"),
            os.path.join(outdir, f"{resn}.frcmod"),
        )

        _post_process_parameterize(cmol, tmpdir, outdir, resn)


def _post_process_parameterize(cmol, tmpdir, outdir, resn):
    from subprocess import call

    # TODO: Move this to parameterize (?)
    mol = Molecule(os.path.join(tmpdir, f"{resn}.mol2"))
    fields = ("element",)
    g1 = makeMolGraph(cmol, "all", fields)
    g2 = makeMolGraph(mol, "all", fields)
    _, _, matching = compareGraphs(
        g1, g2, fields=fields, tolerance=0.5, returnmatching=True
    )
    for pp in matching:  # Rename atoms in reference molecule and copy formal charges
        mol.name[pp[0]] = cmol.name[pp[1]]
        mol.formalcharge[pp[0]] = cmol.formalcharge[pp[1]]

    # Remove the caps
    mask = np.ones(mol.numAtoms, dtype=bool)
    mask[:6] = False
    mask[-6:] = False

    mol.filter(mask, _logger=False)

    # Rename backbone atom types
    backbone_at = {"N": "N", "H": "H", "CA": "CT", "HA": "H1", "C": "C", "O": "O"}
    original = {}
    for key, val in backbone_at.items():
        original[mol.atomtype[mol.name == key][0]] = backbone_at[key]
        mol.atomtype[mol.name == key] = val

    mol2file = os.path.join(tmpdir, f"{resn}.mol2")
    mol.write(mol2file)

    frcmod = os.path.join(outdir, f"{resn}.frcmod")
    prm = parmed.amber.AmberParameterSet(frcmod)

    # Duplicate parameters for renamed atoms which don't exist in the backbone
    _duplicate_parameters(prm, original)

    # Remove unused parameters
    _clean_prm(prm, mol)
    prm.write(frcmod)

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
        ]
    )

    # TODO: Recover the lost formal charge!!!
    mainchain = os.path.join(tmpdir, f"mainchain.{resn.lower()}")
    with open(mainchain, "w") as f:
        f.write("HEAD_NAME N\n")
        f.write("TAIL_NAME C\n")
        f.write("MAIN_CHAIN CA\n")
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


def _clean_prm(prm, mol):
    from moleculekit.util import guessAnglesAndDihedrals

    all_at = np.unique(mol.atomtype)
    angles, dihedrals = guessAnglesAndDihedrals(mol.bonds)
    bond_at = mol.atomtype[mol.bonds].tolist()
    angle_at = mol.atomtype[angles].tolist()
    dihed_at = mol.atomtype[dihedrals].tolist()

    to_delete = []
    for at in prm.atom_types:
        if at not in all_at.tolist():
            to_delete.append(at)
    for at in to_delete:
        del prm.atom_types[at]

    to_delete = []
    for bt in prm.bond_types:
        if not np.any(np.isin(bt, all_at)):
            to_delete.append(bt)
        elif list(bt) not in bond_at and list(bt)[::-1] not in bond_at:
            to_delete.append(bt)
    for bt in to_delete:
        del prm.bond_types[bt]

    to_delete = []
    for bt in prm.angle_types:
        if not np.any(np.isin(bt, all_at)):
            to_delete.append(bt)
        elif list(bt) not in angle_at and list(bt)[::-1] not in angle_at:
            to_delete.append(bt)
    for bt in to_delete:
        del prm.angle_types[bt]

    to_delete = []
    for bt in prm.dihedral_types:
        if not np.any(np.isin(bt, all_at)):
            to_delete.append(bt)
        elif list(bt) not in dihed_at and list(bt)[::-1] not in dihed_at:
            to_delete.append(bt)
    for bt in to_delete:
        del prm.dihedral_types[bt]

    to_delete = []
    for bt in prm.improper_types:
        if not np.any(np.isin(bt, all_at)):
            to_delete.append(bt)
    for bt in to_delete:
        del prm.improper_types[bt]


def _duplicate_parameters(prm, original):
    # Duplicate parameters for renamed atoms which don't exist in the backbone
    new_bt = {}
    for bt in prm.bond_types:
        bt = np.array(bt)
        # If one of the two atom types was replaced duplicate the param
        replaced = np.isin(bt, list(original.keys()))
        if replaced.sum() != 1:
            continue
        bt_new = bt.copy()
        bt_new[replaced] = [original[val] for val in bt[replaced]]
        new_bt[tuple(bt_new)] = prm.bond_types[tuple(bt)]

    for bt in new_bt:
        prm.bond_types[bt] = new_bt[bt]

    new_bt = {}
    for bt in prm.angle_types:
        bt = np.array(bt)
        # If one of the two atom types was replaced duplicate the param
        replaced = np.isin(bt, list(original.keys()))
        if replaced.sum() == 0 or replaced.sum() == 3:
            continue
        bt_new = bt.copy()
        bt_new[replaced] = [original[val] for val in bt[replaced]]
        new_bt[tuple(bt_new)] = prm.angle_types[tuple(bt)]

    for bt in new_bt:
        prm.angle_types[bt] = new_bt[bt]

    new_bt = {}
    for bt in prm.dihedral_types:
        bt = np.array(bt)
        # If one of the two atom types was replaced duplicate the param
        replaced = np.isin(bt, list(original.keys()))
        if replaced.sum() == 0 or replaced.sum() == 4:
            continue
        bt_new = bt.copy()
        bt_new[replaced] = [original[val] for val in bt[replaced]]
        new_bt[tuple(bt_new)] = prm.dihedral_types[tuple(bt)]

    for bt in new_bt:
        prm.dihedral_types[bt] = new_bt[bt]

    new_bt = {}
    for bt in prm.improper_types:
        bt = np.array(bt)
        # If one of the two atom types was replaced duplicate the param
        replaced = np.isin(bt, list(original.keys()))
        if replaced.sum() == 0 or replaced.sum() == 4:
            continue
        bt_new = bt.copy()
        bt_new[replaced] = [original[val] for val in bt[replaced]]
        new_bt[tuple(bt_new)] = prm.improper_types[tuple(bt)]

    for bt in new_bt:
        prm.improper_types[bt] = new_bt[bt]


class _TestCustomResParam(unittest.TestCase):
    def test_custom_residue_parameterization(self):
        from glob import glob

        refdir = home(dataDir=os.path.join("test-custom-residue-param"))
        cifs = glob(os.path.join(refdir, "*.cif"))

        for cif in cifs:
            with tempfile.TemporaryDirectory() as tmpdir, self.subTest(cif):
                parameterizeCustomResidues(cif, tmpdir, method="gaff2")
                params = glob(os.path.join(tmpdir, "*"))
                assert len(params) == 2


if __name__ == "__main__":
    unittest.main(verbosity=2)
