import numpy as np
import pytest
from moleculekit.molecule import Molecule
from htmd.builder.amber import _assign_modrna_terminal_forms
from htmd.builder.builder import BuildError

A6L_CIF = "tests/test_nonstandard_builder/6A6L.cif"


def _rna_5mc_mol():
    # 6A6L chain D is the RNA C-A-U-5MC with 5MC at the 3' terminus
    mol = Molecule(A6L_CIF)
    mol.remove("not chain D", _logger=False)
    mol.remove("water", _logger=False)
    return mol


def test_three_prime_terminal_5mc_renamed_to_5MC3():
    mol = _rna_5mc_mol()
    assert (mol.resname == "5MC").any()
    _assign_modrna_terminal_forms(mol)
    # the 3'-terminal 5MC becomes 5MC3; no bare 5MC remains
    assert (mol.resname == "5MC3").any()
    assert not (mol.resname == "5MC").any()


def test_no_modrna_is_noop():
    mol = Molecule(A6L_CIF)
    mol.remove("nucleic", _logger=False)  # protein only
    before = mol.resname.copy()
    _assign_modrna_terminal_forms(mol)
    assert np.array_equal(before, mol.resname)


def test_lone_modified_residue_raises():
    mol = _rna_5mc_mol()
    # keep only the 5MC residue -> lone one-residue chain
    mol.remove("not resname 5MC", _logger=False)
    with pytest.raises(BuildError):
        _assign_modrna_terminal_forms(mol)
