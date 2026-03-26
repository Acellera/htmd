import numpy as np
import pytest
import scipy.spatial.distance as dist
from moleculekit.molecule import Molecule
from htmd.builder.ionize import ionize, ionizePlace, _ionGetCharge


class _TestIonGetCharge:
    def _test_monovalent_cations(self):
        assert _ionGetCharge("NA") == 1
        assert _ionGetCharge("K") == 1
        assert _ionGetCharge("CS") == 1

    def _test_divalent_cations(self):
        assert _ionGetCharge("MG") == 2
        assert _ionGetCharge("CA") == 2
        assert _ionGetCharge("ZN") == 2

    def _test_anion(self):
        assert _ionGetCharge("CL") == -1

    def _test_unknown_ion_raises(self):
        with pytest.raises(NameError, match="not in the database"):
            _ionGetCharge("XX")


class _TestIonize:
    def _test_neutralize_positive_charge_monovalent(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=5, nwater=1000, neutralize=True
        )
        assert nanion == 5
        assert ncation == 0
        assert anion == "CL"
        assert cation == "NA"

    def _test_neutralize_negative_charge_monovalent(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-3, nwater=1000, neutralize=True
        )
        assert nanion == 0
        assert ncation == 3

    def _test_neutral_system(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=0, nwater=1000, neutralize=True
        )
        assert nanion == 0
        assert ncation == 0

    def _test_neutralize_positive_charge_divalent_anion(self):
        # Even charge: straightforward
        anion, cation, nanion, ncation = ionize(
            None, netcharge=4, nwater=1000, neutralize=True, anion="CL"
        )
        assert nanion * (-1) + ncation * 1 + 4 == 0

    def _test_neutralize_negative_charge_divalent_cation(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-4, nwater=1000, neutralize=True, cation="MG"
        )
        assert nanion * (-1) + ncation * 2 + (-4) == 0

    def _test_neutralize_odd_charge_divalent_cation(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-3, nwater=1000, neutralize=True, cation="MG"
        )
        assert nanion * (-1) + ncation * 2 + (-3) == 0
        assert nanion >= 1  # Need at least one anion to balance

    def _test_saltconc_returns_six_values(self):
        mol = _make_mol_with_atoms(0)
        result = ionize(mol, netcharge=0, nwater=5000, saltconc=0.15)
        assert len(result) == 6

    def _test_saltconc_adds_ions(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=0, nwater=5000, saltconc=0.15
        )
        # num = floor(0.5 + 0.0187 * 0.15 * 5000) = floor(14.525) = 14
        assert nanion == 14
        assert ncation == 14

    def _test_saltconc_zero_neutralize_only(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=2, nwater=5000, saltconc=0
        )
        assert nanion == 2
        assert ncation == 0

    def _test_saltconc_plus_neutralization(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=3, nwater=5000, saltconc=0.15
        )
        # Neutralization: 3 Cl-
        # Salt: 14 NaCl pairs
        assert nanion == 3 + 14
        assert ncation == 0 + 14

    def _test_unknown_cation_raises(self):
        with pytest.raises(NameError):
            ionize(None, netcharge=1, nwater=1000, neutralize=True, cation="XX")

    def _test_unknown_anion_raises(self):
        with pytest.raises(NameError):
            ionize(None, netcharge=1, nwater=1000, neutralize=True, anion="XX")

    def _test_rounding_of_charge(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=2.7, nwater=1000, neutralize=True
        )
        assert nanion == 3
        assert ncation == 0


class _TestIonizePlace:
    def _test_zero_ions_returns_copy(self):
        solvent = _make_water_box(nwaters=50)
        result = ionizePlace(solvent, None, "Cl-", "Na+", "Cl-", "Na+", 0, 0)
        assert result.numAtoms == solvent.numAtoms

    def _test_does_not_modify_input(self):
        solvent, solute = _make_water_and_solute(nwaters=50, solute_atoms=5)
        natoms_before = solvent.numAtoms
        ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert solvent.numAtoms == natoms_before

    def _test_correct_atom_count(self):
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        nanion, ncation = 3, 2
        nions = nanion + ncation
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", nanion, ncation)
        expected = solvent.numAtoms - nions * 3 + nions
        assert result.numAtoms == expected

    def _test_ion_resnames(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 3)
        ion_mask = result.segid == "I"
        assert ion_mask.sum() == 5
        assert np.sum(result.resname[ion_mask] == "Cl-") == 2
        assert np.sum(result.resname[ion_mask] == "Na+") == 3

    def _test_ion_names(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "CLA", "SOD", 1, 1)
        ion_mask = result.segid == "I"
        assert "CLA" in result.name[ion_mask]
        assert "SOD" in result.name[ion_mask]

    def _test_ion_chain_and_segid(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        ion_mask = result.segid == "I"
        assert np.all(result.chain[ion_mask] == "I")
        assert np.all(result.segid[ion_mask] == "I")

    def _test_ions_far_from_solute(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=300, solute_atoms=5, box_size=40)
        dfrom = 5.0
        result = ionizePlace(
            solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2, dfrom=dfrom
        )
        ion_coords = result.coords[result.segid == "I", :, 0]
        solute_coords = solute.coords[:, :, 0]
        dists = dist.cdist(ion_coords, solute_coords)
        assert np.all(dists.min(axis=1) >= dfrom - 0.01)

    def _test_ions_far_from_each_other(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=500, solute_atoms=5, box_size=50)
        dbetween = 5.0
        result = ionizePlace(
            solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 3, 3, dbetween=dbetween
        )
        ion_coords = result.coords[result.segid == "I", :, 0]
        pdists = dist.cdist(ion_coords, ion_coords)
        np.fill_diagonal(pdists, np.inf)
        assert np.all(pdists >= dbetween - 0.01)

    def _test_unique_resids_for_ions(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 3, 3)
        ion_resids = result.resid[result.segid == "I"]
        assert len(np.unique(ion_resids)) == 6

    def _test_no_eligible_water_raises(self):
        solvent, solute = _make_water_and_solute(nwaters=3, solute_atoms=5, box_size=5)
        with pytest.raises(NameError, match="No waters could be found"):
            ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 10, 10, dfrom=100)

    def _test_water_count_decreases(self):
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        nwat_before = solvent.atomselect("water and noh").sum()
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 3)
        nwat_after = result.atomselect("water and noh").sum()
        assert nwat_after == nwat_before - 5

    def _test_no_solute(self):
        solvent = _make_water_box(nwaters=100)
        result = ionizePlace(solvent, None, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert result.numAtoms == solvent.numAtoms - 4 * 3 + 4

    def _test_empty_solute(self):
        solvent = _make_water_box(nwaters=100)
        solute = Molecule()
        solute.empty(0)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert result.numAtoms == solvent.numAtoms - 4 * 3 + 4

    def _test_ion_types_spatially_mixed(self):
        """Anions and cations should be spatially interleaved, not clustered."""
        solvent, solute = _make_water_and_solute(nwaters=500, solute_atoms=5, box_size=50)
        nanion, ncation = 5, 5
        result = ionizePlace(
            solvent, solute, "Cl-", "Na+", "Cl-", "Na+", nanion, ncation, dbetween=5
        )
        ion_mask = result.segid == "I"
        ion_coords = result.coords[ion_mask, :, 0]
        ion_resnames = result.resname[ion_mask]

        # For each ion, find its nearest neighbor and check it's not always the same type
        pdists = dist.cdist(ion_coords, ion_coords)
        np.fill_diagonal(pdists, np.inf)
        nearest = np.argmin(pdists, axis=1)

        same_type_neighbor = sum(
            ion_resnames[i] == ion_resnames[nearest[i]]
            for i in range(len(ion_resnames))
        )
        # With good mixing, fewer than 80% of ions should have a same-type nearest neighbor
        assert same_type_neighbor < 0.8 * (nanion + ncation)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_mol_with_atoms(natoms):
    """Create a minimal molecule with the given number of atoms."""
    mol = Molecule()
    mol.empty(natoms, numFrames=1)
    return mol


def _make_water_and_solute(nwaters=100, solute_atoms=0, box_size=30):
    """Create separate solvent and solute molecules for testing."""
    combined = _make_water_box(nwaters=nwaters, solute_atoms=solute_atoms, box_size=box_size)
    water_sel = combined.atomselect("water")
    solvent = combined.copy()
    solvent.filter(water_sel, _logger=False)
    if solute_atoms > 0:
        solute = combined.copy()
        solute.filter(~water_sel, _logger=False)
    else:
        solute = None
    return solvent, solute


def _make_water_box(nwaters=100, solute_atoms=0, box_size=30):
    """Create a synthetic molecule with water and optional solute atoms.

    Solute atoms are placed at the center of the box.
    Water molecules are distributed on a grid filling the box.
    """
    natoms = solute_atoms + nwaters * 3
    mol = Molecule()
    mol.empty(natoms, numFrames=1)

    idx = 0
    center = box_size / 2

    for i in range(solute_atoms):
        mol.name[idx] = "CA"
        mol.resname[idx] = "ALA"
        mol.resid[idx] = i + 1
        mol.segid[idx] = "P"
        mol.chain[idx] = "A"
        mol.element[idx] = "C"
        mol.coords[idx, :, 0] = [center + i * 0.5, center, center]
        idx += 1

    resid = solute_atoms + 1
    n_per_dim = int(np.ceil(nwaters ** (1 / 3)))
    spacing = box_size / (n_per_dim + 1)
    water_count = 0

    for ix in range(n_per_dim):
        for iy in range(n_per_dim):
            for iz in range(n_per_dim):
                if water_count >= nwaters:
                    break
                ox = spacing * (ix + 0.5)
                oy = spacing * (iy + 0.5)
                oz = spacing * (iz + 0.5)

                for ai, (aname, aelement, dx, dy) in enumerate(
                    [("O", "O", 0, 0), ("H1", "H", 0.96, 0), ("H2", "H", 0, 0.96)]
                ):
                    mol.name[idx] = aname
                    mol.resname[idx] = "WAT"
                    mol.resid[idx] = resid
                    mol.segid[idx] = "W"
                    mol.chain[idx] = "W"
                    mol.element[idx] = aelement
                    mol.coords[idx, :, 0] = [ox + dx, oy + dy, oz]
                    idx += 1

                resid += 1
                water_count += 1
            if water_count >= nwaters:
                break
        if water_count >= nwaters:
            break

    return mol
