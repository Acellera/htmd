import numpy as np
import pytest
import scipy.spatial.distance as dist
from moleculekit.molecule import Molecule
from htmd.builder.ionize import (
    ionize,
    ionizePlace,
    ionizeSalts,
    placeIons,
    resolveSalts,
    _ionGetCharge,
)


class TestIonGetCharge:
    def test_monovalent_cations(self):
        assert _ionGetCharge("NA") == 1
        assert _ionGetCharge("K") == 1
        assert _ionGetCharge("CS") == 1

    def test_divalent_cations(self):
        assert _ionGetCharge("MG") == 2
        assert _ionGetCharge("CA") == 2
        assert _ionGetCharge("ZN") == 2

    def test_anion(self):
        assert _ionGetCharge("CL") == -1

    def test_halides_and_aliases(self):
        for alias in ("F", "fluoride", "BR", "Br-", "IOD", "I", "iodide", "I-"):
            assert _ionGetCharge(alias) == -1

    def test_unknown_ion_raises(self):
        with pytest.raises(NameError, match="not recognized"):
            _ionGetCharge("XX")


class TestIonize:
    def test_neutralize_positive_charge_monovalent(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=5, nwater=1000, neutralize=True
        )
        assert nanion == 5
        assert ncation == 0
        assert anion == "CL"
        assert cation == "NA"

    def test_neutralize_negative_charge_monovalent(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-3, nwater=1000, neutralize=True
        )
        assert nanion == 0
        assert ncation == 3

    def test_neutral_system(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=0, nwater=1000, neutralize=True
        )
        assert nanion == 0
        assert ncation == 0

    def test_neutralize_positive_charge_divalent_anion(self):
        # Even charge: straightforward
        anion, cation, nanion, ncation = ionize(
            None, netcharge=4, nwater=1000, neutralize=True, anion="CL"
        )
        assert nanion * (-1) + ncation * 1 + 4 == 0

    def test_neutralize_negative_charge_divalent_cation(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-4, nwater=1000, neutralize=True, cation="MG"
        )
        assert nanion * (-1) + ncation * 2 + (-4) == 0

    def test_neutralize_odd_charge_divalent_cation(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=-3, nwater=1000, neutralize=True, cation="MG"
        )
        assert nanion * (-1) + ncation * 2 + (-3) == 0
        assert nanion >= 1  # Need at least one anion to balance

    def test_saltconc_returns_six_values(self):
        mol = _make_mol_with_atoms(0)
        result = ionize(mol, netcharge=0, nwater=5000, saltconc=0.15)
        assert len(result) == 6

    def test_saltconc_adds_ions(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=0, nwater=5000, saltconc=0.15
        )
        # num = floor(0.5 + 0.0187 * 0.15 * 5000) = floor(14.525) = 14
        assert nanion == 14
        assert ncation == 14

    def test_saltconc_zero_neutralize_only(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=2, nwater=5000, saltconc=0
        )
        assert nanion == 2
        assert ncation == 0

    def test_saltconc_plus_neutralization(self):
        mol = _make_mol_with_atoms(0)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=3, nwater=5000, saltconc=0.15
        )
        # Neutralization: 3 Cl-
        # Salt: 14 NaCl pairs
        assert nanion == 3 + 14
        assert ncation == 0 + 14

    def test_unknown_cation_raises(self):
        with pytest.raises(NameError):
            ionize(None, netcharge=1, nwater=1000, neutralize=True, cation="XX")

    def test_unknown_anion_raises(self):
        with pytest.raises(NameError):
            ionize(None, netcharge=1, nwater=1000, neutralize=True, anion="XX")

    def test_rounding_of_charge(self):
        anion, cation, nanion, ncation = ionize(
            None, netcharge=2.7, nwater=1000, neutralize=True
        )
        assert nanion == 3
        assert ncation == 0


class TestIonizePlace:
    def test_zero_ions_returns_copy(self):
        solvent = _make_water_box(nwaters=50)
        result = ionizePlace(solvent, None, "Cl-", "Na+", "Cl-", "Na+", 0, 0)
        assert result.numAtoms == solvent.numAtoms

    def test_does_not_modify_input(self):
        solvent, solute = _make_water_and_solute(nwaters=50, solute_atoms=5)
        natoms_before = solvent.numAtoms
        ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert solvent.numAtoms == natoms_before

    def test_correct_atom_count(self):
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        nanion, ncation = 3, 2
        nions = nanion + ncation
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", nanion, ncation)
        expected = solvent.numAtoms - nions * 3 + nions
        assert result.numAtoms == expected

    def test_ion_resnames(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 3)
        ion_mask = result.segid == "I"
        assert ion_mask.sum() == 5
        assert np.sum(result.resname[ion_mask] == "Cl-") == 2
        assert np.sum(result.resname[ion_mask] == "Na+") == 3

    def test_ion_names(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "CLA", "SOD", 1, 1)
        ion_mask = result.segid == "I"
        assert "CLA" in result.name[ion_mask]
        assert "SOD" in result.name[ion_mask]

    def test_ion_chain_and_segid(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        ion_mask = result.segid == "I"
        assert np.all(result.chain[ion_mask] == "I")
        assert np.all(result.segid[ion_mask] == "I")

    def test_ions_far_from_solute(self):
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

    def test_ions_far_from_each_other(self):
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

    def test_unique_resids_for_ions(self):
        np.random.seed(42)
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 3, 3)
        ion_resids = result.resid[result.segid == "I"]
        assert len(np.unique(ion_resids)) == 6

    def test_no_eligible_water_raises(self):
        solvent, solute = _make_water_and_solute(nwaters=3, solute_atoms=5, box_size=5)
        with pytest.raises(NameError, match="No waters could be found"):
            ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 10, 10, dfrom=100)

    def test_water_count_decreases(self):
        solvent, solute = _make_water_and_solute(nwaters=100, solute_atoms=5)
        nwat_before = solvent.atomselect("water and noh").sum()
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 3)
        nwat_after = result.atomselect("water and noh").sum()
        assert nwat_after == nwat_before - 5

    def test_no_solute(self):
        solvent = _make_water_box(nwaters=100)
        result = ionizePlace(solvent, None, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert result.numAtoms == solvent.numAtoms - 4 * 3 + 4

    def test_empty_solute(self):
        solvent = _make_water_box(nwaters=100)
        solute = Molecule()
        solute.empty(0)
        result = ionizePlace(solvent, solute, "Cl-", "Na+", "Cl-", "Na+", 2, 2)
        assert result.numAtoms == solvent.numAtoms - 4 * 3 + 4

    def test_ion_types_spatially_mixed(self):
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

class TestResolveSalts:
    def test_shorthand_defaults_to_nacl(self):
        assert resolveSalts(None, 0.15) == [("NA", "CL", 0.15)]

    def test_shorthand_aliases(self):
        assert resolveSalts(None, 0.1, "POT", "chloride") == [("K", "CL", 0.1)]

    def test_shorthand_ions_are_deprecated(self, caplog):
        with caplog.at_level("WARNING"):
            resolveSalts(None, 0.15)
        assert "deprecated" not in caplog.text
        with caplog.at_level("WARNING"):
            resolveSalts(None, 0.15, "K")
        assert "deprecated" in caplog.text

    def test_shorthand_swapped_ions_raise(self):
        with pytest.raises(ValueError, match="not a cation"):
            resolveSalts(None, 0.1, "CL", "K")

    def test_formulas_in_order_any_case(self):
        salts = resolveSalts([("KCl", 0.15), ("mgcl2", 0.001)])
        assert salts == [("K", "CL", 0.15), ("MG", "CL", 0.001)]

    def test_every_chloride(self):
        formulas = ["NaCl", "KCl", "CsCl", "MgCl2", "CaCl2", "ZnCl2"]
        cations = [c for c, _, _ in resolveSalts([(f, 0.1) for f in formulas])]
        assert cations == ["NA", "K", "CS", "MG", "CA", "ZN"]

    def test_halide_salts(self):
        salts = resolveSalts([("KBr", 0.1), ("NaI", 0.1), ("CaF2", 0.1)])
        assert [s[:2] for s in salts] == [("K", "BR"), ("NA", "IOD"), ("CA", "F")]

    def test_unknown_salt_lists_supported(self):
        with pytest.raises(ValueError, match="Supported: .*MgCl2"):
            resolveSalts([("Na2SO4", 0.15)])

    def test_salt_given_twice_raises(self):
        with pytest.raises(ValueError, match="twice"):
            resolveSalts([("KCl", 0.15), ("kcl", 0.1)])

    def test_list_and_shorthand_raises(self):
        with pytest.raises(ValueError, match="not both"):
            resolveSalts([("KCl", 0.15)], saltconc=0.15)

    def test_negative_concentration_raises(self):
        with pytest.raises(ValueError, match="negative"):
            resolveSalts([("KCl", -0.1)])


class TestIonizeSalts:
    def test_single_salt_matches_ionize(self):
        mol = _make_mol_with_atoms(10)
        _, _, _, _, nanion, ncation = ionize(
            mol, netcharge=-3, nwater=5000, saltconc=0.15, cation="K"
        )
        counts = ionizeSalts(mol.resname, -3, 5000, [("K", "CL", 0.15)])
        assert counts == {"K": ncation, "CL": nanion}

    def test_first_salt_neutralizes(self):
        mol = _make_mol_with_atoms(10)
        # 0.15 M KCl in 5000 waters is 14 units, 0.02 M MgCl2 is 2 units
        counts = ionizeSalts(
            mol.resname, -4, 5000, [("K", "CL", 0.15), ("MG", "CL", 0.02)]
        )
        assert counts == {"K": 4 + 14, "CL": 14 + 2 * 2, "MG": 2}
        charge = sum(_ionGetCharge(ion) * n for ion, n in counts.items())
        assert charge == 4

    def test_divalent_first_salt_neutralizes(self):
        mol = _make_mol_with_atoms(10)
        counts = ionizeSalts(mol.resname, -4, 5000, [("MG", "CL", 0), ("K", "CL", 0)])
        assert counts == {"MG": 2, "CL": 0, "K": 0}

    def test_divalent_cation_halide(self):
        mol = _make_mol_with_atoms(10)
        counts = ionizeSalts(mol.resname, 2, 5000, [("MG", "BR", 0.02)])
        # 2 Br- neutralize, then 2 units of MgBr2
        assert counts == {"MG": 2, "BR": 2 + 4}

    def test_salt_below_one_unit_warns(self, caplog):
        mol = _make_mol_with_atoms(10)
        with caplog.at_level("WARNING"):
            counts = ionizeSalts(
                mol.resname, 0, 5000, [("K", "CL", 0.15), ("MG", "CL", 0.001)]
            )
        assert counts["MG"] == 0
        assert "less than one formula unit" in caplog.text


class TestPlaceIons:
    def test_places_every_type(self):
        solvent, solute = _make_water_and_solute(nwaters=400, solute_atoms=20, box_size=40)
        newmol = placeIons(solvent, solute, {"K": 5, "CL": 7, "MG": 1, "IOD": 2})
        for ion, name, n in [("K", "K", 5), ("CL", "CL", 7), ("MG", "MG", 1), ("IOD", "I", 2)]:
            assert np.sum(newmol.resname == ion) == n
            assert np.sum((newmol.resname == ion) & (newmol.name == name)) == n
        assert np.sum(newmol.atomselect("water and noh")) == 400 - 15

    def test_single_salt_matches_ionizeplace(self):
        solvent, solute = _make_water_and_solute(nwaters=400, solute_atoms=20, box_size=40)
        np.random.seed(3)
        old = ionizePlace(solvent, solute, "CL", "NA", "CL", "NA", 6, 4)
        np.random.seed(3)
        new = placeIons(solvent, solute, {"NA": 4, "CL": 6})
        assert (old.resname == new.resname).all()
        assert np.allclose(old.coords, new.coords)


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
