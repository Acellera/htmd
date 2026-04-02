import numpy as np
import pytest
import shutil
import os
from glob import glob

from moleculekit.molecule import Molecule
from moleculekit.tools.graphalignment import makeMolGraph, compareGraphs
from htmd.builder.noncanonical import (
    _extend_residue,
    parameterizeNonCanonicalResidues,
    triala,
    triala_comp,
    start_map,
    end_map,
)
from htmd.home import home


DATA_DIR = os.path.dirname(os.path.abspath(__file__))
CIF_33X = os.path.join(DATA_DIR, "test-custom-residue-param", "33X.cif")

try:
    from parameterize.main import parameterize_molecule  # noqa: F401

    _PARAMETERIZE_INSTALLED = True
except ImportError:
    _PARAMETERIZE_INSTALLED = False


# ---------------------------------------------------------------------------
# Helpers for comparing parameterization output files
# ---------------------------------------------------------------------------


def _compare_frcmod(
    refFile,
    resFile,
    dihedralForceConstAbsTol=1e-6,
    dihedralPhaseAbsTol=1e-6,
):
    import difflib

    defaultAbsTol = 1e-6
    defaultRelTol = 1e-6

    with open(refFile) as ref, open(resFile) as res:
        refLines, resLines = ref.readlines(), res.readlines()

    # Removes the first line with the HTMD version
    refLines = refLines[1:]
    resLines = resLines[1:]

    if len(refLines) != len(resLines):
        raise AssertionError(
            f"Reference file {refFile} has a different number of lines ({len(refLines)}) "
            f"from file {resFile} ({len(resLines)})"
        )

    for refLine, resLine in zip(refLines, resLines):
        refFields = refLine.split()
        resFields = resLine.split()

        for iField, (refField, resField) in enumerate(zip(refFields, resFields)):
            if len(refFields) == 7 and iField == 2:
                absTol = dihedralForceConstAbsTol
                relTol = 0
            elif len(refFields) == 7 and iField == 3:
                absTol = dihedralPhaseAbsTol
                relTol = 0
            else:
                absTol = defaultAbsTol
                relTol = defaultRelTol

            try:
                refField = float(refField)
                resField = float(resField)
            except ValueError:
                if refField == resField:
                    continue
            else:
                if np.isclose(refField, resField, atol=absTol, rtol=relTol):
                    continue

            print(f"Failed: {refField} == {resField}")
            print(f"Absolute tolerance: {absTol}")
            print(f"Relative tolerance: {relTol}")

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


# ---------------------------------------------------------------------------
# Triala component ordering
# ---------------------------------------------------------------------------


class _TestTrialaComponentOrdering:
    """Verify that the module-level triala decomposition selects the right fragments."""

    def _test_three_components(self):
        assert len(triala_comp) == 3

    def _test_start_component_is_ace_ala(self):
        resnames = set(triala.resname[triala_comp[0]])
        assert resnames == {"ACE", "ALA"}
        assert min(triala_comp[0]) == 0, "Start component should contain atom 0"

    def _test_end_component_is_ala_nme(self):
        resnames = set(triala.resname[triala_comp[2]])
        assert resnames == {"ALA", "NME"}
        assert max(triala_comp[2]) == triala.numAtoms - 1

    def _test_components_are_sorted_by_index(self):
        for i in range(len(triala_comp) - 1):
            assert min(triala_comp[i]) < min(triala_comp[i + 1])

    def _test_start_end_maps_consistent(self):
        assert start_map.sum() == len(triala_comp[0])
        assert end_map.sum() == len(triala_comp[2])
        assert not np.any(start_map & end_map), "start and end maps must not overlap"


# ---------------------------------------------------------------------------
# _extend_residue
# ---------------------------------------------------------------------------


class _TestExtendResidue:
    """Test that _extend_residue correctly flanks a residue with ALA caps."""

    @pytest.fixture
    def ncaa_mol(self):
        return Molecule(CIF_33X)

    def _test_both_caps(self, ncaa_mol):
        xmol = _extend_residue(ncaa_mol, nterm=True, cterm=True)

        resnames = list(np.unique(xmol.resname))
        assert "ACE" in resnames, "N-terminal ACE cap missing"
        assert "NME" in resnames, "C-terminal NME cap missing"
        assert "ALA" in resnames, "Flanking ALA missing"
        assert ncaa_mol.resname[0] in resnames, "Original residue missing"

        assert xmol.numAtoms > ncaa_mol.numAtoms

    def _test_nterm_only(self, ncaa_mol):
        xmol = _extend_residue(ncaa_mol, nterm=True, cterm=False)

        resnames = set(np.unique(xmol.resname))
        assert "ACE" in resnames
        assert "NME" not in resnames

    def _test_cterm_only(self, ncaa_mol):
        xmol = _extend_residue(ncaa_mol, nterm=False, cterm=True)

        resnames = set(np.unique(xmol.resname))
        assert "ACE" not in resnames
        assert "NME" in resnames

    def _test_no_caps(self, ncaa_mol):
        xmol = _extend_residue(ncaa_mol, nterm=False, cterm=False)
        assert xmol.numAtoms == ncaa_mol.numAtoms

    def _test_residue_resid_is_three(self, ncaa_mol):
        resn = ncaa_mol.resname[0]
        xmol = _extend_residue(ncaa_mol, nterm=True, cterm=True)
        assert np.all(xmol.resid[xmol.resname == resn] == 3)

    def _test_resid_ordering(self, ncaa_mol):
        xmol = _extend_residue(ncaa_mol, nterm=True, cterm=True)
        resids = xmol.resid
        assert np.all(np.diff(resids) >= 0), "Residue IDs should be non-decreasing"

    def _test_backbone_connectivity(self, ncaa_mol):
        """The extended molecule should form a single connected graph."""
        import networkx as nx

        xmol = _extend_residue(ncaa_mol, nterm=True, cterm=True)
        g = xmol.toGraph()

        assert nx.is_connected(g), "Extended molecule graph should be fully connected"


# ---------------------------------------------------------------------------
# Padding atoms indexing
# ---------------------------------------------------------------------------


class _TestPaddingAtomsIndexing:
    """Test that padding_atoms is correctly identified via graph matching.

    Simulates the scenario in _post_process_parameterize where the typed
    molecule (from antechamber/parameterize) has lost its original residue
    names and must be re-mapped from a reference molecule.
    """

    @pytest.fixture
    def extended_pair(self):
        """Return (refmol, typed_mol, resn) where typed_mol has flattened resnames."""
        mol = Molecule(CIF_33X)
        resn = mol.resname[0]
        refmol = _extend_residue(mol, nterm=True, cterm=True)

        typed_mol = refmol.copy()
        typed_mol.resname[:] = resn
        typed_mol.resid[:] = 1
        return refmol, typed_mol, resn

    def _run_matching_and_get_padding(self, refmol, typed_mol, resn, use_mol_resname):
        """Run graph matching and compute padding_atoms.

        If use_mol_resname is True, uses mol.resname (the fix).
        If False, uses refmol.resname (the current code).
        """
        fields = ("element",)
        g1 = makeMolGraph(refmol, "all", fields)
        g2 = makeMolGraph(typed_mol, "all", fields)
        _, _, matching = compareGraphs(
            g1, g2, fields=fields, tolerance=0.5, returnmatching=True
        )

        for pp in matching:
            typed_mol.resname[pp[0]] = refmol.resname[pp[1]]

        padding = np.zeros(typed_mol.numAtoms, dtype=bool)
        if use_mol_resname:
            padding[typed_mol.resname != resn] = True
        else:
            padding[refmol.resname != resn] = True
        return padding

    def _test_same_order_both_methods_agree(self, extended_pair):
        """When atom ordering is preserved, both methods give the same result."""
        refmol, typed_mol, resn = extended_pair
        typed_copy = typed_mol.copy()

        pad_old = self._run_matching_and_get_padding(
            refmol, typed_mol, resn, use_mol_resname=False
        )
        typed_mol2 = typed_copy.copy()
        pad_new = self._run_matching_and_get_padding(
            refmol, typed_mol2, resn, use_mol_resname=True
        )

        np.testing.assert_array_equal(pad_old, pad_new)

    def _test_padding_selects_correct_atoms(self, extended_pair):
        refmol, typed_mol, resn = extended_pair

        padding = self._run_matching_and_get_padding(
            refmol, typed_mol, resn, use_mol_resname=True
        )

        expected_padding = refmol.resname != resn
        np.testing.assert_array_equal(padding, expected_padding)

    def _test_padding_count(self, extended_pair):
        refmol, typed_mol, resn = extended_pair

        padding = self._run_matching_and_get_padding(
            refmol, typed_mol, resn, use_mol_resname=True
        )

        n_padding = padding.sum()
        n_ncaa = (~padding).sum()
        ncaa_input = Molecule(CIF_33X)
        expected_ncaa = ncaa_input.numAtoms - 1 - 2  # H removed, OXT+HXT removed
        assert n_ncaa == expected_ncaa
        assert n_padding > 0

    def _test_shuffled_atoms_mol_resname_correct(self, extended_pair):
        """When atoms are reordered, using mol.resname (after matching) is still correct."""
        refmol, typed_mol, resn = extended_pair

        perm = np.random.RandomState(42).permutation(typed_mol.numAtoms)
        inv_perm = np.argsort(perm)

        shuffled = typed_mol.copy()
        for field in shuffled._atom_and_coord_fields:
            arr = getattr(shuffled, field)
            if arr is not None and np.size(arr) > 0:
                if field == "coords":
                    shuffled.__dict__[field] = arr[perm, :, :]
                else:
                    shuffled.__dict__[field] = arr[perm]
        if shuffled.bonds is not None and len(shuffled.bonds):
            shuffled.bonds = inv_perm[shuffled.bonds]

        fields = ("element",)
        g1 = makeMolGraph(refmol, "all", fields)
        g2 = makeMolGraph(shuffled, "all", fields)
        _, _, matching = compareGraphs(
            g1, g2, fields=fields, tolerance=0.5, returnmatching=True
        )

        for pp in matching:
            shuffled.resname[pp[0]] = refmol.resname[pp[1]]

        pad_via_mol = np.zeros(shuffled.numAtoms, dtype=bool)
        pad_via_mol[shuffled.resname != resn] = True

        n_expected_padding = (refmol.resname != resn).sum()
        assert pad_via_mol.sum() == n_expected_padding

        padding_resnames = set(shuffled.resname[pad_via_mol])
        assert (
            resn not in padding_resnames
        ), "Padding should not include the NCAA residue"
        non_padding_resnames = set(shuffled.resname[~pad_via_mol])
        assert non_padding_resnames == {
            resn
        }, "Non-padding should only contain the NCAA residue"


# ---------------------------------------------------------------------------
# End-to-end parameterization tests (native antechamber backend)
# ---------------------------------------------------------------------------

_antechamber_installed = shutil.which("antechamber", mode=os.X_OK) is not None


NCAA_CIF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "noncanonical")


@pytest.mark.skipif(
    not _antechamber_installed, reason="Requires native antechamber (AmberTools)"
)
class _TestNCAAResParamNativeAntechamber:
    def _test_ncaa_native_antechamber_produces_output(self, tmp_path):
        """Smoke test: native antechamber backend produces frcmod and prepi."""

        cif = os.path.join(DATA_DIR, "test-custom-residue-param", "33X.cif")
        parameterizeNonCanonicalResidues(
            cif,
            tmp_path,
            forcefield="GAFF2",
            charge_model="AM1-BCC",
            backend="antechamber_native",
        )
        frcmods = glob(os.path.join(tmp_path, "*.frcmod"))
        prepis = glob(os.path.join(tmp_path, "*.prepi"))
        assert len(frcmods) >= 1, "No frcmod file produced"
        assert len(prepis) == 1, "Expected exactly one prepi file"
        assert os.path.getsize(prepis[0]) > 0, "prepi file is empty"
        assert os.path.getsize(frcmods[0]) > 0, "frcmod file is empty"

    @pytest.mark.parametrize(
        "cif_file",
        sorted(glob(os.path.join(NCAA_CIF_DIR, "*.cif"))),
        ids=lambda p: os.path.splitext(os.path.basename(p))[0],
    )
    def _test_parameterize_ncaa_cifs(self, cif_file, tmp_path):
        """Parameterize each CIF in tests/noncanonical/ with native antechamber."""
        resn = os.path.splitext(os.path.basename(cif_file))[0]
        outdir = str(tmp_path / resn)
        parameterizeNonCanonicalResidues(
            cif_file,
            outdir,
            forcefield="GAFF2",
            charge_model="Gasteiger",
            backend="antechamber_native",
        )
        frcmod = os.path.join(outdir, f"{resn}.frcmod")
        prepi = os.path.join(outdir, f"{resn}.prepi")
        assert os.path.exists(frcmod), f"Missing {resn}.frcmod"
        assert os.path.exists(prepi), f"Missing {resn}.prepi"

        ref_frcmod = os.path.join(NCAA_CIF_DIR, f"{resn}.frcmod")
        ref_prepi = os.path.join(NCAA_CIF_DIR, f"{resn}.prepi")
        _compare_frcmod(ref_frcmod, frcmod)
        _compare_prepis(ref_prepi, prepi)


# ---------------------------------------------------------------------------
# End-to-end parameterization tests (require parameterize backend)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _PARAMETERIZE_INSTALLED, reason="Can only run with parameterize installed"
)
class _TestNCAAResParam:
    @pytest.mark.parametrize(
        "cif_file",
        sorted(glob(os.path.join(DATA_DIR, "test-custom-residue-param", "*.cif"))),
        ids=lambda p: os.path.splitext(os.path.basename(p))[0],
    )
    def _test_ncaa_residue_parameterization(self, cif_file, tmp_path):
        refdir = os.path.join(DATA_DIR, "test-custom-residue-param")
        refresdir = os.path.join(refdir, "gaff2-params")
        parameterizeNonCanonicalResidues(cif_file, tmp_path, forcefield="GAFF2")

        frcmod = glob(os.path.join(tmp_path, "*.frcmod"))[0]
        name = os.path.basename(frcmod)
        _compare_frcmod(os.path.join(refresdir, name), frcmod)

        prepi = glob(os.path.join(tmp_path, "*.prepi"))[0]
        name = os.path.basename(prepi)
        _compare_prepis(os.path.join(refresdir, name), prepi)

    @pytest.mark.parametrize(
        "cif_file",
        sorted(glob(os.path.join(DATA_DIR, "test-custom-residue-param", "*.cif"))),
        ids=lambda p: os.path.splitext(os.path.basename(p))[0],
    )
    def _test_ncaa_residue_parameterization_terminals(self, cif_file, tmp_path):
        refdir = os.path.join(DATA_DIR, "test-custom-residue-param")
        refresdir = os.path.join(refdir, "gaff2-params-terminals")
        parameterizeNonCanonicalResidues(
            cif_file,
            tmp_path,
            forcefield="GAFF2",
            is_cterm=True,
            is_nterm=True,
        )

        frcmod = glob(os.path.join(tmp_path, "*.frcmod"))[0]
        name = os.path.basename(frcmod)
        _compare_frcmod(os.path.join(refresdir, name), frcmod)

        prepi = glob(os.path.join(tmp_path, "*.prepi"))[0]
        name = os.path.basename(prepi)
        _compare_prepis(os.path.join(refresdir, name), prepi)

    # def _test_ncaa_residue_parameterization_xTB(self, tmp_path):
    #     refdir = os.path.join(DATA_DIR, "test-custom-residue-param")
    #     refresdir = os.path.join(refdir, "xtb-params")
    #     cif = os.path.join(refdir, "33X.cif")
    #     parameterizeNonCanonicalResidues(
    #         cif, tmp_path, forcefield="GAFF2", calculator="xTB"
    #     )
    #     frcmod = glob(os.path.join(tmp_path, "*.frcmod"))[0]
    #     name = os.path.basename(frcmod)
    #     _compare_frcmod(os.path.join(refresdir, name), frcmod)

    #     prepi = glob(os.path.join(tmp_path, "*.prepi"))[0]
    #     name = os.path.basename(prepi)
    #     _compare_prepis(os.path.join(refresdir, name), prepi)
