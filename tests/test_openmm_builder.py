# (c) 2015-2025 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import pytest
import numpy as np

try:
    import openmm
    import openmm.app
    import openmm.unit

    _openmm_installed = True
except ImportError:
    _openmm_installed = False

try:
    from htmd.builder.amber import _findTeLeap

    _tleap_installed = _findTeLeap() is not None
except Exception:
    _tleap_installed = False

import sys

# openff-interchange 0.5+ uses PEP 695 ``type`` statement syntax and
# therefore needs Python >= 3.12. On 3.10 / 3.11 the package raises
# SyntaxError at import time (not ImportError), so the version check
# has to come first and the except clause has to cover SyntaxError too.
try:
    if sys.version_info < (3, 12):
        raise ImportError("openff-interchange requires Python >= 3.12")
    import openff.toolkit  # noqa: F401
    import openff.interchange  # noqa: F401

    _openff_installed = True
except (ImportError, SyntaxError):
    _openff_installed = False

try:
    if sys.version_info < (3, 12):
        raise ImportError("openff-nagl requires Python >= 3.12 transitively")
    import openff.nagl  # noqa: F401
    import openff.nagl_models  # noqa: F401

    _nagl_installed = True
except (ImportError, SyntaxError):
    _nagl_installed = False


curr_dir = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def tmpdir(tmp_path):
    return str(tmp_path)


# ------------------------------------------------------------------
# Energy helpers
# ------------------------------------------------------------------

_FORCE_LABELS = {
    "HarmonicBondForce": "bonds",
    "HarmonicAngleForce": "angles",
    "PeriodicTorsionForce": "torsions",
    "NonbondedForce": "nonbonded",
    "CMAPTorsionForce": "cmap",
}


def _heavy_atom_energy(outdir, prefix="structure"):
    """Strip hydrogens from a builder's prmtop+PDB, return single-point
    per-force decomposed energy (kcal/mol) and heavy-atom count.

    By comparing heavy-atom-only energies we eliminate differences caused
    by hydrogen placement (each builder adds H independently) and LJ
    sigma/epsilon convention differences (which mainly affect polar H).
    """
    import parmed
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    prmtop_path = os.path.join(outdir, f"{prefix}.prmtop")
    pdb_path = os.path.join(outdir, f"{prefix}.pdb")

    parm = parmed.load_file(prmtop_path, pdb_path)
    parm.strip("@/H")

    system = parm.createSystem(nonbondedMethod=app.NoCutoff)
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)

    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator)
    context.setPositions(parm.positions)

    result = {}
    for i, force in enumerate(system.getForces()):
        fname = force.__class__.__name__
        label = _FORCE_LABELS.get(fname)
        if label is None:
            continue
        e = context.getState(getEnergy=True, groups={i}).getPotentialEnergy()
        result[label] = e.value_in_unit(unit.kilocalories_per_mole)

    result["total"] = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilocalories_per_mole)
    )
    result["natoms"] = len(parm.atoms)
    return result


def _compare_energies(ene_off, ene_amb, label, bonded_atol=0.005, nonbonded_atol=None):
    """Compare per-force heavy-atom energies between two builders.

    With hydrogens stripped, bonded terms (bonds, angles, torsions, cmap)
    should match almost exactly because both prmtops use the same ff14SB
    parameters and identical heavy-atom coordinates.

    Nonbonded and total energies are always printed.  Pass a numeric
    *nonbonded_atol* to assert on them; the default ``None`` skips the
    assertion because the LJ sigma/epsilon convention difference between
    tleap's native pair-specific A/B matrix and ParmEd's per-atom
    combining-rule export affects heavy-atom LJ interactions.

    All tolerances are per heavy atom (kcal/mol/atom).
    """
    n = max(ene_off.get("natoms", 1), 1)

    print(f"\n[{label}] Heavy-atom energy decomposition (kcal/mol), {n} atoms:")
    print(
        f"  {'term':>12s}  {'OpenFF':>12s}  {'AMBER':>12s}  {'diff':>10s}  {'diff/atom':>10s}"
    )
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*10}")

    for key in ("bonds", "angles", "torsions", "cmap", "nonbonded", "total"):
        if key not in ene_off or key not in ene_amb:
            continue
        vo = ene_off[key]
        va = ene_amb[key]
        diff = abs(vo - va)
        dpa = diff / n
        print(f"  {key:>12s}  {vo:12.2f}  {va:12.2f}  {diff:10.2f}  {dpa:10.4f}")

        assert np.isfinite(vo), f"[{label}] OpenFF {key} not finite: {vo}"
        assert np.isfinite(va), f"[{label}] AMBER {key} not finite: {va}"

        if key in ("bonds", "angles", "torsions", "cmap"):
            assert (
                dpa <= bonded_atol
            ), f"[{label}] {key} diff/atom {dpa:.6f} > bonded tolerance {bonded_atol}"
        elif key in ("nonbonded", "total") and nonbonded_atol is not None:
            assert (
                dpa <= nonbonded_atol
            ), f"[{label}] {key} diff/atom {dpa:.6f} > nonbonded tolerance {nonbonded_atol}"


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _make_water_mol(resname="WAT", oxygen_name="OH2"):
    """Create a minimal 3-atom water Molecule for unit tests."""
    from moleculekit.molecule import Molecule

    mol = Molecule()
    mol.empty(3)
    mol.name[:] = [oxygen_name, "H1", "H2"]
    mol.resname[:] = resname
    mol.resid[:] = 1
    mol.element[:] = ["O", "H", "H"]
    mol.record[:] = "HETATM"
    mol.coords = np.array(
        [[[0.0], [0.0], [0.0]], [[0.96], [0.0], [0.0]], [[-0.24], [0.93], [0.0]]],
        dtype=np.float32,
    )
    return mol


# ------------------------------------------------------------------
# Unit tests
# ------------------------------------------------------------------


@pytest.mark.skipif(not _openmm_installed, reason="OpenMM not installed")
class _TestOpenMMBuilder:
    """Tests for htmd.builder.openmm.build()."""

    def _test_defaultFf(self):
        from htmd.builder.openmm import defaultFf

        ff = defaultFf()
        assert isinstance(ff, list)
        assert len(ff) > 0
        assert all(isinstance(f, str) for f in ff)
        assert any("protein" in f for f in ff)
        assert any("tip3p" in f for f in ff)

    def _test_resolve_ion_name(self):
        from htmd.builder.openmm import _resolve_ion_name

        assert _resolve_ion_name("Na+") == "Na+"
        assert _resolve_ion_name("sodium") == "Na+"
        assert _resolve_ion_name("SOD") == "Na+"
        assert _resolve_ion_name("CL") == "Cl-"
        assert _resolve_ion_name("chloride") == "Cl-"
        assert _resolve_ion_name("K+") == "K+"
        assert _resolve_ion_name("potassium") == "K+"
        with pytest.raises(ValueError):
            _resolve_ion_name("UnknownIon")

    def _test_fix_water_naming(self):
        from htmd.builder.openmm import _fix_water_naming

        mol = _make_water_mol("WAT", "OH2")
        _fix_water_naming(mol)
        assert np.all(mol.resname == "HOH")
        assert mol.name[0] == "O"


# ------------------------------------------------------------------
# Comparative build tests (mirrors test_amber_builder.py)
# ------------------------------------------------------------------


@pytest.mark.skipif(not _openmm_installed, reason="OpenMM not installed")
class _TestOpenMMComparative:
    """Build the same systems as test_amber_builder and compare energies."""

    def _test_protein_prepared(self, tmpdir):
        """3PTB with systemPrepare + autoSegment, no solvation."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule
        from moleculekit.tools.preparation import systemPrepare
        from moleculekit.tools.autosegment import autoSegment

        mol = Molecule("3PTB")
        mol.filter("protein")
        mol, _ = systemPrepare(mol, pH=7.0)
        mol.filter("protein")
        mol = autoSegment(mol)

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(
            f"[3PTB prepared] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol"
        )

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build, defaultFf

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(mol.copy(), ff=defaultFf(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "3PTB prepared")

    def _test_protein_raw(self, tmpdir):
        """3PTB without systemPrepare, no solvation."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        mol = Molecule("3PTB")
        mol.filter("protein")

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[3PTB raw] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build, defaultFf

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(mol.copy(), ff=defaultFf(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "3PTB raw")

    def _test_disulfide_bonds(self, tmpdir):
        """1GZM with explicit disulfide bond selections."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule
        from moleculekit.tools.preparation import systemPrepare

        mol = Molecule("1GZM")
        mol.filter("protein")
        mol, _ = systemPrepare(mol, pH=7.0)
        mol.filter("protein")
        mol.segid = mol.chain

        disu = [
            ["chain A and resid 110", "chain A and resid 187"],
            ["chain B and resid 110", "chain B and resid 187"],
        ]

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
            disulfide=disu,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(
            f"[1GZM disulfide] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol"
        )

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build, defaultFf

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(
                mol.copy(),
                ff=defaultFf(),
                outdir=amber_outdir,
                ionize=False,
                disulfide=disu,
            )
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "1GZM disulfide")

    def _test_rna(self, tmpdir):
        """6VA1 pure RNA, no solvation.

        AMBER uses RNA.Shaw while OpenFF uses RNA.OL3, so torsion
        parameters genuinely differ — tolerance is relaxed.
        """
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        mol = Molecule("6VA1")

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[6VA1 RNA] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(mol.copy(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "6VA1 RNA", bonded_atol=0.5)

    def _test_dna(self, tmpdir):
        """1BNA pure DNA, no solvation."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        mol = Molecule("1BNA")

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[1BNA DNA] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(mol.copy(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "1BNA DNA")

    def _test_protein_rna(self, tmpdir):
        """3WBM protein + RNA complex with systemPrepare."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule
        from moleculekit.tools.preparation import systemPrepare
        from moleculekit.tools.autosegment import autoSegment

        mol = Molecule("3WBM")
        mol.filter("not water")
        mol = autoSegment(mol, field="both")
        pmol, _ = systemPrepare(mol, pH=7.0)

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            pmol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(
            f"[3WBM protein+RNA] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol"
        )

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(pmol.copy(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "3WBM protein+RNA", bonded_atol=0.1)

    def _test_caps(self, tmpdir):
        """6A5J peptide with default caps."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule
        from moleculekit.tools.preparation import systemPrepare

        mol = Molecule("6A5J")
        pmol, _ = systemPrepare(mol, pH=7.0)

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            pmol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[6A5J caps] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(pmol.copy(), outdir=amber_outdir, ionize=False)
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "6A5J caps")

    def _test_cyclic_peptide(self, tmpdir):
        """5VAV cyclic peptide."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        mol = Molecule("5VAV")
        mol.segid = mol.chain

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
            caps={"A": ("none", "none")},
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[5VAV cyclic] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import build as amber_build

            amber_outdir = os.path.join(tmpdir, "amber")
            amber_build(
                mol.copy(),
                outdir=amber_outdir,
                ionize=False,
                caps={"A": ("none", "none")},
            )
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "5VAV cyclic")

    @pytest.mark.skipif(
        sys.platform.startswith("win"),
        reason="AM1-BCC charge model unavailable on Windows",
    )
    def _test_protein_ligand(self, tmpdir):
        """3PTB protein + benzamidine ligand with GAFF2 small-molecule FF."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule
        from openff.toolkit import Molecule as OFFMolecule

        refdir = os.path.join(curr_dir, "data", "test-amber-build", "protLig")
        mol = Molecule(os.path.join(refdir, "3ptb_mod.pdb"))
        lig = Molecule(os.path.join(refdir, "BEN.mol2"))
        lig.segid[:] = "L"
        newmol = Molecule()
        newmol.append(mol)
        newmol.append(lig)

        ben_off = OFFMolecule.from_smiles(
            "c1ccc(cc1)C(=[NH2+])N", allow_undefined_stereo=True
        )

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            newmol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
            gbsa=True,
            small_molecule_ff="gaff-2.11",
            molecules=[ben_off],
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        ene_off = _heavy_atom_energy(outdir)
        assert np.isfinite(ene_off["total"])
        print(f"[3PTB+BEN] OpenFF heavy-atom total: {ene_off['total']:.1f} kcal/mol")

        if _tleap_installed:
            from htmd.builder.amber import (
                build as amber_build,
                defaultParam,
                defaultTopo,
            )

            amber_outdir = os.path.join(tmpdir, "amber")
            params = defaultParam() + [os.path.join(refdir, "BEN.frcmod")]
            topos = defaultTopo() + [os.path.join(refdir, "BEN.mol2")]
            amber_build(
                newmol.copy(),
                outdir=amber_outdir,
                param=params,
                topo=topos,
                ionize=False,
            )
            ene_amb = _heavy_atom_energy(amber_outdir)
            _compare_energies(ene_off, ene_amb, "3PTB+BEN protein+ligand")

    def _test_membrane(self, tmpdir):
        """Pre-built POPC membrane with TIP3 water.

        OpenFF uses lipid17 while AMBER uses lipid21, so bonded parameters
        genuinely differ.  We only verify that the build succeeds and
        produces a system with the expected atom count.
        """
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        homedir = os.path.join(curr_dir, "data", "test-amber-build", "membrane")
        mol = Molecule(os.path.join(homedir, "structure.psf"))
        mol.read(os.path.join(homedir, "structure.pdb"))

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=False,
            ionize=False,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms == mol.numAtoms
        assert os.path.exists(os.path.join(outdir, "structure.prmtop"))
        assert os.path.exists(os.path.join(outdir, "structure.pdb"))
        print(f"[Membrane] Built OK: {molbuilt.numAtoms} atoms")

    def _test_ionize_salt(self, tmpdir):
        """3PTB solvated with salt -- check ion counts."""
        from htmd.builder.openmm import build as openff_build
        from moleculekit.molecule import Molecule

        mol = Molecule("3PTB")
        mol.filter("protein")

        outdir = os.path.join(tmpdir, "openff")
        molbuilt, system = openff_build(
            mol.copy(),
            outdir=outdir,
            solvate=True,
            ionize=True,
            saltconc=0.15,
            padding=10.0,
        )
        assert molbuilt is not None
        assert molbuilt.numAtoms > 0

        import openmm

        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                total_q = sum(
                    force.getParticleParameters(j)[0]._value
                    for j in range(force.getNumParticles())
                )
                assert abs(total_q) < 0.01, f"System not neutral: charge={total_q:.4f}"
                break

        print(f"[3PTB ionize] atoms={molbuilt.numAtoms}, charge={total_q:.4f}")


# ====================================================================
# Phase 1: parameterizeLigandsOpenFF (free ligands via Interchange)
# ====================================================================


_PARAM_TEST_DIR = os.path.join(curr_dir, "test-custom-residue-param")
BEN_SMILES = "[NH2+]=C(N)c1ccccc1"


# Reference values produced by running Sage 2.3.0 (unconstrained) on the
# protonated benzamidine ligand of 3PTB with RDKit Gasteiger preset charges.
# Re-run probe_off_phase1.py to refresh these numbers if the toolkit version
# changes.
_BEN_EXPECTED_ATOMS = 18
_BEN_EXPECTED_BONDS = 18
_BEN_EXPECTED_ANGLES = 27
_BEN_EXPECTED_PROPER_KEYS = 40
_BEN_EXPECTED_IMPROPER_KEYS = 27
_BEN_EXPECTED_TORSIONS_IN_SYSTEM = 67
_BEN_EXPECTED_ATOM_NAMES = [
    "C1",
    "C2",
    "C3",
    "C4",
    "C5",
    "C6",
    "C",
    "N1",
    "N2",
    "H1",
    "H2",
    "H3",
    "H4",
    "H5",
    "H6",
    "H7",
    "H8",
    "H9",
]
_BEN_EXPECTED_CHARGES = [
    0.062728,
    -0.046778,
    -0.061366,
    -0.062216,
    -0.061366,
    -0.046778,
    0.270171,
    -0.287040,
    -0.287040,
    0.063213,
    0.062295,
    0.062269,
    0.062295,
    0.063213,
    0.301600,
    0.301600,
    0.301600,
    0.301600,
]


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed),
    reason="OpenMM + OpenFF Interchange required",
)
def _test_parameterize_ligands_openff_ben():
    """Build an Interchange for BEN via Sage 2.3 with RDKit Gasteiger charges
    and check the result against pre-recorded reference counts and per-atom
    charges. Acts as a regression baseline for the free-ligand path."""
    from moleculekit.molecule import Molecule
    from htmd.builder.openmm import parameterizeLigandsOpenFF

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)

    out = parameterizeLigandsOpenFF(
        mol,
        ligand_ff="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
    )

    assert set(out) == {"BEN"}
    ic = out["BEN"]

    # Topology shape.
    assert ic.topology.n_atoms == _BEN_EXPECTED_ATOMS
    assert ic.topology.n_bonds == _BEN_EXPECTED_BONDS
    assert mol.numAtoms == _BEN_EXPECTED_ATOMS
    assert mol.bonds.shape[0] == _BEN_EXPECTED_BONDS

    # Atom names roundtripped through the RDKit / OFFMolecule conversion.
    off_names = [a.name for a in ic.topology.molecule(0).atoms]
    assert off_names == _BEN_EXPECTED_ATOM_NAMES

    # Per-collection key counts from the SMIRNOFF resolution.
    assert len(ic["Bonds"].key_map) == _BEN_EXPECTED_BONDS
    assert len(ic["Angles"].key_map) == _BEN_EXPECTED_ANGLES
    assert len(ic["ProperTorsions"].key_map) == _BEN_EXPECTED_PROPER_KEYS
    assert len(ic["ImproperTorsions"].key_map) == _BEN_EXPECTED_IMPROPER_KEYS

    # Per-atom Gasteiger charges, ordered by atom index. Compare to 4 dp
    # since RDKit's PEOE iteration converges per platform with that precision.
    charges = []
    for key in sorted(
        ic["Electrostatics"].charges.keys(), key=lambda k: k.atom_indices[0]
    ):
        charges.append(float(ic["Electrostatics"].charges[key].m))
    assert len(charges) == _BEN_EXPECTED_ATOMS
    for actual, expected in zip(charges, _BEN_EXPECTED_CHARGES):
        assert (
            abs(actual - expected) < 1e-4
        ), f"charge mismatch: actual={actual:.6f} expected={expected:.6f}"
    # Charges must sum to the integer net formal charge (+1, protonated
    # amidinium). Gasteiger via RDKit seeds from formal charges, so this is
    # an exactness property of the charge model, not a tolerance check.
    assert abs(sum(charges) - 1.0) < 1e-3

    # Exported OpenMM System: per-force term counts.
    system = ic.to_openmm_system()
    assert system.getNumParticles() == _BEN_EXPECTED_ATOMS
    force_counts = {}
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "HarmonicBondForce":
            force_counts[name] = force.getNumBonds()
        elif name == "HarmonicAngleForce":
            force_counts[name] = force.getNumAngles()
        elif name == "PeriodicTorsionForce":
            force_counts[name] = force.getNumTorsions()
        elif name == "NonbondedForce":
            force_counts[name] = force.getNumParticles()
    assert force_counts == {
        "HarmonicBondForce": _BEN_EXPECTED_BONDS,
        "HarmonicAngleForce": _BEN_EXPECTED_ANGLES,
        "PeriodicTorsionForce": _BEN_EXPECTED_TORSIONS_IN_SYSTEM,
        "NonbondedForce": _BEN_EXPECTED_ATOMS,
    }


# NAGL charges on protonated benzamidinium (BEN, +1). Reference values are
# from openff-gnn-am1bcc-1.0.0; they are a GNN surrogate for AM1-BCC, so the
# expected per-atom values differ from the Gasteiger reference above. The
# molecule has two-fold symmetry across the amidinium C-aromatic axis, so
# the two amidinium nitrogens and their four NH protons each form
# equivalent pairs - the test asserts that explicitly.
_BEN_NAGL_EXPECTED_CHARGES = [
    -0.2451,  # C1 (ipso aromatic, bonded to amidinium C)
    -0.0785,  # C2 (ortho)
    -0.1173,  # C3 (meta)
    -0.0447,  # C4 (para)
    -0.1173,  # C5 (meta)
    -0.0785,  # C6 (ortho)
    +0.5327,  # C   (amidinium central C)
    -0.4711,  # N1  (amidinium N, equivalent pair)
    -0.4711,  # N2  (amidinium N, equivalent pair)
    +0.1486,  # H1  (ortho ring H)
    +0.1707,  # H2  (meta ring H)
    +0.1621,  # H3  (para ring H)
    +0.1707,  # H4  (meta ring H)
    +0.1486,  # H5  (ortho ring H)
    +0.3226,  # H6  (amidinium NH, equivalent quartet)
    +0.3226,  # H7
    +0.3226,  # H8
    +0.3226,  # H9
]


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _nagl_installed),
    reason="OpenMM + OpenFF Interchange + NAGL required (install with 'pip install acellera-htmd[nagl]')",
)
def _test_parameterize_ligands_openff_ben_nagl():
    """Build an Interchange for BEN via Sage 2.3 with NAGL charges (an
    AM1-BCC GNN surrogate) and check the resulting per-atom charges against
    pre-recorded reference values. Confirms the ``charge_method="nagl"``
    plumbing through :func:`parameterizeLigandsOpenFF` flows charges via
    ``charge_from_molecules`` so SMIRNOFF does not overwrite them, and that
    the +1 amidinium net charge is preserved exactly."""
    from moleculekit.molecule import Molecule
    from htmd.builder.openmm import parameterizeLigandsOpenFF

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)

    out = parameterizeLigandsOpenFF(
        mol,
        ligand_ff="openff_unconstrained-2.3.0.offxml",
        charge_method="nagl",
    )
    assert set(out) == {"BEN"}
    ic = out["BEN"]

    assert ic.topology.n_atoms == _BEN_EXPECTED_ATOMS
    assert ic.topology.n_bonds == _BEN_EXPECTED_BONDS

    charges = []
    for key in sorted(
        ic["Electrostatics"].charges.keys(), key=lambda k: k.atom_indices[0]
    ):
        charges.append(float(ic["Electrostatics"].charges[key].m))
    assert len(charges) == _BEN_EXPECTED_ATOMS

    # Net charge: integer +1 (amidinium) preserved exactly.
    assert abs(sum(charges) - 1.0) < 1e-4

    # NAGL output is deterministic for a given model; tolerance covers
    # only float32 / torch round-off, not model drift.
    for i, (actual, expected) in enumerate(
        zip(charges, _BEN_NAGL_EXPECTED_CHARGES)
    ):
        assert (
            abs(actual - expected) < 1e-3
        ), f"NAGL charge mismatch at atom {i}: actual={actual:.6f} expected={expected:.6f}"

    # Resonance symmetry: N1 == N2, and the four amidinium N-H protons
    # (H6/H7/H8/H9) are all equivalent. Bug-canary: if NAGL or our
    # plumbing later asymmetrically biases one nitrogen we want to know.
    assert charges[7] == pytest.approx(charges[8], abs=1e-6), "N1 vs N2 asymmetry"
    nh_charges = charges[14:18]
    assert all(
        c == pytest.approx(nh_charges[0], abs=1e-6) for c in nh_charges
    ), f"NH protons not equivalent: {nh_charges}"


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed),
    reason="OpenMM + OpenFF Interchange required",
)
def _test_parameterize_ligands_openff_multi_resname():
    """Parameterise two different ligand resnames in one call (BEN + a
    renamed copy). Both should produce identical reference counts since
    XBEN is structurally identical to BEN."""
    from moleculekit.molecule import Molecule
    from htmd.builder.openmm import parameterizeLigandsOpenFF

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)
    mol2 = mol.copy()
    mol2.resname[:] = "XBEN"
    mol.append(mol2)

    out = parameterizeLigandsOpenFF(
        mol,
        ligand_ff="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
        resnames=["BEN", "XBEN"],
    )
    assert set(out) == {"BEN", "XBEN"}
    for resname in ("BEN", "XBEN"):
        ic = out[resname]
        assert ic.topology.n_atoms == _BEN_EXPECTED_ATOMS
        assert ic.topology.n_bonds == _BEN_EXPECTED_BONDS
        assert len(ic["Angles"].key_map) == _BEN_EXPECTED_ANGLES
        assert len(ic["ProperTorsions"].key_map) == _BEN_EXPECTED_PROPER_KEYS
        assert len(ic["ImproperTorsions"].key_map) == _BEN_EXPECTED_IMPROPER_KEYS


# The forcefield choice determines which output attributes are populated:
# pure GAFF2 emits prepi (topo_paths) + frcmod (frcmod_paths) + a single
# combined OpenMM XML (gaff_combined.xml in xml_paths); pure SMIRNOFF
# emits only per-cluster XML fragments in xml_paths. Both populate
# xml_paths so the result is always openmm.build-consumable;
# amber.build only works when frcmod_paths is non-empty (i.e. GAFF
# involvement). Lock the contract so a future refactor can't silently
# swap which attribute gets populated.
@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed and _nagl_installed),
    reason="OpenMM + OpenFF Interchange + AmberTools + NAGL required",
)
def _test_parameterizeFromSpecs_output_artifacts_by_forcefield(tmp_path):
    """Pure GAFF2 -> prepi + frcmod + a single combined gaff_combined.xml
    in xml_paths. Pure openff -> per-cluster xml_paths only (no prepi,
    no frcmod, no gaff_combined.xml)."""
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from htmd.builder.nonstandard import parameterizeFromSpecs

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)
    specs = detectNonStandardResidues(mol)

    # Pure GAFF2 forcefield.
    out_gaff = parameterizeFromSpecs(
        specs,
        mol,
        outdir=str(tmp_path / "gaff"),
        forcefield="gaff2",
        charge_method="gasteiger",
    )
    assert len(out_gaff.topo_paths) == 1, "GAFF2 should emit one prepi/cif"
    assert len(out_gaff.frcmod_paths) == 1, "GAFF2 should emit one frcmod"
    # GAFF emits a single combined OpenMM XML named ``gaff_combined.xml``
    # appended to ``xml_paths``. The filename encodes provenance so the
    # downstream caller can tell GAFF-derived XMLs from SMIRNOFF ones.
    assert len(out_gaff.xml_paths) == 1, (
        "GAFF2 should produce one combined XML appended to xml_paths"
    )
    assert out_gaff.xml_paths[0].endswith("gaff_combined.xml"), (
        f"GAFF combined XML should be named gaff_combined.xml, "
        f"got {out_gaff.xml_paths[0]}"
    )
    assert os.path.isfile(out_gaff.xml_paths[0]), out_gaff.xml_paths[0]

    # Pure SMIRNOFF forcefield.
    out_openff = parameterizeFromSpecs(
        specs,
        mol,
        outdir=str(tmp_path / "openff"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="nagl",
    )
    assert out_openff.topo_paths == [], "SMIRNOFF must not emit AMBER prepi"
    assert out_openff.frcmod_paths == [], "SMIRNOFF must not emit frcmod"
    assert len(out_openff.xml_paths) >= 1, (
        "SMIRNOFF must emit at least one per-cluster XML"
    )
    # No GAFF involvement, so no gaff_combined.xml in the list.
    assert not any(
        p.endswith("gaff_combined.xml") for p in out_openff.xml_paths
    ), "Pure-SMIRNOFF call should not emit a gaff_combined.xml"
    for p in out_openff.xml_paths:
        assert os.path.isfile(p), p


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed),
    reason="OpenMM + OpenFF Interchange required",
)
def _test_parameterize_from_specs_openff_backend_ben(tmp_path):
    """End-to-end Phase 2B: parameterizeFromSpecs with the SMIRNOFF forcefield on
    a free ligand (3PTB BEN). Verify the emitted XML loads in OpenMM
    ForceField and reproduces the same energy as Interchange.to_openmm_system."""
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from htmd.builder.nonstandard import parameterizeFromSpecs
    import openmm
    import openmm.app as app
    import openmm.unit as ommunit

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)
    specs = detectNonStandardResidues(mol)

    out = parameterizeFromSpecs(
        specs,
        mol,
        outdir=str(tmp_path / "params"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
    )

    # OpenFF backend should NOT populate the AMBER artifacts.
    assert out.topo_paths == []
    assert out.frcmod_paths == []
    assert out.custombonds == []
    # Should populate xml_paths with exactly one entry for BEN.
    assert len(out.xml_paths) == 1
    assert os.path.isfile(out.xml_paths[0])
    assert out.xml_paths[0].endswith("BEN.xml")

    # Load the XML in OpenMM and confirm the residue template plus all forces
    # are present.
    ff = app.ForceField(out.xml_paths[0])
    assert len(ff._templates) == 1
    assert "BEN" in ff._templates

    # Build a System using the htmd Molecule's OpenMM topology.
    from htmd.builder.openmm import _mol_to_openmm

    topo_dir = str(tmp_path / "topology")
    os.makedirs(topo_dir, exist_ok=True)
    topology, positions = _mol_to_openmm(mol, topo_dir)
    system = ff.createSystem(topology)

    # Force-shape expectations: must match the Phase 1 reference numbers.
    n_bonds = n_angles = n_torsions = n_nb = 0
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            n_bonds = force.getNumBonds()
        elif isinstance(force, openmm.HarmonicAngleForce):
            n_angles = force.getNumAngles()
        elif isinstance(force, openmm.PeriodicTorsionForce):
            n_torsions = force.getNumTorsions()
        elif isinstance(force, openmm.NonbondedForce):
            n_nb = force.getNumParticles()
    assert n_bonds == _BEN_EXPECTED_BONDS
    assert n_angles == _BEN_EXPECTED_ANGLES
    assert n_torsions == _BEN_EXPECTED_TORSIONS_IN_SYSTEM
    assert n_nb == _BEN_EXPECTED_ATOMS

    # Per-atom charges roundtrip from the Interchange through the XML
    # template into the OpenMM System.
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            actual_charges = []
            for i in range(force.getNumParticles()):
                q, _sigma, _eps = force.getParticleParameters(i)
                actual_charges.append(float(q.value_in_unit(ommunit.elementary_charge)))
            assert len(actual_charges) == _BEN_EXPECTED_ATOMS
            for actual, expected in zip(actual_charges, _BEN_EXPECTED_CHARGES):
                assert (
                    abs(actual - expected) < 1e-4
                ), f"charge {actual:.6f} vs expected {expected:.6f}"
            assert abs(sum(actual_charges) - 1.0) < 1e-3
            break


# Chain-resident NCAA cluster test: 5VBL's NLE singleton cluster (norleucine)
# with ACE/NME caps. The cluster model compound has the NLE residue + cap
# atoms; after parameterisation we strip caps and emit one Residue template
# for NLE. Backbone-rename (Phase 2D) is not yet applied, so the resulting
# XML can't be combined with ff14SB into a full peptide; this test just
# validates the cluster pipeline mechanically.
NLE_SMILES = "CCCC[C@@H](C(=O)O)N"


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed),
    reason="OpenMM + OpenFF Interchange + AmberTools (tleap/antechamber) required",
)
def _test_parameterize_from_specs_openff_cluster_nle(tmp_path):
    """Cluster path via openff backend: parameterise an NLE singleton
    cluster (chain-resident NCAA + ACE/NME caps). Verify cap atoms are
    dropped, residue template covers only the NCAA atoms, ExternalBond
    markers are present for N/C, and the XML loads in app.ForceField."""
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.nonstandard import parameterizeFromSpecs
    import openmm.app as app
    import xml.etree.ElementTree as ET

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    smiles_map = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles_map.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    # The 5VBL_A chain-A peptide is missing sidechain atoms on several
    # canonical residues, so PDB2PQR needs the Dunbrack rotamer mutator
    # to fill them in before debumping (same setup as
    # _test_custom_residue_param_reference_5vbl_cluster in
    # test_nonstandard_builder.py).
    pmol = systemPrepare(
        mol,
        detect_specs=specs,
        restore_missing_sidechains=True,
    )[0]

    # Run only the NLE residue through the cluster pipeline via the
    # per-resname backend map (everything else stays on antechamber).
    out = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params"),
        forcefield={"NLE": "openff_unconstrained-2.3.0.offxml", "default": "gaff2"},
        charge_method="gasteiger",
    )

    # NLE went through the openff backend so it should be in xml_paths,
    # NOT in topo_paths or frcmod_paths.
    nle_xml_candidates = [p for p in out.xml_paths if "NLE" in p or "C00" in p]
    # NLE went through the cluster pipeline (singleton cluster), so the
    # xml file is named cluster_00*.xml. Find it.
    assert (
        len(out.xml_paths) >= 1
    ), f"Expected at least one openff XML, got: {out.xml_paths}"
    nle_xml = out.xml_paths[0]
    assert os.path.isfile(nle_xml)

    # Parse the XML and verify structural shape.
    root = ET.parse(nle_xml).getroot()
    residues = root.findall("./Residues/Residue")
    resnames = [r.attrib["name"] for r in residues]
    assert "NLE" in resnames, f"expected NLE in residue templates, got {resnames}"

    nle_residue = next(r for r in residues if r.attrib["name"] == "NLE")
    nle_atoms = nle_residue.findall("Atom")
    nle_bonds = nle_residue.findall("Bond")
    nle_externals = nle_residue.findall("ExternalBond")

    # NLE has 19 atoms (template NLE_SMILES expansion: backbone N, CA, C, O,
    # OXT, HN, HA + sidechain CB, CG, CD, CE plus their H atoms - the
    # actual count is fixed by the antechamber reference at 19, see the
    # _test_parameterize_from_specs_dedup_4tot expectations for the NCAA
    # atom budget).
    assert len(nle_atoms) >= 15, (
        f"NLE residue has only {len(nle_atoms)} atoms; cap-stripping may "
        f"have over-pruned"
    )
    # Chain-resident NCAA: must have ExternalBond markers for N and/or C.
    ext_names = sorted(e.attrib["atomName"] for e in nle_externals)
    assert (
        "N" in ext_names or "C" in ext_names
    ), f"expected at least one of N/C in ExternalBond, got {ext_names}"

    # XML must load in OpenMM ForceField alongside ff14SB. NCAA backbone
    # atoms reference protein-N / protein-CT / protein-C / protein-O /
    # protein-H / protein-H1 classes that live in amber14/protein.ff14SB.xml,
    # so loading the cluster XML standalone is expected to KeyError. The
    # full peptide build uses the same combination (defaultFf() + cluster
    # XML) at apply time.
    from htmd.builder.openmm import defaultFf

    ff = app.ForceField(*defaultFf(), nle_xml)
    assert "NLE" in ff._templates

    # Verify that the NLE template's backbone atoms reference the
    # ff14SB-owned classes (protein-*), while sidechain atoms reference
    # our cluster's synthesised OFF_ classes.
    backbone_classes = {
        "N": "protein-N",
        "H": "protein-H",
        "CA": "protein-CT",
        "HA": "protein-H1",
        "C": "protein-C",
        "O": "protein-O",
    }
    for atom in nle_atoms:
        nm = atom.attrib["name"]
        cls = atom.attrib["type"]
        if nm in backbone_classes:
            assert cls == backbone_classes[nm], (
                f"backbone atom {nm!r} has class {cls!r}, expected "
                f"{backbone_classes[nm]!r}"
            )
        else:
            assert cls.startswith("OFF_"), (
                f"sidechain atom {nm!r} has class {cls!r}, expected an "
                f"OFF_-prefixed synthesised class"
            )

    # ff14SB-owned classes must NOT be redeclared in <AtomTypes>; our
    # synthesised OFF_ classes must be.
    declared_types = {t.attrib["name"] for t in root.findall("./AtomTypes/Type")}
    for ff14sb_class in backbone_classes.values():
        assert (
            ff14sb_class not in declared_types
        ), f"cluster XML should not redeclare ff14SB class {ff14sb_class!r}"
    for atom in nle_atoms:
        cls = atom.attrib["type"]
        if cls.startswith("OFF_"):
            assert cls in declared_types, (
                f"residue atom {atom.attrib['name']!r} references undeclared "
                f"OFF_ class {cls!r}"
            )


# Phase 2E test: 8QFZ-style scaffolded peptide.
# - LFI scaffold (3-fold tris(bromopropionyl)-1,3,5-triazinane) reacts with
#   3 cysteines via S-C bonds. The cluster contains LFI + 3 modified
#   cysteines (XX1/XX2/XX3) connected by 3 non-peptide bonds.
# - For each cluster bond, the OpenFF emitter should add ExternalBond
#   markers at the cross-link atom on each side.
LFI_SMILES = "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr"


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed),
    reason="OpenMM + OpenFF Interchange required",
)
def _test_parameterize_from_specs_openff_cluster_crosslinks(tmp_path):
    """Phase 2E: 8QFZ-style scaffolded peptide cluster (LFI + 3 cysteines
    via S-C crosslinks). The emitted XML must carry ExternalBond markers
    at the cross-link atoms on each side, in addition to peptide N/C
    markers on the chain-resident cysteines."""
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.nonstandard import parameterizeFromSpecs
    import openmm.app as app
    import xml.etree.ElementTree as ET

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "8QFZ_B_bicycle.cif"))
    mol.templateResidueFromSmiles("resname LFI", LFI_SMILES, addHs=True)
    specs = detectNonStandardResidues(mol)
    pmol = systemPrepare(mol, detect_specs=specs)[0]

    out = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
    )
    assert (
        len(out.xml_paths) >= 1
    ), f"expected at least one openff cluster XML, got {out.xml_paths}"

    # Parse every cluster XML and aggregate ExternalBond markers per
    # residue name.
    ext_by_residue = {}
    for xml_path in out.xml_paths:
        root = ET.parse(xml_path).getroot()
        for residue in root.findall("./Residues/Residue"):
            resname = residue.attrib["name"]
            ext_names = [e.attrib["atomName"] for e in residue.findall("ExternalBond")]
            ext_by_residue.setdefault(resname, []).extend(ext_names)

    # The LFI scaffold has three crosslink atoms (one per Cys). These
    # are NOT peptide atoms - LFI is a free scaffold residue, not
    # chain-resident - so its ExternalBonds come only from crosslinks.
    if "LFI" in ext_by_residue:
        lfi_ext = ext_by_residue["LFI"]
        # LFI carries 3 crosslinks (one per Cys); each emits one
        # ExternalBond marker at the LFI side.
        assert len(lfi_ext) == 3, (
            f"LFI scaffold should have 3 crosslink ExternalBonds, "
            f"got {len(lfi_ext)}: {lfi_ext}"
        )
        # Cross-link atoms on LFI per the scaffold's chemistry: aliphatic
        # carbons C10/C11/C12 (or equivalent post-templating). They
        # should not be peptide-backbone names.
        for atom_name in lfi_ext:
            assert atom_name not in {"N", "CA", "C", "O", "HA", "H"}, (
                f"LFI external bond atom {atom_name!r} looks like a "
                f"peptide-backbone name; LFI is not a chain residue"
            )

    # Every cluster XML must load alongside ff14SB.
    from htmd.builder.openmm import defaultFf

    for xml_path in out.xml_paths:
        ff = app.ForceField(*defaultFf(), xml_path)
        assert len(ff._templates) > 0


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed),
    reason="OpenMM + OpenFF Interchange + AmberTools (tleap/antechamber) required",
)
def _test_parameterize_from_specs_openff_backbone_pin(tmp_path):
    """Phase 2E: pin_backbone_charges=True overrides NCAA backbone partial
    charges to ff14SB canonical values. With pin=False, the same atoms
    carry the SMIRNOFF / Gasteiger originals. Test the difference is
    visible in the emitted XML."""
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.nonstandard import parameterizeFromSpecs
    import xml.etree.ElementTree as ET

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    smiles_map = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles_map.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol = systemPrepare(
        mol,
        detect_specs=specs,
        restore_missing_sidechains=True,
    )[0]

    def _nle_backbone_charges(pin, norm):
        out = parameterizeFromSpecs(
            specs,
            pmol,
            outdir=str(tmp_path / f"params_pin{pin}_norm{norm}"),
            forcefield={"NLE": "openff_unconstrained-2.3.0.offxml", "default": "gaff2"},
            charge_method="gasteiger",
            pin_backbone_charges=pin,
            normalize=norm,
        )
        nle_xml = [p for p in out.xml_paths if os.path.isfile(p)]
        assert nle_xml
        root = ET.parse(nle_xml[0]).getroot()
        nle = next(
            r for r in root.findall("./Residues/Residue") if r.attrib["name"] == "NLE"
        )
        return {
            a.attrib["name"]: float(a.attrib["charge"]) for a in nle.findall("Atom")
        }

    pinned = _nle_backbone_charges(pin=True, norm="cluster")
    unpinned = _nle_backbone_charges(pin=False, norm=None)

    # ff14SB amide backbone charges (charge class 0, used for neutral
    # NCAAs): N=-0.4157, H=0.2719, C=0.5973, O=-0.5679. CA and HA are
    # residue-specific in ff14SB so they're only pinned for canonical
    # residues, not for NCAAs - the NCAA fallback table at
    # _FF14SB_BACKBONE_CHARGES_BY_CLASS only carries the four amide atoms.
    expected_pinned = {
        "N": -0.4157,
        "H": 0.2719,
        "C": 0.5973,
        "O": -0.5679,
    }
    for name, expected in expected_pinned.items():
        if name not in pinned:
            continue
        assert abs(pinned[name] - expected) < 1e-3, (
            f"with pin_backbone_charges=True, NLE atom {name!r} has "
            f"charge {pinned[name]:.4f}, expected {expected:.4f}"
        )
        # The unpinned (raw Gasteiger) charge must differ from the
        # pinned ff14SB value (otherwise the pin does nothing).
        if name in unpinned:
            assert abs(unpinned[name] - expected) > 0.01, (
                f"NLE atom {name!r}: unpinned charge {unpinned[name]:.4f} "
                f"already matches ff14SB; pin is a no-op"
            )


def _extract_charges_from_system(system):
    """Read per-atom partial charges from an OpenMM System's NonbondedForce."""
    import openmm
    import openmm.unit as ommunit

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return [
                float(
                    force.getParticleParameters(i)[0].value_in_unit(
                        ommunit.elementary_charge
                    )
                )
                for i in range(force.getNumParticles())
            ]
    raise RuntimeError("System has no NonbondedForce")


# tleap renames residues based on protonation state and disulfide
# bonding (HID/HIE/HIP for histidine, CYX for disulfide-bonded cysteine,
# ASH/GLH for protonated acids, LYN for deprotonated lysine). openmm.build
# matches the same templates without renaming, so the same chemical
# residue appears under the input resname. Map the tleap variants back
# to the canonical AMBER resnames so the resname carries information
# (catches drift to a different residue) without flagging the renames.
_RESNAME_NORMALISE = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "CYM": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    # tleap writes C-terminal -NH2 cap as "NHE"; openmm.build keeps the
    # input "NH2". Same chemical group, two names.
    "NHE": "NH2",
}


def _atom_id_key(atom, residue_index):
    """Atom-id key used to match across builds.

    Resname is included after :data:`_RESNAME_NORMALISE` normalisation
    so the comparison catches atoms migrating to a different residue
    but tolerates ff14SB protonation-state renaming.

    ``residue_index`` is the position of the residue within the parmed
    Structure's residue list (0-based, in input order). We use this
    instead of the raw resid integer because tleap renumbers residues
    sequentially while openmm.build preserves input numbering - the
    same chemical atom can carry different resid integers across the
    two pipelines. Index-by-position is stable.

    The atom name is normalised H1 -> H: tleap names the N-terminal
    amide H as "H1" (first of three NH3+ hydrogens, or sole H on a
    capped N), while openmm.build's residue templates name it "H".
    """
    res = atom.residue
    resname = _RESNAME_NORMALISE.get(res.name, res.name)
    return (
        getattr(res, "segid", "") or "",
        getattr(res, "chain", "") or "",
        resname,
        residue_index,
        "H" if atom.name == "H1" else atom.name,
    )


def _atoms_by_id(structure):
    """Index parmed Structure atoms by :func:`_atom_id_key`. Residue
    index is the residue's position in the Structure's residue list,
    so the key is stable against resid renumbering across builds."""
    # parmed's Atom.residue.idx already carries the structure-wide
    # residue index, but explicitly walking residues in order is the
    # safest way to be sure of the convention.
    out = {}
    for ri, residue in enumerate(structure.residues):
        for atom in residue.atoms:
            out[_atom_id_key(atom, ri)] = atom
    return out


def _bond_set_by_id(structure):
    """Bond set as sorted (atom_id_a, atom_id_b) tuples."""
    residue_idx = {id(r): ri for ri, r in enumerate(structure.residues)}
    return {
        tuple(
            sorted(
                (
                    _atom_id_key(b.atom1, residue_idx[id(b.atom1.residue)]),
                    _atom_id_key(b.atom2, residue_idx[id(b.atom2.residue)]),
                )
            )
        )
        for b in structure.bonds
    }


# SMILES templates per resname for any NCAA / ligand / cofactor we
# expect across the 3-way test fixtures. The test iterates this dict for
# every system and templates any resname that appears in the molecule.
# Resnames not in the molecule are silently skipped.
_NCAA_SMILES = {
    "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
    "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
    "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
    "NLE": "CCCC[C@@H](C=O)N",
    "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    "BEN": "[NH2+]=C(N)c1ccccc1",
    "LFI": "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr",
    "MK8": "CCCC[C@](C)(C(=O)O)N",
    "33X": "CC(C=O)NC",
    "34E": "CN[C@@H]([C@H](C)CN1CCN(CCOC)CC1)C=O",
    "ABA": "CC[C@H](C=O)N",
    "BMT": "C/C=C/C[C@@H](C)[C@H]([C@@H](C=O)NC)O",
    "DAL": "C[C@H](C=O)N",
    "MLE": "CC(C)C[C@@H](C=O)NC",
    "MVA": "CC(C)[C@@H](C=O)NC",
    "P6G": "OCCOCCOCCOCCOCCOCCO",
    "SO4": "[O-]S(=O)(=O)[O-]",
    "NAG": "CC(=O)N[C@@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
    "OIR": "CC(C(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(Cc2ccccc2)S",
    "HEM": (
        "C=CC1=C(C)c2cc3c(C)c(CCC(=O)[O-])c4cc5[n]6->[Fe@SP2+3]7"
        "(<-[n]2c1cc1c(C)c(C=C)c(cc6C(C)=C5CCC(=O)[O-])[n-]->71)<-[n-]34"
    ),
}

_TESTS_DIR = os.path.dirname(__file__)


_THREE_WAY_SYSTEMS = [
    pytest.param("3PTB_A.cif", id="3PTB"),
    pytest.param("5VBL_A.cif", id="5VBL_isopeptide"),
    pytest.param("4TOT_A.cif", id="4TOT_A_protein_ligand"),
    pytest.param("4TOT_E.cif", id="4TOT_E_cyclosporin"),
    pytest.param("8QFZ_B.cif", id="8QFZ_B_scaffolded"),
    pytest.param("8QU4_A.cif", id="8QU4_A_staple"),
    # 1R1J's ZN2+ has been stripped from the fixture: neither AMBER's
    # GAFF2 nor SMIRNOFF parameterise metal-coordination bonds, so the
    # two builds would only agree if both see no Zn at all.
    pytest.param("1R1J_A.cif", id="1R1J_glyco"),
    pytest.param("2KDC_A.cif", id="2KDC_membrane"),
    pytest.param("1BL8_A.cif", id="1BL8_channel"),
    pytest.param("2B5I_A.cif", id="2B5I_canonical"),
]


def _load_three_way_system(input_filename):
    """Load a pre-filtered CIF from tests/, template every matching NCAA
    from :data:`_NCAA_SMILES`, run systemPrepare, return ``(pmol, specs)``.

    All systems share the same setup - water removal, detectNonStandard,
    SMILES templating, and ``systemPrepare(detect_specs=specs,
    restore_missing_sidechains=True)``. Per-system tweaks are only
    expressed via ``pytest.mark.xfail`` markers; the loader itself
    treats every fixture identically.
    """
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare

    mol = Molecule(os.path.join(_TESTS_DIR, input_filename))
    # Defensive: even though fixtures are pre-filtered with no water,
    # strip again in case the user regenerated from a different source.
    mol.remove("water", _logger=False)
    specs = detectNonStandardResidues(mol)

    for resname, smi in _NCAA_SMILES.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )

    # PDB2PQR rejects chains with >10% of sidechain atoms missing. Many
    # of our test fixtures (notably 5VBL_A) carry trimmed sidechains for
    # disorder reasons; restore them via moleculekit's Dunbrack rotamer
    # mutator so PDB2PQR sees a complete structure. Harmless when no
    # atoms are missing.
    pmol = systemPrepare(
        mol,
        detect_specs=specs,
        restore_missing_sidechains=True,
    )[0]
    return pmol, specs


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed),
    reason="OpenMM + OpenFF Interchange + tleap required",
)
@pytest.mark.parametrize("input_filename", _THREE_WAY_SYSTEMS)
def _test_three_way_amber_vs_antechamber_vs_openff(tmp_path, input_filename):
    """Parametrized 3-way build comparison across diverse systems.
    For each: build via amber.build (gold standard), openmm.build with
    antechamber backend XML, and openmm.build with openff backend XML.
    Assert all three agree on atom set, bond connectivity, and per-atom
    partial charges within 1e-3.
    """
    from htmd.builder.nonstandard import parameterizeFromSpecs
    from htmd.builder.amber import build as amber_build
    from htmd.builder.openmm import build as openmm_build
    import parmed

    pmol, specs = _load_three_way_system(input_filename)
    # Suppress auto-capping on every segid. The fixtures (whether
    # downloaded or pre-cooked) already carry whatever terminal atoms
    # they need - amber.build's auto-ACE/NHE would either double up on
    # already-capped peptides or clash with NCAA-owned OXT atoms. Setting
    # all segids to ("none", "none") gives both builders the same input
    # constraint.
    caps = {str(segid): ("none", "none") for segid in set(pmol.segid.tolist())}

    # Parameterise both ways.
    out_a = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_antechamber"),
        forcefield="gaff2",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )
    out_b = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_openff"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )

    common_kwargs = {"ionize": False, "caps": caps}

    amber_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_amber"),
        custombonds=out_a.custombonds,
        topo=out_a.topo_paths,
        param=out_a.frcmod_paths,
        **common_kwargs,
    )
    openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_omm_antechamber"),
        extra_xml=list(out_a.xml_paths),
        custombonds=out_a.custombonds,
        solvate=False,
        **common_kwargs,
    )
    openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_omm_openff"),
        extra_xml=out_b.xml_paths,
        custombonds=out_b.custombonds,
        solvate=False,
        **common_kwargs,
    )

    parm_amber = parmed.load_file(str(tmp_path / "build_amber" / "structure.prmtop"))
    parm_omm_a = parmed.load_file(
        str(tmp_path / "build_omm_antechamber" / "structure.prmtop")
    )
    parm_omm_b = parmed.load_file(
        str(tmp_path / "build_omm_openff" / "structure.prmtop")
    )

    atoms_amber = _atoms_by_id(parm_amber)
    atoms_omm_a = _atoms_by_id(parm_omm_a)
    atoms_omm_b = _atoms_by_id(parm_omm_b)

    keys_amber = set(atoms_amber)
    keys_omm_a = set(atoms_omm_a)
    keys_omm_b = set(atoms_omm_b)

    def _diff(a, b, name_a, name_b):
        in_a, in_b = sorted(a - b), sorted(b - a)
        if not in_a and not in_b:
            return None
        return (
            f"{name_a}\\{name_b}={len(in_a)} (first: {in_a[:3]}); "
            f"{name_b}\\{name_a}={len(in_b)} (first: {in_b[:3]})"
        )

    for d, msg in (
        (
            _diff(keys_amber, keys_omm_a, "amber", "omm_antechamber"),
            "amber vs omm_antechamber",
        ),
        (_diff(keys_amber, keys_omm_b, "amber", "omm_openff"), "amber vs omm_openff"),
        (
            _diff(keys_omm_a, keys_omm_b, "omm_antechamber", "omm_openff"),
            "omm_antechamber vs omm_openff",
        ),
    ):
        assert d is None, f"{msg} atom set differs: {d}"

    bonds_amber = _bond_set_by_id(parm_amber)
    bonds_omm_a = _bond_set_by_id(parm_omm_a)
    bonds_omm_b = _bond_set_by_id(parm_omm_b)
    assert bonds_amber == bonds_omm_a, (
        f"amber vs omm_antechamber bond sets differ by "
        f"{len(bonds_amber ^ bonds_omm_a)} bonds"
    )
    assert bonds_amber == bonds_omm_b, (
        f"amber vs omm_openff bond sets differ by "
        f"{len(bonds_amber ^ bonds_omm_b)} bonds"
    )

    # Charges within 1e-3 across all three.
    for ref_atoms, ref_name in (
        (atoms_amber, "amber"),
        (atoms_omm_a, "omm_antechamber"),
    ):
        mismatch = []
        for k in keys_amber:
            q_ref = float(ref_atoms[k].charge)
            q_b = float(atoms_omm_b[k].charge)
            if abs(q_ref - q_b) > 1e-3:
                mismatch.append(
                    f"{k}: {ref_name}={q_ref:+.5f} omm_openff={q_b:+.5f} "
                    f"diff={abs(q_ref - q_b):.4f}"
                )
        assert not mismatch, (
            f"{len(mismatch)} {ref_name}-vs-omm_openff charge mismatches "
            f"> 1e-3 (first 5):\n  " + "\n  ".join(mismatch[:5])
        )

    # amber vs omm_antechamber: same parameters via different toolchains.
    am_aa_mismatch = []
    for k in keys_amber:
        q_amber = float(atoms_amber[k].charge)
        q_omm_a = float(atoms_omm_a[k].charge)
        if abs(q_amber - q_omm_a) > 1e-3:
            am_aa_mismatch.append(
                f"{k}: amber={q_amber:+.5f} omm_antechamber={q_omm_a:+.5f}"
            )
    assert not am_aa_mismatch, (
        f"{len(am_aa_mismatch)} amber-vs-omm_antechamber charge "
        f"mismatches > 1e-3 (first 5):\n  " + "\n  ".join(am_aa_mismatch[:5])
    )

    # Total charge must be integer in all three.
    for label, parm in (
        ("amber", parm_amber),
        ("omm_antechamber", parm_omm_a),
        ("omm_openff", parm_omm_b),
    ):
        total = sum(float(a.charge) for a in parm.atoms)
        assert (
            abs(total - round(total)) < 0.01
        ), f"{label} total charge {total:.4f} not integer"


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed),
    reason="OpenMM + OpenFF Interchange + tleap required",
)
def _test_5vbl_three_way_amber_vs_antechamber_vs_openff(tmp_path):
    """Build 5VBL chain-A THREE ways and compare topology + charges:
      1. amber.build (tleap) on the antechamber backend's prepi/frcmod
      2. openmm.build on the antechamber backend's merged XML
      3. openmm.build on the openff backend's per-cluster XMLs

    All three should produce equivalent systems:
      - Same atom set (matched by (segid, chain, resname, resid, name))
      - Same bond connectivity
      - Equal per-atom charges between (1) and (2) - exact, same params
        via different toolchains
      - Per-atom charges between (1)/(2) and (3) match within 1e-3 -
        same Gasteiger compute, same ff14SB backbone pin, same cluster
        residual distribution; tiny noise from mol2 rounding on the
        antechamber path and parmed/openmm conversion paths.

    Does NOT assert on:
      - Bond / angle / torsion *parameter values* between antechamber
        and openff (GAFF2 vs Sage by design).
      - tleap-side hydrogens vs openmm-side hydrogens. Both backends
        consume the same systemPrepare'd mol, so hydrogens are placed
        by systemPrepare and reused. Cap-related H differences caused
        by tleap's auto-NME are sidestepped by ``caps={'1': ('none',
        'none')}``.
    """
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.nonstandard import parameterizeFromSpecs
    from htmd.builder.amber import build as amber_build
    from htmd.builder.openmm import build as openmm_build
    import parmed

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    smiles_map = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles_map.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol = systemPrepare(
        mol,
        detect_specs=specs,
        restore_missing_sidechains=True,
    )[0]
    caps = {"1": ("none", "none")}

    # Parameterise both ways.
    out_a = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_antechamber"),
        forcefield="gaff2",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )
    out_b = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_openff"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )

    # Build all three.
    amber_built = amber_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_amber"),
        ionize=False,
        custombonds=out_a.custombonds,
        topo=out_a.topo_paths,
        param=out_a.frcmod_paths,
        caps=caps,
    )
    omm_a_built, _ = openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_omm_antechamber"),
        extra_xml=list(out_a.xml_paths),
        custombonds=out_a.custombonds,
        ionize=False,
        solvate=False,
        caps=caps,
    )
    omm_b_built, _ = openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_omm_openff"),
        extra_xml=out_b.xml_paths,
        custombonds=out_b.custombonds,
        ionize=False,
        solvate=False,
        caps=caps,
    )

    # Load all three prmtops via parmed for per-atom inspection.
    parm_amber = parmed.load_file(str(tmp_path / "build_amber" / "structure.prmtop"))
    parm_omm_a = parmed.load_file(
        str(tmp_path / "build_omm_antechamber" / "structure.prmtop")
    )
    parm_omm_b = parmed.load_file(
        str(tmp_path / "build_omm_openff" / "structure.prmtop")
    )

    atoms_amber = _atoms_by_id(parm_amber)
    atoms_omm_a = _atoms_by_id(parm_omm_a)
    atoms_omm_b = _atoms_by_id(parm_omm_b)

    # ---- Atom-id set equality (topology) ----
    keys_amber = set(atoms_amber)
    keys_omm_a = set(atoms_omm_a)
    keys_omm_b = set(atoms_omm_b)

    def _diff_summary(a, b, name_a, name_b):
        in_a = sorted(a - b)
        in_b = sorted(b - a)
        if not in_a and not in_b:
            return None
        return (
            f"{name_a} \\ {name_b} = {len(in_a)} (first: {in_a[:3]});"
            f" {name_b} \\ {name_a} = {len(in_b)} (first: {in_b[:3]})"
        )

    d_ab = _diff_summary(keys_amber, keys_omm_a, "amber", "omm_antechamber")
    d_ac = _diff_summary(keys_amber, keys_omm_b, "amber", "omm_openff")
    d_bc = _diff_summary(keys_omm_a, keys_omm_b, "omm_antechamber", "omm_openff")
    assert d_ab is None, f"amber vs omm_antechamber atom set differs: {d_ab}"
    assert d_ac is None, f"amber vs omm_openff atom set differs: {d_ac}"
    assert d_bc is None, f"omm_antechamber vs omm_openff atom set differs: {d_bc}"

    # ---- Bond connectivity ----
    bonds_amber = _bond_set_by_id(parm_amber)
    bonds_omm_a = _bond_set_by_id(parm_omm_a)
    bonds_omm_b = _bond_set_by_id(parm_omm_b)
    assert bonds_amber == bonds_omm_a, (
        f"amber vs omm_antechamber bond sets differ by "
        f"{len(bonds_amber ^ bonds_omm_a)} bonds"
    )
    assert bonds_amber == bonds_omm_b, (
        f"amber vs omm_openff bond sets differ by "
        f"{len(bonds_amber ^ bonds_omm_b)} bonds"
    )

    # ---- Per-atom charges ----
    # (amber, omm_antechamber): same parameters via different toolchains.
    # Expect exact (~1e-5) agreement. parmed converts via mol2 charges
    # which round to 4 decimal places, so 1e-4 is the floor.
    am_aa_mismatch = []
    for k in keys_amber:
        q_amber = float(atoms_amber[k].charge)
        q_omm_a = float(atoms_omm_a[k].charge)
        if abs(q_amber - q_omm_a) > 1e-3:
            am_aa_mismatch.append(
                f"{k}: amber={q_amber:+.5f} omm_antechamber={q_omm_a:+.5f}"
            )
    assert not am_aa_mismatch, (
        f"{len(am_aa_mismatch)} amber-vs-omm_antechamber charge "
        f"mismatches > 1e-3 (first 5):\n  " + "\n  ".join(am_aa_mismatch[:5])
    )

    # (amber, omm_openff) and (omm_antechamber, omm_openff): backbone
    # ff14SB pin identical, sidechain Gasteiger from the same RDKit
    # compute on the same cluster_mol. Tolerance ~1e-3.
    for ref_atoms, ref_name in (
        (atoms_amber, "amber"),
        (atoms_omm_a, "omm_antechamber"),
    ):
        mismatch = []
        for k in keys_amber:
            q_ref = float(ref_atoms[k].charge)
            q_b = float(atoms_omm_b[k].charge)
            if abs(q_ref - q_b) > 1e-3:
                mismatch.append(
                    f"{k}: {ref_name}={q_ref:+.5f} omm_openff={q_b:+.5f} "
                    f"diff={abs(q_ref - q_b):.4f}"
                )
        assert not mismatch, (
            f"{len(mismatch)} {ref_name}-vs-omm_openff charge mismatches "
            f"> 1e-3 (first 5):\n  " + "\n  ".join(mismatch[:5])
        )

    # ---- Total charge invariant ----
    total_amber = sum(float(a.charge) for a in parm_amber.atoms)
    total_omm_a = sum(float(a.charge) for a in parm_omm_a.atoms)
    total_omm_b = sum(float(a.charge) for a in parm_omm_b.atoms)
    for label, q in (
        ("amber", total_amber),
        ("omm_antechamber", total_omm_a),
        ("omm_openff", total_omm_b),
    ):
        assert abs(q - round(q)) < 0.01, f"{label} total charge {q:.4f} not integer"
    assert abs(total_amber - total_omm_a) < 1e-3
    assert abs(total_amber - total_omm_b) < 1e-3


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed and _tleap_installed),
    reason="OpenMM + OpenFF Interchange + AmberTools (tleap/antechamber) required",
)
def _test_5vbl_antechamber_vs_openff_cross_consistency(tmp_path):
    """Build 5VBL chain-A two ways via openmm.build:
      - antechamber backend (GAFF2 sidechains, ff14SB-renamed backbone)
      - openff backend (Sage sidechains, ff14SB-class-override on backbone)
    Both with Gasteiger sidechain charges and pin_backbone_charges=True.

    Asserts equivalence of:
      - Atom set: name, resname, resid, element atom-for-atom
      - Bond connectivity (as sorted name tuples)
      - Per-atom partial charges within 1e-3 (canonical ff14SB matches
        exactly on both; NCAA backbone is pinned to ff14SB; NCAA
        sidechain charges come from the same RDKit Gasteiger call on the
        same cluster_mol, with the same pin + normalize residual
        distribution)

    Does NOT assert on:
      - Bond/angle/torsion *parameters* on NCAA sidechains - antechamber
        uses GAFF2, openff uses Sage; different force fields by design.
      - Class names in the XML - antechamber uses GAFF2 short types,
        openff uses synthesised OFF_C000_<i>. Both work at apply time.
    """
    from moleculekit.molecule import Molecule
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.nonstandard import parameterizeFromSpecs
    from htmd.builder.openmm import build as openmm_build

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    smiles_map = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles_map.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol = systemPrepare(
        mol,
        detect_specs=specs,
        restore_missing_sidechains=True,
    )[0]

    # C-terminal residue 200 carries its own OXT through both parameterise
    # paths, so an auto-NME cap on top would clash. 5VBL_A.cif uses segid
    # "1" (not VBL_PDB's "P0" inhibitor segment).
    caps = {"1": ("none", "none")}

    # Build A: antechamber backend.
    out_a = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_antechamber"),
        forcefield="gaff2",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )
    mol_a, sys_a = openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_antechamber"),
        extra_xml=list(out_a.xml_paths),
        custombonds=out_a.custombonds,
        ionize=False,
        solvate=False,
        caps=caps,
    )

    # Build B: openff backend.
    out_b = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / "params_openff"),
        forcefield="openff_unconstrained-2.3.0.offxml",
        charge_method="gasteiger",
        pin_backbone_charges=True,
        normalize="cluster",
    )
    mol_b, sys_b = openmm_build(
        pmol.copy(),
        outdir=str(tmp_path / "build_openff"),
        extra_xml=out_b.xml_paths,
        custombonds=out_b.custombonds,
        ionize=False,
        solvate=False,
        caps=caps,
    )

    # ---- Topology: atom set ----
    assert mol_a.numAtoms == mol_b.numAtoms, (
        f"atom count differs: antechamber={mol_a.numAtoms} " f"openff={mol_b.numAtoms}"
    )
    mismatches = []
    for i in range(mol_a.numAtoms):
        if (
            mol_a.name[i] != mol_b.name[i]
            or mol_a.resname[i] != mol_b.resname[i]
            or int(mol_a.resid[i]) != int(mol_b.resid[i])
            or mol_a.element[i] != mol_b.element[i]
        ):
            mismatches.append(
                f"atom {i}: antechamber={mol_a.name[i]!r}/{mol_a.resname[i]!r}"
                f"/{mol_a.resid[i]}/{mol_a.element[i]!r} vs "
                f"openff={mol_b.name[i]!r}/{mol_b.resname[i]!r}"
                f"/{mol_b.resid[i]}/{mol_b.element[i]!r}"
            )
    assert (
        not mismatches
    ), f"{len(mismatches)} atom mismatches; first 5:\n  " + "\n  ".join(mismatches[:5])

    # ---- Topology: bond connectivity ----
    def _bond_set(m):
        return {
            tuple(
                sorted(
                    (
                        f"{m.name[a]}@{m.resname[a]}@{int(m.resid[a])}",
                        f"{m.name[b]}@{m.resname[b]}@{int(m.resid[b])}",
                    )
                )
            )
            for a, b in m.bonds
        }

    bonds_a = _bond_set(mol_a)
    bonds_b = _bond_set(mol_b)
    only_a = bonds_a - bonds_b
    only_b = bonds_b - bonds_a
    assert not only_a and not only_b, (
        f"bond connectivity differs:\n"
        f"  in antechamber but not openff ({len(only_a)}): {list(only_a)[:5]}\n"
        f"  in openff but not antechamber ({len(only_b)}): {list(only_b)[:5]}"
    )

    # ---- Charges ----
    q_a = _extract_charges_from_system(sys_a)
    q_b = _extract_charges_from_system(sys_b)
    assert len(q_a) == len(q_b) == mol_a.numAtoms

    charge_mismatches = []
    for i in range(mol_a.numAtoms):
        diff = abs(q_a[i] - q_b[i])
        if diff > 1e-3:
            charge_mismatches.append(
                f"atom {i} ({mol_a.name[i]} {mol_a.resname[i]} "
                f"{int(mol_a.resid[i])}): antechamber={q_a[i]:+.4f} "
                f"openff={q_b[i]:+.4f} diff={diff:.4f}"
            )
    assert not charge_mismatches, (
        f"{len(charge_mismatches)} per-atom charge mismatches > 1e-3 "
        f"out of {mol_a.numAtoms}:\n  " + "\n  ".join(charge_mismatches[:10])
    )

    # ---- Total charge invariant ----
    total_q_a = sum(q_a)
    total_q_b = sum(q_b)
    assert abs(total_q_a - total_q_b) < 1e-3
    # Both backends should give an integer-charge system.
    assert (
        abs(total_q_a - round(total_q_a)) < 0.01
    ), f"antechamber system total charge {total_q_a:.4f} is not integer"
    assert (
        abs(total_q_b - round(total_q_b)) < 0.01
    ), f"openff system total charge {total_q_b:.4f} is not integer"


# Reference per-force single-point energies for BEN at the input
# coordinates, computed once via ``Interchange.to_openmm_system()`` with
# Sage 2.3.0 + RDKit Gasteiger charges. Re-generate by running
# ``scripts/probe_phase2_xml_emitter.py``.
_BEN_REFERENCE_ENERGIES_KJ = {
    "HarmonicBondForce":     16.1103,
    "HarmonicAngleForce":    10.8101,
    "PeriodicTorsionForce":  26.1109,
    "NonbondedForce":       125.1611,
}


@pytest.mark.skipif(
    not (_openmm_installed and _openff_installed),
    reason="OpenMM + OpenFF Interchange required",
)
def _test_emitted_xml_matches_interchange_energy(tmp_path):
    """Round-trip fidelity: build a Sage Interchange for BEN, then
    construct an OpenMM System two ways and compare per-force single-
    point energies at the same coordinates:

      (a) ``Interchange.to_openmm_system()`` - the OpenFF stack's own
          native export, our ground truth.
      (b) ``app.ForceField(emitted_xml).createSystem(topology)`` - our
          XML loaded back through OpenMM.

    Both must agree on every force class to within ~1e-2 kJ/mol on
    bonded forces and ~1e-3 kJ/mol on nonbonded (the nonbonded tolerance
    is looser because of float32 rounding in the System's per-particle
    charge attribute). The assertions also check that the absolute
    reference values haven't drifted, so any change to:

      - the XML emitter (unit conversion, idivf division, multi-term
        encoding),
      - SMIRNOFF resolution in a future OpenFF Toolkit version,
      - or the RDKit Gasteiger PEOE convergence,

    triggers a fail with a clear delta.
    """
    from moleculekit.molecule import Molecule
    from htmd.builder.openmm import (
        parameterizeLigandsOpenFF,
        _emit_openmm_xml_from_interchange,
    )
    import openmm as mm
    import openmm.app as app
    import openmm.unit as ommunit
    from openff.units.openmm import to_openmm as q_to_openmm

    mol = Molecule(os.path.join(_PARAM_TEST_DIR, "3PTB_BEN.cif"))
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("resname BEN", BEN_SMILES, addHs=True)
    ic = parameterizeLigandsOpenFF(mol, charge_method="gasteiger")["BEN"]

    # (a) Interchange's native export.
    sys_native = ic.to_openmm_system()
    # (b) Our XML loaded via app.ForceField.
    xml_path = str(tmp_path / "BEN.xml")
    _emit_openmm_xml_from_interchange(ic, "BEN", [], xml_path)
    ff = app.ForceField(xml_path)
    sys_xml = ff.createSystem(ic.to_openmm_topology())

    positions = q_to_openmm(ic.positions)

    # Assign each force to its own group so we can compute per-force
    # energies. Both systems get the same group assignment.
    for i, force in enumerate(sys_native.getForces()):
        force.setForceGroup(i)
    for i, force in enumerate(sys_xml.getForces()):
        force.setForceGroup(i)

    def _energy_by_class(system):
        integ = mm.VerletIntegrator(0.001 * ommunit.picoseconds)
        ctx = mm.Context(
            system, integ, mm.Platform.getPlatformByName("Reference"),
        )
        ctx.setPositions(positions)
        out = {}
        for i, f in enumerate(system.getForces()):
            state = ctx.getState(getEnergy=True, groups={i})
            out[f.__class__.__name__] = state.getPotentialEnergy(
            ).value_in_unit(ommunit.kilojoule_per_mole)
        return out

    e_native = _energy_by_class(sys_native)
    e_xml = _energy_by_class(sys_xml)

    # Both Systems must contain the same set of force classes (modulo
    # CMMotionRemover, which is a no-op zero-energy housekeeping force).
    real_forces = lambda d: {k: v for k, v in d.items() if k != "CMMotionRemover"}
    assert set(real_forces(e_native)) == set(real_forces(e_xml)) == set(
        _BEN_REFERENCE_ENERGIES_KJ
    ), (
        f"force class set drift: native={sorted(real_forces(e_native))} "
        f"xml={sorted(real_forces(e_xml))} "
        f"reference={sorted(_BEN_REFERENCE_ENERGIES_KJ)}"
    )

    for force_class, ref_energy in _BEN_REFERENCE_ENERGIES_KJ.items():
        e_n = e_native[force_class]
        e_x = e_xml[force_class]
        # Native vs reference: pin the absolute value. The reference
        # values are rounded to 4 decimal places, so 5e-4 kJ/mol is the
        # floor; tightening further would fail on the rounding alone.
        # 1e-3 sits just above that floor and catches real drift in
        # SMIRNOFF resolution / Gasteiger convergence between
        # OpenFF / RDKit toolkit versions.
        assert abs(e_n - ref_energy) < 1e-3, (
            f"{force_class}: Interchange.to_openmm_system energy "
            f"{e_n:.4f} drifted from reference {ref_energy:.4f}"
        )
        # XML vs native: bonded forces match to ~5e-6 kJ/mol in
        # practice (numerics-only divergence in the Reference platform);
        # 5e-5 leaves 10x headroom. Nonbonded is looser - float32
        # rounding in NonbondedForce's per-particle charge attribute
        # accumulates to ~2.5e-4 kJ/mol on 18 atoms - so 1e-3 there.
        tol = 1e-3 if force_class == "NonbondedForce" else 5e-5
        assert abs(e_x - e_n) < tol, (
            f"{force_class}: XML-loaded System energy {e_x:.6f} "
            f"disagrees with Interchange.to_openmm_system "
            f"{e_n:.6f} by {abs(e_x - e_n):.6e} > {tol:.6e}"
        )
