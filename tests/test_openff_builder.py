# (c) 2015-2025 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
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
class _TestOpenFFBuilder:
    """Tests for htmd.builder.openff.build()."""

    def _test_defaultFf(self):
        from htmd.builder.openff import defaultFf

        ff = defaultFf()
        assert isinstance(ff, list)
        assert len(ff) > 0
        assert all(isinstance(f, str) for f in ff)
        assert any("protein" in f for f in ff)
        assert any("tip3p" in f for f in ff)

    def _test_resolve_ion_name(self):
        from htmd.builder.openff import _resolve_ion_name

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
        from htmd.builder.openff import _fix_water_naming

        mol = _make_water_mol("WAT", "OH2")
        _fix_water_naming(mol)
        assert np.all(mol.resname == "HOH")
        assert mol.name[0] == "O"


# ------------------------------------------------------------------
# Comparative build tests (mirrors test_amber_builder.py)
# ------------------------------------------------------------------


@pytest.mark.skipif(not _openmm_installed, reason="OpenMM not installed")
class _TestOpenFFComparative:
    """Build the same systems as test_amber_builder and compare energies."""

    def _test_protein_prepared(self, tmpdir):
        """3PTB with systemPrepare + autoSegment, no solvation."""
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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

    def _test_protein_ligand(self, tmpdir):
        """3PTB protein + benzamidine ligand with GAFF2 small-molecule FF."""
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
        from htmd.builder.openff import build as openff_build
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
