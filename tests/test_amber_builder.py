from moleculekit.molecule import Molecule, mol_equal
from htmd.builder.amber import (
    build,
    defaultFf,
    defaultParam,
    defaultTopo,
    defaultAmberHome,
    _defaultAmberSearchPaths,
    _findTeLeap,
)
import numpy as np
import sys
from glob import glob
import pytest
import shutil
import os

reason = "teLeap is not installed. Cannot test amber.build"
tleap_installed = _findTeLeap() is not None


curr_dir = os.path.dirname(os.path.abspath(__file__))


def _standardize_angles_and_dihedrals(mol):
    # 1. STANDARDIZE ANGLES (i - Central - k)
    # Rule: (i, j, k) is the same as (k, j, i).
    # We flip the row if the first index is larger than the last.
    angle_mask = mol.angles[:, 0] > mol.angles[:, 2]
    mol.angles[angle_mask] = mol.angles[angle_mask, ::-1]

    # 2. STANDARDIZE PROPERS (i - j - k - l)
    # Rule: Entire sequence can be reversed (l, k, j, i).
    di_mask = mol.dihedrals[:, 0] > mol.dihedrals[:, 3]
    mol.dihedrals[di_mask] = mol.dihedrals[di_mask, ::-1]

    # 3. STANDARDIZE IMPROPERS (i - j - Central - l)
    # Rule: Central atom (index 2) stays fixed; others are interchangeable.
    periph = [0, 1, 3]
    mol.impropers[:, periph] = np.sort(mol.impropers[:, periph], axis=1)

    # 4. GLOBAL LEXSORT (Deterministic row ordering)
    # We use [:, ::-1].T so lexsort prioritizes Column 0, then 1, etc.
    mol.angles = mol.angles[np.lexsort(mol.angles[:, ::-1].T)]
    mol.dihedrals = mol.dihedrals[np.lexsort(mol.dihedrals[:, ::-1].T)]
    mol.impropers = mol.impropers[np.lexsort(mol.impropers[:, ::-1].T)]


def _compareResultFolders(
    compare, tmpdir, pid, ignore_ftypes=(".log", ".txt", ".frcmod", ".crd")
):
    import filecmp

    def _cutfirstline(infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, "r") as fin:
            data = fin.read().splitlines(True)
        with open(outfile, "w") as fout:
            fout.writelines(data[1:])

    mol2 = Molecule(os.path.join(compare, "structure.prmtop"))
    mol2.read(os.path.join(compare, "structure.pdb"))
    mol = Molecule(os.path.join(tmpdir, "structure.prmtop"))
    mol.read(os.path.join(tmpdir, "structure.pdb"))

    _standardize_angles_and_dihedrals(mol)
    _standardize_angles_and_dihedrals(mol2)
    assert mol_equal(
        mol, mol2, checkFields=Molecule._connectivity_fields, uqBonds=True
    ), f"Bonding structure has changed in {compare} and {tmpdir}"
    assert mol_equal(
        mol, mol2, checkFields=Molecule._atom_fields
    ), f"Atom fields have changed in {compare} and {tmpdir}"
    assert mol_equal(
        mol,
        mol2,
        fieldPrecision={"coords": 2e-3},
        checkFields=Molecule._traj_fields,
        exceptFields=("fileloc"),
    ), f"Traj fields have changed in {compare} and {tmpdir}"

    try:
        from ffevaluation.ffevaluate import FFEvaluate, loadParameters

        # Strip waters to keep energy evaluation tractable for solvated systems.
        water_resnames = ("WAT", "HOH", "TIP3", "TIP4", "SPC")
        mol_e = mol.copy()
        mol_e.filter(~np.isin(mol_e.resname, water_resnames), _logger=False)
        mol2_e = mol2.copy()
        mol2_e.filter(~np.isin(mol2_e.resname, water_resnames), _logger=False)

        prm2 = loadParameters(os.path.join(compare, "structure.prmtop"))
        ffev2 = FFEvaluate(mol2_e, prm2)
        energies2, _, _ = ffev2.calculate(mol2_e.coords)

        prm = loadParameters(os.path.join(tmpdir, "structure.prmtop"))
        ffev = FFEvaluate(mol_e, prm)
        energies, _, _ = ffev.calculate(mol_e.coords)
        ene_diff = np.abs(energies - energies2)
        if any(ene_diff > 1e-3):
            print(f"ENERGY DIFF:\n{ene_diff}")
        else:
            print("SIMILAR ENERGIES")
    except Exception as e:
        print(f"Could not compare energies... {e}")

    files = []
    deletefiles = []
    for f in glob(os.path.join(compare, "*")):
        fname = os.path.basename(f)
        if os.path.splitext(f)[1] in ignore_ftypes:
            continue
        if f.endswith("prmtop"):
            _cutfirstline(f, os.path.join(compare, fname + ".mod"))
            _cutfirstline(
                os.path.join(tmpdir, fname), os.path.join(tmpdir, fname + ".mod")
            )
            files.append(os.path.basename(f) + ".mod")
            deletefiles.append(os.path.join(compare, fname + ".mod"))
        else:
            files.append(os.path.basename(f))

    match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
    if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
        raise RuntimeError(
            f"Different results produced by amber.build for test {pid} between {compare} and {tmpdir} in files {mismatch}."
        )

    for f in deletefiles:
        os.remove(f)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.parametrize("pid", ["3PTB"])
def test_with_protein_prepare(tmp_path, pid):
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.solvate import solvate
    from moleculekit.tools.autosegment import autoSegment

    # Test with systemPrepare
    # pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)

    np.random.seed(1)
    mol = Molecule(pid)
    mol.filter("protein")
    mol, _ = systemPrepare(mol, pH=7.0)
    mol.filter("protein")  # Fix for bad systemPrepare hydrogen placing
    mol = autoSegment(mol)
    smol = solvate(mol)
    ffs = defaultFf()
    _ = build(smol, ff=ffs, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "pp", pid)
    _compareResultFolders(refdir, tmp_path, pid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.parametrize("pid", ["3PTB"])
def test_without_protein_prepare(tmp_path, pid):
    from htmd.builder.solvate import solvate

    # Test without systemPrepare
    # pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
    np.random.seed(1)
    mol = Molecule(pid)
    mol.filter("protein")
    smol = solvate(mol)
    ffs = defaultFf()
    _ = build(smol, ff=ffs, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "nopp", pid)
    _compareResultFolders(refdir, tmp_path, pid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_protein_ligand(tmp_path):
    # Test protein ligand building with parametrized ligand
    refdir = os.path.join(curr_dir, "data", "test-amber-build", "protLig")

    mol = Molecule(os.path.join(refdir, "3ptb_mod.pdb"))
    lig = Molecule(os.path.join(refdir, "BEN.mol2"))
    lig.segid[:] = "L"

    newmol = Molecule()
    newmol.append(mol)
    newmol.append(lig)

    params = defaultParam() + [os.path.join(refdir, "BEN.frcmod")]
    topos = defaultTopo() + [os.path.join(refdir, "BEN.mol2")]

    _ = build(newmol, outdir=tmp_path, param=params, topo=topos, ionize=False)

    resdir = os.path.join(curr_dir, "data", "test-amber-build", "protLig", "results")
    _compareResultFolders(resdir, tmp_path, "3PTB with mol2")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_custom_disulfide_bonds(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)
    pid = "1GZM"
    mol = Molecule(pid)
    mol.filter("protein")
    mol.segid = mol.chain
    smol = solvate(mol)
    ffs = defaultFf()
    disu = [
        ["chain A and resid 110", "chain A and resid 187"],
        ["chain B and resid 110", "chain B and resid 187"],
    ]
    outdir = os.path.join(tmp_path, "with_disulfide_bonds")
    _ = build(smol, ff=ffs, outdir=outdir, disulfide=disu)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "nopp", pid)
    _compareResultFolders(refdir, outdir, pid)

    np.random.seed(1)
    outdir = os.path.join(tmp_path, "without_disulfide_bonds")
    _ = build(smol, ff=ffs, outdir=outdir)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "nopp", pid)
    _compareResultFolders(refdir, outdir, pid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_custom_teleap_imports(tmp_path):
    from htmd.builder.solvate import solvate

    pid = "3PTB"
    np.random.seed(1)
    mol = Molecule(pid)
    mol.filter("protein")
    smol = solvate(mol)
    ffs = defaultFf()

    amberhome = defaultAmberHome()
    teleapimports = [
        os.path.join(amberhome, _defaultAmberSearchPaths["ff"]),
        os.path.join(amberhome, _defaultAmberSearchPaths["lib"]),
        os.path.join(amberhome, _defaultAmberSearchPaths["param"]),
    ]

    _ = build(smol, ff=ffs, outdir=tmp_path, teleapimports=teleapimports)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "nopp", pid)
    _compareResultFolders(refdir, tmp_path, pid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_rna(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)

    mol = Molecule("6VA1")
    smol = solvate(mol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "rna", "6VA1")
    _compareResultFolders(refdir, tmp_path, "6VA1")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_dna(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)

    mol = Molecule("1BNA")
    smol = solvate(mol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "dna", "1BNA")
    _compareResultFolders(refdir, tmp_path, "1BNA")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_protein_rna(tmp_path):
    from htmd.builder.solvate import solvate
    from moleculekit.tools.preparation import systemPrepare
    from moleculekit.tools.autosegment import autoSegment

    np.random.seed(1)

    mol = Molecule("3WBM")
    mol.filter("not water")
    mol = autoSegment(mol, fields=("chain", "segid"))
    pmol, _ = systemPrepare(mol, pH=7.0)
    smol = solvate(pmol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "protein-rna", "3WBM")
    _compareResultFolders(refdir, tmp_path, "3WBM")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_caps(tmp_path):
    from htmd.builder.solvate import solvate
    from moleculekit.tools.preparation import systemPrepare

    np.random.seed(1)

    mol = Molecule("6A5J")
    pmol, _ = systemPrepare(mol, pH=7.0)
    smol = solvate(pmol)

    outdir = os.path.join(tmp_path, "out")
    _ = build(smol, outdir=outdir, ionize=False)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "peptide-cap", "6A5J")
    _compareResultFolders(refdir, outdir, "6A5J")
    shutil.rmtree(outdir)

    np.random.seed(1)

    mol = Molecule("6A5J")
    pmol, _ = systemPrepare(mol, pH=7.0)
    pmol.remove("(resid 1 13 and not backbone) or (resid 13 and name OXT)")
    smol = solvate(pmol)

    _ = build(smol, outdir=outdir, ionize=False)

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "peptide-cap-only-backbone", "6A5J"
    )
    _compareResultFolders(refdir, outdir, "6A5J")
    shutil.rmtree(outdir)

    _ = build(
        smol,
        outdir=outdir,
        ionize=False,
        caps={"resid 1": "ACE", "resid 13": "NME"},
    )

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "peptide-cap-only-backbone", "6A5J"
    )
    _compareResultFolders(refdir, outdir, "6A5J")
    shutil.rmtree(outdir)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.parametrize("pdbid", ["1AWF"])
def test_non_standard_residue_building(tmp_path, pdbid):
    # Builds a protein carrying a non-standard residue (1AWF: mid-chain TYS)
    # from user-supplied prepi/frcmod topology and compares against a reference.
    homedir = os.path.join(curr_dir, "data", "test-amber-build", "non-standard")

    protdir = os.path.join(homedir, pdbid)
    pmol = Molecule(os.path.join(protdir, f"{pdbid}_nolig.pdb"))
    _ = build(
        pmol,
        topo=glob(os.path.join(protdir, "*.prepi")),
        param=glob(os.path.join(protdir, "*.frcmod")),
        ionize=False,
        outdir=tmp_path,
    )
    refdir = os.path.join(homedir, pdbid, "build")
    _compareResultFolders(refdir, tmp_path, pdbid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.skipif(sys.platform.startswith("darwin"), reason="Fails on OSX")
def test_cofactor_building(tmp_path):
    homedir = os.path.join(curr_dir, "data", "test-amber-build", "cofactors")

    mol = Molecule(os.path.join(homedir, "cofactors.pdb"))
    _ = build(mol, ionize=False, outdir=tmp_path)
    refdir = os.path.join(curr_dir, "data", "test-amber-build", "cofactors", "build")

    _compareResultFolders(refdir, tmp_path, "Cofactors")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_post_translational_modifications_building(tmp_path):
    homedir = os.path.join(curr_dir, "data", "test-amber-build", "post-translational")

    mol = Molecule(os.path.join(homedir, "4EFP_nolig.pdb"))
    _ = build(mol, ionize=False, outdir=tmp_path)
    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "post-translational", "build"
    )

    _compareResultFolders(refdir, tmp_path, "post-translational modifications")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_cyclic_peptide_building(tmp_path):
    np.random.seed(1)

    mol = Molecule("5VAV")
    mol.segid = mol.chain
    _ = build(mol, outdir=tmp_path, ionize=False, caps={"A": ("none", "none")})

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "cyclic-peptide", "5VAV"
    )
    _compareResultFolders(refdir, tmp_path, "5VAV")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_membrane(tmp_path):
    np.random.seed(1)

    homedir = os.path.join(curr_dir, "data", "test-amber-build", "membrane")

    mol = Molecule(os.path.join(homedir, "structure.psf"))
    mol.read(os.path.join(homedir, "structure.pdb"))
    _ = build(mol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "membrane", "build")
    _compareResultFolders(refdir, tmp_path, "Membrane")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_cif_building(tmp_path):
    # Tests that mol2 and cif building produce same results
    np.random.seed(1)

    homedir = os.path.join(curr_dir, "data", "test-amber-build", "parameterize-cif")
    ciff = os.path.join(homedir, "MOL-orig.cif")
    fmod = os.path.join(homedir, "MOL.frcmod")

    mol = Molecule("3ptb")
    mol.remove("resname BEN")
    mol2 = Molecule(ciff)
    mol2.segid[:] = "L"
    mol.append(mol2)

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "parameterize-cif", "build"
    )

    build(mol, ionize=False, outdir=tmp_path, topo=[ciff], param=[fmod])
    _compareResultFolders(refdir, tmp_path, "3PTB")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_ionize_salt(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)
    mol = Molecule("3PTB")
    mol.filter("protein")
    smol = solvate(mol, pad=5)
    molbuilt = build(smol, outdir=tmp_path, ionize=True, saltconc=0.15)

    assert molbuilt is not None
    assert molbuilt.numAtoms > 0

    ncl = np.sum(molbuilt.resname == "CL")
    nna = np.sum(molbuilt.resname == "NA")
    assert ncl == 15, f"Expected 15 CL ions, got {ncl}"
    assert nna == 9, f"Expected 9 NA ions, got {nna}"

    assert not os.path.exists(
        os.path.join(tmp_path, "solute_charge.prmtop")
    ), "Temp solute_charge files not cleaned up"
    assert os.path.exists(os.path.join(tmp_path, "tleap_solute.in"))


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_atomtype_decollisioning(tmp_path):
    # Tests that mol2 and cif building produce same results
    np.random.seed(1)

    homedir = os.path.join(
        curr_dir, "data", "test-amber-build", "atomtype-decollisioning"
    )
    mol1f = os.path.join(homedir, "m42.cif")
    mol2f = os.path.join(homedir, "m48.cif")
    fmod1 = os.path.join(homedir, "m42.frcmod")
    fmod2 = os.path.join(homedir, "m48.frcmod")

    mol1 = Molecule(mol1f)
    mol1.resname[:] = "m42"
    mol1.segid[:] = "L1"
    mol2 = Molecule(mol2f)
    mol2.resname[:] = "m48"
    mol2.segid[:] = "L2"
    mol1.append(mol2)

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "atomtype-decollisioning", "build"
    )
    build(
        mol1, ionize=False, outdir=tmp_path, topo=[mol1f, mol2f], param=[fmod1, fmod2]
    )
    _compareResultFolders(
        refdir,
        tmp_path,
        "decollisioning",
        ignore_ftypes=(".log", ".txt", ".frcmod", ".in", ".mol2", ".cif"),
    )


# ----------------------------------------------------------------------
# Regression tests for the bond-directive residue indexing.
#
# tLeap's `<unit>.<N>.<atom>` selector resolves against the residue
# position in the *combined* unit produced by load + combine, NOT the
# PDB resid field. The bond-directive writer in _write_tleap_script must
# therefore compute each atom's position in the load order
# (solute -> water -> cyclic_*, matching the order of the `combine`
# directive it emits). These tests pin that contract: changing the
# combine order without updating the position helper, or vice versa,
# will break them.
# ----------------------------------------------------------------------


def _ala_mol(resid, segid, chain, x0, resname="ALA"):
    """Build a synthetic ALA residue at a chosen position - just enough
    atoms for _tleap_residue_positions to recognise it as a residue.
    """
    mol = Molecule().empty(4)
    mol.name[:] = ["N", "CA", "C", "O"]
    mol.element[:] = ["N", "C", "C", "O"]
    mol.resname[:] = resname
    mol.resid[:] = resid
    mol.chain[:] = chain
    mol.segid[:] = segid
    mol.record[:] = "ATOM"
    mol.coords = np.array(
        [
            [x0 + 0.0, 0.0, 0.0],
            [x0 + 1.5, 0.0, 0.0],
            [x0 + 2.5, 1.0, 0.0],
            [x0 + 2.0, 2.0, 0.0],
        ],
        dtype=np.float32,
    ).reshape(4, 3, 1)
    return mol


def _hoh_mol(resid, segid, x0):
    """Build a synthetic water residue."""
    mol = Molecule().empty(3)
    mol.name[:] = ["O", "H1", "H2"]
    mol.element[:] = ["O", "H", "H"]
    mol.resname[:] = "HOH"
    mol.resid[:] = resid
    mol.chain[:] = "W"
    mol.segid[:] = segid
    mol.record[:] = "HETATM"
    mol.coords = np.array(
        [[x0, 0, 0], [x0 + 1, 0, 0], [x0, 1, 0]], dtype=np.float32
    ).reshape(3, 3, 1)
    return mol


def test_tleap_residue_positions_solute_only():
    """All non-water non-cyclic residues get sequential positions 1..N."""
    from htmd.builder.amber import _tleap_residue_positions

    parts = [_ala_mol(r, "P0", "A", x0=r * 5.0) for r in (1, 2, 3)]
    mol = Molecule()
    for p in parts:
        mol.append(p)

    pos = _tleap_residue_positions(mol, cyc_info=[], include_water=False)
    # Three ALA residues, each with 4 atoms -> positions [1,1,1,1, 2,2,2,2, 3,3,3,3].
    assert list(pos) == [1] * 4 + [2] * 4 + [3] * 4


def test_tleap_residue_positions_solute_then_water():
    """Solute residues get 1..N_solute, then waters get N_solute+1..."""
    from htmd.builder.amber import _tleap_residue_positions

    mol = Molecule()
    mol.append(_ala_mol(1, "P0", "A", x0=0.0))
    mol.append(_hoh_mol(2, "P1", x0=10.0))  # water in the middle of mol order
    mol.append(_hoh_mol(3, "P1", x0=12.0))
    mol.append(_ala_mol(4, "P2", "B", x0=20.0))  # solute again, after waters

    pos = _tleap_residue_positions(mol, cyc_info=[], include_water=True)
    # Two solute residues (4 atoms each) followed by two waters (3 atoms each).
    # Solute first -> positions 1 and 2 - this is the load-order in the
    # final combined unit, NOT the mol order which has waters interleaved.
    ala0_pos = set(pos[mol.atomselect(f"resid 1")].tolist())
    ala4_pos = set(pos[mol.atomselect(f"resid 4")].tolist())
    hoh2_pos = set(pos[mol.atomselect(f"resid 2 and water")].tolist())
    hoh3_pos = set(pos[mol.atomselect(f"resid 3 and water")].tolist())
    assert ala0_pos == {1}, "first solute residue must be at position 1"
    assert ala4_pos == {2}, "second solute residue must be at position 2"
    assert hoh2_pos == {3}, "waters must be appended after all solute"
    assert hoh3_pos == {4}


def test_tleap_residue_positions_cyclic_appended_in_cyc_info_order():
    """Cyclic-segment residues are appended after solute+water, in the
    same order the _write_tleap_script combine emits.
    """
    from htmd.builder.amber import _tleap_residue_positions

    mol = Molecule()
    mol.append(_ala_mol(1, "P0", "A", x0=0.0))  # solute
    mol.append(_hoh_mol(2, "W0", x0=10.0))  # water
    mol.append(_ala_mol(3, "CYC1", "X", x0=20.0))  # first cyclic seg
    mol.append(_ala_mol(4, "CYC2", "Y", x0=30.0))  # second cyclic seg

    # cyc_info entries name the cyclic segs in the order they're combined.
    # See _write_tleap_script: `mol = combine {mol wat cyc_CYC1 cyc_CYC2}`.
    cyc_info = [
        ("cyc_CYC1", "cyclic_CYC1.pdb", 1, 1),
        ("cyc_CYC2", "cyclic_CYC2.pdb", 1, 1),
    ]
    pos = _tleap_residue_positions(mol, cyc_info=cyc_info, include_water=True)

    solute_pos = set(pos[mol.atomselect("resid 1")].tolist())
    water_pos = set(pos[mol.atomselect("water")].tolist())
    cyc1_pos = set(pos[mol.segid == "CYC1"].tolist())
    cyc2_pos = set(pos[mol.segid == "CYC2"].tolist())
    # Order in combined unit: solute (1) -> water (2) -> CYC1 (3) -> CYC2 (4).
    assert solute_pos == {1}
    assert water_pos == {2}
    assert cyc1_pos == {3}
    assert cyc2_pos == {4}


def _cyclic_peptide_mol(closure_dist, add_closure_bond):
    """Synthetic 3-residue protein segment whose first-residue N and
    last-residue C are ``closure_dist`` Angstrom apart. When
    ``add_closure_bond`` is True an explicit bond is added between them
    (as a real cyclic peptide's head-to-tail amide would be in the input).
    """
    mol = Molecule()
    for r in (1, 2, 3):
        mol.append(_ala_mol(r, "P0", "A", x0=r * 5.0))
    mol.guessBonds()
    n_first = int(np.where((mol.resid == 1) & (mol.name == "N"))[0][0])
    c_last = int(np.where((mol.resid == 3) & (mol.name == "C"))[0][0])
    # Place the first N exactly closure_dist away from the last C.
    c_xyz = mol.coords[c_last, :, 0]
    mol.coords[n_first, :, 0] = c_xyz + np.array([closure_dist, 0.0, 0.0], np.float32)
    if add_closure_bond:
        mol.bonds = np.vstack([mol.bonds, [[n_first, c_last]]]).astype(mol.bonds.dtype)
        if mol.bondtype is not None and len(mol.bondtype):
            mol.bondtype = np.append(mol.bondtype, "1")
    return mol


def test_detect_cyclic_segments_honors_explicit_bond():
    """A head-to-tail amide modeled longer than 1.35 A (e.g. 7BTI's
    1.468 A HYP-CYS closure) is still cyclic when the input carries an
    explicit first-N to last-C bond. Distance alone would miss it."""
    from htmd.builder.amber import _detect_cyclic_segments

    mol = _cyclic_peptide_mol(closure_dist=1.468, add_closure_bond=True)
    cyclic = _detect_cyclic_segments(mol)
    assert cyclic == [("P0", 1, 3)]


def test_detect_cyclic_segments_distance_fallback():
    """With no explicit closure bond, a sub-1.35 A first-N/last-C
    distance still flags the segment cyclic (legacy behavior preserved)."""
    from htmd.builder.amber import _detect_cyclic_segments

    mol = _cyclic_peptide_mol(closure_dist=1.32, add_closure_bond=False)
    cyclic = _detect_cyclic_segments(mol)
    assert cyclic == [("P0", 1, 3)]


def test_detect_cyclic_segments_linear_not_cyclic():
    """A linear peptide (no closure bond, ends far apart) is not cyclic -
    guards the explicit-bond addition against false positives."""
    from htmd.builder.amber import _detect_cyclic_segments

    mol = _cyclic_peptide_mol(closure_dist=1.468, add_closure_bond=False)
    cyclic = _detect_cyclic_segments(mol)
    assert cyclic == []


def test_detect_cyclic_segments_ignores_implausibly_long_explicit_bond():
    """A first-N to last-C bond stretched far beyond any real amide (e.g. a
    misassigned LINK/CONECT record) is NOT treated as a cyclic closure -
    honoring it would make tLeap close a ring across a nonsensical gap."""
    from htmd.builder.amber import _detect_cyclic_segments

    mol = _cyclic_peptide_mol(closure_dist=5.0, add_closure_bond=True)
    cyclic = _detect_cyclic_segments(mol)
    assert cyclic == []


def test_detect_modaa_residues_autoloads_modaa_leaprc():
    """A modified amino acid that lives only in AMBER's mod_amino.lib (e.g. MSE,
    ALY) - not the base ff14SB libraries - must auto-load
    leaprc.protein.ff14SB_modAA, otherwise tleap reports 'Unknown residue'."""
    from htmd.builder.amber import _detect_modaa_residues, defaultFf

    mol = Molecule().empty(2)
    mol.resname[:] = "MSE"
    mol.name[:] = ["N", "CA"]
    mol.element[:] = ["N", "C"]
    ff = defaultFf()
    detected = _detect_modaa_residues(mol, ff)
    assert "MSE" in detected
    assert "leaprc.protein.ff14SB_modAA" in ff


def test_detect_modaa_residues_picks_ff19_variant():
    """When the build uses ff19SB, the ff19SB_modAA variant is loaded."""
    from htmd.builder.amber import _detect_modaa_residues

    mol = Molecule().empty(1)
    mol.resname[:] = "ALY"
    mol.name[:] = ["NZ"]
    mol.element[:] = ["N"]
    ff = ["leaprc.protein.ff19SB", "leaprc.water.tip3p"]
    _detect_modaa_residues(mol, ff)
    assert "leaprc.protein.ff19SB_modAA" in ff
    assert "leaprc.protein.ff14SB_modAA" not in ff


def test_detect_modaa_residues_noop_without_modaa():
    """A system with no modAA residue does not pull in the modAA leaprc."""
    from htmd.builder.amber import _detect_modaa_residues, defaultFf

    mol = Molecule().empty(1)
    mol.resname[:] = "ALA"
    mol.name[:] = ["CA"]
    mol.element[:] = ["C"]
    ff = defaultFf()
    detected = _detect_modaa_residues(mol, ff)
    assert detected == []
    assert not any("modAA" in f for f in ff)


def test_glycam_modrna_collision_set_is_1ma_2ma():
    """Pins the exact set of 3-letter codes shared between GLYCAM's
    systematically generated sugar unit names and AMBER's modrna08
    modified-ribonucleotide names. ``_is_glycam_sugar_not_modrna`` tells
    the two apart by checking for a nitrogen atom, which is only valid
    because the current collision set happens to be entirely mannose
    (N-free) codes. If this test fails, a new collision has landed on a
    different code: re-examine ``_is_glycam_sugar_not_modrna`` before
    trusting it, since a GlcNAc/GalNAc/sialic-acid GLYCAM unit (all of
    which DO carry nitrogen) would be misclassified as the modrna
    residue, silently, in all three of its call sites at once."""
    from moleculekit.residues import MODIFIED_NUCLEIC_RESIDUE_NAMES
    from moleculekit.tools.glycans import GLYCAM_UNIT_NAMES

    collision = set(MODIFIED_NUCLEIC_RESIDUE_NAMES) & (
        set(GLYCAM_UNIT_NAMES) | {"ROH"}
    )
    assert collision == {"1MA", "2MA"}


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_detect_glycam_residues_ff_and_hydrogens():
    # needs AMBERHOME to locate the GLYCAM prep files
    from htmd.builder.amber import _detect_glycam_residues, defaultFf

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "3AVE_frag.pdb")
    )
    # simulate prepared naming (subset is enough for this unit test)
    ren = {1: "UYB", 2: "4YB", 3: "VMB", 4: "2MA", 5: "0YB", 6: "2MA", 7: "0YB", 8: "0fA"}
    for resid, code in ren.items():
        mol.resname[(mol.chain == "C") & (mol.resid == resid)] = code
    mol.resname[(mol.chain == "A") & (mol.resid == 297)] = "NLN"
    # 3AVE_frag.pdb is a bare crystallographic fragment with no hydrogens at
    # all; add one synthetic amide hydrogen on the anchor residue, as a real
    # systemPrepare/pdb2pqr protonation would, so the "anchor hydrogens are
    # untouched" assertion below is actually exercised rather than trivially
    # true.
    nd2_idx = np.where(
        (mol.chain == "A") & (mol.resid == 297) & (mol.name == "ND2")
    )[0][0]
    h_atom = Molecule().empty(1)
    h_atom.name[:] = "HD22"
    h_atom.element[:] = "H"
    h_atom.resname[:] = "NLN"
    h_atom.resid[:] = 297
    h_atom.chain[:] = "A"
    h_atom.record[:] = "ATOM"
    h_atom.coords = (
        mol.coords[nd2_idx, :, 0] + np.array([0.5, 0.5, 0.5], dtype=np.float32)
    ).reshape(1, 3, 1)
    mol.append(h_atom)
    mol.segid[:] = np.where(mol.chain == "C", "G1", "P1")
    mol.deleteBonds("all")

    ff = defaultFf()
    bonds = _detect_glycam_residues(mol, ff)
    # leaprc inserted before the water leaprc, not appended
    gidx = ff.index("leaprc.GLYCAM_06j-1")
    widx = [i for i, f in enumerate(ff) if "water" in f][0]
    assert gidx < widx
    # 7 glycosidic + 1 anchor bond
    assert len(bonds) == 8
    # sugar hydrogens stripped, NLN hydrogens kept
    assert not np.any((mol.segid == "G1") & (mol.element == "H"))
    assert np.any((mol.resname == "NLN") & (mol.element == "H"))


def test_detect_glycam_residues_noop_without_sugars():
    from htmd.builder.amber import _detect_glycam_residues, defaultFf

    mol = Molecule("3ptb")
    ff = defaultFf()
    before = list(ff)
    assert _detect_glycam_residues(mol, ff) == []
    assert ff == before


def test_detect_glycam_residues_ignores_tla_ligand_collision():
    """TLA (L-(+)-tartaric acid) is a real PDB Chemical Component Dictionary
    ligand code that also happens to be a valid GLYCAM-06j unit name (a
    fully-branched galactose linkage). Matching purely on resname would
    misdetect a crystallization-buffer tartrate as a glycan and send it
    through glycanBondsFromNames, which fails hard looking for a glycosidic
    partner that does not exist ("Expected exactly one anomeric partner...
    found 0"). glycamUnitMask's sugar-composition gate must keep a real TLA
    ligand out of glycan detection entirely, so it builds like any other
    ordinary HETATM residue."""
    from htmd.builder.amber import _detect_glycam_residues, defaultFf

    mol = Molecule().empty(10)
    mol.resname[:] = "TLA"
    mol.resid[:] = 1
    mol.chain[:] = "L"
    mol.segid[:] = "L1"
    mol.record[:] = "HETATM"
    # Real RCSB Chemical Component Dictionary atom names for TLA - note
    # there is no "O5" atom at all, unlike a genuine GLYCAM sugar unit,
    # which is what tells the two apart.
    mol.name[:] = ["O1", "O11", "C1", "C2", "O2", "C3", "O3", "C4", "O4", "O41"]
    mol.element[:] = ["O", "O", "C", "C", "O", "C", "O", "C", "O", "O"]
    mol.coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [2.5, 1.5, 0.0],
            [3.0, 0.0, 0.0],
            [3.5, 1.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.5, 1.0, 0.0],
            [5.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    ).reshape(10, 3, 1)

    ff = defaultFf()
    before = list(ff)
    assert _detect_glycam_residues(mol, ff) == []
    assert ff == before, "TLA must not trigger the GLYCAM leaprc auto-load"


def test_detect_glycam_residues_raises_on_anchor_without_sugar():
    """A GLYCAM anchor residue (NLN/OLS/OLT/OLP) present without any GLYCAM
    sugar indicates a partially-prepared or hand-edited structure. This must
    raise a clear, actionable error rather than silently building a bare
    protein with a stray GLYCAM-named residue and no glycan at all."""
    from htmd.builder.amber import _detect_glycam_residues, defaultFf

    mol = Molecule().empty(1)
    mol.resname[:] = "NLN"
    mol.name[:] = "ND2"
    mol.element[:] = "N"
    mol.resid[:] = 1
    mol.chain[:] = "A"
    mol.segid[:] = "P1"
    mol.record[:] = "ATOM"
    mol.coords = np.zeros((1, 3, 1), dtype=np.float32)
    ff = defaultFf()
    with pytest.raises(RuntimeError, match="GLYCAM anchor residue"):
        _detect_glycam_residues(mol, ff)


def test_detect_glycam_residues_raises_on_unrenamed_pdb_sugar():
    """A glycan that reaches amber.build still carrying its original PDB
    sugar resnames (NAG/BMA/MAN/FUC/...) instead of GLYCAM unit names -
    because the caller skipped systemPrepare, or bypassed GlycanSpec
    parameterization entirely - must raise a clear, actionable error
    pointing at systemPrepare, instead of dying deep inside tleap with
    "Unknown residue: NAG"."""
    from htmd.builder.amber import _detect_glycam_residues, defaultFf

    # Loaded as-is, with no GLYCAM renaming applied: 3AVE_frag.pdb's chain C
    # still carries the original PDB sugar resnames (NAG, BMA, MAN, FUC).
    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "3AVE_frag.pdb")
    )
    ff = defaultFf()
    with pytest.raises(RuntimeError, match="un-renamed PDB sugar"):
        _detect_glycam_residues(mol, ff)


def test_glycan_break_points():
    from htmd.builder.amber import _glycan_break_points

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "3AVE_frag.pdb")
    )
    ren = {1: "UYB", 2: "4YB", 3: "VMB", 4: "2MA", 5: "0YB", 6: "2MA", 7: "0YB", 8: "0fA"}
    for resid, code in ren.items():
        mol.resname[(mol.chain == "C") & (mol.resid == resid)] = code
    mol.segid[:] = np.where(mol.chain == "C", "G1", "P1")
    points = _glycan_break_points(mol)
    assert points == {("G1", r) for r in range(1, 9)}


def _glycam_atom(resname, resid, segid, chain, name, element, xyz):
    m = Molecule().empty(1)
    m.resname[:] = resname
    m.resid[:] = resid
    m.segid[:] = segid
    m.chain[:] = chain
    m.record[:] = "ATOM"
    m.name[:] = name
    m.element[:] = element
    m.coords = np.array(xyz, dtype=np.float32).reshape(1, 3, 1)
    return m


def _minimal_sugar_mol(resname, resid, segid, chain, x0=0.0):
    """A minimal 2-atom GLYCAM sugar residue (anomeric carbon + ring oxygen,
    at real bonding distance) - just enough for glycamUnitMask's
    composition-and-geometry gate to recognize it."""
    mol = Molecule()
    mol.append(
        _glycam_atom(resname, resid, segid, chain, "C1", "C", [x0, 0.0, 0.0])
    )
    mol.append(
        _glycam_atom(resname, resid, segid, chain, "O5", "O", [x0 + 1.4, 0.0, 0.0])
    )
    return mol


def test_glycan_break_points_isolates_first_glycan_from_preceding_protein():
    """A segid mixing protein and glycan residues is only a warning (see
    _checkMixedSegment), not an error, so amber.build must not rely on the
    caller always splitting them into separate segids. _glycan_break_points
    previously only emitted a break-AFTER point on each glycan residue's own
    resid, which isolates it from the NEXT residue but not from the
    PREVIOUS one - so the first glycan residue of a run landed in the same
    tLeap chain (no TER) as the last protein residue preceding it, and
    tLeap would auto-bond the two together across a bond that must instead
    come exclusively from the explicit `bond` commands."""
    from htmd.builder.amber import _glycan_break_points, _apply_chain_breaks

    mol = Molecule()
    for r in (1, 2, 3):
        mol.append(_ala_mol(r, "P", "A", x0=r * 5.0))
    mol.append(_minimal_sugar_mol("4YB", 4, "P", "A", x0=20.0))
    mol.append(_minimal_sugar_mol("0YB", 5, "P", "A", x0=25.0))

    points = _glycan_break_points(mol)
    assert points == {("P", 3), ("P", 4), ("P", 5)}

    out = _apply_chain_breaks(mol, points, alternate_segids={"P"})
    chain_of = {r: str(out.chain[out.resid == r][0]) for r in (1, 2, 3, 4, 5)}
    assert chain_of[1] == chain_of[2] == chain_of[3]
    assert chain_of[3] != chain_of[4], (
        "the first glycan residue must not share a chain with the "
        "preceding protein residue"
    )
    assert chain_of[4] != chain_of[5]


def test_glycan_break_points_two_adjacent_segid_glycan_trees():
    """Two independent glycan trees under different segids, adjacent in file
    order with no intervening residue, must each keep their own internal
    isolation working (this coverage gap was previously deferred). A segid
    boundary already forces a TER via the PDB writer regardless of chain
    ID, so reusing the same two alternating chain letters across segids
    (see _apply_chain_breaks) must not interfere with either tree's own
    break points."""
    from htmd.builder.amber import _glycan_break_points, _apply_chain_breaks

    mol = Molecule()
    mol.append(_minimal_sugar_mol("4YB", 1, "G1", "C", x0=0.0))
    mol.append(_minimal_sugar_mol("0YB", 2, "G1", "C", x0=5.0))
    mol.append(_minimal_sugar_mol("4YB", 1, "G2", "D", x0=10.0))
    mol.append(_minimal_sugar_mol("0YB", 2, "G2", "D", x0=15.0))

    points = _glycan_break_points(mol)
    assert points == {("G1", 1), ("G1", 2), ("G2", 1), ("G2", 2)}

    out = _apply_chain_breaks(mol, points, alternate_segids={"G1", "G2"})
    c1 = str(out.chain[(out.segid == "G1") & (out.resid == 1)][0])
    c2 = str(out.chain[(out.segid == "G1") & (out.resid == 2)][0])
    d1 = str(out.chain[(out.segid == "G2") & (out.resid == 1)][0])
    d2 = str(out.chain[(out.segid == "G2") & (out.resid == 2)][0])
    assert c1 != c2, "G1's own two glycan residues must still alternate chains"
    assert d1 != d2, "G2's own two glycan residues must still alternate chains"


def test_apply_chain_breaks_alternating():
    from htmd.builder.amber import _apply_chain_breaks

    # 70 single-residue breaks: with alternate_segids only TWO chain IDs are
    # consumed, no matter how many glycan residues exist
    mol = Molecule().empty(70)
    mol.segid[:] = "G1"
    mol.chain[:] = "X"
    mol.resid[:] = np.arange(1, 71)
    mol.resname[:] = "0YB"
    mol.name[:] = "C1"
    mol.element[:] = "C"
    mol.coords = np.zeros((70, 3, 1), dtype=np.float32)
    points = {("G1", int(r)) for r in range(1, 71)}
    out = _apply_chain_breaks(mol, points, alternate_segids={"G1"})
    chains = out.chain
    # every consecutive residue pair differs in chain, using exactly 2 IDs
    assert all(chains[i] != chains[i + 1] for i in range(69))
    assert len(set(chains)) == 2
    assert "X" not in set(chains)  # picked from the free pool


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_glycam_tleap_script(tmp_path):
    # end-to-end script generation without running tleap, but glycan detection
    # still needs AMBERHOME to locate the GLYCAM prep files
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "3AVE_frag.pdb")
    )
    pmol, _specs = systemPrepare(mol)
    build(pmol, outdir=str(tmp_path), ionize=False, execute=False)
    tleap_in = (tmp_path / "tleap.in").read_text()
    source_lines = [ln for ln in tleap_in.splitlines() if ln.startswith("source")]
    assert any("GLYCAM_06j-1" in ln for ln in source_lines)
    # inserted before the water leaprc
    gpos = next(i for i, ln in enumerate(source_lines) if "GLYCAM_06j-1" in ln)
    wpos = next(i for i, ln in enumerate(source_lines) if "water" in ln)
    assert gpos < wpos
    # the anchor bond and at least the branch bonds are explicit
    bond_lines = [ln for ln in tleap_in.splitlines() if ln.startswith("bond ")]
    assert any(".ND2" in ln and ".C1" in ln for ln in bond_lines)
    assert sum(1 for ln in bond_lines if ".O6" in ln or ".O4" in ln or ".O3" in ln or ".O2" in ln) >= 7
    # consecutive glycan residues alternate chain IDs in input.pdb
    input_pdb = (tmp_path / "input.pdb").read_text()
    per_residue = []  # (chain, resid) in file order, one entry per residue
    for ln in input_pdb.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() in (
            "UYB", "4YB", "VMB", "2MA", "0YB", "0fA",
        ):
            entry = (ln[21], ln[22:26].strip())
            if not per_residue or per_residue[-1] != entry:
                per_residue.append(entry)
    assert len(per_residue) == 8
    assert all(
        per_residue[i][0] != per_residue[i + 1][0]
        for i in range(len(per_residue) - 1)
    ), "consecutive glycan residues must differ in chain ID"
    assert len({c for c, _ in per_residue}) == 2, "two alternating chain IDs"


def _mixed_cyclic_mol(closure_dist=1.32):
    """3-residue cyclic segment whose middle residue is non-canonical and
    fails the 'protein' selection (as microcystin's beta-amino-acid Adda does),
    with an explicit first-N to last-C closure bond. The flanking residues are
    canonical."""
    mol = Molecule()
    mol.append(_ala_mol(1, "P0", "A", x0=0.0))
    lig = Molecule().empty(3)
    lig.name[:] = ["C1", "C2", "O1"]
    lig.element[:] = ["C", "C", "O"]
    lig.resname[:] = "LIG"
    lig.resid[:] = 2
    lig.chain[:] = "A"
    lig.segid[:] = "P0"
    lig.record[:] = "HETATM"
    lig.coords = np.array(
        [[5.0, 0, 0], [6.5, 0, 0], [7.0, 1, 0]], np.float32
    ).reshape(3, 3, 1)
    mol.append(lig)
    mol.append(_ala_mol(3, "P0", "A", x0=10.0))
    mol.guessBonds()
    n_first = int(np.where((mol.resid == 1) & (mol.name == "N"))[0][0])
    c_last = int(np.where((mol.resid == 3) & (mol.name == "C"))[0][0])
    mol.coords[n_first, :, 0] = mol.coords[c_last, :, 0] + np.array(
        [closure_dist, 0.0, 0.0], np.float32
    )
    mol.bonds = np.vstack([mol.bonds, [[n_first, c_last]]]).astype(mol.bonds.dtype)
    if mol.bondtype is not None and len(mol.bondtype):
        mol.bondtype = np.append(mol.bondtype, "1")
    return mol


def test_detect_cyclic_segments_allows_noncanonical_residues():
    """A cyclic segment containing a non-canonical residue that fails the
    'protein' selection (e.g. microcystin's beta-amino-acid Adda, 1FJM) must
    still be detected as cyclic. The old gate required the WHOLE segment to be
    canonical protein, which silently dropped such macrocycles."""
    from htmd.builder.amber import _detect_cyclic_segments

    mol = _mixed_cyclic_mol(closure_dist=1.32)
    cyclic = _detect_cyclic_segments(mol)
    assert cyclic == [("P0", 1, 3)]


def test_custombond_break_points_consecutive_only():
    """Only a custombond joining CONSECUTIVE residues is a break-after point.
    A custombond between non-consecutive residues (e.g. a head-to-tail closure
    or a long-range crosslink) does not break auto-sequencing."""
    from htmd.builder.amber import _custombond_break_points
    from moleculekit.molecule import UniqueAtomID

    mol = Molecule()
    for r in (1, 2, 3):
        mol.append(_ala_mol(r, "P0", "A", x0=r * 5.0))
    consec = [
        UniqueAtomID.fromMolecule(mol, "resid 1 and name C"),
        UniqueAtomID.fromMolecule(mol, "resid 2 and name N"),
    ]
    nonconsec = [
        UniqueAtomID.fromMolecule(mol, "resid 1 and name N"),
        UniqueAtomID.fromMolecule(mol, "resid 3 and name C"),
    ]
    pts = _custombond_break_points(mol, [consec, nonconsec])
    assert pts == {("P0", 1)}


def test_apply_chain_breaks_splits_chain_at_break():
    """A break-after point puts the residues before and after it in different
    chains, so tLeap writes a TER and won't head-to-tail auto-sequence them.
    Resids and atom order are unchanged."""
    from htmd.builder.amber import _apply_chain_breaks

    mol = Molecule()
    for r in (1, 2, 3, 4):
        mol.append(_ala_mol(r, "P0", "A", x0=r * 5.0))
    before = mol.resid.copy()
    out = _apply_chain_breaks(mol, {("P0", 2)})
    ch = {r: str(out.chain[out.resid == r][0]) for r in (1, 2, 3, 4)}
    assert ch[1] == ch[2]  # before the break: same chain
    assert ch[3] == ch[4]  # after the break: same chain
    assert ch[2] != ch[3]  # break between 2 and 3 -> distinct chains (TER)
    assert np.array_equal(out.resid, before)  # resids untouched


def test_apply_chain_breaks_noop_without_breaks():
    """No break points -> the molecule is returned unchanged."""
    from htmd.builder.amber import _apply_chain_breaks

    mol = Molecule()
    for r in (1, 2, 3):
        mol.append(_ala_mol(r, "P0", "A", x0=r * 5.0))
    out = _apply_chain_breaks(mol, set())
    assert np.array_equal(out.chain, mol.chain)


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_custombond_directive_matches_combined_unit_position(tmp_path):
    """Regression test: when waters sit between two solute residues in
    the htmd-internal mol order, the bond directive must reference the
    solute residue's position in input.pdb (its load order in the
    combined unit), not its renumbered PDB resid.

    This is the bug that broke the 8QFZ bicyclic-peptide tutorial after
    we stopped stripping waters before the build.
    """
    # Synthesize: ALA(P0/1) + 2 waters (W0/2,3) + ALA(P1/4). The
    # _prepareMolecule renumber will turn these into resids 1, 2, 3, 4
    # in mol order. After splitting waters out to solvent.pdb, the
    # solute-only input.pdb has only 2 ALA residues - so the load-order
    # position of the second ALA is 2, NOT 4.
    mol = Molecule()
    mol.append(_ala_mol(1, "P0", "A", x0=0.0))
    mol.append(_hoh_mol(10, "W0", x0=5.0))
    mol.append(_hoh_mol(11, "W0", x0=7.0))
    mol.append(_ala_mol(50, "P1", "B", x0=20.0))

    custombonds = [
        # User selection strings against the input mol's PDB resids.
        (
            'segid "P0" and chain "A" and resid 1 and name "CA"',
            'segid "P1" and chain "B" and resid 50 and name "CA"',
        ),
    ]

    # execute=False -> write tleap.in without running tleap.
    build(
        mol,
        outdir=str(tmp_path),
        ionize=False,
        custombonds=custombonds,
        execute=False,
    )

    tleap_in = (tmp_path / "tleap.in").read_text()

    # The bond must reference mol.<position_in_combined_unit>, not
    # mol.<htmd-renumbered-resid>. With 2 solute + 2 waters, the second
    # ALA is the 2nd residue of solute -> position 2 in the combined
    # unit. mol.4 would point at the second water, which has no CA atom.
    assert "bond mol.1.CA mol.2.CA" in tleap_in, (
        "custombond directive must use combined-unit residue position 2 "
        "for the post-water solute residue, not its htmd-renumbered resid 4. "
        f"Got:\n{tleap_in}"
    )
    assert "mol.4.CA" not in tleap_in, (
        "Found mol.4.CA in tleap.in - that would point at a water residue, "
        "not the second solute residue. The bond-directive writer is using "
        f"the wrong residue index. tleap.in:\n{tleap_in}"
    )


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_built_molecule_has_box_from_crd(tmp_path):
    """`amber.build()` must propagate the box dimensions tLeap wrote into the
    `.crd` onto the returned Molecule. Without this, downstream consumers
    (e.g. ``acemd.protocols.get_cellular_restraints``) see ``mol.box == 0``
    and fail even though periodicity information was available all along.
    """
    from htmd.builder.solvate import solvate

    np.random.seed(1)
    mol = Molecule("3PTB")
    mol.filter("protein")
    smol = solvate(mol)
    built = build(smol, outdir=str(tmp_path), ionize=False)
    assert not np.all(built.box == 0), (
        f"amber.build() returned a Molecule with box={built.box[:, 0]}; "
        "expected the periodic box tLeap wrote into structure.crd."
    )


def _synthetic_protein(resnames, segid="P0"):
    """Build a single-segment polypeptide with N/CA/C/O backbone atoms per
    residue. After guessBonds it registers as 'protein', and the first/last
    entry of ``resnames`` controls the corresponding terminus."""
    per = ["N", "CA", "C", "O"]
    el = ["N", "C", "C", "O"]
    n = len(resnames) * 4
    mol = Molecule().empty(n)
    names, elems, resids, rn, xyz = [], [], [], [], []
    x = 0.0
    for i, resname in enumerate(resnames):
        for a, e in zip(per, el):
            names.append(a)
            elems.append(e)
            resids.append(i + 1)
            rn.append(resname)
            xyz.append([x, 0.0, 0.0])
            x += 1.4
    mol.name[:] = names
    mol.element[:] = elems
    mol.resid[:] = resids
    mol.resname[:] = rn
    mol.segid[:] = segid
    mol.chain[:] = "A"
    mol.record[:] = "ATOM"
    mol.coords = np.array(xyz, dtype=np.float32).reshape(n, 3, 1)
    mol.guessBonds()
    return mol


def test_default_protein_caps_canonical_termini():
    from htmd.builder.amber import _defaultProteinCaps

    mol = _synthetic_protein(["ALA"] * 12)
    assert _defaultProteinCaps(mol) == {"P0": ["ACE", "NME"]}


def test_default_protein_caps_leaves_noncanonical_terminus_uncapped():
    # A non-canonical terminal residue must not be capped: parameterizeFromSpecs
    # builds it in its terminal form, so tleap must not add ACE/NME there.
    from htmd.builder.amber import _defaultProteinCaps

    mol = _synthetic_protein(["ALA"] * 11 + ["LIG"])
    assert _defaultProteinCaps(mol) == {"P0": ["ACE", "none"]}


def test_default_protein_caps_protonation_variant_terminus_is_capped():
    # HIE / CYX are protonation/naming variants systemPrepare produces; they are
    # standard protein residues and must still be capped.
    from htmd.builder.amber import _defaultProteinCaps

    mol = _synthetic_protein(["CYX"] + ["ALA"] * 10 + ["HIE"])
    assert _defaultProteinCaps(mol) == {"P0": ["ACE", "NME"]}


def test_default_protein_caps_existing_cap_not_recapped():
    from htmd.builder.amber import _defaultProteinCaps

    mol = _synthetic_protein(["ALA"] * 11 + ["NME"])
    assert _defaultProteinCaps(mol) == {"P0": ["ACE", "none"]}


def test_default_protein_caps_is_public():
    from moleculekit.molecule import Molecule
    from moleculekit.tools.autosegment import autoSegment
    from htmd.builder.amber import _defaultProteinCaps, defaultProteinCaps

    mol = Molecule("3ptb")
    mol.filter("protein", _logger=False)
    mol = autoSegment(mol, fields=("segid",), _logger=False)

    caps = defaultProteinCaps(mol)
    assert caps == _defaultProteinCaps(mol)
    assert list(caps.values())[0] == ["ACE", "NME"]


def _assert_glycan_bond(molbuilt, resname_a, name_a, resname_b, name_b):
    idx_a = np.where((molbuilt.resname == resname_a) & (molbuilt.name == name_a))[0]
    idx_b = np.where((molbuilt.resname == resname_b) & (molbuilt.name == name_b))[0]
    bset = {tuple(sorted(b)) for b in molbuilt.bonds.tolist()}
    found = any(
        tuple(sorted((int(a), int(b)))) in bset for a in idx_a for b in idx_b
    )
    assert found, f"missing bond {resname_a}.{name_a} - {resname_b}.{name_b}"


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_glycam_build_3ave_branched_nglycan(tmp_path):
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "3AVE_frag.pdb")
    )
    pmol, _ = systemPrepare(mol)
    molbuilt = build(pmol, outdir=str(tmp_path), ionize=False)
    assert molbuilt is not None
    for code in ("NLN", "UYB", "4YB", "VMB", "2MA", "0YB", "0fA"):
        assert code in molbuilt.resname, code
    _assert_glycan_bond(molbuilt, "NLN", "ND2", "UYB", "C1")
    _assert_glycan_bond(molbuilt, "UYB", "O6", "0fA", "C1")
    _assert_glycan_bond(molbuilt, "VMB", "O3", "2MA", "C1")
    _assert_glycan_bond(molbuilt, "VMB", "O6", "2MA", "C1")
    total_charge = float(np.sum(molbuilt.charge))
    assert abs(total_charge - round(total_charge)) < 1e-3


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_glycam_build_1cvn_free_oligosaccharide(tmp_path):
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "1CVN_frag.pdb")
    )
    pmol, _ = systemPrepare(mol)
    molbuilt = build(pmol, outdir=str(tmp_path), ionize=False)
    assert "ROH" in molbuilt.resname and "VMA" in molbuilt.resname
    _assert_glycan_bond(molbuilt, "ROH", "O1", "VMA", "C1")
    _assert_glycan_bond(molbuilt, "VMA", "O3", "0MA", "C1")
    _assert_glycan_bond(molbuilt, "VMA", "O6", "0MA", "C1")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def test_glycam_build_1g1s_olinked(tmp_path):
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build

    mol = Molecule(
        os.path.join(curr_dir, "data", "test-amber-build", "glycans", "1G1S_frag.pdb")
    )
    pmol, _ = systemPrepare(mol)
    molbuilt = build(pmol, outdir=str(tmp_path), ionize=False)
    for code in ("OLT", "VVB", "WYB", "3LB", "0SA", "0fA", "0LB"):
        assert code in molbuilt.resname, code
    _assert_glycan_bond(molbuilt, "OLT", "OG1", "VVB", "C1")
    _assert_glycan_bond(molbuilt, "3LB", "O3", "0SA", "C2")
    total_charge = float(np.sum(molbuilt.charge))
    assert abs(total_charge - round(total_charge)) < 1e-3  # SIA contributes -1
