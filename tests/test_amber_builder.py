from moleculekit.molecule import Molecule
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
try:
    _findTeLeap()
    tleap_installed = True
except Exception:
    tleap_installed = False


curr_dir = os.path.dirname(os.path.abspath(__file__))


def _compareResultFolders(
    compare, tmpdir, pid, ignore_ftypes=(".log", ".txt", ".frcmod")
):
    import filecmp

    def _cutfirstline(infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, "r") as fin:
            data = fin.read().splitlines(True)
        with open(outfile, "w") as fout:
            fout.writelines(data[1:])

    try:
        from ffevaluation.ffevaluate import FFEvaluate, loadParameters

        prm2 = loadParameters(os.path.join(compare, "structure.prmtop"))
        mol2 = Molecule(os.path.join(compare, "structure.prmtop"))
        mol2.read(os.path.join(compare, "structure.pdb"))
        ffev2 = FFEvaluate(mol2, prm2)
        energies2, _, _ = ffev2.calculate(mol2.coords)

        prm = loadParameters(os.path.join(tmpdir, "structure.prmtop"))
        mol = Molecule(os.path.join(tmpdir, "structure.prmtop"))
        mol.read(os.path.join(tmpdir, "structure.pdb"))
        ffev = FFEvaluate(mol, prm)
        energies, _, _ = ffev.calculate(mol.coords)
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
def _test_with_protein_prepare(tmp_path, pid):
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.solvate import solvate
    from moleculekit.tools.autosegment import autoSegment

    # Test with systemPrepare
    # pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)

    np.random.seed(1)
    mol = Molecule(pid)
    mol.filter("protein")
    mol = systemPrepare(mol, pH=7.0)
    mol.filter("protein")  # Fix for bad systemPrepare hydrogen placing
    mol = autoSegment(mol)
    smol = solvate(mol)
    ffs = defaultFf()
    _ = build(smol, ff=ffs, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "pp", pid)
    _compareResultFolders(refdir, tmp_path, pid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.parametrize("pid", ["3PTB"])
def _test_without_protein_prepare(tmp_path, pid):
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
def _test_protein_ligand(tmp_path):
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
def _test_custom_disulfide_bonds(tmp_path):
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
def _test_custom_teleap_imports(tmp_path):
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
def _test_rna(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)

    mol = Molecule("6VA1")
    smol = solvate(mol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "rna", "6VA1")
    _compareResultFolders(refdir, tmp_path, "6VA1")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_dna(tmp_path):
    from htmd.builder.solvate import solvate

    np.random.seed(1)

    mol = Molecule("1BNA")
    smol = solvate(mol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "dna", "1BNA")
    _compareResultFolders(refdir, tmp_path, "1BNA")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_protein_rna(tmp_path):
    from htmd.builder.solvate import solvate
    from moleculekit.tools.preparation import systemPrepare
    from moleculekit.tools.autosegment import autoSegment

    np.random.seed(1)

    mol = Molecule("3WBM")
    mol.filter("not water")
    mol = autoSegment(mol, field="both")
    pmol = systemPrepare(mol, pH=7.0)
    smol = solvate(pmol)

    _ = build(smol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "protein-rna", "3WBM")
    _compareResultFolders(refdir, tmp_path, "3WBM")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_caps(tmp_path):
    from htmd.builder.solvate import solvate
    from moleculekit.tools.preparation import systemPrepare

    np.random.seed(1)

    mol = Molecule("6A5J")
    pmol = systemPrepare(mol, pH=7.0)
    smol = solvate(pmol)

    outdir = os.path.join(tmp_path, "out")
    _ = build(smol, outdir=outdir, ionize=False)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "peptide-cap", "6A5J")
    _compareResultFolders(refdir, outdir, "6A5J")
    shutil.rmtree(outdir)

    np.random.seed(1)

    mol = Molecule("6A5J")
    pmol = systemPrepare(mol, pH=7.0)
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
@pytest.mark.parametrize("pdbid", ["5VBL", "1AWF"])
def _test_non_standard_residue_building(tmp_path, pdbid):
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
    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "non-standard", pdbid, "build"
    )
    _compareResultFolders(refdir, tmp_path, pdbid)


@pytest.mark.skipif(not tleap_installed, reason=reason)
@pytest.mark.skipif(sys.platform.startswith("darwin"), reason="Fails on OSX")
def _test_cofactor_building(tmp_path):
    homedir = os.path.join(curr_dir, "data", "test-amber-build", "cofactors")

    mol = Molecule(os.path.join(homedir, "cofactors.pdb"))
    _ = build(mol, ionize=False, outdir=tmp_path)
    refdir = os.path.join(curr_dir, "data", "test-amber-build", "cofactors", "build")

    _compareResultFolders(refdir, tmp_path, "Cofactors")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_post_translational_modifications_building(tmp_path):
    homedir = os.path.join(curr_dir, "data", "test-amber-build", "post-translational")

    mol = Molecule(os.path.join(homedir, "4EFP_nolig.pdb"))
    _ = build(mol, ionize=False, outdir=tmp_path)
    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "post-translational", "build"
    )

    _compareResultFolders(refdir, tmp_path, "post-translational modifications")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_cyclic_peptide_building(tmp_path):
    np.random.seed(1)

    mol = Molecule("5VAV")
    mol.segid = mol.chain
    _ = build(mol, outdir=tmp_path, ionize=False, caps={"A": ("none", "none")})

    refdir = os.path.join(
        curr_dir, "data", "test-amber-build", "cyclic-peptide", "5VAV"
    )
    _compareResultFolders(refdir, tmp_path, "5VAV")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_membrane(tmp_path):
    np.random.seed(1)

    homedir = os.path.join(curr_dir, "data", "test-amber-build", "membrane")

    mol = Molecule(os.path.join(homedir, "structure.psf"))
    mol.read(os.path.join(homedir, "structure.pdb"))
    _ = build(mol, outdir=tmp_path)

    refdir = os.path.join(curr_dir, "data", "test-amber-build", "membrane", "build")
    _compareResultFolders(refdir, tmp_path, "Membrane")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_cif_building(tmp_path):
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
def _test_atomtype_decollisioning(tmp_path):
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
