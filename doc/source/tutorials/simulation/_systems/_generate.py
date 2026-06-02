"""Generate vendored pre-built systems for the simulation tutorials.

Run this once whenever the underlying builders change:

    uv --project /home/sdoerr/Work/htmd run python \\
        doc/source/tutorials/simulation/_systems/_generate.py

The output ``trpcage.prmtop`` / ``trpcage.pdb`` and ``apelin-membrane.prmtop`` /
``apelin-membrane.pdb`` pairs are checked into the repo so the simulation
tutorials don't have to repeat the full system-building flow on every docs
build. The membrane build in particular takes ~30 s on CPU.
"""
import shutil
import tempfile
from pathlib import Path

from moleculekit.molecule import Molecule
from moleculekit.opm import get_opm_pdb
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.nonstandard import parameterizeFromSpecs
from htmd.builder.solvate import solvate
from htmd.membranebuilder.build_membrane import buildMembrane


HERE = Path(__file__).resolve().parent


def build_trpcage() -> None:
    print("Building Trp-cage (1L2Y) ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir) / "build"
        mol = Molecule("1L2Y")
        mol = autoSegment(mol, fields=("segid", "chain"))
        prepared, _ = systemPrepare(mol, pH=7.4)
        solvated = solvate(prepared, pad=12)
        built = amber.build(
            solvated, outdir=str(outdir), ionize=True, saltconc=0.15
        )
        built.write(str(HERE / "trpcage.pdb"))
        shutil.copy(outdir / "structure.prmtop", HERE / "trpcage.prmtop")
        print("  wrote trpcage.pdb / trpcage.prmtop")


def build_apelin_membrane() -> None:
    print("Building apelin receptor + membrane (5VBL) ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        builddir = Path(tmpdir) / "build"
        paramsdir = Path(tmpdir) / "params"
        membdir = Path(tmpdir) / "memb"

        mol = Molecule("5VBL")
        mol.remove("water", _logger=False)
        ref, thickness = get_opm_pdb("5VBL", validateElements=False)
        mol.align("protein and name CA", refmol=ref, mode="structure")
        mol = autoSegment(mol, fields=("segid", "chain"))

        specs = detectNonStandardResidues(mol)
        smiles = {
            "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
            "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
            "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
            "NLE": "CCCC[C@@H](C=O)N",
            "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
            "OLC": "CCCCCCCC(=O)OC[C@H](O)CO",
        }
        for resname, smi in smiles.items():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )

        prepared, specs = systemPrepare(mol, pH=7.4, detect_specs=specs)
        out = parameterizeFromSpecs(
            specs, prepared, outdir=str(paramsdir), charge_method="gasteiger"
        )

        memb = buildMembrane(
            [80, 80],
            ratioupper={"popc": 0.7, "chl1": 0.3},
            ratiolower={"popc": 0.7, "chl1": 0.3},
            minimize=0, equilibrate=0, platform="CPU",
            outdir=str(membdir), solute=prepared,
        )
        memb.remove("water", _logger=False)

        system = Molecule()
        system.append(prepared)
        system.append(memb)

        lipid_mask = system.atomselect("lipid")
        xy_min = system.coords[lipid_mask, :2, 0].min(axis=0)
        xy_max = system.coords[lipid_mask, :2, 0].max(axis=0)
        z_min = system.coords[:, 2, 0].min() - 15
        z_max = system.coords[:, 2, 0].max() + 15
        system = solvate(
            system,
            minmax=[
                [xy_min[0], xy_min[1], z_min],
                [xy_max[0], xy_max[1], z_max],
            ],
        )

        built = amber.build(
            system, outdir=str(builddir),
            custombonds=out.custombonds, topo=out.topo_paths,
            param=out.frcmod_paths,
            caps={"P0": ("none", "none")}, ionize=True, saltconc=0.15,
        )
        built.write(str(HERE / "apelin-membrane.pdb"))
        shutil.copy(builddir / "structure.prmtop", HERE / "apelin-membrane.prmtop")
        print("  wrote apelin-membrane.pdb / apelin-membrane.prmtop")


def main() -> None:
    build_trpcage()
    build_apelin_membrane()


if __name__ == "__main__":
    main()
