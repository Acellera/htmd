import os

import numpy as np

from moleculekit.molecule import Molecule

from htmd.builder.solvate import solvate


curr_dir = os.path.dirname(os.path.abspath(__file__))


def _water_triples(mol):
    """Return the (nwat, 3, 3) coordinate array of the solvate-added waters.

    Waters are emitted as consecutive O, H1, H2 triples.
    """
    watidx = np.where(mol.atomselect("water"))[0]
    assert watidx.size % 3 == 0, "water atom count is not a multiple of 3"
    return mol.coords[watidx, :, 0].reshape(-1, 3, 3), watidx


def _solute():
    """A tiny solute so solvate() has something to build a box around."""
    mol = Molecule().empty(2)
    mol.name[:] = ["C1", "C2"]
    mol.element[:] = ["C", "C"]
    mol.resname[:] = "LIG"
    mol.resid[:] = [1, 1]
    mol.segid[:] = "S0"
    mol.coords = np.array(
        [[[0.0], [0.0], [0.0]], [[1.5], [0.0], [0.0]]], dtype=np.float32
    )
    return mol


def test_solvate_returns_bonded_water():
    """Every solvate-added water must carry its two O-H bonds.

    Without them, anything that groups atoms by connectivity (notably
    Molecule.wrap) treats each water atom as an independent molecule.
    """
    smol = solvate(_solute(), pad=8)

    _, watidx = _water_triples(smol)
    nwat = watidx.size // 3
    assert nwat > 0, "solvate added no water"

    watset = set(watidx.tolist())
    watbonds = smol.bonds[np.isin(smol.bonds, list(watset)).any(axis=1)]
    assert watbonds.shape[0] == 2 * nwat, (
        f"expected {2 * nwat} water bonds for {nwat} waters, got {watbonds.shape[0]}"
    )

    # No bond may join two different residues (i.e. two different waters)
    resids = smol.resid[watbonds]
    segids = smol.segid[watbonds]
    assert np.all(resids[:, 0] == resids[:, 1]), "water bond crosses residues"
    assert np.all(segids[:, 0] == segids[:, 1]), "water bond crosses segments"


def test_solvated_water_survives_wrapping():
    """Regression: waters straddling a periodic face must wrap as whole molecules.

    Bond-less water made Molecule.wrap() translate individual O/H atoms across
    the boundary, producing intramolecular "bonds" the length of a box edge.
    This is what corrupted membrane builds.
    """
    smol = solvate(_solute(), pad=8)

    coords = smol.coords[:, :, 0]
    smol.box = (coords.max(axis=0) - coords.min(axis=0)).astype(np.float32).reshape(3, 1)
    smol.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    # Wrap around a corner of the box so plenty of waters straddle a face.
    smol.wrap(wrapcenter=coords.min(axis=0))

    triples, _ = _water_triples(smol)
    oh1 = np.linalg.norm(triples[:, 1] - triples[:, 0], axis=1)
    oh2 = np.linalg.norm(triples[:, 2] - triples[:, 0], axis=1)
    worst = max(oh1.max(), oh2.max())
    assert worst < 1.5, f"water split across the periodic boundary: max O-H = {worst:.2f} A"


def test_custom_bondless_solvent_file_still_gets_bonds(tmp_path):
    """A user-supplied PDB water box has no bonds; solvate must supply them."""
    from htmd.home import home

    ref = Molecule(os.path.join(home(shareDir=True), "solvate", "wat.bcif.gz"))
    assert ref.bonds.shape[0] > 0, "the shipped water box should carry bonds"

    # A hand-made water box PDB carries no CONECT records, which is the case
    # a user hits when passing their own solvent file.
    ref.bonds = np.empty((0, 2), dtype=np.uint32)
    ref.bondtype = np.empty(0, dtype=object)
    custom = str(tmp_path / "mywat.pdb")
    ref.write(custom)
    assert Molecule(custom).bonds.shape[0] == 0, "fixture is expected to be bond-less"

    smol = solvate(_solute(), pad=8, spdb=custom)

    _, watidx = _water_triples(smol)
    nwat = watidx.size // 3
    watbonds = smol.bonds[np.isin(smol.bonds, watidx.tolist()).any(axis=1)]
    assert watbonds.shape[0] == 2 * nwat, (
        f"expected {2 * nwat} water bonds, got {watbonds.shape[0]}"
    )
