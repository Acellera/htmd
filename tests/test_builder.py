import unittest
import numpy as np
from htmd.builder.builder import embed


class TestBuilder(unittest.TestCase):
    def test_embed(self):
        from moleculekit.molecule import Molecule, mol_equal
        from moleculekit.tools.autosegment import autoSegment
        from htmd.home import home
        from os import path

        testdir = path.join(home(), "data", "building-protein-membrane")

        p = Molecule(path.join(testdir, "4dkl.pdb"))
        p.filter("(chain B and protein) or water")
        p = autoSegment(p, "protein", "P")
        m = Molecule(path.join(testdir, "membrane.pdb"))
        a = embed(p, m)

        ref = Molecule(path.join(testdir, "embedded.pdb"))
        assert mol_equal(a, ref, exceptFields=("serial"))

    def test_lipid_inaccessible_points(self):
        """A sealed tube hides space from the bilayer; a slit one does not."""
        from htmd.builder.builder import lipid_inaccessible_points
        from moleculekit.molecule import Molecule

        def tube(open_arc_deg):
            # A cylinder of atoms at radius 10, optionally missing an arc.
            theta = np.linspace(0, 2 * np.pi, 42, endpoint=False)
            theta = theta[theta > np.deg2rad(open_arc_deg)]
            z = np.arange(-18, 18, 1.5)
            t, zz = np.meshgrid(theta, z, indexing="ij")
            xyz = np.stack(
                [10 * np.cos(t).ravel(), 10 * np.sin(t).ravel(), zz.ravel()], axis=1
            )
            mol = Molecule().empty(len(xyz))
            mol.coords = xyz[:, :, None].astype(np.float32)
            mol.element[:] = "C"
            mol.name[:] = "C"
            mol.resname[:] = "TUB"
            mol.resid[:] = np.arange(len(xyz))
            return mol

        sealed = lipid_inaccessible_points(tube(0), thickness=32)
        assert len(sealed), "a sealed tube should hide space from the bilayer"
        # Everything found sits on the tube axis, inside the wall.
        assert np.linalg.norm(sealed[:, :2], axis=1).max() < 10, sealed

        slit = lipid_inaccessible_points(tube(100), thickness=32)
        assert not len(slit), f"a tube open to the bilayer seals nothing, got {len(slit)}"

    def test_autosegment(self):
        from moleculekit.molecule import Molecule
        from moleculekit.tools.autosegment import autoSegment
        from htmd.home import home
        from os import path

        testdir = path.join(home(), "data", "building-protein-membrane")

        mol = Molecule(path.join(testdir, "1ITG_clean.pdb"))
        ref = Molecule(path.join(testdir, "1ITG.pdb"))
        mol = autoSegment(mol, sel="protein")
        assert np.all(mol.segid == ref.segid)

        mol = Molecule(path.join(testdir, "3PTB_clean.pdb"))
        ref = Molecule(path.join(testdir, "3PTB.pdb"))
        mol = autoSegment(mol, sel="protein")
        assert np.all(mol.segid == ref.segid)
