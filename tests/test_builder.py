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
