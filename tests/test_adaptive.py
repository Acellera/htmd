import unittest
import os
from glob import glob

import numpy as np
from moleculekit.molecule import Molecule

from htmd.adaptive.adaptive import _writeInputsFunction


class TestAdaptive(unittest.TestCase):
    def test_input_writer(self):
        from htmd.home import home
        from htmd.simlist import Frame, simlist
        import tempfile

        filedir = home() + "/data/adaptive/"
        sims = simlist(
            glob(os.path.join(filedir, "data", "*", "")),
            glob(os.path.join(filedir, "input", "*", "")),
            glob(os.path.join(filedir, "input", "*", "")),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            f = Frame(sims[0], 0, 5)
            _writeInputsFunction(1, f, 2, tmpdir, "input.coor", "input.xsc")

            ref = Molecule(sims[0])
            ref.dropFrames(keep=5)
            mol = Molecule(sims[0].molfile)
            mol.read(os.path.join(tmpdir, "e2s2_e1s1p0f5", "input.coor"))
            mol.read(os.path.join(tmpdir, "e2s2_e1s1p0f5", "input.xsc"))

            assert np.array_equal(ref.coords, mol.coords)
            assert np.array_equal(ref.box, mol.box)
