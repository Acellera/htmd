import unittest

from htmd.builder.charmm import build, _psfgen_exists


class TestCharmmBuild(unittest.TestCase):
    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
    def test_build(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdbids = ["3PTB", "1A25", "1GZM", "1U5U"]
        for pdb in pdbids:
            with self.subTest(pdb=pdb):
                print("Building {}".format(pdb))
                inFile = os.path.join(
                    preparedInputDir, pdb, "{}-prepared.pdb".format(pdb)
                )
                mol = Molecule(inFile)
                mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

                np.random.seed(1)  # Needed for ions
                smol = solvate(mol)
                topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]
                # 1U5U contains deprotonated arginines (AR0 → patch RN1),
                # which are only defined in the arg0 stream file.
                if pdb == "1U5U":
                    topos += [
                        "str/prot/toppar_all36_prot_arg0.str",
                        "str/misc/toppar_ions_won.str",
                    ]
                tmpdir = tempname()
                _ = build(smol, topo=topos, outdir=tmpdir)

                compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
                assertSameAsReferenceDir(compareDir, tmpdir)

                # shutil.rmtree(tmpdir)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
    def test_customDisulfideBonds(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdb = "1GZM"
        inFile = os.path.join(preparedInputDir, pdb, "{}-prepared.pdb".format(pdb))
        mol = Molecule(inFile)
        mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

        np.random.seed(1)  # Needed for ions
        smol = solvate(mol)
        topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]
        disu = [
            ["segid 1 and resid 110", "segid 1 and resid 187"],
            ["segid 0 and resid 110", "segid 0 and resid 187"],
        ]
        tmpdir = tempname()
        _ = build(smol, topo=topos, outdir=tmpdir, disulfide=disu)

        compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
        assertSameAsReferenceDir(compareDir, tmpdir)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
    def test_disulfideWithInsertion(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdb = "3PTB"

        print("Building {}".format(pdb))
        inFile = os.path.join(preparedInputDir, pdb, "{}-prepared.pdb".format(pdb))
        mol = Molecule(inFile)
        mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

        np.random.seed(1)  # Needed for ions
        smol = solvate(mol)
        topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]

        smol.insertion[smol.resid == 42] = (
            "A"  # Adding an insertion to test that disulfide bonds with insertions work
        )
        tmpdir = tempname()
        _ = build(smol, topo=topos, outdir=tmpdir)
        compareDir = home(dataDir=os.path.join("test-charmm-build", "3PTB_insertion"))
        assertSameAsReferenceDir(compareDir, tmpdir)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
    def test_build_records_the_cell(self):
        """psfgen's PDB has no CRYST1, so charmm.build has to add it."""
        import os

        import numpy as np
        from moleculekit.molecule import Molecule

        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname

        inFile = os.path.join(
            home(dataDir="test-proteinprepare"), "3PTB", "3PTB-prepared.pdb"
        )
        mol = Molecule(inFile)
        mol.filter("protein")

        np.random.seed(1)
        smol = solvate(mol)

        tmpdir = tempname()
        molbuilt = build(
            smol,
            topo=["top/top_all36_prot.rtf", "top/top_water_ions.rtf"],
            outdir=tmpdir,
        )

        assert np.allclose(molbuilt.box.ravel(), smol.box.ravel(), atol=1e-2)
        assert np.allclose(molbuilt.boxangles.ravel(), [90.0, 90.0, 90.0], atol=1e-2)

        pdb = Molecule(os.path.join(tmpdir, "structure.pdb"))
        assert np.allclose(pdb.box.ravel(), smol.box.ravel(), atol=1e-2)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
    def test_build_records_a_non_rectangular_cell(self):
        """_insert_cryst1 must carry a non-rectangular cell, not just 90/90/90."""
        import os

        import numpy as np
        from moleculekit.molecule import Molecule

        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname

        inFile = os.path.join(
            home(dataDir="test-proteinprepare"), "3PTB", "3PTB-prepared.pdb"
        )
        mol = Molecule(inFile)
        mol.filter("protein")

        np.random.seed(1)
        smol = solvate(mol, pad=10, shape="octahedron")

        tmpdir = tempname()
        molbuilt = build(
            smol,
            topo=["top/top_all36_prot.rtf", "top/top_water_ions.rtf"],
            outdir=tmpdir,
        )

        assert np.allclose(molbuilt.box.ravel(), smol.box.ravel(), atol=1e-2)
        assert np.allclose(molbuilt.boxangles.ravel(), [109.4712206] * 3, atol=1e-2)

        pdb = Molecule(os.path.join(tmpdir, "structure.pdb"))
        assert np.allclose(pdb.box.ravel(), smol.box.ravel(), atol=1e-2)
        assert np.allclose(pdb.boxangles.ravel(), [109.4712206] * 3, atol=1e-2)
