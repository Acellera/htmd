import unittest
import os.path as path
from glob import glob as glob
from htmd.simlist import _has_authoritative_bonds, simlist, Sim


class TestSimlist(unittest.TestCase):
    def test_has_authoritative_bonds(self):
        assert _has_authoritative_bonds("structure.psf")
        assert _has_authoritative_bonds("system.prmtop")
        assert _has_authoritative_bonds("system.parm7")
        assert _has_authoritative_bonds("ligand.mol2")
        assert _has_authoritative_bonds(["structure.pdb", "topology.psf"])
        assert _has_authoritative_bonds(["topology.prmtop", "coords.inpcrd"])
        assert not _has_authoritative_bonds("structure.pdb")
        assert not _has_authoritative_bonds("traj.xyz")
        assert not _has_authoritative_bonds(["only.pdb"])
        assert not _has_authoritative_bonds(["a.pdb", "b.coor"])

    def test_simlist_auto_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            glob(path.join(home(dataDir="adaptive"), "input", "*")),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert len(sims[0].molfile) == 2
        assert _singleMolfile(sims)[0]

    def test_simlist_many_structures(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            glob(path.join(home(dataDir="adaptive"), "input", "*", "structure.pdb")),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]

    def test_simlist_single_structure(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile

        sims = simlist(
            glob(path.join(home(dataDir="adaptive"), "data", "*", "")),
            path.join(home(dataDir="adaptive"), "input", "e1s1_1", "structure.pdb"),
        )
        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert not isinstance(sims[0].molfile, list)
        assert _singleMolfile(sims)[0]

    def test_autodetect(self):
        from htmd.home import home
        from htmd.projections.metric import _singleMolfile
        from natsort import natsorted

        datafolders = natsorted(
            glob(path.join(home(dataDir="adaptive"), "data", "*", ""))
        )
        inputfolders = natsorted(
            glob(path.join(home(dataDir="adaptive"), "input", "*"))
        )
        sims = []
        for ff in zip(datafolders, inputfolders):
            sims.append(Sim(ff[0], ff[1], ff[1]))

        x = sims[0].copy()
        assert x == sims[0]
        assert x != sims[1]
        assert len(sims[0].molfile) == 2
        assert _singleMolfile(sims)[0]

        sims2 = []
        for ff in zip(datafolders, inputfolders):
            sims2.append(Sim(ff[0], ff[1]))

        for s1, s2 in zip(sims, sims2):
            assert s1.parent == s2.parent
            assert s1.input == s2.input
            assert s1.trajectory == s2.trajectory
            assert s1.molfile == s2.molfile
            assert s1.numframes == s2.numframes

    def test_sim(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            with open(path.join(tmpdir, "test.pdb"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.psf"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.prmtop"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.mae"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.coor"), "w") as f:
                f.write("fake")
            for i in range(11):
                with open(path.join(tmpdir, f"test_{i}.xtc"), "w") as f:
                    f.write("fake")

            # Passing a single directory
            sim = Sim(tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.prmtop"),
                path.join(tmpdir, "test.pdb"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir

            # Passing same directory to all arguments
            sim = Sim(tmpdir, tmpdir, tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.prmtop"),
                path.join(tmpdir, "test.pdb"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir

            # Passing single topology file and arbitrary input folder
            sim = Sim(tmpdir, path.join(tmpdir, "test.pdb"), "/tmp")
            assert sim.molfile == [path.join(tmpdir, "test.pdb")]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == "/tmp"

        # Coor file reading
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(path.join(tmpdir, "test.psf"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.mae"), "w") as f:
                f.write("fake")
            with open(path.join(tmpdir, "test.coor"), "w") as f:
                f.write("fake")
            for i in range(11):
                with open(path.join(tmpdir, f"test_{i}.xtc"), "w") as f:
                    f.write("fake")

            sim = Sim(tmpdir)
            assert sim.molfile == [
                path.join(tmpdir, "test.psf"),
                path.join(tmpdir, "test.coor"),
            ]
            assert sim.trajectory == [
                path.join(tmpdir, f"test_{i}.xtc") for i in range(11)
            ]
            assert sim.input == tmpdir
