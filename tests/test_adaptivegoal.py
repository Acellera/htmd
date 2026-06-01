import os
import unittest
from htmd.adaptive.adaptivegoal import AdaptiveGoal


class TestAdaptiveGoal(unittest.TestCase):
    def test_adaptive_goal(self):
        from moleculekit.projections.metricdistance import MetricDistance
        from jobqueues.localqueue import LocalCPUQueue
        from htmd.home import home
        import tempfile
        import shutil

        with tempfile.TemporaryDirectory() as tmpdir:
            print(tmpdir)

            gendir = os.path.join(
                home(dataDir="test-adaptive"), "test-ions", "generators"
            )

            shutil.copytree(gendir, os.path.join(tmpdir, "generators"))
            inpdir = os.path.join(tmpdir, "input")
            datdir = os.path.join(tmpdir, "data")

            dist = MetricDistance("name CL", "name NA", periodic="selections")

            app = LocalCPUQueue()
            app.datadir = datdir

            md = AdaptiveGoal()
            md.nmin = 1
            md.nmax = 2
            md.nepochs = 3
            md.updateperiod = 5
            md.projection = dist
            md.goalfunction = lambda mol: -dist.project(mol).flatten()
            md.ticadim = 0
            md.generatorspath = gendir
            md.inputpath = inpdir
            md.datapath = datdir
            md.app = app
            md.savegoal = os.path.join(tmpdir, "goals.pkl")
            # md.run()

            # assert os.path.exists(os.path.join(tmpdir, "goals.pkl"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
