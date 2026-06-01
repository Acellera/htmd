import unittest
import numpy as np
from htmd.metricdata import MetricData


class TestMetricData(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.simlist import simlist, simfilter
        from glob import glob
        from htmd.projections.metric import Metric
        from moleculekit.projections.metricdistance import MetricDistance
        from moleculekit.projections.metricdihedral import MetricDihedral
        from moleculekit.util import tempname
        from htmd.home import home
        from os.path import join

        sims = simlist(
            glob(join(home(dataDir="adaptive"), "data", "*", "")),
            glob(join(home(dataDir="adaptive"), "input", "*")),
        )
        fsims = simfilter(sims, tempname(), "not water")

        metr = Metric(fsims)
        metr.set(
            MetricDistance(
                "protein and resid 10 and name CA",
                "resname BEN and noh",
                periodic="selections",
                metric="contacts",
                groupsel1="residue",
                threshold=4,
            )
        )
        self.data1 = metr.project()

        metr.set(MetricDihedral())
        self.data2 = metr.project()

    def test_combine(self):
        # Testing combining of metrics
        data1 = self.data1.copy()
        data1.combine(self.data2)

        # Testing dimensions
        assert np.array_equal(
            data1.description.shape, (897, 3)
        ), "combine not working correct"
        assert np.array_equal(
            data1.trajectories[0].projection.shape, (6, 897)
        ), "combine not working correct"
        assert np.array_equal(
            np.where(data1.description.type == "contact")[0],
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
        ), "combine not working correct"

    def test_dropping(self):
        # Testing dimension dropping / keeping
        data1 = self.data1.copy()
        assert np.array_equal(data1.description.shape, (9, 3))
        data1.dropDimensions(range(9))
        assert np.array_equal(
            data1.description.shape, (0, 3)
        ), "dropDimensions not working correct"
        assert np.array_equal(
            data1.trajectories[0].projection.shape, (6, 0)
        ), "dropDimensions not working correct"
        assert (
            len(np.where(data1.description.type == "contact")[0]) == 0
        ), "dropDimensions not working correct"

        data2 = self.data2.copy()
        assert np.array_equal(data2.description.shape, (888, 3))
        data2.dropDimensions(keep=range(9))
        assert np.array_equal(
            data2.description.shape, (9, 3)
        ), "dropDimensions not working correct"
        assert np.array_equal(
            data2.trajectories[0].projection.shape, (6, 9)
        ), "dropDimensions not working correct"
        assert (
            len(np.where(data2.description.type == "dihedral")[0]) == 9
        ), "dropDimensions not working correct"

    def test_saving_loading(self):
        from moleculekit.util import tempname

        def checkCorrectness(newdata):
            assert newdata.numTrajectories == 2, "Failed to load trajectories"
            assert newdata.description.shape == (9, 3), "Failed to load pandas data"
            assert newdata.trajectories[0].projection.shape == (6, 9), "Wrong data size"

        savefile = tempname(suffix=".dat")
        self.data1.save(savefile)

        newdata = MetricData(file=savefile)
        checkCorrectness(newdata)

        newdata = MetricData()
        newdata.load(savefile)
        checkCorrectness(newdata)

        # Saving with a parent
        data1 = self.data1.copy()
        data1.parent = self.data2.copy()
        data1.save(savefile)
        newdata = MetricData(file=savefile)
        checkCorrectness(newdata)


if __name__ == "__main__":
    unittest.main(verbosity=2)
