import unittest
import numpy as np
from htmd.model import Model, getStateStatistic, macroAccumulate


class TestModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.simlist import simlist, simfilter
        from glob import glob
        from htmd.projections.metric import Metric
        from moleculekit.projections.metricdistance import MetricDistance
        from moleculekit.util import tempname
        from sklearn.cluster import MiniBatchKMeans
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
        data = metr.project()
        data.cluster(MiniBatchKMeans(n_clusters=4))

        self.model = Model(data)

        self.trans_prob = np.array([[0.6, 0.2, 0.2], [0.3, 0.4, 0.3], [0.2, 0.3, 0.5]])

        self.trans_prob2 = np.array(
            [
                [0.8, 0.15, 0.05, 0.0, 0.0],
                [0.1, 0.75, 0.05, 0.05, 0.05],
                [0.05, 0.1, 0.8, 0.0, 0.05],
                [0.0, 0.2, 0.0, 0.8, 0.0],
                [1e-7, 0.02 - 1e-7, 0.02, 0.0, 0.96],
            ]
        )

    def test_model_saving_loading(self):
        from moleculekit.util import tempname

        modelfile = tempname(suffix=".dat")
        self.model.save(modelfile)

        newmodel = Model(file=modelfile)
        assert newmodel.data.numTrajectories == 2

        # Testing model saving when data has parent
        self.model.data.parent = self.model.data.copy()
        self.model.save(modelfile)
        newmodel = Model(file=modelfile)
        assert newmodel.data.numTrajectories == 2
        assert newmodel.data.parent.numTrajectories == 2
        self.model.data.parent = None

    def test_max_connected_lag(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        maxlag, its = model.maxConnectedLag([100, 500, 999], njobs=1)
        assert maxlag == 999
        assert np.allclose(
            its,
            np.array(
                [
                    [19.86198189, 15.92237368],
                    [81.7366474, 68.67389378],
                    [490.97275421, 330.52311174],
                ]
            ),
        )

    def test_plot_timescales(self):
        from htmd.metricdata import _generate_toy_data
        from matplotlib import pylab as plt
        import tempfile
        import os

        plt.ioff()
        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)

        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            its, lags = model.plotTimescales(
                lags=[1, 50, 100, 500, 800],
                plot=False,
                save=outplot,
                njobs=1,
                results=True,
            )
            assert os.path.exists(outplot)

        assert np.allclose(
            its,
            np.array(
                [
                    [0.98644281, 0.50313423],
                    [8.96353244, 8.72125828],
                    [19.86198189, 15.92237368],
                    [81.7366474, 68.67389378],
                    [183.99300253, 150.92253294],
                ]
            ),
        )
        assert np.array_equal(lags, np.array([1, 50, 100, 500, 800]))

        # Second probability matrix test
        fakedata = _generate_toy_data(self.trans_prob2, n_traj=100, seed=0)
        model = Model(fakedata)

        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            its, lags = model.plotTimescales(
                lags=range(1, 4),
                plot=False,
                save=outplot,
                njobs=1,
                results=True,
            )
            assert os.path.exists(outplot)

        assert np.allclose(
            its,
            np.array(
                [
                    [14.83314609, 4.8958021, 3.53366645, 2.01939226],
                    [14.82089692, 4.86804921, 3.49197211, 2.00651665],
                    [14.81046561, 4.8277952, 3.47217915, 2.02697621],
                ]
            ),
        )

    def test_ck_test(self):
        from htmd.metricdata import _generate_toy_data
        from matplotlib import pylab as plt
        import tempfile
        import os

        plt.ioff()
        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            model.cktest(plot=False, save=outplot)
            assert os.path.exists(outplot)

    def test_msm(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)

        assert np.allclose(self.trans_prob, model.P, atol=1e-1)

        eq_dist = np.array([0.28912558, 0.32754596, 0.38332846])
        assert np.allclose(eq_dist, model.eqDistribution(plot=False), atol=1e-3)

        assert np.array_equal(model.micro_ofcluster, [0, 1, 2])
        assert np.array_equal(model.cluster_ofmicro, [0, 1, 2])
        assert np.array_equal(model.macro_ofmicro, [2, 0, 1])
        assert model.micronum == 3
        assert model.macronum == 3

        # Second probability matrix test
        fakedata = _generate_toy_data(self.trans_prob2, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)

        assert np.allclose(self.trans_prob2, model.P, atol=1e-1)

        eq_dist = np.array([0.14541275, 0.33704423, 0.51754301])
        assert np.allclose(eq_dist, model.eqDistribution(plot=False), atol=1e-3)

        assert np.array_equal(model.micro_ofcluster, [0, 1, 2, 3, 4])
        assert np.array_equal(model.cluster_ofmicro, [0, 1, 2, 3, 4])
        assert np.array_equal(model.macro_ofmicro, [1, 1, 1, 0, 2])
        assert model.micronum == 5
        assert model.macronum == 3

    def test_create_state(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        # Create a new macrostate assigning microstates 0 and 1 to it
        model.createState(microstates=[0, 1])
        assert model.macronum == 2
        assert np.array_equal(model.macro_ofmicro, [1, 1, 0])
        model.markovModel(1, 3)  # Revert changes

        # Create a new cluster out of the first frames of the first trajectory
        model.createState(indexpairs=[[0, 0], [0, 5], [0, 10]])
        model.markovModel(1, 3)
        assert np.array_equal(model.micro_ofcluster, [0, 1, 2, 3])
        assert np.array_equal(model.cluster_ofmicro, [0, 1, 2, 3])
        assert np.array_equal(model.macro_ofmicro, [2, 0, 1, 2])

        # Create a new cluster out of the first frames of the first trajectory without back transitions
        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        model.createState(indexpairs=[[0, 0], [0, 1], [0, 2]])  # No back transitions
        assert np.array_equal(model.micro_ofcluster, [0, 1, 2, -1])
        assert np.array_equal(model.cluster_ofmicro, [0, 1, 2])
        assert np.array_equal(model.macro_ofmicro, [2, 0, 1])

    def test_core_set_model(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.createState(indexpairs=[[0, 0], [0, 5], [0, 10]])  # Create a rare state
        model.markovModel(1, 3)
        assert np.array_equal(model.micro_ofcluster, [0, 1, 2, 3])
        assert np.array_equal(model.cluster_ofmicro, [0, 1, 2, 3])
        assert np.array_equal(model.macro_ofmicro, [2, 0, 1, 2])

        coremodel, _ = model.createCoreSetModel(threshold=0.5)
        coremodel.markovModel(1, 3)
        assert np.array_equal(coremodel.micro_ofcluster, [0, 1, 2])
        assert np.array_equal(coremodel.cluster_ofmicro, [0, 1, 2])
        assert np.array_equal(coremodel.macro_ofmicro, [0, 1, 2])

    def test_plot_fes(self):
        from htmd.metricdata import _generate_toy_data
        from matplotlib import pylab as plt
        import tempfile
        import os

        plt.ioff()

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)

        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            model.plotFES(0, 1, 300, states=True, s=10, plot=False, save=outplot)
            assert os.path.exists(outplot)

        fakedata2 = _generate_toy_data(
            self.trans_prob, n_traj=100, seed=0, cluster=False
        )
        # Plot with second data object
        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            model.plotFES(
                0, 1, 300, states=True, s=10, plot=False, data=fakedata2, save=outplot
            )
            assert os.path.exists(outplot)

    def test_sample_states(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        np.random.seed(0)
        model.sampleStates(states=[0, 1])

    def test_get_state_statistic(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        statestats = getStateStatistic(
            model, fakedata, states=[0, 1], statetype="micro", method=np.mean
        )
        assert np.allclose(statestats[0], [1, 0, 0])
        assert np.allclose(statestats[1], [0, 1, 0])

        statestats = getStateStatistic(
            model, fakedata, states=[0, 1, 2], statetype="cluster", method=np.mean
        )
        assert np.allclose(statestats[0], [1, 0, 0])
        assert np.allclose(statestats[1], [0, 1, 0])
        assert np.allclose(statestats[2], [0, 0, 1])

        model.createState(microstates=[0, 1])
        statestats = getStateStatistic(
            model, fakedata, states=[0, 1], statetype="macro", method=np.mean
        )
        assert np.allclose(statestats[0], [0, 0, 1])
        assert np.allclose(statestats[1], [0.56985649, 0.43014351, 0])

    def test_macro_accumulate(self):
        from htmd.metricdata import _generate_toy_data

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)
        acc = macroAccumulate(model, [1.5, 3.7, 4.5])
        assert np.array_equal(acc, np.array([3.7, 4.5, 1.5]))

        model.createState(microstates=[0, 1])
        acc = macroAccumulate(model, [1.5, 3.7, 4.5])
        assert np.array_equal(acc, np.array([4.5, 5.2]))
