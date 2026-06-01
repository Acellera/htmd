import numpy as np
import unittest
from htmd.kinetics import Kinetics


class TestKinetics(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.trans_prob = np.array([[0.6, 0.2, 0.2], [0.3, 0.4, 0.3], [0.2, 0.3, 0.5]])

    def test_kinetics(self):
        from htmd.metricdata import _generate_toy_data
        from htmd.model import Model
        import tempfile
        import os

        fakedata = _generate_toy_data(self.trans_prob, n_traj=100, seed=0)
        model = Model(fakedata)
        model.markovModel(1, 3)

        kin = Kinetics(model, 300, source=0)
        assert kin.sink == 2
        assert np.array_equal(kin.sinkmicro, [0])
        assert kin.source == 0
        assert np.array_equal(kin.sourcemicro, [1])

        rates = kin.getRates(0, 2, states="macro")
        assert np.allclose(rates.mfpton, 3.8699044076299556)
        assert np.allclose(rates.mfptoff, 4.312080466715925)
        assert np.allclose(rates.kon, 258404315.6281552)
        assert np.allclose(rates.koff, 231906618.56122518)
        assert np.allclose(rates.kdeq, 0.7542502489709302)
        assert np.allclose(rates.g0eq, -0.16813599009752134)

        rates = kin.getRates(0, 2, states="micro")
        assert np.allclose(rates.mfpton, 4.5232360505688645)
        assert np.allclose(rates.mfptoff, 4.17783500152719)
        assert np.allclose(rates.kon, 221080657.48066255)
        assert np.allclose(rates.koff, 239358423.59366855)
        assert np.allclose(rates.kdeq, 1.1703043390973993)
        assert np.allclose(rates.g0eq, 0.09375460063521723)

        with tempfile.TemporaryDirectory() as tmpdir:
            outplot = os.path.join(tmpdir, "test.png")
            kin.plotFluxPathways(
                statetype="macro",
                mode="net_flux",
                fraction=0.1,
                plot=False,
                save=outplot,
            )
            assert os.path.exists(outplot)

            outplot = os.path.join(tmpdir, "mfpt.png")
            kin.plotMarkovModel(plot=False, save=outplot)
            assert os.path.exists(outplot)
