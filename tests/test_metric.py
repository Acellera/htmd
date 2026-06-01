import unittest
from htmd.projections.metric import Metric


class TestMetric(unittest.TestCase):
    def test_set(self):
        from moleculekit.molecule import Molecule
        from moleculekit.projections.metricrmsd import MetricRmsd

        # Testing the set method of Metric
        sims = []
        ref = Molecule("3PTB")
        metr = Metric(sims)
        metr.set(MetricRmsd(ref, "protein and name CA"))
        assert len(metr.projectionlist) == 1

        metr.set(
            [
                MetricRmsd(ref, "protein and name CA"),
                MetricRmsd(ref, "protein and name CA"),
            ]
        )
        assert len(metr.projectionlist) == 2

    def test_function_projections(self):
        from moleculekit.molecule import Molecule

        def foo(mol, ref):
            from moleculekit.util import molRMSD

            mol.wrap("protein")
            mol.align("protein and name CA", refmol=ref)
            return molRMSD(
                mol,
                ref,
                mol.atomselect("protein and name CA"),
                ref.atomselect("protein and name CA"),
            )

        ref = Molecule("3PTB")
        sims = []
        metr = Metric(sims)
        metr.set((foo, (ref,)))
        assert len(metr.projectionlist) == 1
