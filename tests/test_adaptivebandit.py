import unittest


class TestAdaptiveBandit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from htmd.util import tempname
        from htmd.home import home
        import shutil
        import os

        tmpdir = tempname()
        shutil.copytree(home(dataDir="adaptive"), tmpdir)
        os.chdir(tmpdir)
        print(f"Running AdaptiveBandit test in {tmpdir}")

    # def test_adaptive(self):
    #     from sklearn.cluster import MiniBatchKMeans
    #     from jobqueues.localqueue import LocalCPUQueue
    #     from moleculekit.projections.metricdistance import MetricDistance

    #     import numpy as np
    #     import random

    #     np.random.seed(0)  # Needed for the clustering to always give same results
    #     random.seed(0)

    #     md = AdaptiveBandit()
    #     md.app = LocalCPUQueue()
    #     md.generatorspath = "generators"
    #     md.inputpath = "input"
    #     md.datapath = "data"
    #     md.coorname = "input.coor"
    #     md.filter = True
    #     md.filtersel = "all"

    #     md.clustmethod = MiniBatchKMeans
    #     md.projection = MetricDistance(
    #         "protein resid 173 and name CA",
    #         "resname BEN and name C1 C2 C3 C7",
    #         periodic="selections",
    #     )
    #     md.ticadim = 2
    #     md.nmin = 1
    #     md.nmax = 2
    #     md.nepochs = 9999
    #     md.nframes = 1000000

    #     md.reward_method = "mean"
    #     md.exploration = 0.01
    #     md.actionspace = "tica"
    #     md.actionpool = 0
    #     md.recluster = False

    #     md.save = True
    #     md.dryrun = True
    #     md.run()
