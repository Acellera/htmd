import sys
import unittest

from htmd.adaptive.adaptiverun import AdaptiveMD


class TestAdaptiveRun(unittest.TestCase):
    @unittest.skipIf(sys.platform == "win32", "Windows is not supported")
    def test_adaptive_run(self):
        from moleculekit.projections.metricdistance import MetricSelfDistance
        from jobqueues.simqueue import SimQueue
        from htmd.home import home
        from glob import glob
        import re
        import os
        import tempfile
        import shutil

        with tempfile.TemporaryDirectory() as tmpdir:
            source = os.path.join(tmpdir, "adaptivemd")
            shutil.copytree(os.path.join(home(dataDir="adaptivemd")), source)

            class FakeQueue(SimQueue):
                def __init__(self, source, datadir):
                    super().__init__()
                    self._source = source
                    self._datadir = datadir
                    self._dirs = []

                def retrieve(self):
                    basename = os.path.basename(os.path.abspath(self._dirs[0]))
                    regex = re.compile(r"e(\d+)")
                    epoch = int(regex.search(basename).group(1))
                    for ed, dd in zip(
                        glob(os.path.join(self._source, "data", f"e{epoch}s*")),
                        self._dirs,
                    ):
                        basename = os.path.basename(os.path.abspath(dd))
                        shutil.copytree(ed, os.path.join(self._datadir, basename))

                def submit(self, dirs):
                    self._dirs = dirs

                def inprogress(self):
                    return 0

                def memory(self):
                    return 0

                def ncpu(self):
                    return 1

                def ngpu(self):
                    return 0

                def stop(self):
                    pass

            queue = FakeQueue(source, os.path.join(tmpdir, "data"))

            ad = AdaptiveMD()
            ad.app = queue
            ad.nmin = 1
            ad.nmax = 3
            ad.nepochs = 3
            ad.generatorspath = os.path.join(source, "generators")
            ad.inputpath = os.path.join(tmpdir, "input")
            ad.datapath = os.path.join(tmpdir, "data")
            ad.filteredpath = os.path.join(tmpdir, "filtered")
            protsel = "protein and name CA"
            ad.projection = MetricSelfDistance(protsel)
            ad.updateperiod = 1  # execute every 1 sec
            ad.run()

            assert len(glob(os.path.join(tmpdir, "filtered", "e1*", ""))) == 3
            assert len(glob(os.path.join(tmpdir, "filtered", "e2*", ""))) == 2
            assert len(glob(os.path.join(tmpdir, "filtered", "e3*", ""))) == 2


if __name__ == "__main__":
    unittest.main(verbosity=2)
