# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
import unittest
from tempfile import TemporaryDirectory
from subprocess import call


class TestParameterize(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None

        self.dataDir = os.path.join('..', 'data', 'test-param')
        if not os.path.exists(self.dataDir):
            self.dataDir = os.path.join('htmd', 'data', 'test-param')

        suffix = 'travis' if 'TRAVIS' in os.environ else 'local'
        self.dataDir = os.path.abspath(os.path.join(self.dataDir, suffix))

    def _test(self, refDir, resDir):

        shutil.copy(os.path.join(refDir, 'stdin'), resDir)
        molFile = os.path.join(refDir, 'input.mol2')
        if os.path.exists(molFile):
            shutil.copy(molFile, resDir)

        with open(os.path.join(resDir, 'stdin')) as stdin, \
             open(os.path.join(resDir, 'stdout'), 'w') as stdout, \
             open(os.path.join(resDir, 'stderr'), 'w') as stderr:
            arguments = stdin.readline().split()
            returncode = call(arguments, stdout=stdout, stderr=stderr, cwd=resDir)
            self.assertEqual(returncode, 0)

        for directory, _, files in os.walk(refDir):
            for file in files:
                file = os.path.relpath(os.path.join(directory, file), refDir)
                refFile, resFile = os.path.join(refDir, file), os.path.join(resDir, file)

                if file.startswith('minimize') or\
                   file.startswith('esp') or\
                   file.endswith('mol.coor'):
                    continue

                with self.subTest(file=file):
                    self.assertTrue(os.path.exists(refFile))
                    self.assertTrue(os.path.exists(resFile))

                    with open(refFile) as ref, open(resFile) as res:
                        refLines, resLines = ref.readlines(), res.readlines()

                    if file.endswith('frcmod'):
                        # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                        #       Also the first line with the version is removed
                        refLines, resLines = sorted(refLines[1:]), sorted(resLines[1:])

                    self.assertListEqual(refLines, resLines, msg=file)

    def test_doc(self):

        with TemporaryDirectory() as resDir:
            self._test(os.path.join(self.dataDir, 'doc'), resDir)

    def test_h2o2_list(self):

        with TemporaryDirectory() as resDir:
            self._test(os.path.join(self.dataDir, 'list'), resDir)

    def test_h2o2_gaff2(self):

        with TemporaryDirectory() as resDir:
            self._test(os.path.join(self.dataDir, 'gaff2'), resDir)

    def test_h2o2_gaff2_outdir(self):
        pass

    @unittest.skip
    def test_h2o2_gaff2_min(self):

        with TemporaryDirectory() as resDir:
            self._test(os.path.join(self.dataDir, 'gaff2_min'), resDir)

    def test_h2o2_gaff2_min_restart(self):
        pass

    @unittest.skip
    def test_h2o2_gaff2_esp(self):

        with TemporaryDirectory() as resDir:
            #resDir = 'res'
            #os.mkdir(resDir)
            self._test(os.path.join(self.dataDir, 'gaff2_esp'), resDir)

    def test_h2o2_gaff2_esp_restart(self):
        pass

    def test_h2o2_gaff2_dihed_fix(self):
        pass

    def test_h2o2_gaff2_dihed_fix_restart(self):
        pass

    def test_h2o2_gaff2_dihed_opt(self):
        pass

    def test_h2o2_gaff2_dihed_opt_restart(self):
        pass


if __name__ == '__main__':

    os.environ['HTMD_NONINTERACTIVE'] = '1'
    unittest.main()
