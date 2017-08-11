# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
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

        if os.environ.get('TRAVIS_OS_NAME') == 'osx':
            self.skipTest('Mac does not work!')

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

                if file.startswith('minimize') or \
                   file.startswith('esp') or \
                   file.startswith('dihedral-single-point') or \
                   file.startswith('dihedral-opt') or \
                   file.endswith('.coor') or \
                   file.endswith('.svg'):
                    continue

                with self.subTest(file=file):
                    self.assertTrue(os.path.exists(refFile))
                    self.assertTrue(os.path.exists(resFile))

                    with open(refFile) as ref, open(resFile) as res:
                        refLines, resLines = ref.readlines(), res.readlines()

                    # HACK! Before Python 3.6 dict does not preserve order, so the lines are sorted before comparison
                    if file == 'stdout' and sys.version_info.major == 3 and sys.version_info.minor < 6:
                        refLines, resLines = sorted(refLines), sorted(resLines)

                    # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                    #       Also the first line with the version is removed
                    if file.endswith('frcmod'):
                        refLines, resLines = sorted(refLines[1:]), sorted(resLines[1:])

                    self.assertListEqual(refLines, resLines, msg=file)

    def test_doc(self):

        refDir = os.path.join(self.dataDir, 'doc')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_list(self):

        refDir = os.path.join(self.dataDir, 'h2o2_list')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_gaff2(self):

        refDir = os.path.join(self.dataDir, 'h2o2_gaff2')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    @unittest.skip('Finish')
    def test_h2o2_outdir(self):
        pass

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_min(self):

        refDir = os.path.join(self.dataDir, 'h2o2_min')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_min_restart(self):

        refDir = os.path.join(self.dataDir, 'h2o2_min_restart')
        with TemporaryDirectory() as resDir:
            shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
            self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_esp(self):

        refDir = os.path.join(self.dataDir, 'h2o2_esp')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_esp_restart(self):

        refDir = os.path.join(self.dataDir, 'h2o2_esp_restart')
        with TemporaryDirectory() as resDir:
            shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
            self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_fix(self):

        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_dihed_fix_restart(self):

        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix_restart')
        with TemporaryDirectory() as resDir:
            shutil.copytree(os.path.join(refDir, 'dihedral-single-point'),
                            os.path.join(resDir, 'dihedral-single-point'))
            self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_opt(self):

        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt')
        with TemporaryDirectory() as resDir:
            self._test(refDir, resDir)

    def test_h2o2_dihed_opt_restart(self):

        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt_restart')
        with TemporaryDirectory() as resDir:
            shutil.copytree(os.path.join(refDir, 'dihedral-opt'),
                            os.path.join(resDir, 'dihedral-opt'))
            self._test(refDir, resDir)

    def test_benzamidine_gaff(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_gaff')
        with TemporaryDirectory() as resDir:
            #resDir = 'benzamidine_gaff'
            #os.mkdir(resDir)
            self._test(refDir, resDir)

    def test_benzamidine_gaff2(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_gaff2')
        with TemporaryDirectory() as resDir:
            #resDir = 'benzamidine_gaff2'
            #os.mkdir(resDir)
            self._test(refDir, resDir)

    def test_benzamidine_cgenff(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_cgenff')
        with TemporaryDirectory() as resDir:
            #resDir = 'benzamidine_cgenff'
            #os.mkdir(resDir)
            self._test(refDir, resDir)

    @unittest.skip('Finish')
    def test_benzamidine_rtf_prm(self):
        pass

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_full')
        with TemporaryDirectory() as resDir:
            resDir = 'benzamidine_full'
            os.mkdir(resDir)
            self._test(refDir, resDir)

    @unittest.skip('Finish')
    def test_benzamidine_full_restart(self):
        pass

    @unittest.skip('Finish')
    def test_benzamidine_esp_freeze_restart(self):
        pass

    @unittest.skip('Finish')
    def test_benzamidine_dihed_select_restart(self):
        pass


if __name__ == '__main__':

    os.environ['HTMD_NONINTERACTIVE'] = '1'
    unittest.main()
