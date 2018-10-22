# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import shutil
import time
import unittest
from subprocess import call

import numpy as np

from htmd.home import home
from htmd.util import tempname


class TestParameterize(unittest.TestCase):

    def setUp(self):
        if os.environ.get('TRAVIS_OS_NAME') == 'osx':
            self.skipTest('Mac does not work!')

        self.maxDiff = None
        self.dataDir = home(dataDir='test-param')
        self.testDir = os.environ.get('TESTDIR', tempname())

    def _execute(self, refDir, resDir, command):
        os.makedirs(resDir, exist_ok=True)

        molFile = os.path.join(refDir, 'input.mol2')
        self.assertTrue(os.path.exists(molFile))
        shutil.copy(molFile, resDir)

        print('')  # Just for a better readability
        returncode = call(command.split(), cwd=resDir)
        self.assertEqual(returncode, 0)

    def _testFiles(self, refDir, resDir):
        testFiles = []
        exclusions = ('minimize', 'esp', 'dihedral', '.coor', '.svg')
        for root, _, files in os.walk(refDir, followlinks=True):
            for file in files:
                relFile = os.path.relpath(os.path.join(root, file), start=refDir)
                if any([relFile.startswith(exclusion) or relFile.endswith(exclusion) for exclusion in exclusions]):
                    continue
                testFiles.append(os.path.join(root, file))

        print('Compared files:')
        for file in testFiles:
            relFile = os.path.relpath(file, start=refDir)
            refFile = os.path.join(refDir, relFile)
            resFile = os.path.join(resDir, relFile)
            print('  %s' % relFile)

            with self.subTest(refFile=refFile):
                self.assertTrue(os.path.exists(resFile))

                with open(refFile) as ref, open(resFile) as res:
                    refLines, resLines = ref.readlines(), res.readlines()

                # Removes first line with the version
                if file.endswith('frcmod') or file.endswith('rtf') or file.endswith('prm'):
                    refLines, resLines = refLines[1:], resLines[1:]

                refFields = [field for line in refLines for field in line.split()]
                resFields = [field for line in resLines for field in line.split()]
                for refField, resField in zip(refFields, resFields):
                    try:
                        if np.isclose(float(refField), float(resField), rtol=0, atol=1e-5):
                            continue
                    except ValueError:
                        if refField == resField:
                            continue
                    self.assertListEqual(refLines, resLines)  # If there is a mismatch, print a diff of all file

        print('')

    def test_parameterize_speed(self):

        _started_at = time.time()

        command = 'parameterize -h'
        print('')  # Just for a better readability
        returncode = call(command.split())
        self.assertEqual(returncode, 0)

        elapsed = time.time() - _started_at
        self.assertLessEqual(elapsed, 1.5)

    def test_h2o2_list(self):
        refDir = os.path.join(self.dataDir, 'h2o2_list')
        resDir = os.path.join(self.testDir, 'h2o2_list')
        self._execute(refDir, resDir, 'parameterize input.mol2 -l')
        self._testFiles(refDir, resDir)

    def test_h2o2_gaff2(self):
        refDir = os.path.join(self.dataDir, 'h2o2_gaff2')
        resDir = os.path.join(self.testDir, 'h2o2_gaff2')
        self._execute(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_h2o2_outdir(self):
        refDir = os.path.join(self.dataDir, 'h2o2_outdir')
        resDir = os.path.join(self.testDir, 'h2o2_outdir')
        self._execute(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --no-min --no-dihed -o dir')
        self._testFiles(refDir, resDir)

    def test_h2o2_min(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min')
        resDir = os.path.join(self.testDir, 'h2o2_min')
        self._execute(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --no-dihed')
        self._testFiles(refDir, resDir)

    def test_h2o2_min_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min_restart')
        resDir = os.path.join(self.testDir, 'h2o2_min_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-dihed')
        self._testFiles(refDir, resDir)

    def test_h2o2_esp(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp')
        resDir = os.path.join(self.testDir, 'h2o2_esp')
        self._execute(refDir, resDir, 'parameterize input.mol2 --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_h2o2_esp_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp_restart')
        resDir = os.path.join(self.testDir, 'h2o2_esp_restart')
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix')
        self._execute(refDir, resDir, 'parameterize input.mol2 --no-min --charge-type Gasteiger --no-dihed-opt')
        self._testFiles(refDir, resDir)

    def test_h2o2_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_opt(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt')
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min')
        self._testFiles(refDir, resDir)

    def test_h2o2_dihed_opt_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(sys.version_info.major == 3 and sys.version_info.minor > 5, 'Python 3.5 issue')
    def test_h2o2_full_fake(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake')
        self._execute(refDir, resDir, 'parameterize input.mol2 --fake-qm')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(sys.version_info.major == 3 and sys.version_info.minor > 5, 'Python 3.5 issue')
    def test_h2o2_full_fake_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake_restart')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --fake-qm')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethene_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix')
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    def test_ethene_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_glycol_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix')
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    def test_glycol_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    def test_glycol_dihed_fix_restart_2(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix_restart_2')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['O1-C2-C3-O4', 'C2-C3-O4-H10']
        self._execute(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --no-min --no-dihed-opt'.format(' '.join(dihedrals)))
        self._testFiles(refDir, resDir)

    def test_glycol_dihed_select_1_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_1_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_1_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['O1-C2-C3-O4']
        self._execute(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --no-min --no-dihed-opt'.format(' '.join(dihedrals)))
        self._testFiles(refDir, resDir)

    def test_glycol_dihed_select_2_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_2_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_2_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['C2-C3-O4-H10']
        self._execute(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --no-min --no-dihed-opt'.format(' '.join(dihedrals)))
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethanolamine_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix')
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_ethanolamine_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed-opt')
        self._testFiles(refDir, resDir)

    def test_benzamidine_gasteiger(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gasteiger')
        resDir = os.path.join(self.testDir, 'benzamidine_gasteiger')
        self._execute(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_benzamidine_gaff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff')
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_benzamidine_gaff2(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff2')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff2')
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_benzamidine_cgenff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_cgenff')
        resDir = os.path.join(self.testDir, 'benzamidine_cgenff')
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff CGenFF_2b6 --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    def test_benzamidine_rtf_prm(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_rtf_prm')
        resDir = os.path.join(self.testDir, 'benzamidine_rtf_prm')
        os.makedirs(resDir, exist_ok=True)
        shutil.copy(os.path.join(refDir, 'input.rtf'), resDir)
        shutil.copy(os.path.join(refDir, 'input.prm'), resDir)
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff CGenFF_2b6 --rtf-prm input.rtf input.prm --charge-type Gasteiger --no-min --no-dihed')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full(self):
        assert 'HTMD_CONFIG' in os.environ, '"HTMD_CONFIG" environment variable has to be set'
        refDir = os.path.join(self.dataDir, 'benzamidine_full')
        resDir = os.path.join(self.testDir, 'benzamidine_full')
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --basis 6-31G* -q Slurm')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_full_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_full_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --basis 6-31G*')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_esp_freeze_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_esp_freeze_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_esp_freeze_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --fix-charge N14 --basis 6-31G* --no-dihed')
        self._testFiles(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_dihed_select_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_dihed_select_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_dihed_select_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._execute(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 -d C2-C1-C7-N13 --basis 6-31G*')
        self._testFiles(refDir, resDir)


if __name__ == '__main__':
    unittest.main(verbosity=2)
