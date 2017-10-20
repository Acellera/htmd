# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import shutil
import unittest
from subprocess import call

from htmd.home import home
from htmd.util import tempname


class TestParameterize(unittest.TestCase):
    def setUp(self):
        if os.environ.get('TRAVIS_OS_NAME') == 'osx':
            self.skipTest('Mac does not work!')

        self.maxDiff = None
        self.dataDir = home(dataDir='test-param')
        if os.environ.get('TESTDIR') != '':
            self.testDir = os.environ['TESTDIR']
        else:
            self.testDir = tempname()
        print('Test dir: {}'.format(self.testDir))

    def _test(self, refDir, resDir, cmd):

        if not os.path.exists(resDir):
            os.makedirs(resDir)

        molFile = os.path.join(refDir, 'input.mol2')
        if os.path.exists(molFile):
            shutil.copy(molFile, resDir)

        arguments = cmd.split()
        returncode = call(arguments, cwd=resDir)
        self.assertEqual(returncode, 0)

        filestotest = []
        excluded = ('minimize', 'esp', 'dihedral', '.coor', '.svg')
        for root, _, files in os.walk(refDir, followlinks=True):
            for file in files:
                flag = False
                relFile = os.path.relpath(os.path.join(root, file), start=refDir)
                for exc in excluded:
                    if relFile.startswith(exc) or relFile.endswith(exc):
                        flag = True
                if not flag:
                    filestotest.append(os.path.join(root, file))

        print('Compared files:')
        for file in filestotest:
            relFile = os.path.relpath(file, start=refDir)
            print('  %s' % relFile)

            refFile = os.path.join(refDir, relFile)
            resFile = os.path.join(resDir, relFile)

            with self.subTest(refFile=refFile):
                self.assertTrue(os.path.exists(resFile))

                with open(refFile) as ref, open(resFile) as res:
                    refLines, resLines = ref.readlines(), res.readlines()

                # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                #       Also the first line with the version is removed
                if file.endswith('frcmod') or file.endswith('rtf'):
                    refLines, resLines = sorted(refLines[1:]), sorted(resLines[1:])

                # Removes first line with the version
                if file.endswith('prm'):
                    refLines, resLines = refLines[1:], resLines[1:]

                if file.endswith('prm') or file.endswith('frcmod') or \
                        os.path.relpath(file, start=refDir).startswith('energies'):
                    refFields = [field for line in refLines for field in line.split()]
                    resFields = [field for line in resLines for field in line.split()]
                    for refField, resField in zip(refFields, resFields):
                        with self.subTest():
                            try:
                                refFloat = float(refField)
                                resFloat = float(resField)
                                self.assertAlmostEqual(refFloat, resFloat, places=4, msg=refFile)
                            except ValueError:
                                self.assertEqual(refField, resField, msg=refFile)
                else:
                    self.assertListEqual(refLines, resLines, msg=refFile)

        print('')

    def test_h2o2_list(self):
        refDir = os.path.join(self.dataDir, 'h2o2_list')
        resDir = os.path.join(self.testDir, 'h2o2_list')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -l')

    def test_h2o2_gaff2(self):
        refDir = os.path.join(self.dataDir, 'h2o2_gaff2')
        resDir = os.path.join(self.testDir, 'h2o2_gaff2')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-torsions')

    def test_h2o2_outdir(self):
        refDir = os.path.join(self.dataDir, 'h2o2_outdir')
        resDir = os.path.join(self.testDir, 'h2o2_outdir')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-torsions -o dir')

    def test_h2o2_min(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min')
        resDir = os.path.join(self.testDir, 'h2o2_min')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-esp --no-torsions')

    def test_h2o2_min_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min_restart')
        resDir = os.path.join(self.testDir, 'h2o2_min_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-esp --no-torsions')

    def test_h2o2_esp(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp')
        resDir = os.path.join(self.testDir, 'h2o2_esp')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-torsions')

    def test_h2o2_esp_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp_restart')
        resDir = os.path.join(self.testDir, 'h2o2_esp_restart')
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-torsions')

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    def test_h2o2_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_opt(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp')

    def test_h2o2_dihed_opt_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp')

    @unittest.skipUnless(sys.version_info.major == 3 and sys.version_info.minor > 5, 'Python 3.5 issue')
    def test_h2o2_full_fake(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --fake-qm')

    @unittest.skipUnless(sys.version_info.major == 3 and sys.version_info.minor > 5, 'Python 3.5 issue')
    def test_h2o2_full_fake_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake_restart')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --fake-qm')

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethene_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    def test_ethene_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_glycol_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    def test_glycol_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    def test_glycol_dihed_select_1_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_1_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_1_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt --torsion O1-C1-C2-O2')

    def test_glycol_dihed_select_2_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_2_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_2_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt --torsion C1-C2-O2-H6')

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethanolamine_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix')
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_ethanolamine_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 -f GAFF2 --no-min --no-esp --no-geomopt')

    def test_benzamidine_gaff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff')
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --forcefield GAFF --no-min --no-esp --no-torsions')

    def test_benzamidine_gaff2(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff2')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff2')
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --forcefield GAFF2 --no-min --no-esp --no-torsions')

    def test_benzamidine_cgenff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_cgenff')
        resDir = os.path.join(self.testDir, 'benzamidine_cgenff')
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --forcefield CGENFF --no-min --no-esp --no-torsions')

    def test_benzamidine_rtf_prm(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_rtf_prm')
        resDir = os.path.join(self.testDir, 'benzamidine_rtf_prm')
        os.makedirs(resDir)
        shutil.copy(os.path.join(refDir, 'input.rtf'), resDir)
        shutil.copy(os.path.join(refDir, 'input.prm'), resDir)
        self._test(refDir, resDir,
                   'parameterize -m input.mol2 --charge 1 --forcefield CGENFF --rtf input.rtf --prm input.prm --no-min '
                   '--no-esp --no-torsions')

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_full')
        resDir = os.path.join(self.testDir, 'benzamidine_full')
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --basis 6-31G*')

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_full_restart(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_full_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_full_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --basis 6-31G*')

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_esp_freeze_restart(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_esp_freeze_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_esp_freeze_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._test(refDir, resDir,
                   'parameterize -m input.mol2 --charge 1 --basis 6-31G* --no-torsions --freeze-charge N2')

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_dihed_select_restart(self):

        refDir = os.path.join(self.dataDir, 'benzamidine_dihed_select_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_dihed_select_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._test(refDir, resDir, 'parameterize -m input.mol2 --charge 1 --basis 6-31G* --torsion C2-C1-C7-N1')


if __name__ == '__main__':

    unittest.main(verbosity=2)
