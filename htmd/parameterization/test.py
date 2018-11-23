# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
import time
import unittest
from subprocess import call

import numpy as np

from htmd.home import home
from htmd.util import tempname


class TestParameterize(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None  # Make the diff to be complete
        self.dataDir = home(dataDir='test-param')
        self.testDir = os.environ.get('TESTDIR', tempname())

    def _run(self, refDir, resDir, command):

        os.makedirs(resDir, exist_ok=True)

        # Check the existance of the input file and copy it
        molFile = os.path.join(refDir, 'input.mol2')
        self.assertTrue(os.path.exists(molFile))
        shutil.copy(molFile, resDir)

        # Execute the test
        print('', flush=True)
        self.assertEqual(call(command.split(), cwd=resDir), 0)
        print('', flush=True)

    def _test(self, refDir, resDir, energyTermRelTol=1e-6, energyProfileAbsTol=1e-6,
              dihedralForceConstAbsTol=1e-6, dihedralPhaseAbsTol=1e-6):

        assert energyTermRelTol < 1
        assert energyProfileAbsTol < 1
        assert dihedralForceConstAbsTol < 1
        assert dihedralPhaseAbsTol < 10

        # Default tolerances
        defaultAbsTol = 1e-6
        defaultRelTol = 1e-6

        # Find the tested files
        testedFiles = []
        exclusions = ('minimize', 'esp', 'dihedral', '.coor', '.svg', 'random-search.log')
        for root, _, files in os.walk(refDir, followlinks=True):
            for file in files:
                relFile = os.path.relpath(os.path.join(root, file), start=refDir)
                if any([relFile.startswith(exclusion) or relFile.endswith(exclusion) for exclusion in exclusions]):
                    continue
                testedFiles.append(os.path.join(root, file))

        # Test the files
        print('Compared files:', flush=True)
        for file in testedFiles:
            relFile = os.path.relpath(file, start=refDir)
            refFile = os.path.join(refDir, relFile)
            resFile = os.path.join(resDir, relFile)
            print('  {}'.format(relFile))

            # Reset tolerances
            absTol = defaultAbsTol
            relTol = defaultRelTol

            # Set the tolerance for "energies.txt" file
            if relFile.endswith('energies.txt'):
                absTol = 0
                relTol = energyTermRelTol

            # Set the tolerance for "plots/*.dat" files
            if relFile.endswith('.dat'):
                absTol = energyProfileAbsTol
                relTol = 0

            with self.subTest(refFile=refFile):

                self.assertTrue(os.path.exists(resFile),
                                msg='Could not find file {} corresponding to {}'.format(resFile, refFile))

                with open(refFile) as ref, open(resFile) as res:
                    refLines, resLines = ref.readlines(), res.readlines()

                # Removes the first line with the HTMD version
                if file.endswith('frcmod') or file.endswith('rtf') or file.endswith('prm'):
                    refLines = refLines[1:]
                    resLines = resLines[1:]

                # Iterate over lines in the file
                for refLine, resLine in zip(refLines, resLines):
                    refFields = refLine.split()
                    resFields = resLine.split()

                    # Iterate over fields in the line
                    for iField, (refField, resField) in enumerate(zip(refFields, resFields)):

                        # Set tolerance for "mol.frcmod" files
                        if relFile.endswith('.frcmod'):
                            if len(refFields) == 7 and iField == 2:  # Detect dihedral force constant column
                                absTol = dihedralForceConstAbsTol
                                relTol = 0
                            elif len(refFields) == 7 and iField == 3:  # Detect dihedral phase column
                                absTol = dihedralPhaseAbsTol
                                relTol = 0
                            else:
                                absTol = defaultAbsTol
                                relTol = defaultRelTol

                        try:
                            refField = float(refField)
                            resField = float(resField)
                        except ValueError:
                            # The fields cannot be converted to floats, so compare them directly
                            if refField == resField:
                                continue
                        else:
                            # The fields can be converted to floats, so compare them with the tolerances
                            if np.isclose(refField, resField, atol=absTol, rtol=relTol):
                                continue

                        # Print in case of failure
                        print('Failed: {} == {}'.format(refField, resField))
                        print('Absolute tolerance: {}'.format(absTol))
                        print('Relative tolernace: {}'.format(relTol))
                        self.assertListEqual(refLines, resLines)  # If there is a mismatch, print a diff of all file

        print('', flush=True)

    def test_load_time(self):

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
        self._run(refDir, resDir, 'parameterize input.mol2 -l')
        self._test(refDir, resDir)

    def test_h2o2_gaff2(self):
        refDir = os.path.join(self.dataDir, 'h2o2_gaff2')
        resDir = os.path.join(self.testDir, 'h2o2_gaff2')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_h2o2_outdir(self):
        refDir = os.path.join(self.dataDir, 'h2o2_outdir')
        resDir = os.path.join(self.testDir, 'h2o2_outdir')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --min-type None --no-dihed -o dir')
        self._test(refDir, resDir)

    def test_h2o2_min(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min')
        resDir = os.path.join(self.testDir, 'h2o2_min')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --no-dihed')
        self._test(refDir, resDir)

    def test_h2o2_min_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min_restart')
        resDir = os.path.join(self.testDir, 'h2o2_min_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --no-dihed')
        self._test(refDir, resDir)

    def test_h2o2_esp(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp')
        resDir = os.path.join(self.testDir, 'h2o2_esp')
        self._run(refDir, resDir, 'parameterize input.mol2 --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_h2o2_esp_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_esp_restart')
        resDir = os.path.join(self.testDir, 'h2o2_esp_restart')
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._run(refDir, resDir, 'parameterize input.mol2 --min-type None --no-dihed')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix')
        self._run(refDir, resDir, 'parameterize input.mol2 --no-min --charge-type Gasteiger --scan-type None')
        self._test(refDir, resDir)

    def test_h2o2_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir)

    def test_h2o2_zero_searches_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_zero_searches_restart')
        resDir = os.path.join(self.testDir, 'h2o2_zero_searches_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None '
                                  '--scan-type None --dihed-num-searches 0')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_h2o2_dihed_opt(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt')
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None')
        self._test(refDir, resDir)

    def test_h2o2_dihed_opt_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_dihed_opt_restart')
        resDir = os.path.join(self.testDir, 'h2o2_dihed_opt_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None')
        self._test(refDir, resDir, energyTermRelTol=1e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-5, dihedralPhaseAbsTol=1.1e-3)

    def test_h2o2_full_fake(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake')
        self._run(refDir, resDir, 'parameterize input.mol2 --fake-qm')
        self._test(refDir, resDir)

    def test_h2o2_full_fake_restart(self):
        refDir = os.path.join(self.dataDir, 'h2o2_full_fake_restart')
        resDir = os.path.join(self.testDir, 'h2o2_full_fake_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._run(refDir, resDir, 'parameterize input.mol2 --fake-qm')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethene_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix')
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir)

    def test_ethene_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethene_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'ethene_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_glycol_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix')
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir)

    def test_glycol_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir, energyTermRelTol=3e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-4, dihedralPhaseAbsTol=0.6)

    def test_glycol_dihed_fix_restart_2(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_fix_restart_2')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['O1-C2-C3-O4', 'C2-C3-O4-H10']
        self._run(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --min-type None --scan-type None'.format(' '.join(dihedrals)))
        self._test(refDir, resDir, energyTermRelTol=3e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-4, dihedralPhaseAbsTol=0.6)

    def test_glycol_dihed_select_1_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_1_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_1_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['O1-C2-C3-O4']
        self._run(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --min-type None --scan-type None'.format(' '.join(dihedrals)))
        self._test(refDir, resDir, energyTermRelTol=1e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-5, dihedralPhaseAbsTol=0.5)

    def test_glycol_dihed_select_2_restart(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_select_2_restart')
        resDir = os.path.join(self.testDir, 'glycol_dihed_select_2_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        dihedrals = ['C2-C3-O4-H10']
        self._run(refDir, resDir, 'parameterize input.mol2 -d {} --charge-type Gasteiger --min-type None --scan-type None'.format(' '.join(dihedrals)))
        self._test(refDir, resDir, energyTermRelTol=1e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-5, dihedralPhaseAbsTol=0.5)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_ethanolamine_dihed_fix(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix')
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir)

    def test_ethanolamine_dihed_fix_restart(self):
        refDir = os.path.join(self.dataDir, 'ethanolamine_dihed_fix_restart')
        resDir = os.path.join(self.testDir, 'ethanolamine_dihed_fix_restart')
        shutil.copytree(os.path.join(refDir, 'dihedral-single-point'), os.path.join(resDir, 'dihedral-single-point'))
        self._run(refDir, resDir, 'parameterize input.mol2 --charge-type Gasteiger --min-type None --scan-type None')
        self._test(refDir, resDir, energyTermRelTol=5e-5, energyProfileAbsTol=1.1e-3,
                   dihedralForceConstAbsTol=1e-4, dihedralPhaseAbsTol=2.5)

    def test_benzamidine_gasteiger(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gasteiger')
        resDir = os.path.join(self.testDir, 'benzamidine_gasteiger')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_benzamidine_am1_bcc(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_am1_bcc')
        resDir = os.path.join(self.testDir, 'benzamidine_am1_bcc')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 --charge-type AM1-BCC --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_benzamidine_gaff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_benzamidine_gaff2(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_gaff2')
        resDir = os.path.join(self.testDir, 'benzamidine_gaff2')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_benzamidine_cgenff(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_cgenff')
        resDir = os.path.join(self.testDir, 'benzamidine_cgenff')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff CGenFF_2b6 --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    def test_benzamidine_rtf_prm(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_rtf_prm')
        resDir = os.path.join(self.testDir, 'benzamidine_rtf_prm')
        os.makedirs(resDir, exist_ok=True)
        shutil.copy(os.path.join(refDir, 'input.rtf'), resDir)
        shutil.copy(os.path.join(refDir, 'input.prm'), resDir)
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff CGenFF_2b6 --rtf-prm input.rtf input.prm --charge-type Gasteiger --min-type None --no-dihed')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_VERYLONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full(self):
        assert 'HTMD_CONFIG' in os.environ, '"HTMD_CONFIG" environment variable has to be set'
        refDir = os.path.join(self.dataDir, 'benzamidine_full')
        resDir = os.path.join(self.testDir, 'benzamidine_full')
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --basis 6-31G* -q Slurm')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_LONGTESTS') == 'yes', 'Too long')
    def test_benzamidine_full_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_full_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_full_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --basis 6-31G*')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_esp_freeze_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_esp_freeze_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_esp_freeze_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 --fix-charge N14 --basis 6-31G* --no-dihed')
        self._test(refDir, resDir)

    @unittest.skipUnless(os.environ.get('HTMD_UNSTABLETESTS') == 'yes', 'Unstable')
    def test_benzamidine_dihed_select_restart(self):
        refDir = os.path.join(self.dataDir, 'benzamidine_dihed_select_restart')
        resDir = os.path.join(self.testDir, 'benzamidine_dihed_select_restart')
        shutil.copytree(os.path.join(refDir, 'minimize'), os.path.join(resDir, 'minimize'))
        shutil.copytree(os.path.join(refDir, 'esp'), os.path.join(resDir, 'esp'))
        shutil.copytree(os.path.join(refDir, 'dihedral-opt'), os.path.join(resDir, 'dihedral-opt'))
        self._run(refDir, resDir, 'parameterize input.mol2 -c 1 -ff GAFF2 CGenFF_2b6 -d C2-C1-C7-N13 --basis 6-31G*')
        self._test(refDir, resDir)

    def test_h2o2_dihed_opt_mm_fake(self):
        refDir = os.path.join(self.dataDir, 'h2o2_min_mm_dihed_opt_mm_fake')
        resDir = os.path.join(self.testDir, 'h2o2_min_mm_dihed_opt_mm_fake')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --min-type mm --scan-type mm --fake-qm')
        self._test(refDir, resDir)

    def test_glycol_dihed_opt_mm_fake(self):
        refDir = os.path.join(self.dataDir, 'glycol_dihed_opt_mm_fake')
        resDir = os.path.join(self.testDir, 'glycol_dihed_opt_mm_fake')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --min-type None --scan-type mm --fake-qm --dihed-num-searches 0')
        self._test(refDir, resDir, energyProfileAbsTol=0.01)

    def test_glycol_min_mm_fake(self):
        refDir = os.path.join(self.dataDir, 'glycol_min_mm_fake')
        resDir = os.path.join(self.testDir, 'glycol_min_mm_fake')
        self._run(refDir, resDir, 'parameterize input.mol2 -f GAFF2 --charge-type Gasteiger --min-type mm --scan-type None --fake-qm')
        self._test(refDir, resDir, dihedralForceConstAbsTol=1e-4)

    def test_water_full(self):
        refDir = os.path.join(self.dataDir, 'water_full')
        resDir = os.path.join(self.testDir, 'water_full')
        self._run(refDir, resDir, 'parameterize input.mol2')
        self._test(refDir, resDir)


if __name__ == '__main__':
    unittest.main(verbosity=2)
