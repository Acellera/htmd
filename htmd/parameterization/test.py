# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import shutil
import unittest
from tempfile import TemporaryDirectory
from subprocess import call
from filecmp import cmp

os.environ['HTMD_NONINTERACTIVE'] = '1'

class TestParametrize(unittest.TestCase):

    def setUp(self):

        self.dataDir = os.path.join('..', 'data', 'test-param')
        if not os.path.exists(self.dataDir):
            self.dataDir = os.path.join('htmd', 'data', 'test-param')

        suffix = 'travis' if 'TRAVIS' in os.environ else 'local'
        self.dataDir = os.path.abspath(os.path.join(self.dataDir, suffix))

    def _test(self, refDir):

        with TemporaryDirectory() as resDir:
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

                    with self.subTest(file=file):
                        self.assertTrue(os.path.exists(refFile))
                        self.assertTrue(os.path.exists(resFile))

                        if file in ('stdout'):
                            with open(refFile) as ref, open(resFile) as res:
                                refLines, resLines = ref.readlines(), res.readlines()
                            self.assertListEqual(refLines, resLines, msg=file)

                        elif file.endswith('frcmod'):
                            with open(refFile) as ref, open(resFile) as res:
                                # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                                #       Also the first line with the version is removed
                                refLines = sorted(ref.readlines()[1:])
                                resLines = sorted(res.readlines()[1:])
                            self.assertListEqual(refLines, resLines, msg=file)

                        else:
                            self.assertTrue(cmp(refFile, resFile, shallow=False), msg=file)

    def test_doc(self):

        self._test(os.path.join(self.dataDir, 'doc'))

    def test_h2o2_list(self):

        self._test(os.path.join(self.dataDir, 'list'))

    def test_h2o2_gaff2(self):

        self._test(os.path.join(self.dataDir, 'gaff2'))


if __name__ == '__main__':
    unittest.main()
