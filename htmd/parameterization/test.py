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
        self.dataDir = os.path.abspath(self.dataDir)

    def test_doc(self):

        self._test(os.path.join(self.dataDir, 'doc'))

    def test_h2o2_list(self):

        self._test(os.path.join(self.dataDir, 'list'))

    def test_h2o2_gaff2(self):

        self._test(os.path.join(self.dataDir, 'gaff2'))

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

                    with self.subTest(file=file):
                        refFile = os.path.join(refDir, file)
                        resFile = os.path.join(resDir, file)

                        if file in ('stdout'):
                            with open(refFile) as ref:
                                refLines = []
                                for line in ref.readlines():
                                    if re.search('HTMD: Logging setup failed', line):
                                        continue
                                    if re.search('Number of CPUs to use \(default:', line):
                                        continue
                                    if re.search('ncpus:', line):
                                        continue
                                    refLines.append(line)

                            with open(resFile) as res:
                                resLines = []
                                for line in res.readlines():
                                    if re.search('HTMD: Logging setup failed', line):
                                        continue
                                    if re.search('Number of CPUs to use \(default:', line):
                                        continue
                                    if re.search('ncpus:', line):
                                        continue
                                    resLines.append(line)

                            self.assertListEqual(refLines, resLines, msg=file)

                        elif file.endswith('frcmod'):
                            with open(refFile) as ref, open(resFile) as res:
                                # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                                #       Also the first live with the version is removed
                                refLines = sorted(ref.readlines()[1:])
                                resLines = sorted(res.readlines()[1:])

                            self.assertListEqual(refLines, resLines, msg=file)

                        else:
                            self.assertTrue(cmp(refFile, resFile, shallow=False), msg=file)


if __name__ == '__main__':
    unittest.main()
