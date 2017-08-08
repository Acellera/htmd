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


class TestParametrize(unittest.TestCase):

    def tests(self):

        testDir = os.path.join('..', 'data', 'test-param')
        if not os.path.exists(testDir):
            testDir = os.path.join('htmd', 'data', 'test-param')
        testDir = os.path.abspath(testDir)

        refDirs = [os.path.join(testDir, test) for test in os.listdir(testDir)]
        refDirs = [test for test in refDirs if os.path.isdir(test)]

        for refDir in refDirs:
            with self.subTest(refDir=refDir), TemporaryDirectory() as resDir:
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

                        if file in ('stdout'):
                            with open(refFile) as ref, open(resFile) as res:
                                refLines, resLines = ref.readlines(), res.readlines()

                                refLines = [line for line in refLines if not re.search('Deprecation warning:', line)]
                                refLines = [line for line in refLines if not re.search('To import all HTMD shortcuts please from now on use', line)]

                                resLines = [line for line in resLines if not re.search('Deprecation warning:', line)]
                                resLines = [line for line in resLines if not re.search('To import all HTMD shortcuts please from now on use', line)]

                            for refLine, resLine in zip(refLines, resLines):
                                if re.search('HTMD version', refLine):
                                    continue
                                self.assertEqual(refLine, resLine, msg=refFile)

                        elif file.endswith('frcmod'):
                            with open(refFile) as ref, open(resFile) as res:
                                # HACK! FRCMOD file lines are swapped randomly, so first sort them and compare.
                                #       Also the first live with the version is removed
                                refLines, resLines = sorted(ref.readlines()[1:]), sorted(res.readlines()[1:])

                            self.assertEqual(refLines, resLines, msg=refFile)

                        else:
                            self.assertTrue(cmp(refFile, resFile, shallow=False), msg=file)


if __name__ == '__main__':
    unittest.main()
