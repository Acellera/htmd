# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import os
import sys
import unittest

import numpy as np

from htmd.home import home
from htmdmol.molecule import Molecule
from htmd.qm import Psi4
from htmd.charge import fitGasteigerCharges, fitChargesWithAntechamber, fitESPCharges

logger = logging.getLogger(__name__)


@unittest.skipIf('TRAVIS' in os.environ, 'This is a comparison, not a test!')
class ChargeComparison(unittest.TestCase):

    def _run(self, molName):

        self.molName = molName

        logger.info('Molecule: {}'.format(self.molName))

        molFile = os.path.join(home('test-charge'), self.molName + '.mol2')
        self.mol = Molecule(molFile)

        self.new_mols = {}
        self.extras = {}

        self.new_mols['Gasteiger'] = fitGasteigerCharges(self.mol)
        try:
            self.new_mols['AM1-BCC'] = fitChargesWithAntechamber(self.mol, type='bcc')
        except:
            pass

        qm = Psi4()
        qm.theory = 'B3LYP'
        qm.basis = '6-311++G**'

        workDir = os.path.join('tmp', self.molName)
        os.makedirs(workDir, exist_ok=True)

        for factor in [-10, -5, -4, -3, -2, -1]:
            logger.info('Factor: {}'.format(factor))
            key = 'ESP b {}'.format(factor)
            np.random.seed(20181114)  # Make ESP grid generation deterministic
            self.new_mols[key], self.extras[key] = fitESPCharges(self.mol, qm, workDir, restraint_factor=10**factor)

    def _print(self):

        sys.stdout.flush()
        sys.stderr.flush()

        print('Molecule: ' + self.molName)
        print('Atom ' + ' '.join(map(lambda x: '{:9s}'.format(x), self.new_mols.keys())))
        for i, name in enumerate(self.mol.name):
            print('{:4s}'.format(name), end='')
            for mol in self.new_mols.values():
                print('{:10.4f}'.format(mol.charge[i]), end='')
            print()

    def test_00_methanol(self):
        self._run('methanol')
        self._print()

    def test_01_hydroxylamine(self):
        self._run('hydroxylamine')
        self._print()

    def test_02_urea(self):
        self._run('urea')
        self._print()

    def test_03_isobutanol(self):
        self._run('isobutanol')
        self._print()

    def test_04_methylaldoxime(self):
        self._run('methylaldoxime')
        self._print()

    def test_05_dimethylsulphonamide(self):
        self._run('dimethylsulphonamide')
        self._print()

    def test_06_trisulphate(self):
        self._run('trisulphate')
        self._print()

    def test_07_phenyldiphosphate(self):
        self._run('phenyldiphosphate')
        self._print()


if __name__ == '__main__':
    unittest.main(verbosity=2)
