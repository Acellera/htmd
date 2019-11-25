# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import unittest
from tempfile import TemporaryDirectory
import numpy as np

from htmd.home import home
from htmd.qm.base import QMBase
from htmd.qm import Psi4, TeraChem, Gaussian
from htmd.queues.localqueue import LocalCPUQueue
from htmd.queues.slurmqueue import SlurmQueue
from moleculekit.molecule import Molecule
from moleculekit.dihedral import dihedralAngle


# For H2 0.74 A with HF and 3-21G basis
REF_CHARGE_MULTIPLICITY_ENERGIES = {}
REF_CHARGE_MULTIPLICITY_ENERGIES[-1] = {}
REF_CHARGE_MULTIPLICITY_ENERGIES[-1][2] = -542.9961117250292
REF_CHARGE_MULTIPLICITY_ENERGIES[0] = {}
REF_CHARGE_MULTIPLICITY_ENERGIES[0][1] = -704.6590602770854
REF_CHARGE_MULTIPLICITY_ENERGIES[0][3] = -469.197302783436
REF_CHARGE_MULTIPLICITY_ENERGIES[1] = {}
REF_CHARGE_MULTIPLICITY_ENERGIES[1][2] = -347.51346071758957

# For H2 0.74 A with 3-21G basis
REF_THEORY_ENERGIES = {}
REF_THEORY_ENERGIES['HF']     = -704.6590602770854
REF_THEORY_ENERGIES['BLYP']   = -727.8771297555771
REF_THEORY_ENERGIES['PBE']    = -726.1688096094861
REF_THEORY_ENERGIES['B3LYP']  = -734.5593143769217
REF_THEORY_ENERGIES['PBE0']   = -727.6595793319559
REF_THEORY_ENERGIES['B2PLYP'] = -725.3434000281864
REF_THEORY_ENERGIES['wB97']   = -728.9046507790767
REF_THEORY_ENERGIES['wB97X']  = -729.8489202970206
REF_THEORY_ENERGIES['wB97X-D'] = -732.6189385642774

# For H2 0.74 A with BLYP and 3-21G basis
REF_CORRECTION_ENERGIES = {}
REF_CORRECTION_ENERGIES['none'] = -727.8771297555771
REF_CORRECTION_ENERGIES['D']    = -727.8779482287340
REF_CORRECTION_ENERGIES['D3']   = -727.8771485808614

# For H2 0.74 A with HF
REF_BASIS_ENERGIES = {}
REF_BASIS_ENERGIES['3-21G']       = -704.6590602770854

REF_BASIS_ENERGIES['6-31G']       = -707.0503436069422
REF_BASIS_ENERGIES['6-31G*']      = -707.0503436069422
REF_BASIS_ENERGIES['6-31G**']     = -709.8983936135076
REF_BASIS_ENERGIES['6-31+G']      = -707.0503436069422
REF_BASIS_ENERGIES['6-31+G*']     = -707.0503436069422
REF_BASIS_ENERGIES['6-31+G**']    = -709.8983936135076
REF_BASIS_ENERGIES['6-31++G']     = -707.1131282296159
REF_BASIS_ENERGIES['6-31++G*']    = -707.1131282296159
REF_BASIS_ENERGIES['6-31++G**']   = -709.9444092785507

REF_BASIS_ENERGIES['6-311G']      = -707.8237480933639
REF_BASIS_ENERGIES['6-311G*']     = -707.8237480933639
REF_BASIS_ENERGIES['6-311G**']    = -710.6400209800165
REF_BASIS_ENERGIES['6-311+G']     = -707.8237480933639
REF_BASIS_ENERGIES['6-311+G*']    = -707.8237480933639
REF_BASIS_ENERGIES['6-311+G**']   = -710.6400209800165
REF_BASIS_ENERGIES['6-311++G']    = -707.8239227198378
REF_BASIS_ENERGIES['6-311++G*']   = -707.8239227198378
REF_BASIS_ENERGIES['6-311++G**']  = -710.6500303917903

REF_BASIS_ENERGIES['cc-pVDZ']     = -708.2705913710219
REF_BASIS_ENERGIES['cc-pVTZ']     = -710.9482984703429
REF_BASIS_ENERGIES['cc-pVQZ']     = -711.2616851315097

REF_BASIS_ENERGIES['aug-cc-pVDZ'] = -708.3194926157898
REF_BASIS_ENERGIES['aug-cc-pVTZ'] = -710.9898924187445
REF_BASIS_ENERGIES['aug-cc-pVQZ'] = -711.2704032224416

# For H2 0.74 A with HF and 3-21G basis
REF_SOLVET_ENERGIES = {}
REF_SOLVET_ENERGIES['vacuum'] = -704.6590602770854
REF_SOLVET_ENERGIES['PCM']    = -704.7382351554221

# For H2 0.74 A with HF and 3-21G basis
REF_DIPOLE = [0.0, 0.0, 0.0, 0.0]
REF_QUADRUPOLE = [-0.1430, -0.1430, 0.2861, 0.0, 0.0, 0.0]
REF_MULLIKEN = [0.0, 0.0]
REF_ESP_POINTS = np.zeros((10, 3))
REF_ESP_POINTS[:, 2] = np.linspace(-2, 2, 10)
REF_ESP_VALUES = [0.01282272, 0.03737972, 0.17778391, 1.59595817, 4.89089454,
                  4.89089454, 1.59595817, 0.17778391, 0.03737972, 0.01282272]

# For H2 1.00 A with BLYP and 3-21G basis
REF_INIT_COORDS = [[[0], [0], [-0.5]],
                   [[0], [0], [ 0.5]]]
REF_OPT_COORDS  = [[[0], [0], [-0.37537082]],
                   [[0], [0], [ 0.37537082]]]

# For 0.74 and 1.00 A H2 with HF and 3-21G basis
REF_MULTISTRUCTURE_ENERGIES = {}
REF_MULTISTRUCTURE_ENERGIES[0] = -704.6590602770854
REF_MULTISTRUCTURE_ENERGIES[1] = -684.8589617509701

# For Br with HF
# Note: basis sets are substituted
REF_BR_ENERGIES = {}
REF_BR_ENERGIES['3-21G']       = -1614005.9194762693
REF_BR_ENERGIES['6-31+G*']     = -1614007.7511329607
REF_BR_ENERGIES['6-311++G**']  = -1614204.3952243670
REF_BR_ENERGIES['cc-pVDZ']     = -1614005.9194762693
REF_BR_ENERGIES['aug-cc-pVTZ'] = -1614204.3952243670


class _TestBase:

    def assertEqualFloat(self, a, b, tol=1e-10, msg=None):
        message = f'{a} != {b} within rtol = {tol}'
        message = f'{message} : {msg}' if msg else message
        if np.abs(a) < tol:
            self.assertTrue(np.isclose(a, b, atol=tol, rtol=0), msg=message)
        else:
            self.assertTrue(np.isclose(a, b, atol=0, rtol=tol), msg=message)

    def assertEqualFloatList(self, a, b, tol=1e-10, msg=None):

        a = np.array(a).flatten()
        b = np.array(b).flatten()
        self.assertEqual(a.size, b.size, msg=msg)
        for a_, b_ in zip(a, b):
            self.assertEqualFloat(a_, b_, tol=tol, msg=msg)

    def setUp(self):

        self.testDir = None

        molFile = os.path.join(home('test-qm'), 'H2-0.74.mol2')
        self.h2_074 = Molecule(molFile)

        molFile = os.path.join(home('test-qm'), 'H2-1.00.mol2')
        self.h2_100 = Molecule(molFile)

        molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
        self.h2o2_90 = Molecule(molFile)

        molFile = os.path.join(home('test-qm'), 'Br.mol2')
        self.Br = Molecule(molFile)

        self.e_tol = 1e-5 if isinstance(self.qm, TeraChem) else 1e-10

    def test_type(self):

        self.assertIsInstance(self.qm, QMBase)

    def test_defaults(self):

        self.assertEqual(self.qm.molecule, None)
        self.assertEqual(self.qm.multiplicity, 1)
        self.assertEqual(self.qm.theory, 'B3LYP')
        self.assertEqual(self.qm.correction, 'none')
        self.assertEqual(self.qm.basis, '6-31G*')
        self.assertEqual(self.qm.solvent, 'vacuum')
        self.assertEqual(self.qm.esp_points, None)
        self.assertEqual(self.qm.optimize, False)
        self.assertEqual(self.qm.restrained_dihedrals, None)
        self.assertEqual(self.qm.directory, '.')

    def test_queue(self):

        self.assertIsInstance(self.qm.queue, LocalCPUQueue)

    def test_multiplicity(self):

        for charge, multiplicity in ((-1, 2), (0, 1), (0, 3), (1, 2)):
            with self.subTest(charge=charge, multiplicity=multiplicity):

                if isinstance(self.qm, TeraChem):
                    if charge == 0 and multiplicity == 3:
                        self.skipTest('TeraChem bug')
                    if charge == 1 and multiplicity == 2:
                        self.skipTest('TeraChem bug')

                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.h2_074
                    self.qm.multiplicity = multiplicity
                    self.qm.theory = 'HF'
                    self.qm.basis = '3-21G'
                    self.qm.directory = tmpDir
                    self.qm.charge = charge
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored, msg=(charge, multiplicity))
                    self.assertEqualFloat(REF_CHARGE_MULTIPLICITY_ENERGIES[charge][multiplicity],
                                          result.energy, tol=self.e_tol, msg=(charge, multiplicity))

    def test_theories(self):

        for theory in self.qm.THEORIES:
            with self.subTest(theory=theory):

                if isinstance(self.qm, TeraChem) and theory in ('B2PLYP', 'wB97X-D'):
                    self.skipTest(f'TeraChem does not support: {theory}')

                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.h2_074
                    self.qm.theory = theory
                    self.qm.basis = '3-21G'
                    self.qm.directory = tmpDir
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored, msg=theory)
                    self.e_tol = 1e-4 if isinstance(self.qm, TeraChem) and theory in ('wB97', 'wB97X') else self.e_tol
                    self.assertEqualFloat(REF_THEORY_ENERGIES[theory], result.energy, tol=self.e_tol, msg=theory)

    def test_corrections(self):

        for correction in self.qm.CORRECTIONS:
            with self.subTest(correction=correction):
                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.h2_074
                    self.qm.theory = 'BLYP' # Using BLYP as HF-D isn't available
                    self.qm.correction = correction
                    self.qm.basis = '3-21G'
                    self.qm.directory = tmpDir
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored, msg=correction)
                    self.e_tol = 2e-4 if isinstance(self.qm, TeraChem) and correction == 'D3' else self.e_tol
                    self.assertEqualFloat(REF_CORRECTION_ENERGIES[correction], result.energy, tol=self.e_tol, msg=correction)

    def test_basis_sets(self):

        for basis in self.qm.BASIS_SETS:
            with self.subTest(basis=basis):

                if isinstance(self.qm, TeraChem) and basis in ('cc-pVTZ', 'aug-cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVQZ'):
                    self.skipTest(f'TeraChem does not support: {basis}')

                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.h2_074
                    self.qm.theory = 'HF'
                    self.qm.basis = basis
                    self.qm.directory = tmpDir
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored, msg=basis)
                    # Large basis sets are unstable
                    self.e_tol = 100 * self.e_tol if basis in ('cc-pVTZ', 'aug-cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVQZ') else self.e_tol
                    self.assertEqualFloat(REF_BASIS_ENERGIES[basis], result.energy, tol=self.e_tol, msg=basis)

    def test_solvents(self):

        for solvent in self.qm.SOLVENTS:
            with self.subTest(solvent=solvent):
                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.h2_074
                    self.qm.theory = 'HF'
                    self.qm.basis = '3-21G'
                    self.qm.solvent = solvent
                    self.qm.directory = tmpDir
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored, msg=solvent)
                    self.e_tol = 5e-5 if isinstance(self.qm, TeraChem) and solvent == 'PCM' else self.e_tol
                    self.assertEqualFloat(REF_SOLVET_ENERGIES[solvent], result.energy, tol=self.e_tol, msg=solvent)

    def test_properties(self):

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.molecule = self.h2_074
            self.qm.theory = 'HF'
            self.qm.basis = '3-21G'
            if isinstance(self.qm, Psi4):
                self.qm.esp_points = REF_ESP_POINTS
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            tol = 1e-5 if isinstance(self.qm, TeraChem) else 1e-10
            self.assertEqualFloatList(REF_DIPOLE, result.dipole, tol=tol)
            self.assertEqualFloatList(REF_MULLIKEN, result.mulliken, tol=tol)
            if isinstance(self.qm, Psi4):
                self.assertEqualFloatList(REF_QUADRUPOLE, result.quadrupole)
                self.assertEqualFloatList(REF_ESP_VALUES, result.esp_values)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.esp_points = None
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            self.assertEqual(None, result.esp_values)

    def test_optimization(self):

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.molecule = self.h2_100
            self.qm.theory = 'BLYP' # HF fails
            self.qm.basis = '3-21G'
            self.qm.optimize = True
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            tol = 1e-4 if isinstance(self.qm, TeraChem) else 1e-5
            self.assertEqualFloatList(REF_OPT_COORDS, result.coords, tol=tol)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.optimize = False
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            self.assertEqualFloatList(REF_INIT_COORDS, result.coords)

    def test_restrained_dihedrals(self):

        quad = [2, 0, 1, 3]
        angle = np.rad2deg(dihedralAngle(self.h2o2_90.coords[quad, :, 0]))
        self.assertEqualFloat(89.999544881803772, angle, tol=1e-7)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.molecule = self.h2o2_90
            self.qm.theory = 'BLYP' # HF fails
            self.qm.basis = '3-21G'
            self.qm.optimize = True
            self.qm.restrained_dihedrals = np.array([quad])
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            angle = np.rad2deg(dihedralAngle(result.coords[quad, :, 0]))
            self.assertEqualFloat(89.999541178019271, angle, tol=1e-7)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.restrained_dihedrals = None
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            angle = np.rad2deg(dihedralAngle(result.coords[quad, :, 0]))
            if isinstance(self.qm, Psi4):
                self.assertEqualFloat(179.51690845119924, angle, tol=1e-6) # Unstable results
            else:
                self.assertEqualFloat(-168.9488713666722, angle, tol=1e-5) # Unstable results

    def test_directory(self):

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            tmpDir2 = os.path.join(tmpDir, 'test')
            self.qm.molecule = self.h2_074
            self.qm.theory = 'HF'
            self.qm.basis = '3-21G'
            self.qm.directory = tmpDir2
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            self.assertTrue(os.path.exists(os.path.join(tmpDir2, '00000', 'run.sh')))

    def test_multistructure(self):

        mol = self.h2_074
        mol.appendFrames(self.h2_100)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.molecule = mol
            self.qm.theory = 'HF'
            self.qm.basis = '3-21G'
            self.qm.directory = tmpDir
            results = self.qm.run()
            self.assertEqual(2, len(results))
            for ires, result in enumerate(results):
                self.assertFalse(result.errored, msg=ires)
                self.assertEqualFloat(REF_MULTISTRUCTURE_ENERGIES[ires], result.energy, tol=self.e_tol, msg=ires)

    def test_basis_set_substitution(self):

        if isinstance(self.qm, TeraChem):
            self.skipTest('Not implemented')

        for basis in ('3-21G', '6-31+G*', '6-311++G**', 'cc-pVDZ', 'aug-cc-pVTZ'):
            with self.subTest(basis=basis):
                with TemporaryDirectory(dir=self.testDir) as tmpDir:
                    self.qm.molecule = self.Br
                    self.qm.multiplicity = 2
                    self.qm.theory = 'HF'
                    self.qm.basis = basis
                    self.qm.directory = tmpDir
                    result = self.qm.run()[0]
                    self.assertFalse(result.errored)
                    self.assertEqualFloat(REF_BR_ENERGIES[basis], result.energy, tol=self.e_tol)


class _TestPsi4Local(_TestBase, unittest.TestCase):

    def setUp(self):

        self.qm = Psi4()
        super().setUp()


class _TestPsi4Slurm(_TestBase, unittest.TestCase):

    def setUp(self):

        self.skipTest('No Slurm tests')

        if 'TRAVIS' in os.environ:
           self.skipTest('No Psi4 Slurm tests on Travis')

        # For Slurm, the test directory has to be on a shared filesystem
        self.testDir = os.getcwd()

        self.qm = Psi4()
        self.qm.queue = SlurmQueue()
        self.qm.queue.partition = 'normalGPU' # TODO use CPU partition when available
        self.qm.queue.jobname = 'Psi4_test'
        super().setUp()

    def test_queue(self):

        self.assertIsInstance(self.qm.queue, SlurmQueue)


class _TestTeraChemLocal(_TestBase, unittest.TestCase):

    def setUp(self):

        if 'TRAVIS' in os.environ:
           self.skipTest('No TeraChem tests on Travis')

        self.qm = TeraChem()
        super().setUp()


class _TestGaussian(_TestBase, unittest.TestCase):

    def setUp(self):

        # TODO finish
        self.skipTest('No Gaussian, folk!')

        self.qm = Gaussian()
        super().setUp()


if __name__ == '__main__':
    unittest.main(verbosity=2)
