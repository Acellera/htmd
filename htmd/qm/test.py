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
from htmd.qm import Psi4, Gaussian
from htmd.queues.localqueue import LocalCPUQueue
from htmd.queues.slurmqueue import SlurmQueue
from htmd.queues.acecloudqueue import AceCloudQueue
from moleculekit.molecule import Molecule
from htmd.numbautil import dihedralAngle


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


class TestBase:

    def setUp(self):

        self.testDir = None

        molFile = os.path.join(home('test-qm'), 'H2-0.74.mol2')
        self.h2_074 = Molecule(molFile)

        molFile = os.path.join(home('test-qm'), 'H2-1.00.mol2')
        self.h2_100 = Molecule(molFile)

        molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
        self.h2o2_90 = Molecule(molFile)

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
            with TemporaryDirectory(dir=self.testDir) as tmpDir:
                self.qm.molecule = self.h2_074
                self.qm.multiplicity = multiplicity
                self.qm.theory = 'HF'
                self.qm.basis = '3-21G'
                self.qm.directory = tmpDir
                self.qm.charge = charge
                result = self.qm.run()[0]
                self.assertFalse(result.errored, msg=(charge, multiplicity))
                self.assertAlmostEqual(REF_CHARGE_MULTIPLICITY_ENERGIES[charge][multiplicity],
                                       result.energy, msg=(charge, multiplicity))

    def test_theories(self):

        for theory in self.qm.THEORIES:
            with TemporaryDirectory(dir=self.testDir) as tmpDir:
                self.qm.molecule = self.h2_074
                self.qm.theory = theory
                self.qm.basis = '3-21G'
                self.qm.directory = tmpDir
                result = self.qm.run()[0]
                self.assertFalse(result.errored, msg=theory)
                self.assertAlmostEqual(REF_THEORY_ENERGIES[theory], result.energy, msg=theory)

    def test_corrections(self):

        for correction in self.qm.CORRECTIONS:
            with TemporaryDirectory(dir=self.testDir) as tmpDir:
                self.qm.molecule = self.h2_074
                self.qm.theory = 'BLYP' # Using BLYP as HF-D isn't available
                self.qm.correction = correction
                self.qm.basis = '3-21G'
                self.qm.directory = tmpDir
                result = self.qm.run()[0]
                self.assertFalse(result.errored, msg=correction)
                self.assertAlmostEqual(REF_CORRECTION_ENERGIES[correction], result.energy, msg=correction)

    def test_basis_sets(self):

        for basis in self.qm.BASIS_SETS:
            with TemporaryDirectory(dir=self.testDir) as tmpDir:
                self.qm.molecule = self.h2_074
                self.qm.theory = 'HF'
                self.qm.basis = basis
                self.qm.directory = tmpDir
                result = self.qm.run()[0]
                self.assertFalse(result.errored, msg=basis)
                places = 5 if basis in ('cc-pVTZ', 'aug-cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVQZ') else 7 # Large basis sets are unstable
                self.assertAlmostEqual(REF_BASIS_ENERGIES[basis], result.energy, places=places, msg=basis)

    def test_solvents(self):

        for solvent in self.qm.SOLVENTS:
            with TemporaryDirectory(dir=self.testDir) as tmpDir:
                self.qm.molecule = self.h2_074
                self.qm.theory = 'HF'
                self.qm.basis = '3-21G'
                self.qm.solvent = solvent
                self.qm.directory = tmpDir
                result = self.qm.run()[0]
                self.assertFalse(result.errored, msg=solvent)
                self.assertAlmostEqual(REF_SOLVET_ENERGIES[solvent], result.energy, msg=solvent)

    @unittest.skip(reason='joblib 0.11 breaks it on travis?')  # TODO: bring back when joblib back to 0.12
    def test_properties(self):

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.molecule = self.h2_074
            self.qm.theory = 'HF'
            self.qm.basis = '3-21G'
            self.qm.esp_points = REF_ESP_POINTS
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            self.assertTrue(np.all(np.isclose(REF_DIPOLE, result.dipole)))
            self.assertTrue(np.all(np.isclose(REF_QUADRUPOLE, result.quadrupole)))
            self.assertTrue(np.all(np.isclose(REF_MULLIKEN, result.mulliken)))
            self.assertTrue(np.all(np.isclose(REF_ESP_VALUES, result.esp_values)))

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
            self.assertTrue(np.all(np.isclose(REF_OPT_COORDS, result.coords)))

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.optimize = False
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            self.assertTrue(np.all(np.isclose(REF_INIT_COORDS, result.coords)))

    def test_restrained_dihedrals(self):

        quad = [2, 0, 1, 3]
        angle = np.rad2deg(dihedralAngle(self.h2o2_90.coords[quad, :, 0]))
        self.assertAlmostEqual(89.999544881803772, angle, places=5)

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
            self.assertAlmostEqual(89.999541178019271, angle, places=5)

        with TemporaryDirectory(dir=self.testDir) as tmpDir:
            self.qm.restrained_dihedrals = None
            self.qm.directory = tmpDir
            result = self.qm.run()[0]
            self.assertFalse(result.errored)
            angle = np.rad2deg(dihedralAngle(result.coords[quad, :, 0]))
            self.assertAlmostEqual(179.51690845119924, angle, places=3) # Unstable results

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
                self.assertAlmostEqual(REF_MULTISTRUCTURE_ENERGIES[ires], result.energy, msg=ires)


class TestPsi4Local(TestBase, unittest.TestCase):

    def setUp(self):

        if os.environ.get('TRAVIS_OS_NAME') == 'osx':
            self.skipTest('Psi4 does not work on Mac')  # TODO fix!

        self.qm = Psi4()
        super().setUp()


class TestPsi4Slurm(TestBase, unittest.TestCase):

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


class TestPsi4AceCloud(TestBase, unittest.TestCase):

    def setUp(self):

        # TODO finish
        self.skipTest('No AceCloud tests')

        self.qm = Psi4()
        self.qm.queue = AceCloudQueue()
        super().setUp()


class TestGaussian(TestBase, unittest.TestCase):

    def setUp(self):

        # TODO finish
        self.skipTest('No Gaussian, folk!')

        self.qm = Gaussian()
        super().setUp()


if __name__ == '__main__':
    unittest.main(verbosity=2)
