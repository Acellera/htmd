# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import time
import logging
import unittest
import numpy as np
import nlopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from htmd.numbautil import dihedralAngle
from htmd.ffevaluation.ffevaluate import FFEvaluate

logger = logging.getLogger(__name__)


class DihedralFitting:
    """
    Dihedral parameter fitting

    Capabilities
    ------------
    - Fit parameters from QM energies
    - Fit multiple dihedral angles simultaneously
    - Global parameter optimization

    Attributes
    ----------
    molecule : FFMolecule
        Molecule object
    dihedrals : list of lists
        List of dihedral angles. Each angle is define by 4 atom indices.
    qm_results : List of QMResult
        Reference QM results. The number of results has to be the same to the number of dihedrals.
    result_directory : str
        The directory to save plots
    zerod_parameters : bool
        If set to True, the initial parameter values are set to zeros, else the existing parameters are used as a guess.
    """

    MAX_DIHEDRAL_MULTIPLICITY = 6

    def __init__(self):

        self.molecule = None
        self.dihedrals = []
        self.qm_results = []
        self.result_directory = None
        self.zeroed_parameters = False

        self.parameters = None
        self.loss = None

        self._names = None
        self._equivalent_indices = None

        self._valid_qm_results = None
        self._reference_energies = None
        self._coords = None
        self._angle_values = None

        self._initial_energies = None
        self._target_energies = None
        self._fitted_energies = None

        self._parameterizable_dihedrals = None
        self._parameterizable_dihedral_atomtypes = None

    @property
    def numDihedrals(self):
        """Number of dihedral angles"""
        return len(self.dihedrals)

    def _getValidQMResults(self):
        """
        Get a set of valid QM results
        """

        all_valid_results = []
        for results in self.qm_results:

            # TODO: The test with fake QMResult does not work with this, because there is no molecule
            # dihedral_atomnames = tuple(self.molecule.name[self.dihedrals[self.qm_results.index(results)]])

            # Remove failed QM results
            # TODO print removed QM jobs
            valid_results = [result for result in results if not result.errored]

            # Remove QM results with too high QM energies (>20 kcal/mol above the minimum)
            # TODO print removed QM jobs
            if valid_results:
                qm_min = np.min([result.energy for result in valid_results])
                valid_results = [result for result in valid_results if (result.energy - qm_min) < 20]
            else:
                raise RuntimeError('No valid results.')

            if len(valid_results) < 13:
                raise RuntimeError('Fewer than 13 valid QM points. Not enough to fit.')

            all_valid_results.append(valid_results)

        return all_valid_results

    def _setup(self):

        if len(self.dihedrals) != len(self.qm_results):
            raise ValueError('The number of dihedral and QM result sets has to be the same!')

        # Get dihedral names
        self._names = ['-'.join(self.molecule.name[dihedral]) for dihedral in self.dihedrals]

        # Get equivalent dihedral atom indices
        self._equivalent_indices = []
        for idihed, dihedral in enumerate(self.dihedrals):
            found = False
            for parameterizableDihedral in self._parameterizable_dihedrals:
                if np.all(list(parameterizableDihedral[0]) == dihedral):
                    self._equivalent_indices.append(parameterizableDihedral)
                    found = True
                    break
            if not found:
                raise ValueError('%s is not recognized as a parameterizable dihedral\n' % self._names[idihed])

        # Get reference QM energies and rotamer coordinates
        self._valid_qm_results = self._getValidQMResults()
        self._reference_energies = []
        self._coords = []
        for results in self._valid_qm_results:
            self._reference_energies.append(np.array([result.energy for result in results]))
            self._coords.append([result.coords for result in results])

        # Calculate dihedral angle values
        # [# of dihedrals, # of conformations, # of equivalents]
        self._angle_values = []
        for rotamer_coords, equivalent_indices in zip(self._coords, self._equivalent_indices):
            angle_values = []
            for coords in rotamer_coords:
                angle_values.append([dihedralAngle(coords[indices, :, 0]) for indices in equivalent_indices])
            self._angle_values.append(np.array(angle_values))
        self._angle_values_rad = [angle_values[:, :, None] for angle_values in self._angle_values]

        self._parameterizable_dihedral_atomtypes = [tuple(self.molecule.atomtype[idx]) for idx in self.dihedrals]

        # Calculated initial MM energies
        ff = FFEvaluate(self.molecule, self.parameters)
        self._initial_energies = []
        for rotamer_coords in self._coords:
            self._initial_energies.append(np.array([ff.calculateEnergies(coords[:, :, 0])['total'] for coords in rotamer_coords]))

    def _getBounds(self):
        """
        Get parameter bounds
        """

        nterms = self.MAX_DIHEDRAL_MULTIPLICITY * self.numDihedrals
        lower_bounds = np.zeros(2 * nterms + self.numDihedrals)
        upper_bounds = np.empty_like(lower_bounds)

        # Set force constant and phase bounds
        upper_bounds[:nterms] = 10
        upper_bounds[nterms:2*nterms] = 2 * np.pi

        # Set offset bounds
        lower_bounds[-self.numDihedrals:] = -10
        upper_bounds[-self.numDihedrals:] = 10

        return lower_bounds, upper_bounds

    def _objective(self, x, grad):
        """
        Objective function for the parameter fitting.
        """

        k0, phi0 = np.reshape(x[:-self.numDihedrals], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))
        offset = x[-self.numDihedrals:]

        n = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)

        actual_energies = []
        for i in range(self.numDihedrals):
            phis = self._angle_values[i][:, :, None]
            energies = np.sum(k0[i] * (1 + np.cos(n * phis - phi0[i])), axis=(1, 2)) + offset[i]
            actual_energies.append(energies)

        all_actual_energies = np.concatenate(actual_energies)
        all_target_energies = np.concatenate(self._target_energies)
        rmsd = np.sqrt(np.mean((all_actual_energies - all_target_energies)**2))

        if grad is not None:
            if grad.size > 0:

                grad_k0 = []
                grad_phi0 = []
                grad_offset = []

                for i in range(self.numDihedrals):

                    # Compute partial derivatives
                    phis = self._angle_values[i][:, :, None]
                    dL_dV = (actual_energies[i] - self._target_energies[i])/(rmsd*all_actual_energies.size)
                    dV_dk0 = np.sum(1 + np.cos(n * phis - phi0[i]), axis=1)
                    dV_dphi0 = np.sum(k0[i] * np.sin(n * phis - phi0[i]), axis=1)

                    # Compute gradients with the chain rule
                    grad_k0.append(dL_dV @ dV_dk0)
                    grad_phi0.append(dL_dV @ dV_dphi0)
                    grad_offset.append(np.sum(dL_dV, keepdims=True))

                # Pack gradients
                grad[:] = np.concatenate(grad_k0 + grad_phi0 + grad_offset)

        return rmsd

    def _paramsToVector(self, params, dihedral_atomtypes):
        """
        Convert the parameter objects to a vector.
        """
        vector = []
        for k in dihedral_atomtypes:
            assert len(params.dihedral_types[k]) == self.MAX_DIHEDRAL_MULTIPLICITY
            for term in params.dihedral_types[k]:
                vector.append(term.phi_k)
        for k in dihedral_atomtypes:
            for term in params.dihedral_types[k]:
                vector.append(np.deg2rad(term.phase))
        for i in range(self.numDihedrals):
            vector.append(0)  # The offset

        vector = np.array(vector)

        return vector

    def _optimize_CRS2_LM(self, vector):
        """
        Controlled random search with local mutations
        """

        # Create a global optimizer
        opt = nlopt.opt(nlopt.GN_CRS2_LM, vector.size)
        opt.set_min_objective(self._objective)
        lower_bounds, upper_bounds = self._getBounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        neval = 10000 * opt.get_dimension()  # TODO allow to tune this parameter
        opt.set_maxeval(neval)

        # Optimize parameters
        vector = opt.optimize(vector)  # TODO check optimizer status
        self.loss = opt.last_optimum_value()
        assert self._objective(vector, None) == self.loss

        # Create a local optimizer
        opt = nlopt.opt(nlopt.LN_BOBYQA, opt.get_dimension())
        opt.set_min_objective(self._objective)
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.set_xtol_rel(1e-3)
        opt.set_maxeval(neval)
        opt.set_initial_step(1e-3 * (upper_bounds-lower_bounds))

        # Optimize parameters
        vector = opt.optimize(vector)  # TODO check optimizer status
        self.loss = opt.last_optimum_value()
        assert self._objective(vector, None) == self.loss

        return vector

    def _optimize_random_search(self, vector):
        """
        Naive random search
        """

        # Create a local optimizer
        opt = nlopt.opt(nlopt.LD_LBFGS, vector.size)
        logger.info('Local optimizer: {}'.format(opt.get_algorithm_name()))
        opt.set_min_objective(self._objective)

        # Set bounds
        lower_bounds, upper_bounds = self._getBounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)

        # Set convergence criteria
        opt.set_xtol_rel(1e-3)
        opt.set_maxeval(100 * opt.get_dimension())

        # Initialize
        best_loss = self._objective(vector, None)
        best_vector = vector
        logger.info('Initial RMSD: {:.6f} kcal/mol'.format(best_loss))

        # Naive random search
        niter = 10 * opt.get_dimension()  # TODO allow to tune this parameter
        logger.info('Number of random searches: {}'.format(niter))
        for i in range(niter):

            try:
                vector = opt.optimize(vector)  # TODO check optimizer status
                loss = opt.last_optimum_value()

            except RuntimeError:
                pass

            else:
                if loss < best_loss:
                    best_loss = loss
                    best_vector = vector
                    logger.info('Current RMSD: {:.6f} kcal/mol'.format(best_loss))

            vector = np.random.uniform(low=lower_bounds, high=upper_bounds)

        self.loss = best_loss
        logger.info('Final RMSD: {:.6f} kcal/mol'.format(best_loss))

        return best_vector

    def _vectorToParams(self, parameters, dihedral_atomtypes, vector):
        nparams = len(dihedral_atomtypes) * self.MAX_DIHEDRAL_MULTIPLICITY
        assert vector.size == 2 * nparams + self.numDihedrals

        for i, k in enumerate(dihedral_atomtypes):
            for j, t in enumerate(parameters.dihedral_types[k]):
                t.phi_k = vector[i*self.MAX_DIHEDRAL_MULTIPLICITY+j]
                t.phase = np.rad2deg(vector[i*self.MAX_DIHEDRAL_MULTIPLICITY+j+nparams])

    def _fit(self):

        # Save the initial parameters
        vector = self._paramsToVector(self.parameters, self._parameterizable_dihedral_atomtypes)
        if self.zeroed_parameters:
            vector[:] = 0

        # Evaluate the MM potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # the QM potential and this modified MM potential
        for key in self._parameterizable_dihedral_atomtypes:
            for term in self.parameters.dihedral_types[key]:
                term.phi_k = 0

        # Now evaluate the FF without the dihedral being fitted
        self._target_energies = []
        ff = FFEvaluate(self.molecule, self.parameters)
        for rotamer_coords, ref_energies in zip(self._coords, self._reference_energies):
            energies = ref_energies - np.array([ff.calculateEnergies(coords[:, :, 0])['total'] for coords in rotamer_coords])
            energies -= np.min(energies)
            self._target_energies.append(energies)

        # Optimize the parameters
        logger.info('Start parameter optimization')
        start = time.clock()
        # vector = self._optimize_CRS2_LM(vector)  # TODO this should work better, but it doesn't
        vector = self._optimize_random_search(vector)
        finish = time.clock()
        logger.info('Finished parameter optimization after %f s' % (finish-start))

        # Update the target dihedral with the optimized parameters
        self._vectorToParams(self.parameters, self._parameterizable_dihedral_atomtypes, vector)

        return self.loss

    def _check(self):

        # Evaluate the fitted energies
        self._fitted_energies = []
        ffeval = FFEvaluate(self.molecule, self.parameters)
        for rotamer_coords in self._coords:
            self._fitted_energies.append(np.array([ffeval.calculateEnergies(coords[:, :, 0])['total'] for coords in rotamer_coords]))

        # TODO make the self-consistency test numerically robust
        #reference_energies = np.concatenate([energies - np.mean(energies) for energies in self._reference_energies])
        #fitted_energies = np.concatenate([energies - np.mean(energies) for energies in self._fitted_energies])
        #check_loss = np.sqrt(np.mean((fitted_energies - reference_energies)**2))
        #assert np.isclose(self.loss, check_loss)

        if self.result_directory:
            os.makedirs(self.result_directory, exist_ok=True)
            self.plotConformerEnergies()
            for idihed in range(len(self.dihedrals)):
                self.plotDihedralEnergies(idihed)

    def run(self):

        self._setup()
        self._fit()
        self._check()

        return self.loss

    def plotDihedralEnergies(self, idihed, write_data=True):
        """
        Plot conformer energies for a specific dihedral angle, including QM, original and fitted MM energies.
        """

        angle = np.rad2deg(self._angle_values[idihed][:, 0])
        reference_energy = self._reference_energies[idihed] - np.min(self._reference_energies[idihed])
        initial_energy = self._initial_energies[idihed] - np.min(self._initial_energies[idihed])
        fitted_energy = self._fitted_energies[idihed] - np.min(self._fitted_energies[idihed])
        indices = np.argsort(angle)

        path = os.path.join(self.result_directory, self._names[idihed])

        if write_data:
            fmtsz = 8
            header = ''.join('{:{size}}'.format(s, size=fmtsz) for s in ['# angle', 'QM_ref', 'MM_init', 'MM_fit'])
            data = np.column_stack((angle[indices], reference_energy[indices], initial_energy[indices],
                                    fitted_energy[indices]))
            np.savetxt(path + '.dat', data, fmt='%{size}.3f'.format(size=fmtsz), header=header, comments='')

        plt.figure()
        plt.title(self._names[idihed])
        plt.xlabel('Dihedral angle, deg')
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        plt.ylabel('Energy, kcal/mol')
        plt.plot(angle[indices], reference_energy[indices], 'r-', marker='o', lw=3, label='QM')
        plt.plot(angle[indices], initial_energy[indices], 'g-', marker='o', label='MM initial')
        plt.plot(angle[indices], fitted_energy[indices], 'b-', marker='o', label='MM fitted')
        plt.legend()
        plt.savefig(path + '.svg')
        plt.close()

    def plotConformerEnergies(self, write_data=True):
        """
        Plot all conformer QM energies versus MM energies with the fitted parameters
        """

        qm_energy = np.concatenate(self._reference_energies)[:, None]
        mm_energy = np.concatenate(self._fitted_energies)[:, None]
        qm_energy -= np.min(qm_energy)
        mm_energy -= np.min(mm_energy)

        regression = LinearRegression(fit_intercept=False)
        regression.fit(qm_energy, mm_energy)
        prediction = regression.predict(qm_energy)

        path = os.path.join(self.result_directory, 'conformer-energies')

        if write_data:
            fmtsz = 8
            header = ''.join('{:{size}}'.format(s, size=fmtsz) for s in ['# QM', 'MM'])
            data = np.column_stack((qm_energy, mm_energy))
            np.savetxt(path + '.dat', data, fmt='%{size}.3f'.format(size=fmtsz), header=header, comments='')

        plt.figure()
        plt.title('Conformer Energies MM vs QM')
        plt.xlabel('QM energy, kcal/mol')
        plt.ylabel('MM energy, kcal/mol')
        plt.plot(qm_energy, mm_energy, 'ko')
        plt.plot(qm_energy, prediction, 'r-', lw=2)
        plt.savefig(path + '.svg')
        plt.close()


class TestDihedralFitting(unittest.TestCase):

    def setUp(self):
        self.df = DihedralFitting()

    def test_numDihedrals(self):
        self.df.dihedrals = [[0, 1, 2, 3]]
        self.assertEqual(self.df.numDihedrals, 1)

    def test_getValidQMResults(self):
        from htmd.qm import QMResult

        results = [QMResult() for _ in range(20)]
        for result in results:
            result.energy = 0.
        self.df.qm_results = [results]
        self.assertEqual(len(self.df._getValidQMResults()[0]), 20)

        results[1].errored = True
        results[19].errored = True
        self.assertEqual(len(self.df._getValidQMResults()[0]), 18)

        results[10].energy = -5
        results[12].energy = 12
        results[15].energy = 17
        self.assertEqual(len(self.df._getValidQMResults()[0]), 17)

    def test_getBounds(self):

        for ndihed in range(1, 3):
            with self.subTest(ndihed=ndihed):
                nterm = DihedralFitting.MAX_DIHEDRAL_MULTIPLICITY * ndihed
                self.df.dihedrals = [[0, 0, 0, 0]] * ndihed
                self.assertEqual(ndihed, self.df.numDihedrals)
                lower_bounds, upper_bounds = self.df._getBounds()
                self.assertListEqual(list(lower_bounds), [0] * 2 * nterm + [-10] * ndihed)
                self.assertListEqual(list(upper_bounds), [10] * nterm + [2*np.pi] * nterm + [10] * ndihed)

    def test_paramsToVector(self):
        from parmed.parameters import ParameterSet
        from parmed.topologyobjects import DihedralTypeList, DihedralType

        params = ParameterSet()
        dihlist = DihedralTypeList()
        for i in range(6):
            dihtype = DihedralType(float(i)+10, i+1, float(i)+20)
            dihlist.append(dihtype)
        params.dihedral_types[('x', 'x', 'x', 'x')] = dihlist

        self.df.dihedrals = [(0, 0, 0, 0),]
        vector = self.df._paramsToVector(params, [('x', 'x', 'x', 'x'),])
        self.assertListEqual(list(vector), [10, 11, 12, 13, 14, 15,
                                            np.deg2rad(20), np.deg2rad(21), np.deg2rad(22),
                                            np.deg2rad(23), np.deg2rad(24), np.deg2rad(25), 0.])

    def test_vectorToParams(self):
        from parmed.parameters import ParameterSet
        from parmed.topologyobjects import DihedralTypeList, DihedralType

        params = ParameterSet()
        dihlist = DihedralTypeList()
        for i in range(6):
            dihtype = DihedralType(float(i)+10, i+1, float(i)+20)
            dihlist.append(dihtype)
        params.dihedral_types[('x', 'x', 'x', 'x')] = dihlist

        self.df.dihedrals = [[0, 1, 2, 3]]
        vector = np.array([30, 31, 32, 33, 34, 35, 40, 41, 42, 43, 44, 45, 50])
        self.df._vectorToParams(params, [('x', 'x', 'x', 'x'),], vector)

        self.assertEqual(len(params.dihedral_types[('x', 'x', 'x', 'x')]), 6)
        for i, param in enumerate(params.dihedral_types[('x', 'x', 'x', 'x')]):
            self.assertEqual(param.phi_k, i+30)
            self.assertEqual(param.per, i+1)
            self.assertAlmostEqual(np.deg2rad(param.phase), i+40)

    def test_objective(self):

        from scipy.misc import derivative

        np.random.seed(20181010)

        for ndihed, nequiv, nconf, ref_value in [(1, 1, 1, 372.32948041618585),
                                                 (1, 1, 5, 308.20314159433246),
                                                 (1, 3, 1, 745.73230831710710),
                                                 (2, 1, 1, 368.27031744452563),
                                                 (2, 3, 5, 832.25851847289550)]:
            with self.subTest(ndihed=ndihed, nequiv=nequiv, nconf=nconf):

                self.df.dihedrals = [[0]*4]*ndihed
                self.df._angle_values = 100*np.random.random((ndihed, nconf, nequiv))
                self.df._target_energies = 100*np.random.random((ndihed, nconf))

                vector = 100*np.random.random(13*ndihed)
                grad = np.zeros_like(vector)
                value = self.df._objective(vector, grad)
                self.assertAlmostEqual(ref_value, value)

                for i in range(vector.size):

                    def func(x):
                        v = vector.copy()
                        v[i] = x
                        return self.df._objective(v, None)

                    # Compute gradient numerically
                    ref_grad = derivative(func, vector[i], dx=1e-3, order=5)

                    self.assertAlmostEqual(ref_grad, grad[i])

    # Note: the rest methods are tested indirectly via the "parameterize" tests in test.py


if __name__ == '__main__':

    unittest.main(verbosity=2)
