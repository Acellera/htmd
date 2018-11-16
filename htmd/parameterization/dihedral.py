# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import copy
import logging
import os
import time
import unittest

import matplotlib
matplotlib.use('Agg') # Use Agg to work on Travis

import matplotlib.pyplot as plt
import nlopt
import numpy as np

from sklearn.linear_model import LinearRegression

from htmd.numbautil import dihedralAngle
from htmd.ffevaluation.ffevaluate import FFEvaluate
from htmd.parameterization.parameterset import findDihedralType
from htmd.parameterization.util import detectChiralCenters

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
        self._equivalent_dihedrals = None

        self._valid_qm_results = None
        self._reference_energies = None
        self._coords = None
        self._angle_values = None

        self._initial_energies = None
        self._const_energies = None
        self._target_energies = None
        self._actual_energies = None
        self._fitted_energies = None

        self._parameterizable_dihedrals = None
        self._dihedral_atomtypes = None

        self._plot_directory = None

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

            # Remove results with wrong chiral centers
            new_results = []
            for result in valid_results:

                # Convert QM result into Molecule
                mol = self.molecule.copy()
                mol.coords = result.coords

                # Detect changes of chiral centers
                chiral_centers = detectChiralCenters(mol)[0]
                if self.molecule.chiral_centers == chiral_centers:
                    new_results.append(result)
                else:
                    logger.warning('Rotamer is removed due to a change of chiral centers: '
                                   '{} --> {}'.format(self.molecule.chiral_centers, chiral_centers))
            valid_results = new_results

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
        self._equivalent_dihedrals = []
        for idihed, dihedral in enumerate(self.dihedrals):
            found = False
            for parameterizableDihedral in self._parameterizable_dihedrals:
                if np.all(list(parameterizableDihedral[0]) == dihedral):
                    self._equivalent_dihedrals.append(parameterizableDihedral)
                    found = True
                    break
            if not found:
                raise ValueError('%s is not recognized as a parameterizable dihedral\n' % self._names[idihed])

        # Get dihedral atom types
        self._dihedral_atomtypes = [findDihedralType(tuple(self.molecule.atomtype[dihedral]), self.parameters) for dihedral in self.dihedrals]

        # Get reference QM energies and rotamer coordinates
        self._valid_qm_results = self._getValidQMResults()
        self._reference_energies = []
        self._coords = []
        for results in self._valid_qm_results:
            self._reference_energies.append(np.array([result.energy for result in results]))
            self._coords.append([result.coords for result in results])

        # Calculate dihedral angle values
        # [# of scans, # of dihedrals, # of conformations, # of equivalents]
        self._angle_values = []
        for scan_coords in self._coords:
            scan_angle_values = []
            for equivalent_indices in self._equivalent_dihedrals:
                angle_values = []
                for coords in scan_coords:
                    angle_values.append([dihedralAngle(coords[indices, :, 0]) for indices in equivalent_indices])
                scan_angle_values.append(np.array(angle_values))
            self._angle_values.append(scan_angle_values)

        # Calculated initial MM energies
        ff = FFEvaluate(self.molecule, self.parameters)
        self._initial_energies = []
        for scan_coords in self._coords:
            energies = [ff.calculateEnergies(coords[:, :, 0])['total'] for coords in scan_coords]
            self._initial_energies.append(np.array(energies))

        # Make result directories
        self._plot_directory = os.path.join(self.result_directory, 'plots')
        os.makedirs(self._plot_directory, exist_ok=True)

    def _getBounds(self):
        """
        Get parameter bounds
        """

        nterms = self.MAX_DIHEDRAL_MULTIPLICITY * self.numDihedrals
        lower_bounds = np.zeros(2 * nterms + 1)
        upper_bounds = np.empty_like(lower_bounds)

        # Set force constant and phase bounds
        upper_bounds[:nterms] = 10
        upper_bounds[nterms:2*nterms] = 2 * np.pi

        # Set offset bounds
        lower_bounds[-1] = -20
        upper_bounds[-1] = 20

        return lower_bounds, upper_bounds

    def _objective(self, x, grad=None):
        """
        Objective function for the parameter fitting.
        """

        k0, phi0 = np.reshape(x[:-1], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))
        offset = x[-1]

        n = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)

        self._actual_energies = [0] * self.numDihedrals
        for iscan in range(self.numDihedrals):
            for idihed in range(self.numDihedrals):
                phis = self._angle_values[iscan][idihed][:, :, None]
                energies = np.sum(k0[idihed] * (1 + np.cos(n * phis - phi0[idihed])), axis=(1, 2)) + offset
                self._actual_energies[iscan] += energies

        all_actual_energies = np.concatenate(self._actual_energies)
        all_target_energies = np.concatenate(self._target_energies)
        rmsd = np.sqrt(np.mean((all_actual_energies - all_target_energies)**2))

        if grad is not None:
            if grad.size > 0:

                grad_k0 = [0] * self.numDihedrals
                grad_phi0 = [0] * self.numDihedrals
                grad_offset = 0

                for iscan in range(self.numDihedrals):

                    # Compute partial derivatives
                    dL_dV = self._actual_energies[iscan] - self._target_energies[iscan]
                    dL_dV /= rmsd*all_actual_energies.size

                    for idihed in range(self.numDihedrals):

                        # Compute partial derivatives
                        phis = self._angle_values[iscan][idihed][:, :, None]
                        dV_dk0 = np.sum(1 + np.cos(n * phis - phi0[idihed]), axis=1)
                        dV_dphi0 = np.sum(k0[idihed] * np.sin(n * phis - phi0[idihed]), axis=1)

                        # Compute gradients with the chain rule
                        grad_k0[idihed] += dL_dV @ dV_dk0
                        grad_phi0[idihed] += dL_dV @ dV_dphi0

                    grad_offset += np.sum(dL_dV)

                grad_offset = [[self.numDihedrals * grad_offset]]

                # Pack gradients
                grad[:] = np.concatenate(grad_k0 + grad_phi0 + grad_offset)

        return rmsd

    def _paramsToVector(self, parameters):
        """
        Convert the parameter objects to a vector.
        """
        phi_ks = []
        phases = []

        for atomtypes in self._dihedral_atomtypes:
            assert len(parameters.dihedral_types[atomtypes]) == self.MAX_DIHEDRAL_MULTIPLICITY
            for i, term in enumerate(parameters.dihedral_types[atomtypes]):
                phi_ks.append(term.phi_k)
                phases.append(np.deg2rad(term.phase))
                assert term.per == i + 1  # Check if the periodicity is correct

        return np.array(phi_ks + phases + [0])

    def _vectorToParams(self, vector):
        """
        Convert a vector to a parameter object
        """

        assert vector.size == 2 * len(self._dihedral_atomtypes) * self.MAX_DIHEDRAL_MULTIPLICITY + 1
        phi_k, phase = np.reshape(vector[:-1], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))

        parameters = copy.deepcopy(self.parameters)
        for i, atomtypes in enumerate(self._dihedral_atomtypes):
            assert len(parameters.dihedral_types[atomtypes]) == self.MAX_DIHEDRAL_MULTIPLICITY
            for j, term in enumerate(parameters.dihedral_types[atomtypes]):
                term.phi_k = phi_k[i, j]
                term.phase = np.rad2deg(phase[i, j])
                assert term.per == j + 1  # Check if the periodicity is still correct

        return parameters

    def _optimizeWithRandomSearch(self, vector):
        """
        Naive random search
        """

        # Create a local optimizer
        opt = nlopt.opt(nlopt.LD_LBFGS, vector.size)
        opt.set_vector_storage(opt.get_dimension())
        logger.info('Local optimizer: {} with {} vector storage'.format(opt.get_algorithm_name(),
                                                                        opt.get_vector_storage()))

        # Set bounds
        lower_bounds, upper_bounds = self._getBounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)

        # Set convergence criteria
        opt.set_xtol_rel(1e-3)
        opt.set_maxeval(100 * opt.get_dimension())

        # Initialize
        best_loss = self._objective(vector)
        best_vector = vector
        logger.info('Initial RMSD: {:.6f} kcal/mol'.format(best_loss))
        opt.set_min_objective(self._objective)

        # Naive random search
        niter = 10 * opt.get_dimension()  # TODO allow to tune this parameter
        logger.info('Number of random searches: {}'.format(niter))
        with open(os.path.join(self.result_directory, 'random-search.log'), 'w') as log:
            log.write('{:6s} {:6s} {:10s} {}\n'.format('# Step', 'Status', 'Loss', 'Vector'))
            for i in range(niter):

                try:
                    vector = opt.optimize(vector)
                    loss = opt.last_optimum_value()
                    status = opt.last_optimize_result()

                except RuntimeError:
                    loss = -1
                    status = -1

                else:
                    if loss < best_loss:
                        best_loss = loss
                        best_vector = vector
                        logger.info('Current RMSD: {:.6f} kcal/mol'.format(best_loss))

                string = ' '.join(['{:10.6f}'.format(value) for value in vector])
                log.write('{:6d} {:6d} {:10.6f} {}\n'.format(i, status, loss, string))

                vector = np.random.uniform(low=lower_bounds, high=upper_bounds)

                if best_loss < 1e-3:
                    logger.info('Terminate optization: small RMSD reached!')
                    break

        self.loss = best_loss
        logger.info('Final RMSD: {:.6f} kcal/mol'.format(best_loss))

        return best_vector

    def _evaluateConstTerms(self):
        """
        Evalutate constant MM terms
        """

        parameters = copy.deepcopy(self.parameters)

        # Disable parameterizable (i.e. non-constant) terms
        for atomtypes in self._dihedral_atomtypes:
            for term in parameters.dihedral_types[atomtypes]:
                term.phi_k = 0
                assert term.per > 0 # Guard from messing up with improper dihedrals

        # Evaluate MM energies
        const_energies = []
        ff = FFEvaluate(self.molecule, parameters)
        for scan_coords in self._coords:
            energies = [ff.calculateEnergies(coords[:, :, 0])['total'] for coords in scan_coords]
            const_energies.append(np.array(energies))

        return const_energies

    def _fit(self):

        # Save the initial parameters
        vector = self._paramsToVector(self.parameters)

        # Evaluate constant terms
        self._const_energies = self._evaluateConstTerms()

        # Evaluate target energies for fitting
        self._target_energies = []
        for ref_energies, const_energies in zip(self._reference_energies, self._const_energies):
            energies = ref_energies - const_energies
            self._target_energies.append(energies)

        shift = np.min(np.concatenate(self._target_energies))
        self._target_energies = [energies - shift for energies in self._target_energies]

        # Check self-consistency of computed energies
        self._objective(vector)
        test_energies = zip(self._initial_energies, self._const_energies, self._actual_energies)
        for initial_energies, const_energies, actual_enegies in test_energies:
            assert np.allclose(initial_energies, const_energies + actual_enegies, rtol=0, atol=1e-5) # TODO debug

        # Zero the initial parameters, so they are not used to start the parameter fitting
        if self.zeroed_parameters:
            vector[:] = 0

        # Optimize the parameters
        logger.info('Start parameter optimization')
        start = time.clock()
        vector = self._optimizeWithRandomSearch(vector)
        finish = time.clock()
        logger.info('Finished parameter optimization after %f s' % (finish-start))

        upper_bounds, lower_bounds = self._getBounds()
        if np.isclose(vector[-1], upper_bounds[-1], atol=0.01) or np.isclose(vector[-1], lower_bounds[-1], atol=0.01):
            raise AssertionError('Fitting hit upper/lower bound of the offset. Please report this issue.')

        # Update parameters
        self.parameters = self._vectorToParams(vector)

        return self.loss

    def _check(self):

        # Evaluate the fitted energies
        self._fitted_energies = []
        ffeval = FFEvaluate(self.molecule, self.parameters)
        for scan_coords in self._coords:
            energies = [ffeval.calculateEnergies(coords[:, :, 0])['total'] for coords in scan_coords]
            self._fitted_energies.append(np.array(energies))

        # Check the self-consistency of fitting
        reference_energies = np.concatenate(self._reference_energies)
        reference_energies -= np.mean(reference_energies)
        fitted_energies = np.concatenate(self._fitted_energies)
        fitted_energies -= np.mean(fitted_energies)
        loss = np.sqrt(np.mean((fitted_energies - reference_energies)**2))
        assert np.isclose(self.loss, loss, rtol=0, atol=1e-5)

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

        angle = np.rad2deg(self._angle_values[idihed][idihed][:, 0])
        reference_energy = self._reference_energies[idihed] - np.min(self._reference_energies[idihed])
        initial_energy = self._initial_energies[idihed] - np.min(self._initial_energies[idihed])
        fitted_energy = self._fitted_energies[idihed] - np.min(self._fitted_energies[idihed])
        indices = np.argsort(angle)

        path = os.path.join(self._plot_directory, self._names[idihed])

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

        path = os.path.join(self._plot_directory, 'conformer-energies')

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
                self.assertListEqual(list(lower_bounds), [0] * 2 * nterm + [-20])
                self.assertListEqual(list(upper_bounds), [10] * nterm + [2*np.pi] * nterm + [20])

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
        self.df._dihedral_atomtypes = [('x', 'x', 'x', 'x')]
        vector = self.df._paramsToVector(params)
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
        self.df.parameters = params
        self.df._dihedral_atomtypes = [('x', 'x', 'x', 'x')]
        vector = np.array([30, 31, 32, 33, 34, 35, 40, 41, 42, 43, 44, 45, 50])

        new_params = self.df._vectorToParams(vector)

        self.assertFalse(params is new_params)
        self.assertEqual(len(new_params.dihedral_types[('x', 'x', 'x', 'x')]), 6)
        for i, param in enumerate(new_params.dihedral_types[('x', 'x', 'x', 'x')]):
            self.assertEqual(param.phi_k, i+30)
            self.assertEqual(param.per, i+1)
            self.assertAlmostEqual(np.deg2rad(param.phase), i+40)

    def test_objective(self):

        from scipy.misc import derivative

        np.random.seed(20181010)

        for ndihed, nequiv, nconf, ref_value in [(1, 1, 1, 372.32948041618585),
                                                 (1, 1, 5, 308.20314159433246),
                                                 (1, 3, 1, 745.73230831710710),
                                                 (2, 1, 1, 614.93885227468590),
                                                 (2, 3, 5, 1972.6738129847004)]:
            with self.subTest(ndihed=ndihed, nequiv=nequiv, nconf=nconf):

                self.df.dihedrals = [[0]*4]*ndihed
                self.df._angle_values = 100*np.random.random((ndihed, ndihed, nconf, nequiv))
                self.df._target_energies = 100*np.random.random((ndihed, nconf))

                vector = 100*np.random.random(12*ndihed+1)
                grad = np.zeros_like(vector)
                value = self.df._objective(vector, grad)
                self.assertAlmostEqual(ref_value, value)

                for i in range(vector.size):

                    def func(x):
                        v = vector.copy()
                        v[i] = x
                        return self.df._objective(v)

                    # Compute gradient numerically
                    ref_grad = derivative(func, vector[i], dx=1e-3, order=5)

                    self.assertAlmostEqual(ref_grad, grad[i])

    # Note: the rest methods are tested indirectly via the "parameterize" tests in test.py


if __name__ == '__main__':

    unittest.main(verbosity=2)
