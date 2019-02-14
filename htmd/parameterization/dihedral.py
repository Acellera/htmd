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

from htmd.ffevaluation.ffevaluate import FFEvaluate
from htmd.numbautil import dihedralAngle
from htmd.parameterization.detect import detectParameterizableDihedrals
from htmd.parameterization.parameterset import findDihedralType

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
        self.num_iterations = 3
        self.fit_type = 'Iterative'

        self.parameters = None
        self.loss = None

        self._names = None

        self._reference_energies = None
        self._coords = None
        self._angle_values = None

        self._initial_energies = None
        self._fitted_energies = None

        self._dihedral_atomtypes = None

    @property
    def numDihedrals(self):
        """Number of dihedral angles"""
        return len(self.dihedrals)

    def _setup(self):

        if len(self.dihedrals) != len(self.qm_results):
            raise ValueError('The number of dihedral and QM result sets has to be the same!')

        # Get dihedral names
        self._names = ['-'.join(self.molecule.name[dihedral]) for dihedral in self.dihedrals]

        # Get all equivalent dihedrals
        all_equivalent_dihedrals = detectParameterizableDihedrals(self.molecule)
        all_equivalent_dihedrals = {tuple(dihedrals[0]): dihedrals for dihedrals in all_equivalent_dihedrals}

        # Choose the selected dihedrals
        _equivalent_dihedrals = []
        for dihedral, name in zip(self.dihedrals, self._names):
            if tuple(dihedral) not in all_equivalent_dihedrals:
                raise ValueError('{} is not a parameterizable dihedral!'.format(name))
            _equivalent_dihedrals.append(all_equivalent_dihedrals[tuple(dihedral)])

        # Get dihedral atom types
        self._dihedral_atomtypes = [findDihedralType(tuple(self.molecule.atomtype[dihedral]), self.parameters) for dihedral in self.dihedrals]

        # Get reference QM energies and rotamer coordinates
        self._reference_energies = []
        self._coords = []
        for results in self.qm_results:
            self._reference_energies.append(np.array([result.energy for result in results]))
            self._coords.append([result.coords for result in results])

        # Calculate dihedral angle values
        # [# of scans, # of dihedrals, # of conformations, # of equivalents]
        self._angle_values = []
        for scan_coords in self._coords:
            scan_angle_values = []
            for equivalent_indices in _equivalent_dihedrals:
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
        os.makedirs(self.result_directory, exist_ok=True)

    def _getBounds(self, numDihedrals):
        """
        Get parameter bounds
        """

        nterms = self.MAX_DIHEDRAL_MULTIPLICITY * numDihedrals
        lower_bounds = np.zeros(2 * nterms + 1)
        upper_bounds = np.empty_like(lower_bounds)

        # Set force constant and phase bounds
        upper_bounds[:nterms] = 10
        upper_bounds[nterms:2*nterms] = 2 * np.pi

        # Set offset bounds
        lower_bounds[-1] = -25
        upper_bounds[-1] = 25

        return lower_bounds, upper_bounds

    def _unpackVector(self, vector):
        k0, phi0 = np.reshape(vector[:-1], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))
        offset = vector[-1]
        n = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)
        return k0, phi0, offset, n

    @staticmethod
    def _evaluateScanEnergies(k0, phi0, offset, n, angle_values, numScans, dihedrals):
        actual_energies = []

        for iscan in range(numScans):
            current_energy = 0
            for i, idihed in enumerate(dihedrals):
                phis = angle_values[iscan][idihed][:, :, None]
                energies = np.sum(k0[i] * (1 + np.cos(n * phis - phi0[i])), axis=(1, 2)) + offset
                current_energy += energies
            actual_energies.append(current_energy)

        return actual_energies

    def _objective(self, x, grad, angle_values, target_energies, idihed=None):
        """
        Objective function for the parameter fitting.
        """
        dihedrals = np.arange(self.numDihedrals) if idihed is None else [idihed,]

        k0, phi0, offset, n = self._unpackVector(x)
        actual_energies = DihedralFitting._evaluateScanEnergies(k0, phi0, offset, n, angle_values, self.numDihedrals, dihedrals)

        all_actual_energies = np.concatenate([actual_energies[d] for d in dihedrals])
        all_target_energies = np.concatenate([target_energies[d] for d in dihedrals])
        rmsd = np.sqrt(np.mean((all_actual_energies - all_target_energies)**2))

        if grad is not None:
            if grad.size > 0:

                grad_k0 = [0] * len(dihedrals)
                grad_phi0 = [0] * len(dihedrals)
                grad_offset = 0

                for iscan in dihedrals:

                    # Compute partial derivatives
                    dL_dV = actual_energies[iscan] - target_energies[iscan]
                    dL_dV /= rmsd*all_actual_energies.size

                    for k, idihed in enumerate(dihedrals):
                        # Compute partial derivatives
                        phis = angle_values[iscan][idihed][:, :, None]
                        dV_dk0 = np.sum(1 + np.cos(n * phis - phi0[k]), axis=1)
                        dV_dphi0 = np.sum(k0[k] * np.sin(n * phis - phi0[k]), axis=1)

                        # Compute gradients with the chain rule
                        grad_k0[k] += dL_dV @ dV_dk0
                        grad_phi0[k] += dL_dV @ dV_dphi0

                    grad_offset += np.sum(dL_dV)

                grad_offset = [[len(dihedrals) * grad_offset]]

                # Pack gradients
                grad[:] = np.concatenate(grad_k0 + grad_phi0 + grad_offset)

        return rmsd

    def _paramsToVector(self, parameters, idihed=None):
        """
        Convert the parameter objects to a vector.
        """
        phi_ks = []
        phases = []

        dihedrals = np.arange(self.numDihedrals) if idihed is None else [idihed,]

        for idih in dihedrals:
            atomtypes = self._dihedral_atomtypes[idih]
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

    def _optimizeWithIterativeScheme(self, vector, target_energies):
        best_loss = None
        best_vector = None
        for i in range(self.num_iterations):
            logger.info('Dihedral fitting iteration {}'.format(i+1))

            # Evaluate constant terms for single dihedral fitting
            current_params = self._vectorToParams(vector)
            const_energies_dihedral = self._evaluateConstTermsPerDihedral(current_params)
            target_energies_dihedral = self._calculateTargetEnergies(const_energies_dihedral)

            # Fit the individual dihedrals
            dihedral_vectors = []
            for idihed in range(self.numDihedrals):
                logger.info('Optimizing dihedral {}/{}'.format(idihed+1, self.numDihedrals))
                subvector = self._paramsToVector(current_params, idihed)
                subvector = self._optimizeWithRandomSearch(subvector, target_energies_dihedral, num_searches=20, idihed=idihed, _logger=False)
                dihedral_vectors.append(subvector)

            # Pack the individually fitted parameters into a vector
            vector = DihedralFitting._subvectorsToVector(dihedral_vectors, vector)

            # Perform global optimization of all dihedrals
            logger.info('Globally optimizing all dihedrals')
            vector = self._optimizeWithRandomSearch(vector, target_energies, best_vector=best_vector, best_loss=best_loss, num_searches=1)
            best_vector = vector.copy()
            best_loss = self._objective(vector, None, self._angle_values, target_energies)

        return vector

    def _optimizeWithRandomSearch(self, vector, target_energies, idihed=None, convergence_rmsd=1e-3,
                                  num_searches=None, best_vector=None, best_loss=None, _logger=True):
        """
        Naive random search
        """
        numDihedrals = self.numDihedrals if idihed is None else 1

        # Create a local optimizer
        opt = nlopt.opt(nlopt.LD_LBFGS, vector.size)
        opt.set_vector_storage(opt.get_dimension())

        # Set bounds
        lower_bounds, upper_bounds = self._getBounds(numDihedrals)
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)

        # Set convergence criteria
        opt.set_xtol_rel(1e-4)
        opt.set_maxeval(100 * opt.get_dimension())

        # Initialize
        if best_loss is None:
            best_loss = self._objective(vector, None, self._angle_values, target_energies, idihed)
        if best_vector is None:
            best_vector = vector
        if _logger:
            logger.info('Initial RMSD: {:.6f} kcal/mol'.format(best_loss))

        opt.set_min_objective(lambda x, grad: self._objective(x, grad, self._angle_values, target_energies, idihed))

        # Decide the number of the random searches
        num_searches = 10 * opt.get_dimension() if num_searches is None else int(num_searches)
        if num_searches < 0:
            raise ValueError('The number of random searches has to be possive, but it is {}'.format(num_searches))

        with open(os.path.join(self.result_directory, 'random-search.log'), 'w') as log:
            log.write('{:6s} {:6s} {:10s} {}\n'.format('# Step', 'Status', 'Loss', 'Vector'))
            for i in range(num_searches):

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
                        if _logger:
                            logger.info('Current RMSD: {:.6f} kcal/mol'.format(best_loss))

                string = ' '.join(['{:10.6f}'.format(value) for value in vector])
                log.write('{:6d} {:6d} {:10.6f} {}\n'.format(i, status, loss, string))

                vector = np.random.uniform(low=lower_bounds, high=upper_bounds)

                if best_loss < convergence_rmsd:
                    if _logger:
                        logger.info('Terminate optimization: small RMSD reached!')
                    break

        self.loss = best_loss
        if _logger:
            logger.info('Best RMSD: {:.6f} kcal/mol'.format(best_loss))

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

    def _evaluateConstTermsPerDihedral(self, parameters):
        # for key in parameters.dihedral_types: print(key, [term.phi_k for term in parameters.dihedral_types[key]])
        # Evaluate MM energies
        const_energies = []
        for idihed in range(self.numDihedrals):
            # Zero the current dihedral values to calculate all other "constant" terms
            modified_parameters = copy.deepcopy(parameters)
            atomtypes = self._dihedral_atomtypes[idihed]
            for term in modified_parameters.dihedral_types[atomtypes]:
                term.phi_k = 0
                assert term.per > 0 # Guard from messing up with improp
                    
            ff = FFEvaluate(self.molecule, modified_parameters)
            scan_coords = self._coords[idihed]
            energies = [ff.calculateEnergies(coords[:, :, 0])['total'] for coords in scan_coords]
            const_energies.append(np.array(energies))

        return const_energies

    @staticmethod
    def _subvectorsToVector(dihedral_vectors, prev_vector):
        halfsize = int((len(dihedral_vectors[0])-1)/2)
        vector = np.hstack([vv[:halfsize] for vv in dihedral_vectors] +
                        [vv[halfsize:-1] for vv in dihedral_vectors] +
                        [prev_vector[-1]]) # Keep last offset from previous iteration
        return vector

    def _calculateTargetEnergies(self, const_energies):
        target_energies = []
        for ref_energies, const_energies_dih in zip(self._reference_energies, const_energies):
            target_energies.append(ref_energies - const_energies_dih)
        shift = np.min(np.concatenate(target_energies))
        target_energies = [energies - shift for energies in target_energies]
        return target_energies

    def _fit(self):
        # Save the initial parameters
        vector = self._paramsToVector(self.parameters)

        # Evaluate constant terms
        const_energies = self._evaluateConstTerms()

        # Evaluate target energies for fitting
        target_energies = self._calculateTargetEnergies(const_energies)

        # Check self-consistency of computed energies
        k0, phi0, offset, n = self._unpackVector(vector)
        actual_energies = DihedralFitting._evaluateScanEnergies(k0, phi0, offset, n, self._angle_values, self.numDihedrals, np.arange(self.numDihedrals))
        test_energies = zip(self._initial_energies, const_energies, actual_energies)
        for initial_energies, const_energies, actual_enegies in test_energies:
            assert np.allclose(initial_energies, const_energies + actual_enegies, rtol=0, atol=1e-5) # TODO debug

        # Zero the initial parameters, so they are not used to start the parameter fitting
        if self.zeroed_parameters:
            vector[:] = 0

        # Optimize each dihedral individually
        logger.info('Start parameter optimization with {} iterations'.format(self.num_iterations))
        start = time.clock()
        if self.fit_type.lower() == 'iterative':
            vector = self._optimizeWithIterativeScheme(vector, target_energies)
        elif self.fit_type.lower() == 'nrs':
            vector = self._optimizeWithRandomSearch(vector, target_energies)                            
        finish = time.clock()
        logger.info('Finished parameter optimization after %f s' % (finish-start))

        upper_bounds, lower_bounds = self._getBounds(self.numDihedrals)
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
        # HACK: without searches, the offset is not computed. So the test will not pass!
        if self.num_iterations != 0:
            assert np.isclose(self.loss, loss, rtol=0, atol=1e-2)

    def run(self):

        self._setup()
        self._fit()
        self._check()

        return self.parameters

    def plotDihedralEnergies(self, idihed, directory='.', ref_name = 'Ref', write_data=True):
        """
        Plot conformer energies for a specific dihedral angle, including QM, original and fitted MM energies.
        """

        path = os.path.join(directory, self._names[idihed])

        # Get data
        angle_values = self._angle_values[idihed][idihed][:, 0]
        reference_energies = self._reference_energies[idihed]
        initial_energies = self._initial_energies[idihed]
        fitted_energies = self._fitted_energies[idihed]

        # Convert and offset data
        angle_values = np.rad2deg(angle_values)
        reference_energies -= np.min(reference_energies)
        initial_energies -= np.min(initial_energies)
        fitted_energies -= np.min(fitted_energies)

        # Sort data
        indices = np.argsort(angle_values)
        angle_values = angle_values[indices]
        reference_energies = reference_energies[indices]
        initial_energies = initial_energies[indices]
        fitted_energies = fitted_energies[indices]

        if write_data:
            fmtsz = 8
            header = ''.join('{:{size}}'.format(s, size=fmtsz) for s in ['# angle', 'ref', 'MM_init', 'MM_fit'])
            data = np.column_stack((angle_values, reference_energies, initial_energies, fitted_energies))
            np.savetxt(path + '.dat', data, fmt='%{size}.3f'.format(size=fmtsz), header=header, comments='')

        # Impose periodic boundaries
        angle_values = np.concatenate([[angle_values[-1] - 360], angle_values, [angle_values[0] + 360]])
        reference_energies = np.concatenate([[reference_energies[-1]], reference_energies, [reference_energies[0]]])
        initial_energies = np.concatenate([[initial_energies[-1]], initial_energies, [initial_energies[0]]])
        fitted_energies = np.concatenate([[fitted_energies[-1]], fitted_energies, [fitted_energies[0]]])

        plt.figure()
        plt.title('Dihedral angle: {}'.format(self._names[idihed]))
        plt.xlabel('Dihedral angle [deg]')
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        plt.ylabel('Energy [kcal/mol]')
        plt.plot(angle_values, reference_energies, 'r-', marker='o', lw=3, label=ref_name)
        plt.plot(angle_values, initial_energies, 'g-', marker='o', lw=2, label='MM initial')
        plt.plot(angle_values, fitted_energies, 'b-', marker='o', lw=2, label='MM fitted')
        plt.legend()
        plt.savefig(path + '.svg')
        plt.close()

    def plotConformerEnergies(self, directory='.', ref_name='Ref', write_data=True):
        """
        Plot all conformer QM energies versus MM energies with the fitted parameters
        """

        path = os.path.join(directory, 'conformer-energies')

        # Get data
        qm_energy = np.concatenate(self._reference_energies)[:, None]
        mm_energy = np.concatenate(self._fitted_energies)[:, None]

        # Offset data
        qm_energy -= np.min(qm_energy)
        mm_energy -= np.min(mm_energy)

        # Fit a linear regression
        regression = LinearRegression(fit_intercept=False)
        regression.fit(qm_energy, mm_energy)
        prediction = regression.predict(qm_energy)

        if write_data:
            fmtsz = 8
            header = ''.join('{:{size}}'.format(s, size=fmtsz) for s in ['# ref', 'MM'])
            data = np.column_stack((qm_energy, mm_energy))
            np.savetxt(path + '.dat', data, fmt='%{size}.3f'.format(size=fmtsz), header=header, comments='')

        plt.figure()
        plt.title('{} vs MM energies'.format(ref_name))
        plt.xlabel('{} energy [kcal/mol]'.format(ref_name))
        plt.ylabel('MM energy [kcal/mol]')
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

    def test_getBounds(self):

        for ndihed in range(1, 3):
            with self.subTest(ndihed=ndihed):
                nterm = DihedralFitting.MAX_DIHEDRAL_MULTIPLICITY * ndihed
                self.df.dihedrals = [[0, 0, 0, 0]] * ndihed
                self.assertEqual(ndihed, self.df.numDihedrals)
                lower_bounds, upper_bounds = self.df._getBounds(ndihed)
                self.assertListEqual(list(lower_bounds), [0] * 2 * nterm + [-25])
                self.assertListEqual(list(upper_bounds), [10] * nterm + [2*np.pi] * nterm + [25])

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
                angle_values = 100*np.random.random((ndihed, ndihed, nconf, nequiv))
                target_energies = 100*np.random.random((ndihed, nconf))

                vector = 100*np.random.random(12*ndihed+1)
                grad = np.zeros_like(vector)
                value = self.df._objective(vector, grad, angle_values, target_energies)
                self.assertAlmostEqual(ref_value, value)

                for i in range(vector.size):

                    def func(x):
                        v = vector.copy()
                        v[i] = x
                        return self.df._objective(v, None, angle_values, target_energies)

                    # Compute gradient numerically
                    ref_grad = derivative(func, vector[i], dx=1e-3, order=5)

                    self.assertAlmostEqual(ref_grad, grad[i])

    # Note: the rest methods are tested indirectly via the "parameterize" tests in test.py


if __name__ == '__main__':

    unittest.main(verbosity=2)
