# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import logging
import unittest
import numpy as np
import nlopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from htmd.numbautil import dihedralAngle
from htmd.parameterization.ffevaluate import FFEvaluate

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

        self._all_target_energies = None
        self._angle_values_rad = None

    @property
    def numDihedrals(self):
        """Number of dihedral angles"""

        return len(self.dihedrals)

    def _getEquivalentDihedrals(self, dihedral):
        """
        Find equivalent dihedral angles to the specificied one.
        """
        #TODO maybe this should be moved to FFMolecule

        types = [self.molecule._rtf.type_by_index[index] for index in dihedral]

        all_dihedrals = []
        for dihedral_indices in self.molecule.dihedrals:
            dihedral_types = [self.molecule._rtf.type_by_index[index] for index in dihedral_indices]
            if types == dihedral_types or types == dihedral_types[::-1]:
                all_dihedrals.append(dihedral_indices)

        # Now for each of the uses, remove any which are equivalent
        unique_dihedrals = [dihedral]
        groups = [self.molecule._equivalent_group_by_atom[index] for index in dihedral]
        for dihed in all_dihedrals:
            dihedral_groups = [self.molecule._equivalent_group_by_atom[index] for index in dihed]
            if groups != dihedral_groups and groups != dihedral_groups[::-1]:
                unique_dihedrals.append(dihed)

        return unique_dihedrals

    def _makeDihedralUnique(self, dihedral):
        """
        Duplicate atom types of the dihedral, so its parameters are unique.
        """
        # TODO check symmetry

        # Duplicate the atom types of the dihedral
        for i in range(4):
            self.molecule._duplicateAtomType(dihedral[i])

        equivalent_dihedrals = self._getEquivalentDihedrals(dihedral)
        if len(equivalent_dihedrals) > 1:
            print(dihedral)
            print(len(equivalent_dihedrals))
            print(equivalent_dihedrals)
            raise ValueError("Dihedral term still not unique after duplication")

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
            for rotatableDihedral in self.molecule._rotatable_dihedrals:
                if np.all(rotatableDihedral.atoms == dihedral):
                    self._equivalent_indices.append([dihedral] + rotatableDihedral.equivalents)
                    break
            else:
                raise ValueError('%s is not recognized as a rotable dihedral\n' % self._names[idihed])

        # Get dihedral parameters
        for dihedral in self.dihedrals:
            self._makeDihedralUnique(dihedral)
        all_types = [tuple([self.molecule._rtf.type_by_index[index] for index in dihedral]) for dihedral in self.dihedrals]
        self.parameters = sum([self.molecule._prm.dihedralParam(*types) for types in all_types], [])

        # Get reference QM energies and rotamer coordinates
        self._valid_qm_results = self._getValidQMResults()
        self._reference_energies = []
        self._coords = []
        for results in self._valid_qm_results:
            self._reference_energies.append(np.array([result.energy for result in results]))
            self._coords.append([result.coords for result in results])

        # Calculate dihedral angle values for the fitted equivalent dihedral
        self._angle_values = []
        for rotamer_coords, equivalent_indices in zip(self._coords, self._equivalent_indices):
            angle_values = []
            for coords in rotamer_coords:
                angle_values.append([np.rad2deg(dihedralAngle(coords[indices, :, 0])) for indices in equivalent_indices])
            self._angle_values.append(np.array(angle_values))
        self._angle_values_rad = [np.deg2rad(angle_values)[:, :, None] for angle_values in self._angle_values]

        # Calculated initial MM energies
        ff = FFEvaluate(self.molecule)
        self._initial_energies = []
        for rotamer_coords in self._coords:
            self._initial_energies.append(np.array([ff.run(coords[:, :, 0])['total'] for coords in rotamer_coords]))

    def _getBounds(self):
        """
        Get parameter bounds
        """

        nterms = self.MAX_DIHEDRAL_MULTIPLICITY * self.numDihedrals
        lower_bounds = np.zeros(2 * nterms + self.numDihedrals)
        upper_bounds = np.empty_like(lower_bounds)

        # Set force constant and phase bounds
        upper_bounds[:nterms] = 10
        upper_bounds[nterms:2*nterms] = 360

        # Set offset bounds
        lower_bounds[-self.numDihedrals:] = -10
        upper_bounds[-self.numDihedrals:] = 10

        return lower_bounds, upper_bounds

    def _objective(self, x, _):
        """
        Objective function for the parameter fitting.
        """

        k0, phi0 = np.reshape(x[:-self.numDihedrals], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))
        offset = x[-self.numDihedrals:]

        n = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)
        phi0 = np.deg2rad(phi0)

        all_energies = []
        for i in range(self.numDihedrals):
            phis = self._angle_values_rad[i]
            energies = np.sum(k0[i] * (1 + np.cos(n * phis - phi0[i])), axis=(1, 2)) + offset[i]
            all_energies.append(energies)

        all_energies = np.concatenate(all_energies)
        rmsd = np.sqrt(np.mean((all_energies - self._all_target_energies)**2))

        return rmsd

    def _paramsToVector(self, params):
        """
        Convert the parameter objects to a vector.
        """

        assert [param.n for param in params] == list(range(1, self.MAX_DIHEDRAL_MULTIPLICITY+1)) * self.numDihedrals

        vector = [param.k0 for param in params] + [param.phi0 for param in params] + [0] * self.numDihedrals
        vector = np.array(vector)

        return vector

    def _vectorToParams(self, vector, params):
        """
        Copy results from the vector to the parameter objects.
        """

        nparams = len(params)
        assert vector.size == 2 * nparams + self.numDihedrals

        for i, param in enumerate(params):
            param.k0 = float(vector[i])
            param.phi0 = float(vector[i + nparams])

        return params

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
        opt = nlopt.opt(nlopt.LN_BOBYQA, vector.size)
        opt.set_min_objective(self._objective)
        lower_bounds, upper_bounds = self._getBounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.set_xtol_rel(1e-3)
        opt.set_maxeval(10000 * opt.get_dimension())
        opt.set_initial_step(1e-3 * (upper_bounds - lower_bounds))

        # Optimize the initial vector
        logger.info('Initial RMSD: %f kcal/mol' % self._objective(vector, None))
        best_vector = opt.optimize(vector)  # TODO check optimizer status
        best_loss = opt.last_optimum_value()
        assert self._objective(best_vector, None) == best_loss
        logger.info('Current RMSD: %f kcal/mol' % best_loss)

        # Naive random search
        for i in range(opt.get_dimension()):  # TODO allow to tune this parameter

            # Get random vector and optimize it
            random_vector = np.random.uniform(low=lower_bounds, high=upper_bounds)
            vector = opt.optimize(random_vector)  # TODO check optimizer status

            if opt.last_optimum_value() < best_loss:
                best_loss = opt.last_optimum_value()
                best_vector = vector
                logger.info('Current RMSD: %f kcal/mol' % best_loss)

        self.loss = best_loss

        return best_vector

    def _fit(self):

        # Save the initial parameters
        vector = self._paramsToVector(self.parameters)
        if self.zeroed_parameters:
            vector[:] = 0

        # Evaluate the MM potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # the QM potential and this modified MM potential
        for param in self.parameters:
            param.k0 = 0
        self.molecule._prm.updateDihedral(self.parameters)

        # Now evaluate the FF without the dihedral being fitted
        self._target_energies = []
        ff = FFEvaluate(self.molecule)
        for rotamer_coords, ref_energies in zip(self._coords, self._reference_energies):
            energies = ref_energies - np.array([ff.run(coords[:, :, 0])['total'] for coords in rotamer_coords])
            energies -= np.min(energies)
            self._target_energies.append(energies)
        self._all_target_energies = np.concatenate(self._target_energies)

        # Optimize the parameters
        logger.info('Start parameter optimization')
        # vector = self._optimize_CRS2_LM(vector)  # TODO this should work better, but it doesn't
        vector = self._optimize_random_search(vector)
        logger.info('Final RMSD: %f kcal/mol' % self._objective(vector, None))
        logger.info('Finished parameter optimization')

        # Update the target dihedral with the optimized parameters
        self.parameters = self._vectorToParams(vector, self.parameters)
        self.molecule._prm.updateDihedral(self.parameters)

        return self.loss

    def _check(self):

        # Evaluate the fitted energies
        self._fitted_energies = []
        ffeval = FFEvaluate(self.molecule)
        for rotamer_coords in self._coords:
            self._fitted_energies.append(np.array([ffeval.run(coords[:, :, 0])['total'] for coords in rotamer_coords]))

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

        angle = self._angle_values[idihed][:, 0]
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

    def test_getEquivalentDihedrals(self):

        from htmd.home import home
        from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod

        molFile = os.path.join(home('test-param'), 'glycol.mol2')
        self.df.molecule = FFMolecule(molFile, method=FFTypeMethod.GAFF2)

        self.assertListEqual(self.df._getEquivalentDihedrals([0, 1, 2, 3]), [[0, 1, 2, 3]])
        self.assertListEqual(self.df._getEquivalentDihedrals([4, 0, 1, 2]), [[4, 0, 1, 2]])
        self.assertListEqual(self.df._getEquivalentDihedrals([5, 1, 2, 7]), [[5, 1, 2, 7]])

    def test_makeDihedralUnique(self):

        from htmd.home import home
        from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod

        molFile = os.path.join(home('test-param'), 'glycol.mol2')
        self.df.molecule = FFMolecule(molFile, method=FFTypeMethod.GAFF2)
        types = [self.df.molecule._rtf.type_by_index[i] for i in range(self.df.molecule.numAtoms)]
        self.assertListEqual(types, ['oh', 'c3', 'c3', 'oh', 'ho', 'h1', 'h1', 'h1', 'h1', 'ho'])

        self.df.molecule = FFMolecule(molFile, method=FFTypeMethod.GAFF2)
        self.df._makeDihedralUnique([0, 1, 2, 3])
        types = [self.df.molecule._rtf.type_by_index[i] for i in range(self.df.molecule.numAtoms)]
        self.assertListEqual(types, ['ohx0', 'c3x0', 'c3x0', 'ohx0', 'ho', 'h1', 'h1', 'h1', 'h1', 'ho'])

        self.df.molecule = FFMolecule(molFile, method=FFTypeMethod.GAFF2)
        self.df._makeDihedralUnique([4, 0, 1, 2])
        types = [self.df.molecule._rtf.type_by_index[i] for i in range(self.df.molecule.numAtoms)]
        self.assertListEqual(types, ['ohx0', 'c3x0', 'c3x0', 'ohx0', 'hox0', 'h1', 'h1', 'h1', 'h1', 'hox0'])

        self.df.molecule = FFMolecule(molFile, method=FFTypeMethod.GAFF2)
        self.df._makeDihedralUnique([5, 1, 2, 7])
        types = [self.df.molecule._rtf.type_by_index[i] for i in range(self.df.molecule.numAtoms)]
        self.assertListEqual(types, ['oh', 'c3x0', 'c3x0', 'oh', 'ho', 'h1x0', 'h1x0', 'h1x0', 'h1x0', 'ho'])

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

        self.df.dihedrals = [[0, 1, 2, 3]]
        lower_bounds, upper_bounds = self.df._getBounds()
        self.assertListEqual(list(lower_bounds), [0.,  0.,  0.,  0.,  0.,  0.,   0.,   0.,   0.,   0.,   0.,   0., -10.])
        self.assertListEqual(list(upper_bounds), [10., 10., 10., 10., 10., 10., 360., 360., 360., 360., 360., 360.,  10.])

    def test_paramsToVector(self):

        from htmd.parameterization.ff import TorsPrm

        self.df.dihedrals = [[0, 1, 2, 3]]
        params = [TorsPrm(['x', 'x', 'x', 'x'], k0=float(i)+10, n=i+1, phi0=float(i)+20) for i in range(6)]
        vector = self.df._paramsToVector(params)
        self.assertListEqual(list(vector), [10., 11., 12., 13., 14., 15., 20., 21., 22., 23., 24., 25., 0.])

    def test_vectorToParams(self):

        from htmd.parameterization.ff import TorsPrm

        self.df.dihedrals = [[0, 1, 2, 3]]
        params = [TorsPrm(['x', 'x', 'x', 'x'], k0=float(i)+10, n=i+1, phi0=float(i)+20) for i in range(6)]
        vector = np.array([30., 31., 32., 33., 34., 35., 40., 41., 42., 43., 44., 45., 50.])
        params = self.df._vectorToParams(vector, params)

        self.assertEqual(len(params), 6)
        for i, param in enumerate(params):
            self.assertEqual(param.k0, i+30)
            self.assertEqual(param.n, i+1)
            self.assertEqual(param.phi0, i+40)

    # Note: the rest methods are tested indirectly via the "parameterize" tests in test.py


if __name__ == '__main__':

    unittest.main(verbosity=2)
