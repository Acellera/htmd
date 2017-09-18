# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import numpy as np
import nlopt
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from htmd.molecule.util import dihedralAngle
from htmd.parameterization.ffevaluate import FFEvaluate


class DihedralFitting:
    """
    Dihedral parameter fitting from QM data
    """

    MAX_DIHEDRAL_MULTIPLICITY = 6

    def __init__(self):

        self.molecule = None
        self.dihedrals = []
        self.qm_results = []
        self.result_directory = None

        self.paramters = None
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
    def num_dihedrals(self):
        return len(self.dihedrals)

    def _countUsesOfDihedral(self, dihedral):
        """
        Return the number of uses of the dihedral specified by the types of the 4 atom indices
        """

        types = [self.molecule._rtf.type_by_index[index] for index in dihedral]

        all_uses = []
        for dihedral_indices in self.molecule.dihedrals:
            dihedral_types = [self.molecule._rtf.type_by_index[index] for index in dihedral_indices]
            if types == dihedral_types or types == dihedral_types[::-1]:
                all_uses.append(dihedral_indices)

        # Now for each of the uses, remove any which are equivalent
        unique_uses = [dihedral]
        groups = [self.molecule._equivalent_group_by_atom[index] for index in dihedral]
        for dihed in all_uses:
            dihedral_groups = [self.molecule._equivalent_group_by_atom[index] for index in dihed]
            if groups != dihedral_groups and groups != dihedral_groups[::-1]:
                unique_uses.append(dihed)

        return len(unique_uses), unique_uses

    def _makeDihedralUnique(self, dihedral):
        """
        Create a new type for (arbitrarily) a middle atom of the dihedral, so that the dihedral we are going to modify
        is unique
        """
        # TODO check symmetry

        # Duplicate the dihedrals types so this modified term is unique
        for i in range(4):
            if not ("x" in self.molecule._rtf.type_by_index[dihedral[i]]):
                self.molecule.duplicateTypeOfAtom(dihedral[i])

        number_of_uses, uses = self._countUsesOfDihedral(dihedral)
        if number_of_uses > 1:
            print(dihedral)
            print(number_of_uses)
            print(uses)
            raise ValueError("Dihedral term still not unique after duplication")

    def _get_valid_qm_results(self):

        all_valid_results = []

        for results in self.qm_results:

            # Remove failed QM results
            # TODO print removed QM jobs
            valid_results = [result for result in results if not result.errored]

            # Remove QM results with too high QM energies (>20 kcal/mol above the minimum)
            # TODO print removed QM jobs
            qm_min = np.min([result.energy for result in valid_results])
            valid_results = [result for result in valid_results if (result.energy - qm_min) < 20]

            if len(valid_results) < 13:
                raise RuntimeError("Fewer than 13 valid QM points. Not enough to fit!")

            all_valid_results.append(valid_results)

        return all_valid_results

    def _setup(self):

        if len(self.dihedrals) != len(self.qm_results):
            raise ValueError('The number of dihedral and QM result sets has to be the same!')

        # Get dihedral names
        self._names = ['%s-%s-%s-%s' % tuple(self.molecule.name[dihedral]) for dihedral in self.dihedrals]

        # Get equivalent dihedral atom indices
        self._equivalent_indices = []
        for idihed, dihedral in enumerate(self.dihedrals):
            for rotableDihedral in self.molecule._soft_dihedrals:
                if np.all(rotableDihedral.atoms == dihedral):
                    self._equivalent_indices.append([dihedral] + rotableDihedral.equivalents)
                    break
            else:
                raise ValueError('%s is not recognized as a rotable dihedral\n' % self._names[idihed])

        # Get dihedral parameters
        for dihedral in self.dihedrals:
            self._makeDihedralUnique(dihedral)
        all_types = [tuple([self.molecule._rtf.type_by_index[index] for index in dihedral]) for dihedral in self.dihedrals]
        self.parameters = sum([self.molecule._prm.dihedralParam(*types) for types in all_types], [])

        # Get reference QM energies and rotamer coordinates
        self._valid_qm_results = self._get_valid_qm_results()
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
                angle_values.append([dihedralAngle(coords[indices, :, 0]) for indices in equivalent_indices])
            self._angle_values.append(np.array(angle_values))
        self._angle_values_rad = [np.deg2rad(angle_values)[:, :, None] for angle_values in self._angle_values]

        # Calculated initial MM energies
        ff = FFEvaluate(self.molecule)
        self._initial_energies = []
        for rotamer_coords in self._coords:
            self._initial_energies.append(np.array([ff.run(coords[:, :, 0])['total'] for coords in rotamer_coords]))

    def _get_bounds(self):

        nterms = self.MAX_DIHEDRAL_MULTIPLICITY*self.num_dihedrals
        lower_bounds = np.zeros(2*nterms + self.num_dihedrals)
        upper_bounds = np.empty_like(lower_bounds)

        # Set force constant and phase bounds
        upper_bounds[:nterms] = 10
        multiplicities = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)
        lower_bounds[nterms:2*nterms] = -180/np.tile(multiplicities, self.num_dihedrals)
        upper_bounds[nterms:2*nterms] =  180/np.tile(multiplicities, self.num_dihedrals)

        # Set offset bounds
        lower_bounds[-self.num_dihedrals:] = -10
        upper_bounds[-self.num_dihedrals:] = 10

        return lower_bounds, upper_bounds

    def _objective(self, x, _):
        """
        Evaluate the objective function of the torsion with the input params for each of the phi's poses
        """

        k0, phi0 = np.reshape(x[:-self.num_dihedrals], (2, -1, self.MAX_DIHEDRAL_MULTIPLICITY))
        offset = x[-self.num_dihedrals:]

        n = np.arange(1, self.MAX_DIHEDRAL_MULTIPLICITY + 1)
        phi0 = np.deg2rad(phi0)

        all_energies = []
        for i in range(self.num_dihedrals):
            phis = self._angle_values_rad[i]
            energies = np.sum(k0[i]*(1 + np.cos(n*phis - phi0[i])), axis=(1, 2)) + offset[i]
            all_energies.append(energies)

        all_energies = np.concatenate(all_energies)
        rmsd = np.sqrt(np.mean((all_energies - self._all_target_energies)**2))

        return rmsd

    def _params_to_vector(self, params):

        assert [param.n for param in params] == list(range(1, self.MAX_DIHEDRAL_MULTIPLICITY+1))*self.num_dihedrals

        vector = [param.k0 for param in params] + [param.phi0 for param in params] + [0]*self.num_dihedrals
        vector = np.array(vector)

        return vector

    def _vector_to_params(self, vector, params):

        nparams = len(params)
        assert vector.size == 2*nparams + self.num_dihedrals

        for i, param in enumerate(params):
            param.k0 = float(vector[i])
            param.phi0 = float(vector[i + nparams])

        return params

    def _fit(self):

        # TODO do not modify the molecule

        # Save the initial parameters as the best ones
        vector = self._params_to_vector(self.parameters)

        # Evalaute the mm potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # The QM potential and the this modified mm potential
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

        # Create a global optimizer
        opt = nlopt.opt(nlopt.GN_CRS2_LM, vector.size)
        opt.set_min_objective(self._objective)
        lower_bounds, upper_bounds = self._get_bounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        neval = 10000*opt.get_dimension()  # TODO allow to tune this parameter
        opt.set_maxeval(neval)

        # Optimize parameters
        vector = opt.optimize(vector.copy())
        self.loss = opt.last_optimum_value()
        assert self._objective(vector, None) == self.loss
        # TODO check optimizer status

        opt = nlopt.opt(nlopt.LN_BOBYQA, opt.get_dimension())
        opt.set_min_objective(self._objective)
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.set_xtol_rel(1e-4)
        opt.set_maxeval(neval)
        opt.set_initial_step(1e-4*(upper_bounds-lower_bounds))

        vector = opt.optimize(vector)
        self.loss = opt.last_optimum_value()
        assert self._objective(vector, None) == self.loss
        # TODO check optimizer status

        # Update the target dihedral with the optimized parameters
        self.parameters = self._vector_to_params(vector, self.parameters)
        self.molecule._prm.updateDihedral(self.parameters)

        return self.loss

    def _check(self):

        # Evaluate the fitted energies
        self._fitted_energies = []
        ffeval = FFEvaluate(self.molecule)
        for rotamer_coords in self._coords:
            self._fitted_energies.append(np.array([ffeval.run(coords[:, :, 0])['total'] for coords in rotamer_coords]))

        reference_energies = np.concatenate([energies - np.mean(energies) for energies in self._reference_energies])
        fitted_energies = np.concatenate([energies - np.mean(energies) for energies in self._fitted_energies])
        check_loss = np.sqrt(np.mean((fitted_energies - reference_energies)**2))
        assert np.isclose(self.loss, check_loss)

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

    def plotDihedralEnergies(self, idihed):

        angle = self._angle_values[idihed][:, 0]
        reference_energy = self._reference_energies[idihed] - np.min(self._reference_energies[idihed])
        initial_energy = self._initial_energies[idihed] - np.min(self._initial_energies[idihed])
        fitted_energy = self._fitted_energies[idihed] - np.min(self._fitted_energies[idihed])
        indices = np.argsort(angle)

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
        plt.savefig(os.path.join(self.result_directory, self._names[idihed] + '.svg'))
        plt.close()

    def plotConformerEnergies(self):

        qm_energy = np.concatenate(self._reference_energies)[:, None]
        mm_energy = np.concatenate(self._fitted_energies)[:, None]
        qm_energy -= np.min(qm_energy)
        mm_energy -= np.min(mm_energy)

        regression = LinearRegression(fit_intercept=False)
        regression.fit(qm_energy, mm_energy)
        prediction = regression.predict(qm_energy)

        plt.figure()
        plt.title('Conformer Energies MM vs QM')
        plt.xlabel('QM energy, kcal/mol')
        plt.ylabel('MM energy, kcal/mol')
        plt.plot(qm_energy, mm_energy, 'ko')
        plt.plot(qm_energy, prediction, 'r-', lw=2)
        plt.savefig(os.path.join(self.result_directory, 'conformer-energies.svg'))
        plt.close()


if __name__ == '__main__':
    pass