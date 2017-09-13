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

    def __init__(self):

        self.molecule = None
        self.dihedrals = []
        self.qm_results = []
        self.result_directory = None

        self._names = None
        self._types = None
        self._parameters = None
        self._offsets = None
        self._equivalent_indices = None

        self._valid_qm_results = None
        self._reference_energies = None

        self._coords = None
        self._angle_values = None

        self._initial_energies = None
        self._target_energies = None
        self._fitted_energies = None

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

        assert len(self.dihedrals) == len(self.qm_results)

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

        # Get dihedral atom types
        for dihedral in self.dihedrals:
            self._makeDihedralUnique(dihedral)
        self._types = [tuple([self.molecule._rtf.type_by_index[index] for index in dihedral]) for dihedral in self.dihedrals]

        # Get dihedral parameters
        self._parameters = [self.molecule._prm.dihedralParam(*types) for types in self._types]
        self._offsets = np.zeros(len(self.dihedrals))

        # Get valid QM results
        self._valid_qm_results = self._get_valid_qm_results()

        # Get reference QM energies
        self._reference_energies = []
        for results in self._valid_qm_results:
            energies = np.array([result.energy for result in results])
            energies -= np.mean(energies)
            self._reference_energies.append(energies)

        # Get rotamer coordinates
        self._coords = []
        for results in self._valid_qm_results:
            self._coords.append([result.coords for result in results])

        # Calculate dihedral angle values
        self._angle_values = []
        for idihed, rotamer_coords in enumerate(self._coords):
            angle_values = []
            for coords in rotamer_coords:
                angles = [dihedralAngle(coords[indices, :, 0]) for indices in self._equivalent_indices[idihed]]
                angle_values.append(angles)
            self._angle_values.append(np.array(angle_values))

        # Get initial MM energies
        ff = FFEvaluate(self.molecule)
        self._initial_energies = []
        for rotamer_coords in self._coords:
            energies = np.array([ff.run(coords[:, :, 0])['total'] for coords in rotamer_coords])
            energies -= np.mean(energies)
            self._initial_energies.append(energies)

        # Create a directory for results, i.e. plots
        if self.result_directory:
            os.makedirs(self.result_directory, exist_ok=True)

    def _get_bounds(self):

        lower_bounds = np.zeros(13)
        lower_bounds[-1] = -10

        upper_bounds = np.zeros(13)
        upper_bounds[0:6] = 20
        upper_bounds[6:12] = 180/np.arange(1, 7)
        upper_bounds[-1] = 10

        return lower_bounds, upper_bounds

    def _objective(self, x, idihed):
        """
        Evaluate the objective function of the torsion with the input params for each of the phi's poses
        """

        k0 = x[0:6]
        phi0 = np.deg2rad(x[6:12])
        offset = x[12]

        n = np.arange(1, 7)
        phis = np.deg2rad(self._angle_values[idihed])[:, :, None]  # rotamers x equivalent dihedral values

        energies = np.sum(k0 * (1. + np.cos(n * phis - phi0)), axis=(1, 2)) + offset
        chisq = np.sum((energies - self._target_energies[idihed])**2)

        return chisq

    @staticmethod
    def _params_to_vector(params, offset):

        vector = [param.k0 for param in params]
        vector += [param.phi0 for param in params]
        vector += [offset]
        vector = np.array(vector)

        return vector

    @staticmethod
    def _vector_to_params(vector, params):

        nparams = len(params)
        assert vector.size == 2*nparams + 1

        for i, param in enumerate(params):
            param.k0 = float(vector[i])
            param.phi0 = float(vector[i + nparams])

        offset = float(vector[-1])

        return params, offset

    def _fitDihedral(self, idihed):

        # Save the initial parameters as the best ones
        best_params = self._params_to_vector(self._parameters[idihed], self._offsets[idihed])

        # Evalaute the mm potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # The QM potential and the this modified mm potential
        for param in self._parameters[idihed]:
            param.k0 = 0
        self.molecule._prm.updateDihedral(self._parameters[idihed])

        # Now evaluate the FF without the dihedral being fitted
        ff = FFEvaluate(self.molecule)
        self._target_energies[idihed] = -np.array([ff.run(coords[:, :, 0])['total'] for coords in self._coords[idihed]])
        self._target_energies[idihed] += self._reference_energies[idihed]
        self._target_energies[idihed] -= np.min(self._target_energies[idihed])

        # Create global optimizer
        opt = nlopt.opt(nlopt.G_MLSL_LDS, best_params.size)
        opt.set_min_objective(lambda x, _: self._objective(x, idihed))
        lower_bounds, upper_bounds = self._get_bounds()
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.set_maxeval(int(1e4)*opt.get_dimension())
        opt.set_population(opt.get_dimension())

        # Create local optimizer
        local_opt = nlopt.opt(nlopt.LN_BOBYQA, opt.get_dimension())
        local_opt.set_xtol_rel(1e-6)
        local_opt.set_initial_step(1e-1)
        opt.set_local_optimizer(local_opt)

        # Calculate the best loss
        best_chisq = self._objective(best_params, idihed)

        # Optimize parameters
        initial = np.random.uniform(low=lower_bounds, high=upper_bounds)
        xopt = opt.optimize(initial)
        chisq = opt.last_optimum_value()

        if chisq < best_chisq:
            best_chisq = chisq
            best_params = xopt.copy()

        # Update the target dihedral with the optimized parameters
        params, self._offsets[idihed] = self._vector_to_params(best_params, self._parameters[idihed])
        self.molecule._prm.updateDihedral(params)

        # Finally evaluate the fitted potential
        ffeval = FFEvaluate(self.molecule)
        energies = np.array([ffeval.run(coords[:, :, 0])['total'] for coords in self._coords[idihed]])
        energies -= np.mean(energies)
        final_chisq = np.sum((energies - self._reference_energies[idihed])**2)

        assert np.isclose(best_chisq, final_chisq)

        return final_chisq, energies

    def run(self):

        self._setup()

        self._target_energies = [None]*len(self.dihedrals)
        self._fitted_energies = [None]*len(self.dihedrals)
        scores = np.ones(len(self.dihedrals))
        converged = False
        iteration = 1

        while not converged:
            print("\nIteration %d" % iteration)

            last_scores = scores
            scores = np.zeros(len(self.dihedrals))

            for i, name in enumerate(self._names):
                print('\n == Fitting torsion %s ==\n' % name)

                chisq, self._fitted_energies[i] = self._fitDihedral(i)
                scores[i] = chisq

                rating = 'GOOD'
                if chisq > 10:
                    rating = 'CHECK'
                if chisq > 100:
                    rating = 'BAD'
                print('Chi^2 score : %f : %s' % (chisq, rating))
                sys.stdout.flush()

            if iteration > 1:
                converged = True
                for j in range(len(scores)):
                    # Check convergence
                    try:
                        relerr = (scores[j] - last_scores[j]) / last_scores[j]
                    except:
                        relerr = 0.
                    if np.isnan(relerr):
                        relerr = 0.
                    convstr = "- converged"
                    if np.fabs(relerr) > 1.e-2:
                        convstr = ""
                        converged = False
                    print(" Dihedral %d relative error : %f %s" % (j, relerr, convstr))

            iteration += 1

        print(" Fitting converged at iteration %d" % (iteration - 1))

        if self.result_directory:
            self.plotConformerEnergies()
            for idihed in range(len(self._names)):
                self.plotDihedralEnergies(idihed)

    def plotDihedralEnergies(self, idihed):

        plt.figure()
        plt.title(self._names[idihed])
        plt.xlabel('Dihedral angle, deg')
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        plt.ylabel('Energy, kcal/mol')
        angles = self._angle_values[idihed][:, 0]
        refs = self._reference_energies[idihed] - np.min(self._reference_energies[idihed])
        inits = self._initial_energies[idihed] - np.min(self._initial_energies[idihed])
        fits = self._fitted_energies[idihed] - np.min(self._fitted_energies[idihed])
        plt.plot(angles, refs, 'r-', marker='o', lw=3, label='QM')
        plt.plot(angles, inits, 'g-', marker='o', label='MM initial')
        plt.plot(angles, fits, 'b-', marker='o', label='MM fitted', )
        plt.legend()
        plt.savefig(os.path.join(self.result_directory, self._names[idihed] + '.svg'))
        plt.close()

    def plotConformerEnergies(self):

        qm_energy = np.concatenate(self._reference_energies)[:, None]
        mm_energy = np.concatenate(self._fitted_energies)[:, None]
        qm_energy -= np.min(qm_energy)
        mm_energy -= np.min(mm_energy)

        regr = LinearRegression(fit_intercept=False)
        regr.fit(qm_energy, mm_energy)
        prediction = regr.predict(qm_energy)

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