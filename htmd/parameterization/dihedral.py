# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from htmd.molecule.util import dihedralAngle
from htmd.parameterization.ffevaluate import FFEvaluate
from htmd.progress.progress import ProgressBar


class DihedralFittingData:

    def __init__(self):

        self.name = None
        self.atom_indices = None
        self.equivalents = None

        self.qm_energies = []
        self.mm_initial_energies = []
        self.mm_delta_energies = []

        self.coords = []
        self.anlge_values = []

        self.mm_fitted_energies = []


class DihedralFitting:
    """
    Dihedral parameter fitting from QM data
    """

    def __init__(self):

        self.molecule = None
        self.dihedrals = []
        self.qm_resutls = []
        self.result_directory = None

        self._data = []

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

    def _makeFittingData(self, dihedral, qm_results):
        """
        Extract the valid QM poses and energies from the QM result set.
        Evaluate the MM on those poses.
        """

        # Remove failed QM results
        # TODO print failed QM jobs
        results = [result for result in qm_results if not result.errored]

        # Remove QM results with too high QM energies (>20 kcal/mol above the minimum)
        # TODO print removed QM jobs
        qm_min = np.min([result.energy for result in results])
        results = [result for result in results if (result.energy - qm_min) < 20]

        if len(results) < 13:
            raise RuntimeError("Fewer than 13 valid QM points. Not enough to fit!")

        fittingData = DihedralFittingData()
        fittingData.name = '%s-%s-%s-%s' % tuple(self.molecule.name[dihedral])
        fittingData.atom_indices = dihedral
        for rotableDihedral in self.molecule._soft_dihedrals:
                if np.all(rotableDihedral.atoms == dihedral):
                    fittingData.equivalents = rotableDihedral.equivalents
                    break
        else:
            raise ValueError('%s are not recognized as the rotable dihedrals\n' % fittingData.name)


        fittingData.coords = [result.coords for result in results]

        # Calculate angle values for all equivaltent dihedrals
        equivalentDihedrals = [fittingData.atom_indices] + fittingData.equivalents
        fittingData.anlge_values = []
        for coords in fittingData.coords:
            angles = [dihedralAngle(coords[indices, :, 0]) for indices in equivalentDihedrals]
            fittingData.anlge_values.append(angles)
        fittingData.anlge_values = np.array(fittingData.anlge_values)

        fittingData.qm_energies = np.array([result.energy for result in results])
        fittingData.qm_energies -= np.min(fittingData.qm_energies)

        ff = FFEvaluate(self.molecule)
        fittingData.mm_initial_energies = np.array([ff.run(result.coords[:, :, 0])['total'] for result in results])
        fittingData.mm_initial_energies -= np.min(fittingData.mm_initial_energies)

        return fittingData

    @staticmethod
    def _makeBounds(i):

        start = np.zeros(13)
        bounds = []

        for j in range(6):
            bounds.append((-20., 20.))

        for j in range(6):
            if i & (2 ** j):
                bounds.append((180., 180.))
                start[6 + j] = 180.
            else:
                bounds.append((0., 0.))

        bounds.append((-10., 10.))

        return bounds, start

    @staticmethod
    def _objective(x, data):
        """
        Evaluate the torsion with the input params for each of the phi's poses
        """

        k0 = x[0:6]
        phi0 = np.deg2rad(x[6:12])
        offset = x[12]

        n = np.arange(6) + 1
        phis = np.deg2rad(data.anlge_values)[:, :, None]  # rotamers x equivalent dihedral values

        energies = np.sum(k0 * (1. + np.cos(n * phis - phi0)), axis=(1, 2)) + offset
        chisq = np.sum((energies - data.mm_delta_energies)**2)

        return chisq

    def _fitDihedral(self, fittingData):

        # Get the initial parameters of the dihedral we are going to fit
        types = tuple([self.molecule._rtf.type_by_index[index] for index in fittingData.atom_indices])
        param = self.molecule._prm.dihedralParam(*types)

        # Save these parameters as the best fit (fit to beat)
        best_param = np.zeros(13)
        for i, term in enumerate(param):
            best_param[i] = term.k0
            best_param[i + 6] = term.phi0
        best_param[12] = 0.

        # Evalaute the mm potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # The QM potential and the this modified mm potential
        for term in param:
            term.k0 = term.phi0 = 0.
        self.molecule._prm.updateDihedral(param)

        # Now evaluate the ff without the dihedral being fitted
        ffeval = FFEvaluate(self.molecule)
        mm_zeroed = np.array([ffeval.run(coords[:, :, 0])['total'] for coords in fittingData.coords])
        fittingData.mm_delta_energies = fittingData.qm_energies - mm_zeroed
        fittingData.mm_delta_energies -= np.min(fittingData.mm_delta_energies)

        # Optimize parameters
        best_chisq = DihedralFitting._objective(best_param, fittingData)
        bar = ProgressBar(64, description="Fitting")
        for i in range(64):
            bar.progress()
            bounds, start = DihedralFitting._makeBounds(i)
            xopt = optimize.minimize(self._objective, start, args=fittingData, method="L-BFGS-B", bounds=bounds,
                                     options={'disp': False})
            chisq = DihedralFitting._objective(xopt.x, fittingData)
            if chisq < best_chisq:
                best_chisq = chisq
                best_param = xopt.x
        bar.stop()

        # Update the target dihedral with the optimized parameters
        for i in range(6):
            param[i].k0 = best_param[i]
            param[i].phi0 = best_param[i + 6]
        self.molecule._prm.updateDihedral(param)

        # Finally evaluate the fitted potential
        ffeval = FFEvaluate(self.molecule)
        fittingData.mm_fitted_energies = np.array([ffeval.run(coords[:, :, 0])['total'] for coords in fittingData.coords])
        fittingData.mm_fitted_energies -= np.min(fittingData.mm_fitted_energies)
        chisq = np.sum((fittingData.mm_fitted_energies - fittingData.qm_energies)**2)

        return chisq

    def plotDihedralEnergies(self, fittingData):

        plt.figure()
        plt.title(fittingData.name)
        plt.xlabel('Dihedral angle, deg')
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        plt.ylabel('Energy, kcal/mol')
        plt.plot(fittingData.anlge_values[:, 0], fittingData.qm_energies, 'r-', marker='o', lw=3, label='QM')
        plt.plot(fittingData.anlge_values[:, 0], fittingData.mm_initial_energies, 'g-', marker='o', label='MM original')
        plt.plot(fittingData.anlge_values[:, 0], fittingData.mm_fitted_energies, 'b-', marker='o', label='MM fitted', )
        plt.legend()
        plt.savefig(os.path.join(self.result_directory, fittingData.name + '.svg'))
        plt.close()

    def plotConformerEnergies(self, fits):

        qm_energy = np.concatenate([fit.qm_energies for fit in fits])[:, None]
        mm_energy = np.concatenate([fit.mm_fitted_energies for fit in fits])[:, None]
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

    def run(self):

        if self.result_directory:
            os.makedirs(self.result_directory, exist_ok=True)

        assert len(self.dihedrals) == len(self.qm_resutls)

        for dihedral in self.dihedrals:
            self._makeDihedralUnique(dihedral)

        self._data = [self._makeFittingData(dihedral, qm_results)
                      for dihedral, qm_results in zip(self.dihedrals, self.qm_resutls)]

        scores = np.ones(len(self.dihedrals))
        converged = False
        iteration = 1

        while not converged:
            print("\nIteration %d" % iteration)

            last_scores = scores
            scores = np.zeros(len(self.dihedrals))

            for i, fittingData in enumerate(self._data):
                print('\n == Fitting torsion %s ==\n' % fittingData.name)

                chisq = self._fitDihedral(fittingData)
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
            self.plotConformerEnergies(self._data)
            for fittingData in self._data:
                self.plotDihedralEnergies(fittingData)


if __name__ == '__main__':
    pass