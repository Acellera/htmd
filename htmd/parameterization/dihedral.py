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
        self.phi = []
        self.qm = []
        self.mm_original = []
        self.mm_zeroed = []
        self.mm_delta = []
        self.mm_fitted = []
        self.coords = []

        self.phis = []  # TODO factor out
        self.chisq = None


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

        data = DihedralFittingData()
        data.name = '%s-%s-%s-%s' % tuple(self.molecule.name[dihedral])
        data.phi = np.array([dihedralAngle(result.coords[dihedral, :, 0]) for result in results])

        data.qm = np.array([result.energy for result in results])
        data.qm -= np.min(data.qm)

        ff = FFEvaluate(self.molecule)
        data.mm_original = np.array([ff.run(result.coords[:, :, 0])['total'] for result in results])
        data.mm_original -= np.min(data.mm_original)

        data.coords = [result.coords for result in results]

        return data

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
    def _objective(params, data):
        """
        Evaluate the torsion with the input params for each of the phi's poses
        """

        chisq = 0.
        for phis, mm_delta in zip(np.deg2rad(data.phis), data.mm_delta):
            energy = 0.
            for phi in phis:
                for j in range(6):
                    n = j + 1
                    phi0 = np.deg2rad(params[6 + j])
                    energy += params[j] * (1. + np.cos(n * phi - phi0))
            energy += params[12]
            chisq += (mm_delta - energy)**2

        return chisq

    def _fitDihedral(self, rotableDihedral, fittingData):

        atoms = rotableDihedral.atoms
        equivs = rotableDihedral.equivalents

        # Get the initial parameters of the dihedral we are going to fit
        types = tuple([self.molecule._rtf.type_by_index[atom] for atom in atoms])
        param = self.molecule._prm.dihedralParam(*types)

        # Save these parameters as the best fit (fit to beat)
        best_param = np.zeros(13)
        for t in range(6):
            best_param[t] = param[t].k0
            best_param[t + 6] = param[t].phi0
        best_param[12] = 0.

        # Evalaute the mm potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # The QM potential and the this modified mm potential
        for t in param:
            t.k0 = t.phi0 = 0.
        self.molecule._prm.updateDihedral(param)

        # Now evaluate the ff without the dihedral being fitted
        ffeval = FFEvaluate(self.molecule)
        fittingData.mm_zeroed = np.array([ffeval.run(coords[:, :, 0])['total'] for coords in fittingData.coords])
        fittingData.mm_delta = fittingData.qm - fittingData.mm_zeroed
        fittingData.mm_zeroed -= np.min(fittingData.mm_zeroed)
        fittingData.mm_delta -= np.min(fittingData.mm_delta)

        # Now measure all of the soft dihedrals phis that are mapped to this dihedral
        # TODO get rid of this insanity
        fittingData.phis = []
        for iframe in range(len(fittingData.phi)):
            fittingData.phis.append([fittingData.phi[iframe]])
            for atoms in equivs:
                fittingData.phis[iframe].append(dihedralAngle(fittingData.coords[iframe][atoms, :, 0]))

        # Optimize parameters
        best_chisq = DihedralFitting._objective(best_param, fittingData)
        bar = ProgressBar(64, description="Fitting")
        for i in range(64):
            bar.progress()
            bounds, start = DihedralFitting._makeBounds(i)
            xopt = optimize.minimize(DihedralFitting._objective, start, args=(fittingData,),
                                     method="L-BFGS-B", bounds=bounds, options={'disp': False})
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
        fittingData.mm_fitted = np.array([ffeval.run(coords[:, :, 0])['total'] for coords in fittingData.coords])
        fittingData.mm_fitted -= np.min(fittingData.mm_fitted)
        fittingData.chisq = np.sum((fittingData.mm_fitted - fittingData.qm)**2)

        return fittingData

    def plotTorsionFit(self, fit):

        plt.figure()
        plt.title(fit.name)
        plt.xlabel('Dihedral angle, deg')
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        plt.ylabel('Energy, kcal/mol')
        plt.plot(fit.phi, fit.qm, 'r-', marker='o', lw=3, label='QM')
        plt.plot(fit.phi, fit.mm_original, 'g-', marker='o', label='MM original')
        plt.plot(fit.phi, fit.mm_fitted, 'b-', marker='o', label='MM fitted',)
        plt.legend()
        plt.savefig(os.path.join(self.result_directory, fit.name + '.svg'))
        plt.close()

    def plotConformerEnergies(self, fits):

        qm_energy = np.concatenate([fit.qm for fit in fits])[:, None]
        mm_energy = np.concatenate([fit.mm_fitted for fit in fits])[:, None]
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

        rotableDihedrals = []
        for dihedral in self.dihedrals:
            for dihed in self.molecule._soft_dihedrals:
                if np.all(dihed.atoms == dihedral):
                    rotableDihedrals.append(dihed)
        if len(rotableDihedrals) != len(self.dihedrals):
            raise ValueError('Some dihedrals are not recognized as the rotable dihedrals\n')

        for dihedral in self.dihedrals:
            self._makeDihedralUnique(dihedral)

        self._data = [self._makeFittingData(dihedral, qm_results)
                      for dihedral, qm_results in zip(self.dihedrals, self.qm_resutls)]

        scores = np.ones(len(self.dihedrals))
        converged = False
        iteration = 1

        while not converged:
            rets = []

            print("\nIteration %d" % iteration)

            last_scores = scores
            scores = np.zeros(len(self.dihedrals))

            for idx, rotableDihedral in enumerate(rotableDihedrals):
                name = '%s-%s-%s-%s' % tuple(self.molecule.name[rotableDihedral.atoms])
                print('\n == Fitting torsion %s ==\n' % name)
                try:
                    ret = self._fitDihedral(rotableDihedral, self._data[idx])
                    scores[idx] = ret.chisq
                    rets.append(ret)

                    rating = 'GOOD'
                    if ret.chisq > 10:
                        rating = 'CHECK'
                    if ret.chisq > 100:
                        rating = 'BAD'
                    print('Chi^2 score : %f : %s' % (ret.chisq, rating))
                    sys.stdout.flush()

                    if self.result_directory:
                        self.plotTorsionFit(ret)

                except Exception as e:
                    print("Error in fitting")
                    print(str(e))
                    raise

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
            self.plotConformerEnergies(rets)


if __name__ == '__main__':
    pass