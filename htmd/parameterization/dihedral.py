# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import sys
import numpy as np

from htmd.molecule.util import dihedralAngle
from htmd.parameterization.ffevaluate import FFEvaluate


class QMFittingSet:

    def __init__(self):
        self.name = None
        self.phi = []
        self.qm = []
        self.mm_original = []
        self.mm_zeroed = []
        self.mm_delta = []
        self.mm_fitted = []
        self.coords = []


class DihedralFitting:
    """
    Dihedral parameter fitting from QM data
    """

    def __init__(self):

        self.molecule = None
        self.dihedrals = []
        self.qm_resutls = []

    @staticmethod
    def _makeFittingData(molecule, dihedral, qm_results):
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

        data = QMFittingSet()
        data.name = '%s-%s-%s-%s' % tuple(molecule.name[dihedral])
        data.phi = np.array([dihedralAngle(result.coords[dihedral, :, 0]) for result in results])

        data.qm = np.array([result.energy for result in results])
        data.qm -= np.min(data.qm)

        ff = FFEvaluate(molecule)
        data.mm_original = np.array([ff.run(result.coords[:, :, 0])['total'] for result in results])
        data.mm_original -= np.min(data.mm_original)

        data.coords = [result.coords for result in results]

        return data

    @staticmethod
    def _make_bounds(i):

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
        for phis, mm_delta in zip(data.phis, data.mm_delta):
            energy = 0.
            for phi in phis:
                phi = np.deg2rad(phi)
                for j in range(6):
                    n = j + 1
                    phi0 = np.deg2rad(params[6 + j])
                    energy += params[j] * (1. + np.cos(n * phi - phi0))
            energy += params[12]
            chisq += (mm_delta - energy)**2

        return chisq

    @staticmethod
    def fitTorsionParameters(molecule, dihedrals, geomopt):

        from htmd.parameterization.ffmolecule import FFMolecule

        if not isinstance(molecule, FFMolecule):
            raise TypeError('molecule')

        scores = np.ones(len(dihedrals))
        converged = False
        iteration = 1
        ref_mm = dict()

        while not converged:
            rets = []

            print("\nIteration %d" % iteration)

            last_scores = scores
            scores = np.zeros(len(dihedrals))

            for idx, dihedral in enumerate(dihedrals):
                name = '%s-%s-%s-%s' % tuple(molecule.name[dihedral])

                print("\n == Fitting torsion {} ==\n".format(name))
                try:
                    ret = molecule.fitSoftTorsion(dihedral, geomopt=geomopt)
                    rets.append(ret)

                    if iteration == 1:
                        ref_mm[name] = ret

                    rating = "GOOD"
                    if ret.chisq > 10:
                        rating = "CHECK"
                    if ret.chisq > 100:
                        rating = "BAD"
                    print("Torsion %s Chi^2 score : %f : %s" % (name, ret.chisq, rating))

                    sys.stdout.flush()
                    scores[idx] = ret.chisq

                    # Always use the mm_orig from first iteration (unmodified)
                    ret.mm_original = ref_mm[name].mm_original
                    molecule.plotTorsionFit(ret)

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

        if len(rets):
            fit = molecule.plotConformerEnergies(rets)
            print("\n Fit of conformer energies: RMS %f Variance %f" % (fit[0], fit[1]))


if __name__ == '__main__':
    pass