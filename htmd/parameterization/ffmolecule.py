# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import numpy as np
import scipy.optimize as optimize

from htmd.molecule.molecule import Molecule
from htmd.molecule.util import dihedralAngle
from htmd.molecule.vdw import VDW
from htmd.molecule.vmdparser import guessAnglesAndDihedrals
from htmd.parameterization.detectsoftdihedrals import detectSoftDihedrals
from htmd.parameterization.detectequivalents import detectEquivalents
from htmd.parameterization.fftype import FFTypeMethod, FFType
from htmd.parameterization.ff import RTF, PRM
from htmd.parameterization.ffevaluate import FFEvaluate
from htmd.parameterization.esp import ESP
from htmd.qm import Psi4
from htmd.progress.progress import ProgressBar

class QMFittingSet:

    def __init__(self):
        self.phi = []
        self.qm = []
        self.mm_original = []
        self.mm_zeroed = []
        self.mm_delta = []
        self.mm_fitted = []
        self.coords = []


class FFMolecule(Molecule):
    """
    filename -- a mol2 format input geometry
    rtf, prm -- rtf, prm files
    method  -- if rtf, prm == None, guess atom types according to this method ( of enum FFTypeMethod )
    """

    def __init__(self, filename=None, name=None, rtf=None, prm=None, netcharge=None, method=FFTypeMethod.CGenFF_2b6,
                 qm=None, outdir="./"):

        self.method = method
        self.outdir = outdir

        if not filename.endswith(".mol2"):
            raise ValueError("Input file must be mol2 format")

        super().__init__(filename=filename, name=name)

        self.natoms = self.numAtoms # TODO remove

        # Guess bonds
        if len(self.bonds)==0:
           print("No bonds found. Guessing them")
           self.bonds =  self._guessBonds()

        # Guess angles and dihedrals
        self.angles, self.dihedrals = guessAnglesAndDihedrals(self.bonds)

        # Detect equivalent atoms
        equivalents = detectEquivalents(self)
        self._equivalent_atom_groups = equivalents[0]  # list of groups of equivalent atoms
        self._equivalent_atoms = equivalents[1]  # list of equivalent atoms, indexed by atom
        self._equivalent_group_by_atom = equivalents[2]  # mapping from atom index to equivalent atom group

        # Detect rotable dihedrals
        self._soft_dihedrals = detectSoftDihedrals(self, equivalents)

        # Set total charge
        if netcharge is None:
            self.netcharge = int(round(np.sum(self.charge)))
        else:
            self.netcharge = int(round(netcharge))

        # Canonicalise the atom naming.
        self._rename_mol()

        if rtf and prm:
            # If the user has specified explicit RTF and PRM files go ahead and load those
            self._rtf = RTF(rtf)
            self._prm = PRM(prm)
        elif method == FFTypeMethod.NONE:
            # Don't assign any atom types
            pass
        else:
            # Otherwise make atom types using the specified method
            fftype = FFType(self, method=self.method)
            self._rtf = fftype._rtf
            self._prm = fftype._prm

        self.qm = qm if qm else Psi4()

        self.theory_name = self.qm.theory
        self.basis_name = self.qm.basis
        self.basis_name = re.sub('\+', 'plus', self.basis_name)  # Replace '+' with 'plus'
        self.basis_name = re.sub('\*', 'star', self.basis_name)  # Replace '*' with 'star'
        self.solvent_name = self.qm.solvent.lower()

        if hasattr(self, '_rtf'):
            self.impropers = np.array(self._rtf.impropers)

        # Set atom masses
        if self.masses.size == 0:
            if hasattr(self, '_rtf'):
                self.masses[:] = [self._rtf.mass_by_type[self._rtf.type_by_index[i]] for i in range(self.numAtoms)]
            else:
                self.masses[:] = [VDW.masses[VDW.elements.index(element)] for element in self.element]

        self.report()

    def copy(self):

        # HACK! Circumvent 'qm' coping problem
        qm, self.qm = self.qm, None
        copy = super().copy()
        self.qm = copy.qm = qm

        return copy

    def report(self):
        print("Net Charge: {}".format(self.netcharge))
        print("Equivalent atom groups:")
        for i in self._equivalent_atom_groups:
            for j in i:
                print(" {}".format(self.name[j]), end="")
            print("")

        print("Soft torsions:")
        for i in self._soft_dihedrals:
            for j in i.atoms:
                print(" {}".format(self.name[j]), end="")
            print("")

    def _rename_mol(self):
        """
        This fixes up the atom naming and reside name to be consistent.
        NB this scheme matches what MATCH does.
        Don't change it or the naming will be inconsistent with the RTF.
        """

        sufices = dict()

        print('\nRename atoms:')
        for i in range(len(self.name)):
            name = self.name[i].upper()

            # This fixes the specific case where a name is 3 or 4 characters, as X-TOOL seems to make
            if re.match('^[A-Z]{3,4}$', name):
               name = name[:-2] # Remove the last 2 characters

            # Remove any character that isn't alpha
            name = re.sub('[^A-Z]*', '', name)

            sufices[name] = sufices.get(name, 0) + 1

            name += str(sufices[name])
            print(' %-4s --> %-4s' % (self.name[i], name))

            self.name[i] = name
            self.resname[i] = "MOL"

        print()

    def output_directory_name(self):
        return self.theory_name + "-" + self.basis_name + "-" + self.solvent_name 

    def minimize(self):

        assert self.numFrames == 1

        mindir = os.path.join(self.outdir, "minimize", self.output_directory_name())
        os.makedirs(mindir, exist_ok=True)

        self.qm.molecule = self
        self.qm.esp_points = None
        self.qm.optimize = True
        self.qm.restrained_dihedrals = None
        self.qm.directory = mindir
        results = self.qm.run()
        if results[0].errored:
            raise RuntimeError("QM Optimization failed")

        # Replace coordinates with the minimized set
        self.coords = results[0].coords

    def _fitDihedral_objective(self, params):
        """
        Evaluate the torsion with the input params for each of the phi's poses
        """

        chisq = 0.
        for t in range(self._fitDihedral_results.N):
            energy = 0.
            for s in range(len(self._fitDihedral_results.phis[t])):
                phi = np.deg2rad(self._fitDihedral_results.phis[t][s])
                for j in range(6):
                    n = j + 1
                    phi0 = np.deg2rad(params[6 + j])
                    energy += params[j] * (1. + np.cos(n * phi - phi0))
            energy += params[12]

            diff = self._fitDihedral_results.mm_delta[t] - energy
            chisq += diff**2

        return chisq

    @property
    def centreOfMass(self):
        return np.dot(self.masses, self.coords[:, :, self.frame])/np.sum(self.masses)

    def removeCOM(self):
        """Relocate centre of mass to the origin"""

        for frame in range(self.numFrames):
            self.frame = frame
            self.coords[:, :, frame] -= self.centreOfMass

    def fitCharges(self, fixed=[]):

        # Cereate an ESP directory
        espDir = os.path.join(self.outdir, "esp", self.output_directory_name() )
        os.makedirs(espDir, exist_ok=True)

        # Get ESP points
        point_file = os.path.join(espDir, "00000", "grid.dat")
        if os.path.exists(point_file):
            # Load a point file if one exists from a previous job
            esp_points = np.loadtxt(point_file)
            print("Reusing previously-generated point cloud")
        else:
            # Generate ESP points
            esp_points = ESP.generate_points(self)[0]

        # Run QM simulation
        self.qm.molecule = self
        self.qm.esp_points = esp_points
        self.qm.optimize = False
        self.qm.restrained_dihedrals = None
        self.qm.directory = espDir
        qm_results = self.qm.run()
        if qm_results[0].errored:
            raise RuntimeError("QM Calculation failed")

        # Safeguard QM code from changing coordinates :)
        assert np.all(np.isclose(self.coords, qm_results[0].coords, atol=1e-6))

        # Fit ESP charges
        self.esp = ESP()
        self.esp.molecule = self
        self.esp.qm_results = qm_results
        self.esp.fixed = fixed
        esp_result = self.esp.run()
        esp_charges, esp_loss = esp_result['charges'], esp_result['loss']

        # Update the charges
        self.charge[:] = esp_charges
        self._rtf.updateCharges(esp_charges)

        return esp_loss, qm_results[0].dipole

    def getDipole(self):
        """Calculate the dipole moment (in Debyes) of the molecule"""

        coords = self.coords[:, :, self.frame] - self.centreOfMass

        dipole = np.zeros(4)
        dipole[:3] = np.dot(self.charge, coords)
        dipole[3] = np.linalg.norm(dipole[:3]) # Total dipole moment
        dipole /= 0.20819434 # e * Ang --> Debye (https://en.wikipedia.org/wiki/Debye)

        return dipole

    def getSoftTorsions(self):

        return [dihedral.atoms.copy() for dihedral in self._soft_dihedrals]

    def fitSoftTorsion(self, dihedral, geomopt=True):

        bkp_coords = self.coords.copy()

        phi_to_fit = None
        frozens = []

        for d in self._soft_dihedrals:
            if (d.atoms == dihedral).all():
                phi_to_fit = d
                frozens.append(d.atoms)
            else:
                if not geomopt:
                    frozens.append(d.atoms)

        if not phi_to_fit:
            raise ValueError("specified phi is not a recognised soft dihedral")

        self._makeDihedralUnique(phi_to_fit)

        atoms = phi_to_fit.atoms
        equivs = phi_to_fit.equivalents

        # Number of rotamers for each dihedral to compute
        nrotamer = 36

        # Create a copy of molecule with nrotamer frames
        mol = self.copy()
        for _ in range(nrotamer-1):
            mol.appendFrames(self)
        assert mol.numFrames == nrotamer

        # Set rotamer coordinates
        angles = np.linspace(-np.pi, np.pi, num=nrotamer, endpoint=False)
        for frame, dihedral in enumerate(angles):
            mol.frame = frame
            mol.setDihedral(atoms, dihedral, bonds=mol.bonds)

        dirname = 'dihedral-opt' if geomopt else 'dihedral-single-point'
        dihedral_name = "%s-%s-%s-%s" % tuple(self.name[atoms])
        fitdir = os.path.join(self.outdir, dirname, dihedral_name, self.output_directory_name())
        os.makedirs(fitdir, exist_ok=True)

        self.qm.molecule = mol
        self.qm.esp_points = None
        self.qm.optimize = geomopt
        self.qm.restrained_dihedrals = np.array(frozens)
        self.qm.directory = fitdir
        results = self.qm.run()

        ret = self._makeDihedralFittingSetFromQMResults(atoms, results)

        # Get the initial parameters of the dihedral we are going to fit
        types = tuple([self._rtf.type_by_index[atom] for atom in atoms])
        param = self._prm.dihedralParam(*types)

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
        self._prm.updateDihedral(param)

        ffeval = FFEvaluate(self)

        # Now evaluate the ff without the dihedral being fitted
        for t in range(ret.N):
            mm_zeroed = ffeval.run(ret.coords[t][:, :, 0])['total']
            ret.mm_delta.append(ret.qm[t] - mm_zeroed)
            ret.mm_zeroed.append(mm_zeroed)

        mmin1 = min(ret.mm_zeroed)
        mmin2 = min(ret.mm_delta)
        for t in range(ret.N):
            ret.mm_zeroed[t] = ret.mm_zeroed[t] - mmin1
            ret.mm_delta[t] = ret.mm_delta[t] - mmin2

        self._fitDihedral_results = ret
        self._fitDihedral_phi = param

        # Now measure all of the soft dihedrals phis that are mapped to this dihedral
        ret.phis = []
        for iframe in range(ret.N):
            ret.phis.append([ret.phi[iframe]])
            for atoms in equivs:
                dihedral = dihedralAngle(ret.coords[iframe][atoms, :, 0])
                ret.phis[iframe].append(dihedral)

        best_chisq = self._fitDihedral_objective(best_param)

        bar = ProgressBar(64, description="Fitting")
        for iframe in range(64):

            bounds, start = self._fitDihedral_make_bounds(iframe)

            xopt = optimize.minimize(self._fitDihedral_objective, start, method="L-BFGS-B", bounds=bounds,
                                     options={'disp': False})

            chisq = self._fitDihedral_objective(xopt.x)
            if chisq < best_chisq:
                best_chisq = chisq
                best_param = xopt.x
            bar.progress()
        bar.stop()

        # Update the target dihedral with the optimized parameters
        for iframe in range(6):
            param[iframe].k0 = best_param[0 + iframe]
            param[iframe].phi0 = best_param[6 + iframe]
        self._prm.updateDihedral(param)

        # Finally evaluate the fitted potential
        ffeval = FFEvaluate(self)

        for t in range(ret.N):
            ret.mm_fitted.append(ffeval.run(ret.coords[t][:, :, 0])['total'])
        mmin = min(ret.mm_fitted)

        chisq = 0.
        for t in range(ret.N):
            ret.mm_fitted[t] = ret.mm_fitted[t] - mmin
            delta = ret.mm_fitted[t] - ret.qm[t]
            chisq += delta**2
        ret.chisq = chisq

        # TODO Score it
        self.coords = bkp_coords

        return ret

    @staticmethod
    def _fitDihedral_make_bounds(i):

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

    def _makeDihedralFittingSetFromQMResults(self, atoms, results):
        # Extract the valid QM poses and energies from the QM result set
        # Evaluate the MM on those poses
        ffeval = FFEvaluate(self)

        ret = QMFittingSet()
        ret.name = "%s-%s-%s-%s" % (
            self._rtf.names[atoms[0]], self._rtf.names[atoms[1]], self._rtf.names[atoms[2]], self._rtf.names[atoms[3]])

        qmin = 1.e100
        for q in results:
            if not q.errored:
                if q.energy < qmin:
                    qmin = q.energy

        completed = 0
        for q in results:
            if not q.errored:
                if (q.energy - qmin) < 20.:  # Only fit against QM points < 20 kcal above the minimum
                    mmeval = ffeval.run(q.coords[:, :, 0])
                    angle = dihedralAngle(q.coords[atoms, :, 0])
                    if mmeval["vdw"] < 200:
                        completed += 1
                        ret.qm.append(q.energy - qmin)
                        ret.mm_original.append(mmeval['total'])
                        ret.coords.append(q.coords)
                        ret.phi.append(angle)
                    else:
                        print("Omitting optimised pose for phi=%f (MM VDW too high)" % angle)
                else:
                    print("Omitting optimised QM pose (QM energy too high %f)" % q.energy)

        mmin = min(ret.mm_original)
        # roughly align the qm with the mm
        for q in range(completed):
            ret.mm_original[q] = ret.mm_original[q] - mmin
        ret.N = completed

        if completed < 13:
            raise RuntimeError("Fewer than 13 valid QM points. Not enough to fit!")

        return ret

    def _makeDihedralUnique(self, phi_to_fit):
        #    (number_of_uses, uses) = self._countUsesOfDihedral( phi_to_fit.atoms )
        #    if( number_of_uses > 1 ):
        # Create a new type for (arbitrarily) a middle atom of the dihedral
        # So that the dihedral we are going to modify is unique
        # TODO -- check symmetry
        #    print( "Dihedral term is not unique. Copying type.." ) # Used %d times, by:" % ( number_of_uses ) )
        # print( uses )

        # Duplicate the dihedrals types so this modified term is unique
#        print("Duplicating types..")
        for i in range(4):
            if not ("x" in self._rtf.type_by_index[phi_to_fit.atoms[i]]):
                self._duplicateTypeOfAtom(phi_to_fit.atoms[i])

        number_of_uses, uses = self._countUsesOfDihedral(phi_to_fit.atoms)
        if number_of_uses > 1:
            print(phi_to_fit.atoms)
            print(number_of_uses)
            print(uses)
            raise ValueError("Dihedral term still not unique after duplication")

    def _countUsesOfDihedral(self, indices):
        """
        Return the number of uses of the dihedral specified by the types of the 4 atom indices
        """

        types = [self._rtf.type_by_index[index] for index in indices]

        all_uses = []
        for dihedral_indices in self.dihedrals:
            dihedral_types = [self._rtf.type_by_index[index] for index in dihedral_indices]
            if types == dihedral_types or types == dihedral_types[::-1]:
                all_uses.append(dihedral_indices)

        # Now for each of the uses, remove any which are equivalent
        unique_uses = [indices]
        groups = [self._equivalent_group_by_atom[index] for index in indices]
        for dihedral in all_uses:
            dihedral_groups = [self._equivalent_group_by_atom[index] for index in dihedral]
            if groups != dihedral_groups and groups != dihedral_groups[::-1]:
                unique_uses.append(dihedral)

        return len(unique_uses), unique_uses

    def _duplicateTypeOfAtom(self, atom_index):
        # This duplicates the type of the specified atom
        # First get the type
        type_ = self._rtf.type_by_index[atom_index]

        # perhaps the type is already a duplicate? if so
        # remove the duplicated suffix
        type_ = re.sub('x\d+$', '', type_)

        # make the new type name
        i = 0
        while ('%sx%d' % (type_, i)) in self._rtf.types:
            i += 1
        newtype = '%sx%d' % (type_, i)
        print("Creating new type %s from %s for atom %s" % (newtype, type_, self._rtf.names[atom_index]))

        # TODO: move to RTF class
        # duplicate the type in the fields RTF --
        self._rtf.type_by_index[atom_index] = newtype
        self._rtf.mass_by_type[newtype] = self._rtf.mass_by_type[type_]
        self._rtf.types.append(newtype)
        self._rtf.type_by_name[self._rtf.names[atom_index]] = newtype
        self._rtf.type_by_index[atom_index] = newtype
        self._rtf.typeindex_by_type[newtype] = self._rtf.typeindex_by_type[type_] + 1000
        self._rtf.element_by_type[newtype] = self._rtf.element_by_type[type_]

        # Now also reset the type of  any atoms that share equivalency
        for index in self._equivalent_atoms[atom_index]:
            if atom_index != index:
                if "x" in self._rtf.type_by_index[index]:
                    raise RuntimeError(
                        "Equivalent atom already has a duplicated type: {} {}".format(index, self._rtf.type_by_index[index]))
                self._rtf.type_by_index[index] = newtype
                self._rtf.type_by_name[self._rtf.names[index]] = newtype

        # the PRM parameters will be automatically duplicated by forcing an ff evaluation
        FFEvaluate(self).run(self.coords[:, :, 0])

    def plotConformerEnergies(self, fits):

        import matplotlib.pyplot as plt
        from sklearn.linear_model import LinearRegression

        qm_energy = np.concatenate([fit.qm for fit in fits])[:, None]
        mm_energy = np.concatenate([fit.mm_fitted for fit in fits])[:, None]
        qm_energy -= np.min(qm_energy)
        mm_energy -= np.min(mm_energy)

        regr = LinearRegression(fit_intercept=False)
        regr.fit(qm_energy, mm_energy)
        prediction = regr.predict(qm_energy)
        rms = np.sqrt(np.mean((prediction - mm_energy)**2))
        score = regr.score(qm_energy, mm_energy)

        plotdir = os.path.join(self.outdir, 'parameters', self.method.name, self.output_directory_name(), 'plots')
        os.makedirs(plotdir, exist_ok=True)

        plt.figure()
        plt.title('Conformer Energies MM vs QM')
        plt.xlabel('QM energy, kcal/mol')
        plt.ylabel('MM energy, kcal/mol')
        plt.plot(qm_energy, mm_energy, 'ko')
        plt.plot(qm_energy, prediction, 'r-', lw=2)
        plt.savefig(os.path.join(plotdir, 'conformer-energies.svg'))
        plt.close()

        return rms, score, regr.coef_

    def plotTorsionFit(self, fit):

        import matplotlib.pyplot as plt

        plotdir = os.path.join(self.outdir, 'parameters', self.method.name, self.output_directory_name(), 'plots')
        os.makedirs(plotdir, exist_ok=True)

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
        plt.savefig(os.path.join(plotdir, fit.name + '.svg'))
        plt.close()

    def write(self, filename, sel=None, type=None, typemap=None):

        if hasattr(self, "_rtf"):  # Update base Molecule's attributes so write() works correctly
            self.segid[:] = 'L'
            self.charge[:] = [self._rtf.charge_by_name[name] for name in self.name]
            self.atomtype[:] = [self._rtf.type_by_name[name] for name in self.name]

        atomtype_backup = np.copy(self.atomtype)
        if typemap:
            self.atomtype[:] = [typemap[atomtype] for atomtype in self.atomtype]

        super().write(filename, sel=sel, type=type)

        self.atomtype[:] = atomtype_backup
