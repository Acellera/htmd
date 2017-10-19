# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import numpy as np
from scipy import constants as const

from htmd.molecule.molecule import Molecule
from htmd.molecule import vdw
from htmd.molecule.util import dihedralAngle, guessAnglesAndDihedrals
from htmd.parameterization.detectsoftdihedrals import detectSoftDihedrals
from htmd.parameterization.detectequivalents import detectEquivalents
from htmd.parameterization.fftype import FFTypeMethod, FFType
from htmd.parameterization.ff import RTF, PRM
from htmd.parameterization.ffevaluate import FFEvaluate
from htmd.parameterization.esp import ESP
from htmd.parameterization.dihedral import DihedralFitting
from htmd.qm import Psi4, FakeQM
from copy import deepcopy


class FFMolecule(Molecule):
    """
    filename -- a mol2 format input geometry
    rtf, prm -- rtf, prm files
    method  -- if rtf, prm == None, guess atom types according to this method ( of enum FFTypeMethod )
    """

    def __init__(self, filename=None, name=None, rtf=None, prm=None, netcharge=None, method=FFTypeMethod.CGenFF_2b6,
                 qm=None, outdir="./", mol=None):

        self.method = method
        self.outdir = outdir

        if filename is not None and not filename.endswith(".mol2"):
            raise ValueError("Input file must be mol2 format")

        if mol is None:
            super().__init__(filename=filename, name=name)
        else:
            for v in mol.__dict__:
                self.__dict__[v] = deepcopy(mol.__dict__[v])

        self.natoms = self.numAtoms # TODO remove

        # Guess bonds
        if len(self.bonds)==0:
           print("No bonds found. Guessing them")
           self.bonds =  self._guessBonds()

        # Guess angles and dihedrals
        self.angles, self.dihedrals = guessAnglesAndDihedrals(self.bonds, cyclicdih=True)

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

        if hasattr(self, '_rtf'):
            self.impropers = np.array(self._rtf.impropers)

        # Set atom masses
        if self.masses.size == 0:
            if hasattr(self, '_rtf'):
                self.masses[:] = [self._rtf.mass_by_type[self._rtf.type_by_index[i]] for i in range(self.numAtoms)]
            else:
                self.masses[:] = [vdw.massByElement(element) for element in self.element]

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

        basis = self.qm.basis
        basis = re.sub('\+', 'plus', basis)  # Replace '+' with 'plus'
        basis = re.sub('\*', 'star', basis)  # Replace '*' with 'star'

        name = self.qm.theory + '-' + basis + '-' + self.qm.solvent

        return name

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
        dipole *= 1e11*const.elementary_charge*const.speed_of_light # e * Ang --> Debye (https://en.wikipedia.org/wiki/Debye)

        return dipole

    def getSoftTorsions(self):

        return [dihedral.atoms.copy() for dihedral in self._soft_dihedrals]

    def fitSoftTorsions(self, dihedrals, geomopt=True):
        """
        Dihedrals passed as 4 atom indices
        """

        # Create molecules with rotamers
        molecules = []
        for dihedral in dihedrals:

            nrotamers = 36  # Number of rotamers for each dihedral to compute

            # Create a copy of molecule with "nrotamers" frames
            mol = self.copy()
            while mol.numFrames < nrotamers:
                mol.appendFrames(self)
            assert mol.numFrames == nrotamers

            # Set rotamer coordinates
            angles = np.linspace(-np.pi, np.pi, num=nrotamers, endpoint=False)
            for frame, angle in enumerate(angles):
                mol.frame = frame
                mol.setDihedral(dihedral, angle, bonds=mol.bonds)

            molecules.append(mol)

        # Run QM calculation of the rotamers
        dirname = 'dihedral-opt' if geomopt else 'dihedral-single-point'
        qm_results = []
        for dihedral, molecule in zip(dihedrals, molecules):
            name = "%s-%s-%s-%s" % tuple(self.name[dihedral])
            fitdir = os.path.join(self.outdir, dirname, name, self.output_directory_name())
            os.makedirs(fitdir, exist_ok=True)

            self.qm.molecule = molecule
            self.qm.esp_points = None
            self.qm.optimize = geomopt
            self.qm.restrained_dihedrals = np.array([dihedral])
            self.qm.directory = fitdir
            qm_results.append(self.qm.run())  # TODO submit all jobs at once

        # Fit the dihedral parameters
        df = DihedralFitting()
        df.molecule = self
        df.dihedrals = dihedrals
        df.qm_results = qm_results
        df.result_directory = os.path.join(self.outdir, 'parameters', self.method.name,
                                           self.output_directory_name(), 'plots')

        # In case of FakeQM, the initial parameters are set to zeros.
        # It prevents DihedralFitting class from cheating :D
        if isinstance(self.qm, FakeQM):
            df.zeroed_parameters = True

        df.run()
        # TODO explicit parameter update

        print(' RMSD: %s kcal/mol\n' % df.loss)

    def duplicateTypeOfAtom(self, atom_index):
        """Duplicate the type of the specified atom"""

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
                        "Equivalent atom already has a duplicated type: {} {}".format(index,
                                                                                      self._rtf.type_by_index[index]))
                self._rtf.type_by_index[index] = newtype
                self._rtf.type_by_name[self._rtf.names[index]] = newtype

        # the PRM parameters will be automatically duplicated by forcing an ff evaluation
        FFEvaluate(self).run(self.coords[:, :, 0])

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
