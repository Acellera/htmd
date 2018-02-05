# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import logging
from copy import deepcopy
import numpy as np
from scipy import constants as const
import periodictable

from htmd.molecule.molecule import Molecule
from htmd.molecule import vdw
from htmd.molecule.util import guessAnglesAndDihedrals
from htmd.parameterization.detectsoftdihedrals import detectSoftDihedrals
from htmd.parameterization.detectequivalents import detectEquivalents
from htmd.parameterization.fftype import FFTypeMethod, FFType
from htmd.parameterization.ff import RTF, PRM
from htmd.parameterization.ffevaluate import FFEvaluate
from htmd.parameterization.esp import ESP
from htmd.parameterization.dihedral import DihedralFitting
from htmd.qm import Psi4, FakeQM2

logger = logging.getLogger(__name__)


class FFMolecule(Molecule):
    """
    filename -- a mol2 format input geometry
    rtf, prm -- rtf, prm files
    method  -- if rtf, prm == None, guess atom types according to this method ( of enum FFTypeMethod )
    """

    _ATOM_TYPE_REG_EX = re.compile('^\S+x\d+$')

    def __init__(self, filename=None, name=None, rtf=None, prm=None, netcharge=None, method=FFTypeMethod.CGenFF_2b6,
                 qm=None, outdir="./", mol=None, acCharges=None):

        if filename is not None and not filename.endswith('.mol2'):
            raise ValueError('Input file must be mol2 format')

        if mol is None:
            super().__init__(filename=filename, name=name)
        else:
            for v in mol.__dict__:
                self.__dict__[v] = deepcopy(mol.__dict__[v])

        # Guess bonds
        if len(self.bonds) == 0:
            logger.warning('No bonds found! Guessing them...')
            self.bonds = self._guessBonds()

        # Guess angles and dihedrals
        self.angles, self.dihedrals = guessAnglesAndDihedrals(self.bonds, cyclicdih=True)

        # Detect equivalent atoms
        equivalents = detectEquivalents(self)
        self._equivalent_atom_groups = equivalents[0]  # List of groups of equivalent atoms
        self._equivalent_atoms = equivalents[1]  # List of equivalent atoms, indexed by atom
        self._equivalent_group_by_atom = equivalents[2]  # Mapping from atom index to equivalent atom group

        # Detect rotatable dihedrals
        self._rotatable_dihedrals = detectSoftDihedrals(self, equivalents)

        # Set total charge
        if netcharge is None:
            self.netcharge = int(round(np.sum(self.charge)))
        else:
            self.netcharge = int(round(netcharge))

        # Canonicalise the names
        self._rename()

        # Assign atom types, charges, and initial parameters
        self.method = method
        if rtf and prm:
            # If the user has specified explicit RTF and PRM files go ahead and load those
            self._rtf = RTF(rtf)
            self._prm = PRM(prm)
            logger.info('Reading FF parameters from %s and %s' % (rtf, prm))
        elif method == FFTypeMethod.NONE:
            pass  # Don't assign any atom types
        else:
            # Otherwise make atom types using the specified method
            fftype = FFType(self, method=self.method, acCharges=acCharges)
            logger.info('Assigned atom types with %s' % self.method.name)
            self._rtf = fftype._rtf
            self._prm = fftype._prm

        if hasattr(self, '_rtf'):
            self.atomtype[:] = [self._rtf.type_by_name[name] for name in self.name]
            self.charge[:] = [self._rtf.charge_by_name[name] for name in self.name]
            self.impropers = np.array(self._rtf.impropers)

            # Check if atom type names are compatible
            for type_ in self._rtf.types:
                if re.match(FFMolecule._ATOM_TYPE_REG_EX, type_):
                    raise ValueError('Atom type %s is incompatable. It cannot finish with "x" + number!' % type_)

        # Set atom masses
        # TODO: maybe move to molecule
        if self.masses.size == 0:
            if hasattr(self, '_rtf'):
                self.masses[:] = [self._rtf.mass_by_type[self._rtf.type_by_index[i]] for i in range(self.numAtoms)]
            else:
                self.masses[:] = [vdw.massByElement(element) for element in self.element]

        self.qm = qm if qm else Psi4()
        self.outdir = outdir

    def copy(self):

        # HACK! Circumvent 'qm' coping problem
        qm, self.qm = self.qm, None
        copy = super().copy()
        self.qm = copy.qm = qm

        return copy

    def printReport(self):

        print('\n == Molecule report ==\n')

        print('Total number of atoms: %d' % self.numAtoms)
        print('Total charge: %d' % self.netcharge)

        print('Equivalent atom groups:')
        for atom_group in self._equivalent_atom_groups:
            print('  ' + ', '.join(self.name[atom_group]))

        print('Rotatable dihedral angles:')
        for dihedral in self._rotatable_dihedrals:
            print('  ' + '-'.join(self.name[dihedral.atoms]))
            if dihedral.equivalents:
                print('    Equivalents:')
            for equivalent_dihedral in dihedral.equivalents:
                print('      ' + '-'.join(self.name[equivalent_dihedral]))

    @staticmethod
    def guessElementFromName(name):
        '''
        Guess element from an atom name

        >>> from htmd.parameterization.ffmolecule import FFMolecule
        >>> FFMolecule.guessElementFromName('C')
        'C'
        >>> FFMolecule.guessElementFromName('C1')
        'C'
        >>> FFMolecule.guessElementFromName('C42')
        'C'
        >>> FFMolecule.guessElementFromName('C7S')
        'C'
        >>> FFMolecule.guessElementFromName('HN1')
        'H'
        >>> FFMolecule.guessElementFromName('CL')
        'Cl'
        >>> FFMolecule.guessElementFromName('CA1')
        'Ca'
        '''

        symbol = name.capitalize()

        while symbol:
            try:
                element = periodictable.elements.symbol(symbol)
            except ValueError:
                symbol = symbol[:-1]
            else:
                return element.symbol

        raise ValueError('Cannot guess element from atom name: {}'.format(name))


    def _rename(self):
        """
        This fixes up the atom naming and reside name to be consistent.
        NB this scheme matches what MATCH does.
        Don't change it or the naming will be inconsistent with the RTF.
        """

        self.segid[:] = 'L'
        logger.info('Rename segment to %s' % self.segid[0])
        self.resname[:] = 'MOL'
        logger.info('Rename residue to %s' % self.resname[0])

        sufices = {}
        for i in range(self.numAtoms):
            name = self.guessElementFromName(self.name[i]).upper()

            sufices[name] = sufices.get(name, 0) + 1
            name += str(sufices[name])

            logger.info('Rename atom %d: %-4s --> %-4s' % (i, self.name[i], name))
            self.name[i] = name

    def qm_method_name(self):

        basis = self.qm.basis
        basis = re.sub('\+', 'plus', basis)  # Replace '+' with 'plus'
        basis = re.sub('\*', 'star', basis)  # Replace '*' with 'star'

        name = self.qm.theory + '-' + basis + '-' + self.qm.solvent

        return name

    def minimize(self):

        assert self.numFrames == 1

        mindir = os.path.join(self.outdir, "minimize", self.qm_method_name())
        os.makedirs(mindir, exist_ok=True)

        self.qm.molecule = self
        self.qm.esp_points = None
        self.qm.optimize = True
        self.qm.restrained_dihedrals = None
        self.qm.directory = mindir
        results = self.qm.run()
        if results[0].errored:
            raise RuntimeError('\nQM minimization failed! Check logs at %s\n' % mindir)

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
        espDir = os.path.join(self.outdir, "esp", self.qm_method_name())
        os.makedirs(espDir, exist_ok=True)

        # Get ESP points
        point_file = os.path.join(espDir, "00000", "grid.dat")
        if os.path.exists(point_file):
            # Load a point file if one exists from a previous job
            esp_points = np.loadtxt(point_file)
            logger.info('Reusing ESP grid from %s' % point_file)
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
            raise RuntimeError('\nQM calculation failed! Check logs at %s\n' % espDir)

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
        for name, charge in zip(self.name, self.charge):
            logger.info('Set charge %4s: %6.3f' % (name, charge))

        return esp_loss, qm_results[0].dipole

    def getDipole(self):
        """Calculate the dipole moment (in Debyes) of the molecule"""

        coords = self.coords[:, :, self.frame] - self.centreOfMass

        dipole = np.zeros(4)
        dipole[:3] = np.dot(self.charge, coords)
        dipole[3] = np.linalg.norm(dipole[:3]) # Total dipole moment
        dipole *= 1e11*const.elementary_charge*const.speed_of_light # e * Ang --> Debye (https://en.wikipedia.org/wiki/Debye)

        return dipole

    def getRotatableDihedrals(self):

        return [dihedral.atoms.copy() for dihedral in self._rotatable_dihedrals]

    def fitDihedrals(self, dihedrals, geomopt=True):
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

        # Create directories for QM data
        directories = []
        dihedral_directory = 'dihedral-opt' if geomopt else 'dihedral-single-point'
        for dihedral in dihedrals:
            dihedral_name = '-'.join(self.name[dihedral])
            directory = os.path.join(self.outdir, dihedral_directory, dihedral_name, self.qm_method_name())
            os.makedirs(directory, exist_ok=True)
            directories.append(directory)

        # Setup and submit QM calculations
        for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
            self.qm.molecule = molecule
            self.qm.esp_points = None
            self.qm.optimize = geomopt
            self.qm.restrained_dihedrals = np.array([dihedral])
            self.qm.directory = directory
            self.qm.setup()
            self.qm.submit()

        # Wait and retrieve QM calculation data
        qm_results = []
        for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
            self.qm.molecule = molecule
            self.qm.esp_points = None
            self.qm.optimize = geomopt
            self.qm.restrained_dihedrals = np.array([dihedral])
            self.qm.directory = directory
            self.qm.setup() # QM object is reused, so it has to be properly set up before retrieving.
            qm_results.append(self.qm.retrieve())

        # Fit the dihedral parameters
        df = DihedralFitting()
        df.molecule = self
        df.dihedrals = dihedrals
        df.qm_results = qm_results
        df.result_directory = os.path.join(self.outdir, 'parameters', self.method.name,
                                           self.qm_method_name(), 'plots')

        # In case of FakeQM, the initial parameters are set to zeros.
        # It prevents DihedralFitting class from cheating :D
        if isinstance(self.qm, FakeQM2):
            df.zeroed_parameters = True

        # Fit dihedral parameters
        df.run()

        # Update atom types
        self.atomtype[:] = [self._rtf.type_by_name[name] for name in self.name]

    def _duplicateAtomType(self, atom_index):
        """Duplicate the type of the specified atom

           Duplicated types are named: original_name + "x" + number, e.g. ca --> cax0
        """

        # Get a type name
        type_ = self._rtf.type_by_index[atom_index]

        # if the type is already duplicated
        if re.match(FFMolecule._ATOM_TYPE_REG_EX, type_):
            return

        # Create a new atom type name
        i = 0
        while ('%sx%d' % (type_, i)) in self._rtf.types:
            i += 1
        newtype = '%sx%d' % (type_, i)
        logger.info('Create a new atom type %s from %s' % (newtype, type_))

        # Duplicate the type in RTF
        # TODO: move to RTF class
        self._rtf.type_by_index[atom_index] = newtype
        self._rtf.mass_by_type[newtype] = self._rtf.mass_by_type[type_]
        self._rtf.types.append(newtype)
        self._rtf.type_by_name[self._rtf.names[atom_index]] = newtype
        self._rtf.type_by_index[atom_index] = newtype
        self._rtf.typeindex_by_type[newtype] = self._rtf.typeindex_by_type[type_] + 1000
        self._rtf.element_by_type[newtype] = self._rtf.element_by_type[type_]

        # Rename the atom types of the equivalent atoms
        for index in self._equivalent_atoms[atom_index]:
            if atom_index != index:
                assert not re.match(FFMolecule._ATOM_TYPE_REG_EX, self._rtf.type_by_index[index])
                self._rtf.type_by_index[index] = newtype
                self._rtf.type_by_name[self._rtf.names[index]] = newtype

        # PRM parameters are duplicated during FF evaluation
        # TODO: move to PRM class
        FFEvaluate(self).run(self.coords[:, :, 0])

    def write(self, filename, sel=None, type=None, typemap=None):

        # TODO: remove type mapping
        if typemap:
            mol = self.copy()
            mol.atomtype[:] = [typemap[atomtype] for atomtype in self.atomtype]
            mol.write(filename, sel=sel, type=type)
        else:
            if filename.endswith('.rtf'):
                self._rtf.write(filename)
            elif filename.endswith('.prm'):
                self._prm.write(filename)
            else:
                super().write(filename, sel=sel, type=type)

    def writeParameters(self, original_molecule=None):

        paramDir = os.path.join(self.outdir, 'parameters', self.method.name, self.qm_method_name())
        os.makedirs(paramDir, exist_ok=True)

        typemap = None
        extensions = ('mol2', 'pdb', 'coor')

        if self.method == FFTypeMethod.CGenFF_2b6:
            extensions += ('psf', 'rtf', 'prm')

            # TODO: remove?
            f = open(os.path.join(paramDir, "input.namd"), "w")
            tmp = '''parameters mol.prm
paraTypeCharmm on
coordinates mol.pdb
bincoordinates mol.coor
temperature 0
timestep 0
1-4scaling 1.0
exclude scaled1-4
outputname .out
outputenergies 1
structure mol.psf
cutoff 20.
switching off
stepsPerCycle 1
rigidbonds none
cellBasisVector1 50. 0. 0.
cellBasisVector2 0. 50. 0.
cellBasisVector3 0. 0. 50.
run 0'''
            print(tmp, file=f)
            f.close()

        elif self.method in (FFTypeMethod.GAFF, FFTypeMethod.GAFF2):
            # types need to be remapped because Amber FRCMOD format limits the type to characters
            # writeFrcmod does this on the fly and returns a mapping that needs to be applied to the mol
            # TODO: get rid of this mapping
            frcFile = os.path.join(paramDir, 'mol.frcmod')
            typemap = self._prm.writeFrcmod(self._rtf, frcFile)  # TODO move to FFMolecule.write
            logger.info('Write FRCMOD file: %s' % frcFile)

            tleapFile = os.path.join(paramDir, 'tleap.in')
            with open(tleapFile, 'w') as file_:
                file_.write('loadAmberParams mol.frcmod\n')
                file_.write('A = loadMol2 mol.mol2\n')
                file_.write('saveAmberParm A structure.prmtop mol.crd\n')
                file_.write('quit\n')
            logger.info('Write tleap input file: %s' % tleapFile)

            # TODO: remove?
            f = open(os.path.join(paramDir, "input.namd"), "w")
            tmp = '''parmfile structure.prmtop
amber on
coordinates mol.pdb
bincoordinates mol.coor
temperature 0
timestep 0
1-4scaling 0.83333333
exclude scaled1-4
outputname .out
outputenergies 1
cutoff 20.
switching off
stepsPerCycle 1
rigidbonds none
cellBasisVector1 50. 0. 0.
cellBasisVector2 0. 50. 0.
cellBasisVector3 0. 0. 50.
run 0'''
            print(tmp, file=f)
            f.close()

        else:
            raise NotImplementedError

        for ext in extensions:
            file_ = os.path.join(paramDir, "mol." + ext)
            self.write(file_, typemap=typemap)
            logger.info('Write %s file: %s' % (ext.upper(), file_))

        if original_molecule:
            molFile = os.path.join(paramDir, 'mol-orig.mol2')
            original_molecule.write(molFile, typemap=typemap)
            logger.info('Write MOL2 file (with original coordinates): %s' % molFile)

if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
