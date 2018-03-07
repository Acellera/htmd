# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import scipy as sp
from scipy import constants as const
from htmd.numbautil import dihedralAngle


class FFEvaluate:
    """
    Compute potential energy of a molecule

    Support CHARMM and AMBER force fields.

    Parameters
    ----------
    ffmol : FFMolecule
        molecule object containing topology and force field parameters

    Examples
    --------

    # Create a FFMolecule object of benzamidine and assign AMBER FF parameters with GAFF2
    >>> import os
    >>> from htmd.home import home
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2)

    # Create FFEvaluate object of benzamidine
    >>> from htmd.parameterization.ffevaluate import FFEvaluate
    >>> mm = FFEvaluate(mol)
    >>> mm # doctest: +ELLIPSIS
    <htmd.parameterization.ffevaluate.FFEvaluate object at 0x...>

    >>> energies = mm.run(mol.coords[:, :, 0])
    >>> energies['bond'] # doctest: +ELLIPSIS
    5.91752585...
    >>> energies['angle'] # doctest: +ELLIPSIS
    2.961617...
    >>> energies['dihedral'] # doctest: +ELLIPSIS
    2.67383...
    >>> energies['improper'] # doctest: +ELLIPSIS
    0.006973...
    >>> energies['vdw'] # doctest: +ELLIPSIS
    4.629441...
    >>> energies['elec'] # doctest: +ELLIPSIS
    0.0...
    >>> energies['total'] # doctest: +ELLIPSIS
    16.18939...
    """

    ELEC_FACTOR = 1/(4*const.pi*const.epsilon_0) # Coulomb's constant
    ELEC_FACTOR *= const.elementary_charge**2 # Convert elementary charges to Coulombs
    ELEC_FACTOR /= const.angstrom # Convert Angstroms to meters
    ELEC_FACTOR *= const.Avogadro/(const.kilo*const.calorie) # Convert J to kcal/mol

    def __init__(self, molecule):

        self.mol = molecule
        self.natoms = self.mol.numAtoms
        self.rtf = self.mol._rtf
        self.prm = self.mol._prm

        # 1-2 and 1-3 exclusion matrix
        self.excl = sp.sparse.lil_matrix((self.natoms, self.natoms))
        for bond in self.mol.bonds:
            self.excl[bond[0], bond[1]] = self.excl[bond[1], bond[0]] = 1
        for angle in self.mol.angles:
            self.excl[angle[0], angle[2]] = self.excl[angle[2], angle[0]] = 1

        # 1-4 van der Waals scaling matrix
        self.s14 = sp.sparse.lil_matrix((self.natoms, self.natoms))
        for dihed in self.mol.dihedrals:
            self.s14[dihed[0], dihed[3]] = self.s14[dihed[3], dihed[0]] = 1

        # 1-4 electrostatic scaling matrix
        self.e14 = sp.sparse.lil_matrix((self.natoms, self.natoms))
        for dihed in self.mol.dihedrals:
            types = tuple([self.rtf.type_by_index[atom] for atom in dihed])
            parameters = self.prm.dihedralParam(*types)
            # Increament by 1, so intentionally to 0-scaled dihedrals can be distinguished
            self.e14[dihed[0], dihed[3]] = self.e14[dihed[3], dihed[0]] = parameters[0].e14 + 1

    def _evaluate_elec(self, coords):

        energy = 0.
        for i in range(self.natoms):
            qi = self.mol.charge[i]

            for j in range(i + 1, self.natoms):
                if self.excl[i, j]:
                    continue

                qj = self.mol.charge[j]
                scale = self.e14[i, j] - 1 if self.e14[i, j] != 0 else 1  # If 0, assume it's not a 1-4 term
                dist = np.linalg.norm(coords[j, :] - coords[i, :])
                energy += self.ELEC_FACTOR * scale * qi * qj / dist

        return energy

    def _evaluate_vdw(self, coords):

        energy = 0.
        for i in range(self.natoms):
            for j in range(i + 1, self.natoms):
                if self.excl[i, j]:
                    continue

                A, B = self.prm.vdwParam(self.rtf.type_by_index[i], self.rtf.type_by_index[j], self.s14[i, j])
                dist = np.linalg.norm(coords[j, :] - coords[i, :])
                energy += A/dist**12 - B/dist**6

        return energy

    def _evaluate_bonds(self, coords):

        energy = 0.
        for bond in self.mol.bonds:
            types = tuple([self.rtf.type_by_index[atom] for atom in bond])
            parameters = self.prm.bondParam(*types)
            dist = np.linalg.norm(coords[bond[0], :] - coords[bond[1], :])
            energy += parameters.k0 * (dist - parameters.r0)**2

        return energy

    def _evaluate_angles(self, coords):

        energy = 0.
        for angle in self.mol.angles:

            # TODO: move htmd.molecule.util
            r23 = coords[angle[2], :] - coords[angle[1], :]
            r21 = coords[angle[0], :] - coords[angle[1], :]
            cos_theta = np.dot(r21, r23) / (np.linalg.norm(r21) * np.linalg.norm(r23))
            cos_theta = np.clip(cos_theta, -1.0, 1.0)
            theta = np.arccos(cos_theta)

            types = tuple([self.rtf.type_by_index[atom] for atom in angle])
            parameters = self.prm.angleParam(*types)
            theta0 = np.deg2rad(parameters.theta0)
            energy += parameters.k0 * (theta - theta0)**2

            # TODO: no idea what is happening here. Dead code!
            if (parameters.kUB is None) and (parameters.kUB != 0.):
                r13 = r23 - r21
                r12len = np.linalg.norm(r13)
                dist = r12len / parameters.rUB
                energy += parameters.kUB * dist / r12len

        return energy

    def _evaluate_dihedrals(self, coords):

        energy = 0.
        for dihedral in self.mol.dihedrals:
            types = tuple([self.rtf.type_by_index[atom] for atom in dihedral])
            parameters = self.prm.dihedralParam(*types)
            energy += self._evaluateTorsion(coords[dihedral, :], parameters)

        return energy

    def _evaluate_impropers(self, coords):

        energy = 0.
        for improper in self.mol.impropers:
            types = tuple([self.rtf.type_by_index[atom] for atom in improper])
            parameters = self.prm.improperParam(*types)
            energy += self._evaluateTorsion(coords[improper, :], parameters)

        return energy

    @staticmethod
    def _evaluateTorsion(coords, torsions):

        phi = dihedralAngle(coords)

        energy = 0.
        for torsion in torsions:
            k = torsion.k0
            n = torsion.n
            phi0 = np.deg2rad(torsion.phi0)

            if n > 0:
                energy += k * (1. + np.cos(n * phi - phi0))  # This is also AMBER improper
            else:
                energy += k * (phi - phi0)**2  # This is a CHARMM improper

        return energy

    def run(self, coords):
        """
        Compute potential energy of the molecule with given atomic coordinates

        Parameters
        ----------
        coords : numpy.ndarray
            Coordinates of the molecule

        Return
        ------
        energies : dict
            Dictionary containing potential energy and its components
        """

        assert coords.ndim == 2
        assert coords.shape[1] == 3

        energy = dict()
        energy['elec'] = self._evaluate_elec(coords)
        energy['vdw'] = self._evaluate_vdw(coords)
        energy['bond'] = self._evaluate_bonds(coords)
        energy['angle'] = self._evaluate_angles(coords)
        energy['dihedral'] = self._evaluate_dihedrals(coords)
        energy['improper'] = self._evaluate_impropers(coords)
        energy['total'] = sum(energy.values())

        return energy


def _openmm_energy_charmm(psfFile, rtfFile, prmFile, coords):

    import parmed
    from simtk import unit
    from simtk import openmm

    # Read PSF and PRM files
    psf = parmed.charmm.CharmmPsfFile(psfFile)
    prm = parmed.charmm.CharmmParameterSet(rtfFile, prmFile)

    # Create OpenMM
    system = psf.createSystem(prm)
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 2 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')
    context = openmm.Context(system, integrator, platform)

    # Run OpenMM with given coordinates
    context.setPositions(coords * unit.angstrom)
    energies = parmed.openmm.energy_decomposition(psf, context)

    return energies


def _openmm_energy_amber(mol2File, frcmodFile, coords, tempDir=None):

    from tempfile import TemporaryDirectory
    from subprocess import call
    import parmed
    from simtk import openmm
    from simtk import unit

    with TemporaryDirectory() as tmpDir:
        tmpDir = tempDir if tempDir else tmpDir

        # Create "tleap" input
        with open(os.path.join(tmpDir, 'tleap.inp'), 'w') as file:
            file.writelines(('loadAmberParams %s\n' % frcmodFile,
                             'MOL = loadMol2 %s\n' % mol2File,
                             'saveAmberParm MOL mol.prmtop mol.inpcrd\n',
                             'quit'))

        # Run "tleap" to generate mol.prmtop
        with open(os.path.join(tmpDir, 'tleap.out'), 'w') as out:
            call(('tleap', '-f', 'tleap.inp'), cwd=tmpDir, stdout=out)

        # Read PRMTOP file
        prmtop = parmed.amber.LoadParm(os.path.join(tmpDir, 'mol.prmtop'))

    # Create OpenMM
    system = prmtop.createSystem()
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 2 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')
    context = openmm.Context(system, integrator, platform)

    # Run OpenMM with given coordinates
    context.setPositions(coords*unit.angstrom)
    energies = parmed.openmm.energy_decomposition(prmtop, context)

    return energies

if __name__ == '__main__':

    import os
    from tempfile import TemporaryDirectory
    from htmd.home import home
    from htmd.parameterization.ffmolecule import FFMolecule
    from htmd.parameterization.fftype import FFTypeMethod

    np.random.seed(20170801)  # Make the tests deterministic

    molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
    methods = (FFTypeMethod.CGenFF_2b6, FFTypeMethod.GAFF, FFTypeMethod.GAFF2)

    # TODO: remove then MATCH is fixed on Mac
    methods = methods[1:] if os.environ.get('TRAVIS_OS_NAME') == 'osx' else methods

    for method in methods:
        mol = FFMolecule(molFile, method=method)

        # Generate random charges
        for name in mol._rtf.charge_by_name:
            mol._rtf.charge_by_name[name] = 0.1*np.random.randn()

        # Generate a list of original and randomly distorted coordinates
        coords = mol.coords[:, :, 0]
        coordsList = [coords] + [coords + 0.01*np.random.randn(*coords.shape) for _ in range(9)]

        for coords in coordsList:

            with TemporaryDirectory() as tmpDir:

                if method == FFTypeMethod.CGenFF_2b6:
                    psfFile = os.path.join(tmpDir, 'mol.psf')
                    rtfFile = os.path.join(tmpDir, 'mol.rtf')
                    prmFile = os.path.join(tmpDir, 'mol.prm')
                    mol.write(psfFile)
                    mol._rtf.write(rtfFile)
                    mol._prm.write(prmFile)
                    reference = _openmm_energy_charmm(psfFile, rtfFile, prmFile, coords)

                elif method in (FFTypeMethod.GAFF, FFTypeMethod.GAFF2):
                    mol2File = os.path.join(tmpDir, 'mol.mol2')
                    frcmodFile = os.path.join(tmpDir, 'mol.frcmod')
                    map = mol._prm.writeFrcmod(mol._rtf, frcmodFile)
                    mol.write(mol2File, typemap=map)
                    reference = _openmm_energy_amber(mol2File, frcmodFile, coords)

                else:
                    assert False

            ff = FFEvaluate(mol)
            result = ff.run(coords)

            if not np.isclose(reference['total'], result['total'], atol=5e-5):
                print('\nReference:')
                for term in reference:
                    print(term, reference[term])
                print('\nResult:')
                for term in result:
                    print(term, result[term])
                assert False

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
