# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import logging
import re
import os
import subprocess
from tempfile import TemporaryDirectory

logger = logging.getLogger(__name__)


def getEquivalentsAndDihedrals(mol):

    from htmd.molecule.util import guessAnglesAndDihedrals
    from htmd.parameterization.detect import detectParameterizableDihedrals, detectEquivalentAtoms

    mol = mol.copy()

    # Guess bonds
    if len(mol.bonds) == 0:
        logger.warning('No bonds found! Guessing them...')
        mol.bonds = mol._guessBonds()

    mol.angles, mol.dihedrals = guessAnglesAndDihedrals(mol.bonds)
    equivalents = detectEquivalentAtoms(mol)
    all_dihedrals = detectParameterizableDihedrals(mol)
    return mol, equivalents, all_dihedrals


def guessElements(mol, fftypemethod):
    """
    Guess element from an atom name
    """

    from htmd.parameterization.fftype import fftypemethods

    if fftypemethod not in fftypemethods:
        raise ValueError('Invalid "fftypemethod": {}. Valid methods: {}'
                         ''.format(fftypemethod, ','.join(fftypemethods)))

    elements = {}
    elements['CGenFF_2b6'] = ['H', 'C', 'N', 'O', 'F', 'S', 'P', 'Cl', 'Br', 'I']
    elements['GAFF']       = ['H', 'C', 'N', 'O', 'F', 'S', 'P', 'Cl', 'Br', 'I']
    elements['GAFF2']      = ['H', 'C', 'N', 'O', 'F', 'S', 'P', 'Cl', 'Br', 'I']

    mol = mol.copy()

    for i, name in enumerate(mol.name):

        candidates = [element for element in elements[fftypemethod] if name.capitalize().startswith(element)]

        if len(candidates) == 1:
            mol.element[i] = candidates[0]
            continue

        if candidates == ['C', 'Cl']:

            if len(mol.bonds) == 0:
                raise RuntimeError('No chemical bonds found in the molecule')

            # Create a molecular graph
            import networkx as nx
            graph = nx.Graph()
            graph.add_edges_from(mol.bonds)

            if len(graph[i]) in (2, 3, 4):
                mol.element[i] = 'C'
                continue

            if len(graph[i]) == 1:
                mol.element[i] = 'Cl'
                continue

        raise ValueError('Cannot guess element from atom name: {}. '
                         'It does not match any of the expected elements ({}) for {}.'
                         ''.format(name, elements[fftypemethod], fftypemethod))

    return mol


def centreOfMass(mol):
    return np.dot(mol.masses, mol.coords[:, :, mol.frame]) / np.sum(mol.masses)


def getDipole(mol):
    """Calculate the dipole moment (in Debyes) of the molecule"""
    from scipy import constants as const

    if mol.masses.sum() == 0:
        logger.warning('No masses found in Molecule. Cannot calculate dipole.')
        return np.zeros(4)
    else:
        coords = mol.coords[:, :, mol.frame] - centreOfMass(mol)

    dipole = np.zeros(4)
    dipole[:3] = np.dot(mol.charge, coords)
    dipole[3] = np.linalg.norm(dipole[:3]) # Total dipole moment
    dipole *= 1e11*const.elementary_charge*const.speed_of_light # e * Ang --> Debye (https://en.wikipedia.org/wiki/Debye)

    return dipole


def _qm_method_name(qm):
    basis = qm.basis
    basis = re.sub('\+', 'plus', basis)  # Replace '+' with 'plus'
    basis = re.sub('\*', 'star', basis)  # Replace '*' with 'star'
    name = qm.theory + '-' + basis + '-' + qm.solvent
    return name


def getFixedChargeAtomIndices(mol, fix_charge):
    fixed_atom_indices = []
    for fixed_atom_name in fix_charge:

        if fixed_atom_name not in mol.name:
            raise ValueError('Atom {} is not found. Check --fix-charge arguments'.format(fixed_atom_name))

        for aton_index in range(mol.numAtoms):
            if mol.name[aton_index] == fixed_atom_name:
                fixed_atom_indices.append(aton_index)
                logger.info('Charge of atom {} is fixed to {}'.format(fixed_atom_name, mol.charge[aton_index]))
    return fixed_atom_indices


def minimize(mol, qm, outdir, min_type='qm', mm_minimizer=None):
    assert mol.numFrames == 1

    if min_type == 'qm':
        mindir = os.path.join(outdir, "minimize", _qm_method_name(qm))
        os.makedirs(mindir, exist_ok=True)

        qm.molecule = mol
        qm.esp_points = None
        qm.optimize = True
        qm.restrained_dihedrals = None
        qm.directory = mindir
        results = qm.run()
        if results[0].errored:
            raise RuntimeError('\nQM minimization failed! Check logs at %s\n' % mindir)

        mol = mol.copy()
        # Replace coordinates with the minimized set
        mol.coords = np.atleast_3d(np.array(results[0].coords, dtype=np.float32))
    elif min_type == 'mm':
        mol = mol.copy()
        mol.coords[:, :, 0] = mm_minimizer.minimize(mol.coords)
    elif min_type == 'None':
        pass
    else:
        raise RuntimeError('Invalid minimization mode {}. Check parameterize help.'.format(min_type))

    return mol


def fitDihedrals(mol, qm, method, prm, all_dihedrals, dihedrals, outdir, dihed_opt_type='qm', mm_minimizer=None):
    """
    Dihedrals passed as 4 atom indices
    """
    from htmd.parameterization.dihedral import DihedralFitting
    from htmd.qm import FakeQM2

    # Create molecules with rotamers
    molecules = []
    for dihedral in dihedrals:
        nrotamers = 36  # Number of rotamers for each dihedral to compute

        # Create a copy of molecule with "nrotamers" frames
        tmpmol = mol.copy()
        tmpmol.coords = np.tile(tmpmol.coords, (1, 1, nrotamers))
        tmpmol.box = np.zeros((3, tmpmol.numFrames), dtype=np.float32)

        # Set rotamer coordinates
        angles = np.linspace(-np.pi, np.pi, num=nrotamers, endpoint=False)
        for frame, angle in enumerate(angles):
            tmpmol.frame = frame
            tmpmol.setDihedral(dihedral, angle, bonds=tmpmol.bonds)

        molecules.append(tmpmol)

    # Minimize with MM if requested
    if dihed_opt_type == 'mm':
        for dihedral, tmpmol in zip(dihedrals, molecules):
            for f in range(tmpmol.numFrames):
                tmpmol.coords[:, :, f] = mm_minimizer.minimize(tmpmol.coords[:, :, f], restrained_dihedrals=[dihedral,])

    # Create directories for QM data
    directories = []
    dihedral_directory = 'dihedral-opt' if dihed_opt_type == 'qm' else 'dihedral-single-point'
    for dihedral in dihedrals:
        dihedral_name = '-'.join(mol.name[dihedral])
        directory = os.path.join(outdir, dihedral_directory, dihedral_name, _qm_method_name(qm))
        os.makedirs(directory, exist_ok=True)
        directories.append(directory)

    # Setup and submit QM calculations
    for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
        qm.molecule = molecule
        qm.esp_points = None
        qm.optimize = (dihed_opt_type == 'qm')
        qm.restrained_dihedrals = np.array([dihedral])
        qm.directory = directory
        qm.setup()
        qm.submit()

    # Wait and retrieve QM calculation data
    qm_results = []
    for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
        qm.molecule = molecule
        qm.esp_points = None
        qm.optimize = (dihed_opt_type == 'qm')
        qm.restrained_dihedrals = np.array([dihedral])
        qm.directory = directory
        qm.setup()  # QM object is reused, so it has to be properly set up before retrieving.
        qm_results.append(qm.retrieve())

    # Fit the dihedral parameters
    df = DihedralFitting()
    df.parmedMode = True
    df.parameters = prm
    df._parameterizable_dihedrals = all_dihedrals
    df.molecule = mol
    df.dihedrals = dihedrals
    df.qm_results = qm_results
    df.result_directory = os.path.join(outdir, 'parameters', method, _qm_method_name(qm))

    # In case of FakeQM, the initial parameters are set to zeros.
    # It prevents DihedralFitting class from cheating :D
    if isinstance(qm, FakeQM2):
        df.zeroed_parameters = True

    # Fit dihedral parameters
    df.run()

    return df.parameters

def guessBondType(mol):

    """
    Guess bond type with Antechamber

    Parameters
    ----------
    mol: Molecule
        Molecule to guess the bond types

    Return
    ------
    results: Molecule
        Copy of the molecule with the bond type set

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmd.molecule.molecule import Molecule

    >>> molFile = os.path.join(home('test-qm'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.bondtype[:] = "un"

    >>> new_mol = guessBondType(mol)
    >>> assert new_mol is not mol
    >>> new_mol.bondtype
    array(['1', '1'], dtype=object)

    >>> molFile = os.path.join(home('test-param'), 'H2O2.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.bondtype[:] = "un"

    >>> new_mol = guessBondType(mol)
    >>> assert new_mol is not mol
    >>> new_mol.bondtype
    array(['1', '1', '1'], dtype=object)

    >>> molFile = os.path.join(home('test-param'), 'benzamidine.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.bondtype[:] = "un"

    >>> new_mol = guessBondType(mol)
    >>> assert new_mol is not mol
    >>> new_mol.bondtype
    array(['ar', 'ar', '1', 'ar', '1', 'ar', '1', 'ar', '1', 'ar', '1', '1',
           '2', '1', '1', '1', '1', '1'], dtype=object)

    """

    from htmd.molecule.molecule import Molecule

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    mol = mol.copy()

    with TemporaryDirectory() as tmpDir:
        old_name = os.path.join(tmpDir, 'old.mol2')
        new_name = os.path.join(tmpDir, 'new.mol2')

        mol.write(old_name)

        cmd = ['antechamber',
               '-fi', 'mol2', '-i', old_name,
               '-fo', 'mol2', '-o', new_name]
        if subprocess.call(cmd, cwd=tmpDir) != 0:
            raise RuntimeError('"antechamber" failed')

        mol.bondtype[:] = Molecule(new_name).bondtype

    return mol

def makeAtomNamesUnique(mol):
    """
    Make atom names unique by appending/incrementing terminal digits.
    Already unique names are preserved.

    Parameters
    ----------
    mol: Molecule
        Molecule to make atom name unique

    Return
    ------
    results: Molecule
        Copy of the molecule with the atom names set

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmd.molecule.molecule import Molecule
    >>> molFile = os.path.join(home('test-param'), 'H2O2.mol2')
    >>> mol = Molecule(molFile)

    >>> mol.name[:] = ['A', 'A', 'A', 'A']
    >>> new_mol = makeAtomNamesUnique(mol)
    >>> assert new_mol is not mol
    >>> new_mol.name
    array(['A', 'A0', 'A1', 'A2'], dtype=object)

    >>> mol.name[:] = ['A', 'A', 'A', 'A0']
    >>> new_mol = makeAtomNamesUnique(mol)
    >>> assert new_mol is not mol
    >>> new_mol.name
    array(['A', 'A1', 'A2', 'A0'], dtype=object)

    >>> mol.name[:] = ['A', 'B', 'A', 'B']
    >>> new_mol = makeAtomNamesUnique(mol)
    >>> assert new_mol is not mol
    >>> new_mol.name
    array(['A', 'B', 'A0', 'B0'], dtype=object)

    >>> mol.name[:] = ['A', 'B', 'C', 'D']
    >>> new_mol = makeAtomNamesUnique(mol)
    >>> assert new_mol is not mol
    >>> new_mol.name
    array(['A', 'B', 'C', 'D'], dtype=object)

    >>> mol.name[:] = ['1A', '1A', 'A1B1', 'A1B1']
    >>> new_mol = makeAtomNamesUnique(mol)
    >>> assert new_mol is not mol
    >>> new_mol.name
    array(['1A', '1A0', 'A1B1', 'A1B2'], dtype=object)
    """

    from htmd.molecule.molecule import Molecule

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be an instance of {}'.format(Molecule))

    mol = mol.copy()

    for i, name in enumerate(mol.name):
        while np.sum(name == mol.name) > 1: # Check for identical names
            j = np.flatnonzero(name == mol.name)[1] # Get the second identical name index
            prefix, sufix = re.match('(.*?\D*)(\d*)$', mol.name[j]).groups()
            sufix = 0 if sufix == '' else int(sufix)
            while prefix + str(sufix) in mol.name: # Search for a unique name
                sufix += 1
            mol.name[j] = prefix + str(sufix)

    return mol

def detectChiralCenters(mol):
    """
    Detect chiral centers

    Parameters
    ----------
    mol: Molecule
        Molecule to detect chiral centers

    Return
    ------
    results: List of tuples
        List of chircal centers, where the chiral centers are tuples made of an atom index and a label ('S', 'R', '?').

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmd.molecule.molecule import Molecule

    >>> molFile = os.path.join(home('test-param'), 'H2O2.mol2')
    >>> mol = Molecule(molFile)
    >>> detectChiralCenters(mol)
    []

    >>> molFile = os.path.join(home('test-param'), 'fluorchlorcyclopronol.mol2')
    >>> mol = Molecule(molFile)
    >>> detectChiralCenters(mol)
    [(0, 'R'), (2, 'S'), (4, 'R')]
    """

    from htmd.molecule.molecule import Molecule
    from rdkit.Chem import MolFromMol2File, AssignAtomChiralTagsFromStructure, FindMolChiralCenters

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    # Set atom types to elements, overwise rdkit refuse to read a MOL2 file
    htmd_mol = mol.copy()
    htmd_mol.atomtype = htmd_mol.element

    # Convert Molecule to rdkit Mol
    with TemporaryDirectory() as tmpDir:
        filename = os.path.join(tmpDir, 'mol.mol2')
        htmd_mol.write(filename)
        rdkit_mol = MolFromMol2File(filename, removeHs=False)
    assert mol.numAtoms == rdkit_mol.GetNumAtoms()

    # Detect chiral centers and assign their labels
    AssignAtomChiralTagsFromStructure(rdkit_mol)
    chiral_centers = FindMolChiralCenters(rdkit_mol, includeUnassigned=True)

    return chiral_centers

if __name__ == '__main__':

    import sys
    import doctest

    sys.exit(doctest.testmod().failed)
