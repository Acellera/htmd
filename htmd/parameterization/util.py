import numpy as np
import logging
import re
import os

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


def minimize(mol, qm, outdir):
    assert mol.numFrames == 1

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
    return mol


def fitDihedrals(mol, qm, method, prm, all_dihedrals, dihedrals, outdir, geomopt=True):
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

        # Set rotamer coordinates
        angles = np.linspace(-np.pi, np.pi, num=nrotamers, endpoint=False)
        for frame, angle in enumerate(angles):
            tmpmol.frame = frame
            tmpmol.setDihedral(dihedral, angle, bonds=tmpmol.bonds)

        molecules.append(tmpmol)

    # Create directories for QM data
    directories = []
    dihedral_directory = 'dihedral-opt' if geomopt else 'dihedral-single-point'
    for dihedral in dihedrals:
        dihedral_name = '-'.join(mol.name[dihedral])
        directory = os.path.join(outdir, dihedral_directory, dihedral_name, _qm_method_name(qm))
        os.makedirs(directory, exist_ok=True)
        directories.append(directory)

    # Setup and submit QM calculations
    for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
        qm.molecule = molecule
        qm.esp_points = None
        qm.optimize = geomopt
        qm.restrained_dihedrals = np.array([dihedral])
        qm.directory = directory
        qm.setup()
        qm.submit()

    # Wait and retrieve QM calculation data
    qm_results = []
    for molecule, dihedral, directory in zip(molecules, dihedrals, directories):
        qm.molecule = molecule
        qm.esp_points = None
        qm.optimize = geomopt
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
    df.result_directory = os.path.join(outdir, 'parameters', method, _qm_method_name(qm), 'plots')

    # In case of FakeQM, the initial parameters are set to zeros.
    # It prevents DihedralFitting class from cheating :D
    if isinstance(qm, FakeQM2):
        df.zeroed_parameters = True

    # Fit dihedral parameters
    df.run()

    return prm


