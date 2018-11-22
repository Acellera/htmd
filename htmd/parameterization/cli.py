# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import argparse
import logging
import numpy as np

from htmd.version import version
from htmd.parameterization.fftype import fftypemethods

logger = logging.getLogger(__name__)


def getArgumentParser():

    parser = argparse.ArgumentParser(description='Acellera small molecule parameterization tool')

    parser.add_argument('filename', help='Molecule file in MOL2 format')
    parser.add_argument('-c', '--charge', type=int,
                        help='Total charge of the molecule (default: sum of partial charges)')
    parser.add_argument('-l', '--list', action='store_true', help='List parameterizable dihedral angles')
    parser.add_argument('--rtf-prm', nargs=2, metavar='<filename>', help='CHARMM RTF and PRM files')
    parser.add_argument('-ff', '--forcefield', default='GAFF2', choices=fftypemethods,
                        help='Initial atom type and parameter assignment (default: %(default)s)')
    parser.add_argument('--fix-charge', nargs='+', default=[], metavar='<atom name>',
                        help='Fix atomic charge during charge fitting (default: none)')
    parser.add_argument('-d', '--dihedral', nargs='+', default=[], metavar='A1-A2-A3-A4',
                        help='Select dihedral angle to parameterize (default: all parameterizable dihedral angles)')
    parser.add_argument('--code', default='Psi4', choices=['Psi4', 'Gaussian'],
                        help='QM code (default: %(default)s)')
    parser.add_argument('--theory', default='B3LYP', choices=['HF', 'B3LYP', 'wB97X-D'],
                        help='QM level of theory (default: %(default)s)')
    parser.add_argument('--basis', default='cc-pVDZ', choices=['6-31G*', '6-31+G*', '6-311G**', '6-311++G**', 'cc-pVDZ',
                                                               'aug-cc-pVDZ'],
                        help='QM basis set (default: %(default)s)')
    parser.add_argument('--environment', default='vacuum', choices=['vacuum', 'PCM'],
                        help='QM environment (default: %(default)s)')
    parser.add_argument('--no-min', action='store_false', dest='minimize',
                        help='DDEPRECATED: use `--min-type` instead')
    parser.add_argument('--min-type', default='qm', dest='min_type', choices=['None', 'qm', 'mm'],
                        help='Type of initial structure optimization (default: %(default)s)')
    parser.add_argument('--charge-type', default='ESP', choices=['None', 'Gasteiger', 'AM1-BCC', 'ESP'],
                        help='Partial atomic charge type (default: %(default)s)')
    parser.add_argument('--no-dihed', action='store_false', dest='fit_dihedral',
                        help='Do not perform QM scanning of dihedral angles')
    parser.add_argument('--no-dihed-opt', action='store_false', dest='optimize_dihedral',
                        help='DEPRECATED: use `--scan-type` instead')
    parser.add_argument('--scan-type', default='qm', dest='dihed_opt_type', choices=['None', 'qm', 'mm'],
                        help='Type of structure optimization when scanning dihedral angles (default: %(default)s)')
    parser.add_argument('--dihed-num-searches', default=None, type=int,
                        help='Number of random searches during the dihedral parameter fitting')
    parser.add_argument('-q', '--queue', default='local', choices=['local', 'Slurm', 'LSF', 'AceCloud'],
                        help='QM queue (default: %(default)s)')
    parser.add_argument('-n', '--ncpus', default=None, type=int, help='Number of CPU per QM job (default: queue '
                                                                      'defaults)')
    parser.add_argument('-m', '--memory', default=None, type=int, help='Maximum amount of memory in MB to use.')
    parser.add_argument('--groupname', default=None, help=argparse.SUPPRESS)
    parser.add_argument('-o', '--outdir', default='./', help='Output directory (default: %(default)s)')
    parser.add_argument('--seed', default=20170920, type=int,
                        help='Random number generator seed (default: %(default)s)')
    parser.add_argument('--version', action='version', version=version())

    # Enable replacement of any real QM class with FakeQM.
    # This is intedended for debugging only and should be kept hidden.
    parser.add_argument('--fake-qm', action='store_true', default=False, dest='fake_qm', help=argparse.SUPPRESS)

    # NNP module name
    parser.add_argument('--nnp', help=argparse.SUPPRESS)

    # Debug mode
    parser.add_argument('--debug', action='store_true', default=False, dest='debug', help=argparse.SUPPRESS)

    return parser


def _prepare_molecule(args):

    from htmd.molecule.molecule import Molecule
    from htmd.parameterization.util import makeAtomNamesUnique, guessElements, detectChiralCenters

    logger.info('=== Molecule ===')

    # Check if the molecule faile exist
    if not os.path.exists(args.filename):
        raise ValueError('File {} cannot be found'.format(args.filename))

    # Check the file extension
    if os.path.splitext(args.filename)[1] != '.mol2':
        raise RuntimeError('{} has to be in MOL2 format'.format(args.filename))

    # Read the molecule file
    mol = Molecule(args.filename, guessNE='bonds', guess=('angles', 'dihedrals'))
    logger.info('Read a molecule from file: {}'.format(args.filename))

    # Check if the file contain just one conformation
    if mol.numFrames != 1:
        raise RuntimeError('{} has to contain only one conformation, but found {}'.format(args.filename, mol.numFrames))

    # Check the number of atoms
    if mol.numAtoms < 2:
        raise RuntimeError('Molecule has to contain more than one atom')
    logger.info('Number of atoms: {}'.format(mol.numAtoms))

    # Check the number of bonds
    if mol.bonds.size < 1:
        raise RuntimeError('Molecule has to have at least one bond')
    logger.info('Number of bonds: {}'.format(mol.bonds.size))

    # Set the molecular charge
    charge = int(round(np.sum(mol.charge)))
    if args.charge is None:
        args.charge = charge
        logger.info('Molecular charge is set to {} by adding up the atomic charges in {}'
                    ''.format(args.charge, args.filename))
    else:
        logger.info('Molecular charge is set to {}'.format(args.charge))
        if args.charge_type == 'None' and args.charge != charge:
            raise ValueError('The molecular charge is set to {}, but the partial atomic charges in {} '
                             'add up to {}'.format(args.charge, args.filename, charge))

    # Make atom names unique if needed
    if np.unique(mol.name).size != mol.numAtoms:
        logger.warning('Atom names in the molecule are not unique!')
        new_mol = makeAtomNamesUnique(mol)
        for i, (old_name, new_name) in enumerate(zip(mol.name, new_mol.name)):
            if old_name != new_name:
                logger.warning('Renamed atom {:3d}: {:4s} --> {:4s}'.format(i, old_name, new_name))
        mol = new_mol

    # Guess elements
    # TODO: it should not depend on FF
    mol = guessElements(mol, args.forcefield)
    logger.info('Elements detected:')
    for name, element in zip(mol.name, mol.element):
        logger.info('    {:6s}: {:2s}'.format(name, element))

    # Check residue names
    if not np.all(mol.resname == mol.resname[0]):
        raise RuntimeError('All atoms of the molecule need to have the same residue name')
    logger.info('Residue name: {}'.format(mol.resname[0]))

    # Set segment ID
    mol.segid[:] = 'L' # Note: it is need to write complete PDB files
    logger.info('Segment ID: {}'.format(mol.segid[0]))

    # Detect chiral centers
    chiral_centers = detectChiralCenters(mol)
    if len(chiral_centers) > 0:
        logger.info('Chiral centers:')
        for atom_index, chiral_label in chiral_centers:
            logger.info(' {:4} {}'.format(mol.name[atom_index], chiral_label))

    return mol


def _select_dihedrals(mol, args):

    from htmd.parameterization.detect import detectParameterizableDihedrals

    logger.info('=== Dihedral angles ===')

    # Detect parameterizable dihedral angles
    parameterizable_dihedrals = detectParameterizableDihedrals(mol)
    logger.info('Parameterizable dihedral angles (and their equivalents):')
    for i, equivalent_dihedrals in enumerate(parameterizable_dihedrals):
        names = ['-'.join(mol.name[list(dihedral)]) for dihedral in equivalent_dihedrals]
        line = '    {:2d}: {}'.format(i+1, names[0])
        if names[1:]:
            line += ' ({})'.format(', '.join(names[1:]))
        logger.info(line)

    # Print parameterize dihedral list to a file and exit
    # TODO factor this out
    if args.list:
        dihedral_file = 'torsions.txt'
        with open(dihedral_file, 'w') as file:
            for dihedrals in parameterizable_dihedrals:
                name = '-'.join(mol.name[list(dihedrals[0])])
                file.write('{}\n'.format(name))
            logger.info('Write the list of the parameterizable dihedral angles to {}'.format(dihedral_file))
        sys.exit()

    # Select dihedrals to parameterize
    selected_dihedrals = []
    if args.fit_dihedral:
        if args.dihedral:
            name2dihedral = {'-'.join(mol.name[list(dihedrals[0])]): dihedrals[0] for dihedrals in parameterizable_dihedrals}
            for dihedral_name in args.dihedral:
                if dihedral_name not in name2dihedral.keys():
                    raise ValueError('%s is not recognized as a rotatable dihedral angle' % dihedral_name)
                selected_dihedrals.append(list(name2dihedral[dihedral_name]))
        else:
            # By default parameterize all the dihedrals
            selected_dihedrals = [list(dihedrals[0]) for dihedrals in parameterizable_dihedrals]

    if len(selected_dihedrals) > 0:
        logger.info('Selected dihedral angles:')
        for i, dihedral in enumerate(selected_dihedrals):
            name = '-'.join(mol.name[list(dihedral)])
            logger.info('    {:2d}: {}'.format(i+1, name))
    else:
        logger.info('No dihedral angles selected')

    return selected_dihedrals


def _get_queue(args):

    logger.info('=== Computation queue ===')

    # Create a queue
    logger.info('Queue type: {}'.format(args.queue))
    if args.queue == 'local':
        from htmd.queues.localqueue import LocalCPUQueue
        queue = LocalCPUQueue()
    elif args.queue == 'Slurm':
        from htmd.queues.slurmqueue import SlurmQueue
        queue = SlurmQueue(_configapp=args.code.lower())
    elif args.queue == 'LSF':
        from htmd.queues.lsfqueue import LsfQueue
        queue = LsfQueue(_configapp=args.code.lower())
    elif args.queue == 'PBS':
        from htmd.queues.pbsqueue import PBSQueue
        queue = PBSQueue()  # TODO: configure
    elif args.queue == 'AceCloud':
        from htmd.queues.acecloudqueue import AceCloudQueue
        queue = AceCloudQueue()  # TODO: configure
        queue.groupname = args.groupname
        queue.hashnames = True
    else:
        raise AssertionError()

    # Configure the queue
    if args.ncpus:
        logger.info('Overriding ncpus to {}'.format(args.ncpus))
        queue.ncpu = args.ncpus
    if args.memory:
        logger.info('Overriding memory to {}'.format(args.memory))
        queue.memory = args.memory

    return queue


def _get_qm_calculator(args, queue):

    from htmd.qm import Psi4, Gaussian, FakeQM2

    logger.info('=== QM configuration ===')

    # Create a QM object
    logger.info('Code: {}'.format(args.code))
    if args.code == 'Psi4':
        qm = Psi4()
    elif args.code == 'Gaussian':
        qm = Gaussian()
    else:
        raise AssertionError()

    # Override with a FakeQM object
    if args.fake_qm:
        qm = FakeQM2()
        logger.warning('Using FakeQM')

    # Configure the QM object
    qm.theory = args.theory
    logger.info('Theory: {}'.format(qm.theory))
    qm.basis = args.basis
    logger.info('Basis sets: {}'.format(qm.basis))
    qm.solvent = args.environment
    logger.info('Environment: {}'.format(qm.solvent))
    qm.queue = queue
    qm.charge = args.charge

    # Update B3LYP to B3LYP-D3
    # TODO: this is silent and not documented stuff
    if qm.theory == 'B3LYP':
        qm.correction = 'D3'
        logger.warning('The dispersion correction is set to {}'.format(qm.correction))

    # Update basis sets
    # TODO: this is silent and not documented stuff
    if args.charge < 0 and qm.solvent == 'vacuum':
        if qm.basis == '6-31G*':
            qm.basis = '6-31+G*'
        if qm.basis == '6-311G**':
            qm.basis = '6-311++G**'
        if qm.basis == 'cc-pVDZ':
            qm.basis = 'aug-cc-pVDZ'
        logger.warning('Basis sets are changed to {}'.format(qm.basis))

    return qm

def _get_nnp_calculator(args, queue):

    if args.nnp:
        import importlib
        from htmd.qm.custom import CustomQM

        logger.info('=== NNP configuration ===')

        # Get NNP calculator
        nnp_module = importlib.import_module(args.nnp)
        logger.info('NNP module: {}'.format(nnp_module))
        nnp_calculator = nnp_module.get_calculator()
        logger.info('NNP calculator: {}'.format(nnp_calculator))

        # Create a custom "QM" object
        nnp = CustomQM(verbose=False)
        nnp.theory = args.nnp
        nnp.basis = ''
        nnp.solvent = ''
        nnp.queue = queue
        nnp.charge = args.charge
        nnp.calculator = nnp_calculator

    else:
        nnp = None

    return nnp


def _get_initial_parameters(mol, args):

    from htmd.parameterization.fftype import fftype

    logger.info('=== Atom type and initial parameter assignment ===')
    logger.info('Method: {}'.format(args.forcefield))

    # Get RTF and PRM file names
    rtfFile, prmFile = args.rtf_prm if args.rtf_prm else None, None

    # Assing atom types and initial force field parameters
    _charge = mol.charge.copy()
    parameters, mol = fftype(mol, method=args.forcefield, rtfFile=rtfFile, prmFile=prmFile, netcharge=args.charge)
    assert np.all(mol.charge == _charge), 'fftype is meddling with charges!'

    logger.info('Atom types:')
    for name, type in zip(mol.name, mol.atomtype):
        logger.info('   {:4s} : {}'.format(name, type))

    # TODO write initial parameter to a file

    return mol, parameters


def _fit_initial_charges(mol, args):

    from htmd.charge import fitGasteigerCharges
    from htmd.parameterization.util import guessBondType

    logger.info('=== Initial atomic charge fitting ===')

    if args.charge_type == 'None':
        logger.info('Initial atomic charges are taken from {}'.format(args.filename))

    elif args.charge_type in ('Gasteiger', 'AM1-BCC', 'ESP'):
        if args.min_type == 'mm':
            logger.info('Method: Gasteiger')

            # TODO move to _prepare_molecule
            if np.any(mol.bondtype == "un"):
                logger.info('Guessing bond types')
                mol = guessBondType(mol)

            mol = fitGasteigerCharges(mol)

            # Print the initial charges
            logger.info('Initial atomic charges:')
            for name, charge in zip(mol.name, mol.charge):
                logger.info('   {:4s}: {:6.3f}'.format(name, charge))
            logger.info('Molecular charge: {:6.3f}'.format(np.sum(mol.charge)))

        elif args.min_type in ('None', 'qm'):
            logger.info('Initial atomic charges are not required')

        else:
            raise AssertionError()

    else:
        raise AssertionError()

    return mol


def _fit_charges(mol, args, qm):

    from htmd.charge import fitGasteigerCharges, fitChargesWithAntechamber, fitESPCharges, symmetrizeCharges
    from htmd.parameterization.util import guessBondType, getFixedChargeAtomIndices, getDipole, _qm_method_name
    from htmd.parameterization.detect import detectEquivalentAtoms

    logger.info('=== Atomic charge fitting ===')
    logger.info('Method: {}'.format(args.charge_type))

    if args.charge_type == 'None':

        # TODO move to argument validation
        if len(args.fix_charge) > 0:
            logger.warning('Flag --fix-charge does not have effect!')

        logger.info('Atomic charges are taken from {}'.format(args.filename))

    elif args.charge_type == 'Gasteiger':

        # TODO move to argument validation
        if len(args.fix_charge) > 0:
            logger.warning('Flag --fix-charge does not have effect!')

        # TODO move to _prepare_molecule
        if np.any(mol.bondtype == "un"):
            logger.info('Guessing bond types')
            mol = guessBondType(mol)

        mol = fitGasteigerCharges(mol)

        charge = int(round(np.sum(mol.charge)))
        if args.charge != charge:
            raise RuntimeError('Molecular charge is set to {}, but Gasteiger atomic charges add up to {}.'.format(
                args.charge, charge))

    elif args.charge_type == 'AM1-BCC':

        # TODO move to argument validation
        if len(args.fix_charge) > 0:
            logger.warning('Flag --fix-charge does not have effect!')

        mol = fitChargesWithAntechamber(mol, type='bcc', molCharge=args.charge)
        mol = symmetrizeCharges(mol)

    elif args.charge_type == 'ESP':

        # Detect equivalent atom groups
        logger.info('Equivalent atom groups:')
        atom_groups = [group for group in detectEquivalentAtoms(mol)[0] if len(group) > 1]
        for atom_group in atom_groups:
            logger.info('    {}'.format(', '.join(mol.name[list(atom_group)])))

        # Select the atoms with fixed charges
        fixed_atom_indices = getFixedChargeAtomIndices(mol, args.fix_charge)

        # Create an ESP directory
        espDir = os.path.join(args.outdir, "esp", _qm_method_name(qm))
        os.makedirs(espDir, exist_ok=True)

        charge = int(round(np.sum(mol.charge)))
        if args.charge != charge:
            logger.warning('Molecular charge is set to {}, but atomic charges add up to {}'
                           ''.format(args.charge, charge))
            if len(args.fix_charge) > 0:
                raise RuntimeError('Flag --fix-charge cannot be used when atomic charges are inconsistent with passed '
                                   'molecular charge {}'.format(args.charge))
            mol.charge[:] = args.charge/mol.numAtoms

        # Set random number generator seed
        if args.seed:
            np.random.seed(args.seed)

        # Fit ESP charges
        mol, extra = fitESPCharges(mol, qm, espDir, fixed=fixed_atom_indices)

        # Print QM dipole
        logger.info('QM dipole: {:6.3f} {:6.3f} {:6.3f}; total: {:6.3f}'.format(*extra['qm_dipole']))

    else:
        raise ValueError()

    # Print MM dipole
    mm_dipole = getDipole(mol)
    logger.info('MM dipole: {:6.3f} {:6.3f} {:6.3f}; total: {:6.3f}'.format(*mm_dipole))

    # Print the new charges
    logger.info('Atomic charges:')
    for name, charge in zip(mol.name, mol.charge):
        logger.info('   {:4s}: {:6.3f}'.format(name, charge))
    logger.info('Molecular charge: {:6.3f}'.format(np.sum(mol.charge)))

    return mol


def _printEnergies(molecule, parameters, filename):

    from htmd.ffevaluation.ffevaluate import FFEvaluate

    energies = FFEvaluate(molecule, parameters).calculateEnergies(molecule.coords[:, :, 0])

    string = '''
== Diagnostic Energies ==

Bond     : {BOND_ENERGY}
Angle    : {ANGLE_ENERGY}
Dihedral : {DIHEDRAL_ENERGY}
Improper : {IMPROPER_ENERGY}
Electro  : {ELEC_ENERGY}
VdW      : {VDW_ENERGY}

'''.format(BOND_ENERGY=energies['bond'],
           ANGLE_ENERGY=energies['angle'],
           DIHEDRAL_ENERGY=energies['dihedral'],
           IMPROPER_ENERGY=energies['improper'],
           ELEC_ENERGY=energies['elec'],
           VDW_ENERGY=energies['vdw'])

    with open(filename, 'w') as file_:
        file_.write(string)
    logger.info('Write energy file: {}'.format(filename))

    logger.info('Diagnostic energies:')
    for name in ('bond', 'angle', 'dihedral', 'improper', 'elec', 'vdw'):
        logger.info('   {:8s} : {:10.3f} kcal/mol'.format(name, energies[name]))


def _output_results(mol, parameters, original_coords, qm, args):

    from htmd.parameterization.util import _qm_method_name
    from htmd.parameterization.writers import writeParameters

    logger.info('=== Results ===')

    # Output the FF parameters and other files
    # TODO get rid of QM
    paramoutdir = os.path.join(args.outdir, 'parameters', args.forcefield, _qm_method_name(qm))
    # TODO split into separate writer
    writeParameters(paramoutdir, mol, parameters, args.forcefield, args.charge, original_coords=original_coords)

    # Write energy file
    energyFile = os.path.join(paramoutdir, 'energies.txt')
    _printEnergies(mol, parameters, energyFile)


def main_parameterize(arguments=None):

    from htmd.parameterization.parameterset import recreateParameters, createMultitermDihedralTypes, inventAtomTypes
    from htmd.parameterization.util import minimize, fitDihedrals, detectChiralCenters

    logger.info('===== Parameterize =====')

    # Parse arguments
    parser = getArgumentParser()
    args = parser.parse_args(args=arguments)

    # Configure loggers
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug(sys.argv[1:])

    # Deprecation warnings
    # TODO remove at some point of time
    if args.minimize is not parser.get_default('minimize'):
        raise DeprecationWarning('Use `--min-type` instead.')
    if args.optimize_dihedral is not parser.get_default('optimize_dihedral'):
        raise DeprecationWarning('Use `--scan-type` instead.')

    # Print arguments
    logger.info('=== Arguments ===')
    for key, value in vars(args).items():
        if key in ('fake_qm',):  # Hidden
            continue
        logger.info('{:>20s}: {:s}'.format(key, str(value)))        

    # Get a molecule and check its validity
    mol = _prepare_molecule(args)

    # Copy the molecule to preserve initial coordinates
    orig_coor = mol.coords.copy()

    # Select dihedral angles to parameterize
    selected_dihedrals = _select_dihedrals(mol, args)

    # Get a queue
    queue = _get_queue(args)

    # Get a QM calculator
    qm = _get_qm_calculator(args, queue)

    # Get a NNP calculator
    nnp = _get_nnp_calculator(args, queue)

    # Assign atom types and initial force field parameters
    mol, parameters = _get_initial_parameters(mol, args)

    # Assign initial atomic charges, if needed
    mol = _fit_initial_charges(mol, args)

    # Get a MM calculator
    # TODO refactor
    mm_minimizer = None
    if args.min_type == 'mm' or args.dihed_opt_type == 'mm':
        from htmd.qm.custom import OMMMinimizer
        mm_minimizer = OMMMinimizer(mol, parameters)

    # Geometry minimization
    # TODO refactor
    if args.min_type != 'None':
        logger.info(' === Geometry minimization ===')

        if args.min_type == 'mm':
            logger.info('Model: MM with the initial force field parameters')
        elif args.min_type == 'qm':
            logger.info('Model: reference method')
        else:
            raise ValueError()

        # Set parameters for the fake QM
        if args.fake_qm:
            qm._parameters = parameters

        # Detect chiral centers
        initial_chiral_centers = detectChiralCenters(mol)

        # Minimize molecule
        mol = minimize(mol, qm, args.outdir, min_type=args.min_type, mm_minimizer=mm_minimizer)

        # TODO print minimization status

        # Check if the chiral center hasn't changed during the minimization
        chiral_centers = detectChiralCenters(mol)
        if initial_chiral_centers != chiral_centers:
            raise RuntimeError('Chiral centers have changed during the minimization: '
                               '{} --> {}'.format(initial_chiral_centers, chiral_centers))

    # Fit charges
    mol = _fit_charges(mol, args, qm)

    # Scan dihedrals and fit parameters
    # TODO refactor
    if len(selected_dihedrals) > 0:

        logger.info('=== Dihedral angle scanning and parameter fitting ===')

        if args.dihed_opt_type == 'None':
            logger.info('Dihedral scanning: static')
        elif args.dihed_opt_type == 'mm':
            logger.info('Dihedral scanning: mimimized with MM (using the initial force field parameters)')
        elif args.dihed_opt_type == 'qm':
            logger.info('Dihedral scanning: mimimized with the reference method')
        else:
            raise ValueError()

        # Invent new atom types for dihedral atoms
        old_types = mol.atomtype
        mol, initial_types = inventAtomTypes(mol, selected_dihedrals)
        parameters = recreateParameters(mol, initial_types, parameters)
        parameters = createMultitermDihedralTypes(parameters)
        logger.info('Assign atom with new atom types:')
        for name, old_type, new_type in zip(mol.name, old_types, mol.atomtype):
            if old_type != new_type:
                logger.info('   {:4s} : {:6s} --> {:6s}'.format(name, old_type, new_type))

        # Set parameters for the fake QM
        if args.fake_qm:
            qm._parameters = parameters

        # Set random number generator seed
        if args.seed:
            np.random.seed(args.seed)

        # Fit the parameters
        # TODO separate scanning and fitting
        parameters = fitDihedrals(mol, qm, args.forcefield, parameters, selected_dihedrals, args.outdir,
                                  dihed_opt_type=args.dihed_opt_type, mm_minimizer=mm_minimizer,
                                  num_searches=args.dihed_num_searches)

    # Output the parameters and other results
    _output_results(mol, parameters, orig_coor, qm, args)


if __name__ == '__main__':

    arguments = sys.argv[1:] if len(sys.argv) > 1 else ['-h']
    main_parameterize(arguments=arguments)
