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
    parser.add_argument('--dihed-num-iterations', default=3, type=int,
                        help='Number of iterations during the dihedral parameter fitting')
    parser.add_argument('--dihed-fit-type', default='NRS', choices=['iterative', 'NRS'],
                        help='Dihedral fitting method. Can be either iterative or naive random search (NRS).')
    parser.add_argument('-q', '--queue', default='local', choices=['local', 'Slurm', 'LSF'],
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

    # PlayQueue arguments
    parser.add_argument('--pm-token', help=argparse.SUPPRESS)
    parser.add_argument('--max-jobs', type=int, default=sys.maxsize, help=argparse.SUPPRESS)

    # Debug mode
    parser.add_argument('--debug', action='store_true', default=False, dest='debug', help=argparse.SUPPRESS)

    return parser


def _printArguments(args, filename=None):

    if filename:
        logger.propagate = False  # Turn off logging to stdout
        fh = logging.FileHandler(filename, mode='w')
        logger.addHandler(fh)

    logger.info('=== Arguments ===')
    for key, value in sorted(vars(args).items()):
        if key in ('fake_qm', 'max_jobs', 'pm_token'):  # Hidden
            continue
        logger.info('{:>20s}: {:s}'.format(key, str(value)))

    if filename:
        logger.propagate = True  # Turn on logging to stdout
        logger.removeHandler(fh)


def _prepare_molecule(args):

    from moleculekit.molecule import Molecule
    from htmd.parameterization.fixes import fixPhosphateTypes
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

    # Fix atom and bond types
    mol = fixPhosphateTypes(mol)

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
    elif args.queue == 'PlayQueue':
        from htmd.queues.playqueue import PlayQueue
        queue = PlayQueue()  # TODO: configure
        queue.token = args.pm_token
        queue.app = 'Psi4'
        queue.max_jobs = args.max_jobs
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
    nnp.queue = queue
    nnp.charge = args.charge
    nnp.calculator = nnp_calculator

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


def _fit_initial_charges(mol, args, atom_types):

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

            mol = fitGasteigerCharges(mol, atom_types=atom_types)

            charge = int(round(np.sum(mol.charge)))
            if args.charge != charge:
                logger.warning(f'Molecular charge is {args.charge}, but Gasteiger atomic charges add up to {charge}!')
                args.charge = charge

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


def _fit_charges(mol, args, qm, atom_types):

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

        mol = fitGasteigerCharges(mol, atom_types=atom_types)

        charge = int(round(np.sum(mol.charge)))
        if args.charge != charge:
            logger.warning(f'Molecular charge is {args.charge}, but Gasteiger atomic charges add up to {charge}!')
            args.charge = charge

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

Bond     : {BOND_ENERGY:12.5f} kcal/mol
Angle    : {ANGLE_ENERGY:12.5f} kcal/mol
Dihedral : {DIHEDRAL_ENERGY:12.5f} kcal/mol
Improper : {IMPROPER_ENERGY:12.5f} kcal/mol
Electro  : {ELEC_ENERGY:12.5f} kcal/mol
VdW      : {VDW_ENERGY:12.5f} kcal/mol

'''.format(BOND_ENERGY=energies['bond'],
           ANGLE_ENERGY=energies['angle'],
           DIHEDRAL_ENERGY=energies['dihedral'],
           IMPROPER_ENERGY=energies['improper'],
           ELEC_ENERGY=energies['elec'],
           VDW_ENERGY=energies['vdw'])

    for line in string.split('\n'):
        logger.info(line)
    with open(filename, 'w') as file_:
        file_.write(string)
    logger.info('Write energy file: {}'.format(filename))


def _output_results(mol, parameters, original_coords, args):

    from htmd.parameterization.writers import writeParameters

    logger.info('=== Results ===')

    paramDir = os.path.join(args.outdir, 'parameters', args.forcefield)
    os.makedirs(paramDir, exist_ok=True)

    # Write arguments
    argumentsFile = os.path.join(paramDir, 'arguments.txt')
    _printArguments(args, filename=argumentsFile)
    logger.info('Write the list of  to {}'.format(argumentsFile))

    # Output the FF parameters and other files
    # TODO split into separate writer
    writeParameters(paramDir, mol, parameters, args.forcefield, args.charge, original_coords=original_coords)

    # Write energy file
    energyFile = os.path.join(paramDir, 'energies.txt')
    _printEnergies(mol, parameters, energyFile)


def main_parameterize(arguments=None, progress=None):

    from htmd.parameterization.parameterset import recreateParameters, createMultitermDihedralTypes, inventAtomTypes
    from htmd.parameterization.util import detectChiralCenters, scanDihedrals, filterQMResults, minimize

    progress = progress if callable(progress) else lambda x, **kwds: None

    logger.info('===== Parameterize =====')

    # Parse arguments
    parser = getArgumentParser()
    args = parser.parse_args(args=arguments)
    args.queue = 'PlayQueue' if args.pm_token else args.queue
    _printArguments(args)

    # Validate arguments
    if args.fake_qm and args.nnp:
        raise ValueError('FakeQM and NNP are not compatible')

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

    # Get a molecule and check its validity
    progress('Preparing the molecule')
    mol = _prepare_molecule(args)

    # Preserve the initial molecule
    initial_mol = mol.copy()

    # Select dihedral angles to parameterize
    progress('Detecting dihedral angles')
    selected_dihedrals = _select_dihedrals(mol, args)

    # Get a queue
    queue = _get_queue(args)

    # Get QM calculators
    qm_calculator = None
    qm_name = ''
    if args.min_type == 'qm' or args.charge_type == 'ESP' or (args.fit_dihedral and not args.nnp):
        qm_calculator = _get_qm_calculator(args, queue)
        qm_name = '{}/{}'.format(qm_calculator.theory, qm_calculator.basis)
        if args.fake_qm:
            qm_name = 'Fake QM'

    # Get NNP calculators
    nnp_calculator = None
    nnp_name = ''
    if args.nnp:
        nnp_calculator = _get_nnp_calculator(args, queue)
        nnp_name = args.nnp

    # Set the reference calculator
    if args.nnp:
        ref_calculator = nnp_calculator
        ref_name = nnp_name
    else:
        ref_calculator = qm_calculator
        ref_name = qm_name
    logger.info('Reference method: {}'.format(ref_name))

    # Assign atom types and initial force field parameters
    progress('Assigning atom types and initial parameters')
    mol, parameters = _get_initial_parameters(mol, args)

    # Assign initial atomic charges, if needed
    progress('Assigning initial atomic charges')
    mol = _fit_initial_charges(mol, args, initial_mol.atomtype)

    # Geometry minimization
    # TODO refactor
    if args.min_type != 'None':
        progress('Optimizing geometry')
        logger.info(' === Geometry minimization ===')

        if args.min_type == 'mm':
            logger.info('Model: MM with the initial force field parameters')
        elif args.min_type == 'qm':
            logger.info('Model: the reference method')
        else:
            raise ValueError()

        # Set parameters for the fake QM
        if args.fake_qm:
            assert not args.nnp
            ref_calculator._parameters = parameters

        # Detect chiral centers
        initial_chiral_centers = detectChiralCenters(mol, atom_types=initial_mol.atomtype)

        mm_minimizer = None
        if args.min_type == 'mm':
            from htmd.qm.custom import OMMMinimizer
            mm_minimizer = OMMMinimizer(mol, parameters)

        # Minimize molecule
        mol = minimize(mol, ref_calculator, args.outdir, min_type=args.min_type, mm_minimizer=mm_minimizer)

        # TODO print minimization status
        # Check if the chiral center hasn't changed during the minimization
        chiral_centers = detectChiralCenters(mol, atom_types=initial_mol.atomtype)
        if initial_chiral_centers != chiral_centers:
            raise RuntimeError('Chiral centers have changed during the minimization: '
                               '{} --> {}'.format(initial_chiral_centers, chiral_centers))

    # Fit charges
    progress('Assigning atomic charges')
    mol = _fit_charges(mol, args, qm_calculator, initial_mol.atomtype)

    # Scan dihedrals and fit parameters
    # TODO refactor
    if len(selected_dihedrals) > 0:

        from htmd.parameterization.dihedral import DihedralFitting  # Slow import

        num_jobs = 0 if args.nnp else 36 * len(selected_dihedrals)
        progress('Scanning dihedral angles', num_jobs=num_jobs)
        logger.info('=== Dihedral angle scanning ===')

        if args.dihed_opt_type == 'None':
            logger.info('Dihedral scanning: static')
        elif args.dihed_opt_type == 'mm':
            logger.info('Dihedral scanning: minimized with MM (using the initial force field parameters)')
        elif args.dihed_opt_type == 'qm':
            logger.info('Dihedral scanning: minimized with the reference method')
        else:
            raise ValueError()

        # Recreate MM minimizer now that charges have been fitted
        mm_minimizer = None
        if args.dihed_opt_type == 'mm':
            from htmd.qm.custom import OMMMinimizer
            mm_minimizer = OMMMinimizer(mol, parameters)

        # Set parameters for the fake QM
        if args.fake_qm:
            assert not args.nnp
            ref_calculator._parameters = parameters

        # Scan dihedral angles
        scan_results = scanDihedrals(mol, ref_calculator, selected_dihedrals, args.outdir,
                                     scan_type=args.dihed_opt_type, mm_minimizer=mm_minimizer)

        # Filter scan results
        scan_results = filterQMResults(scan_results, mol=initial_mol)
        logger.info('Valid rotamers:')
        for idihed, (dihedral, results) in enumerate(zip(selected_dihedrals, scan_results)):
            dihed_name = '-'.join(mol.name[list(dihedral)])
            logger.info('  {:2d}: {}: {}'.format(idihed, dihed_name, len(results)))

            if len(results) < 13:
                raise RuntimeError('Less than 13 valid rotamers for {} dihedral. '
                                   'Not enough for fitting!'.format(dihed_name))

        logger.info('=== Dihedral parameter fitting ===')

        # Invent new atom types for dihedral atoms
        progress('Creating new atom types')
        old_types = mol.atomtype
        mol, initial_types = inventAtomTypes(mol, selected_dihedrals)
        parameters = recreateParameters(mol, initial_types, parameters)
        parameters = createMultitermDihedralTypes(parameters)
        logger.info('Assign atom with new atom types:')
        for name, old_type, new_type in zip(mol.name, old_types, mol.atomtype):
            if old_type != new_type:
                logger.info('   {:4s} : {:6s} --> {:6s}'.format(name, old_type, new_type))

        # Set random number generator seed
        if args.seed:
            np.random.seed(args.seed)

        # Fit the dihedral parameters
        progress('Fitting dihedral angle parameters')
        df = DihedralFitting()
        df.parameters = parameters
        df.molecule = mol
        df.dihedrals = selected_dihedrals
        df.qm_results = scan_results
        df.num_iterations = args.dihed_num_iterations
        df.fit_type = args.dihed_fit_type
        df.result_directory = os.path.join(args.outdir, 'parameters', args.forcefield)

        # In case of FakeQM, the initial parameters are set to zeros.
        # It prevents DihedralFitting class from cheating :D
        if args.fake_qm:
            df.zeroed_parameters = True

        # Fit dihedral parameters
        parameters = df.run()

        # Plot dihedral profiles
        plot_dir = os.path.join(args.outdir, 'parameters', args.forcefield, 'plots')
        os.makedirs(plot_dir, exist_ok=True)
        df.plotConformerEnergies(plot_dir, ref_name=ref_name)
        for idihed in range(len(df.dihedrals)):
            df.plotDihedralEnergies(idihed, plot_dir, ref_name=ref_name)

    # Output the parameters and other results
    progress('Outputing results')
    _output_results(mol, parameters, initial_mol.coords, args)


if __name__ == '__main__':

    arguments = sys.argv[1:] if len(sys.argv) > 1 else ['-h']
    main_parameterize(arguments)
