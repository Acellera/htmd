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
    parser.add_argument('-ff', '--forcefield', nargs='+', default=['GAFF2'], choices=fftypemethods,
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
                        help='Do not perform QM structure minimization')
    parser.add_argument('--charge-type', default='ESP', choices=['None', 'Gasteiger', 'ESP'],
                        help='Partial atomic charge type (default: %(default)s)')
    parser.add_argument('--no-dihed', action='store_false', dest='fit_dihedral',
                        help='Do not perform QM scanning of dihedral angles')
    parser.add_argument('--no-dihed-opt', action='store_false', dest='optimize_dihedral',
                        help='Do not perform QM structure optimisation when scanning dihedral angles')
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

    # QMML module name
    parser.add_argument('--qmml', help=argparse.SUPPRESS)

    return parser

def _prepare_molecule(args):

    from htmd.molecule.molecule import Molecule
    from htmd.parameterization.util import makeAtomNamesUnique, guessElements

    mol = Molecule(args.filename, guessNE=['bonds'], guess=[])

    # Check if the file contain just one conformation
    if mol.numFrames != 1:
        raise RuntimeError('{} has to contain only one molecule, but found {}'.format(args.filename, mol.numFrames))

    # Make atom names unique if needed
    if np.unique(mol.name).size != mol.numAtoms:
        logger.warning('Atom names in the molecule are not unique!')
        new_mol = makeAtomNamesUnique(mol)
        for i, (old_name, new_name) in enumerate(zip(mol.name, new_mol.name)):
            if old_name != new_name:
                logger.warning('Rename atom {:3d}: {:4s} --> {:4s}'.format(i, old_name, new_name))
        mol = new_mol

    # Guess elements
    # TODO: it should not depend on FF
    mol = guessElements(mol, args.forcefield[0])

    # Set segment ID
    # Note: it is need to write complete PDB files
    mol.segid[:] = 'L'

    # TODO: check charge

    # TODO: check bonds

    return mol


def printEnergies(molecule, parameters, filename):
    from htmd.ffevaluation.ffevaluate import FFEvaluate
    assert molecule.numFrames == 1
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

    sys.stdout.write(string)
    with open(filename, 'w') as file_:
        file_.write(string)


def printReport(mol, netcharge, equivalents, all_dihedrals):

    print('\n == Molecule report ==\n')

    print('Total number of atoms: %d' % mol.numAtoms)
    print('Total charge: %d' % netcharge)

    print('Equivalent atom groups:')
    for atom_group in equivalents[0]:
        if len(atom_group) > 1:
            print('  ' + ', '.join(mol.name[list(atom_group)]))

    print('Parameterizable dihedral angles:')
    for equivalent_dihedrals in all_dihedrals:
        dihedral, equivalent_dihedrals = equivalent_dihedrals[0], equivalent_dihedrals[1:]
        print('  ' + '-'.join(mol.name[list(dihedral)]))
        if equivalent_dihedrals:
            print('    Equivalents:')
            for dihedral in equivalent_dihedrals:
                print('      ' + '-'.join(mol.name[list(dihedral)]))


def _fit_charges(mol, args, qm):

    from htmd.charge import fitGasteigerCharges, fitESPCharges
    from htmd.parameterization.util import guessBondType, getFixedChargeAtomIndices, getDipole, _qm_method_name

    logger.info('=== Fitting atomic charges ===')

    if args.charge_type == 'None':

        if len(args.fix_charge) > 0:
            logger.warning('Flag --fix-charge does not have effect!')

        logger.info('Atomic charges are taken from {}'.format(args.filename))

    elif args.charge_type == 'Gasteiger':

        if len(args.fix_charge) > 0:
            logger.warning('Flag --fix-charge does not have effect!')

        if np.any(mol.bondtype == "un"):
            logger.info('Guessing bond types')
            mol = guessBondType(mol)

        mol = fitGasteigerCharges(mol)

        charge = int(round(np.sum(mol.charge)))
        if args.charge != charge:
            raise RuntimeError('Molecular charge is set to {}, but Gasteiger atomic charges add up to {}.'.format(
                args.charge, charge))

    elif args.charge_type == 'ESP':

        # Set random number generator seed
        if args.seed:
            np.random.seed(args.seed)

        # Select the atoms with fixed charges
        fixed_atom_indices = getFixedChargeAtomIndices(mol, args.fix_charge)

        # Create an ESP directory
        espDir = os.path.join(args.outdir, "esp", _qm_method_name(qm))
        os.makedirs(espDir, exist_ok=True)

        # Fit ESP charges
        mol, extra = fitESPCharges(mol, qm, espDir, fixed=fixed_atom_indices)
        logger.info('QM dipole: %f %f %f; %f' % tuple(extra['qm_dipole']))

    else:
        raise ValueError()

    mm_dipole = getDipole(mol)
    if np.all(np.isfinite(mm_dipole)):
        logger.info('MM dipole: %f %f %f; %f' % tuple(mm_dipole))
    else:  # TODO fix
        logger.warning('MM dipole cannot be computed. Check if elements are detected correctly.')

    # Print the new charges
    logger.info('Atomic charges:')
    for name, charge in zip(mol.name, mol.charge):
        logger.info('    {}: {:6.3f}'.format(name, charge))
    logger.info('Molecular charge: {:6.3f}'.format(np.sum(mol.charge)))

    return mol


def main_parameterize(arguments=None):

    args = getArgumentParser().parse_args(args=arguments)

    # Get a molecule and check its validity
    mol = _prepare_molecule(args)

    # Get RTF and PRM file names
    rtfFile, prmFile = None, None
    if args.rtf_prm:
        rtfFile, prmFile = args.rtf_prm

    # Create a queue for QM
    from htmd.queues.localqueue import LocalCPUQueue
    from htmd.queues.slurmqueue import SlurmQueue
    from htmd.queues.lsfqueue import LsfQueue
    from htmd.queues.pbsqueue import PBSQueue
    from htmd.queues.acecloudqueue import AceCloudQueue

    if args.queue == 'local':
        queue = LocalCPUQueue()
    elif args.queue == 'Slurm':
        queue = SlurmQueue(_configapp=args.code.lower())
    elif args.queue == 'LSF':
        queue = LsfQueue(_configapp=args.code.lower())
    elif args.queue == 'PBS':
        queue = PBSQueue()  # TODO: configure
    elif args.queue == 'AceCloud':
        queue = AceCloudQueue()  # TODO: configure
        queue.groupname = args.groupname
        queue.hashnames = True
    else:
        raise NotImplementedError

    # Override default ncpus
    if args.ncpus:
        logger.info('Overriding ncpus to {}'.format(args.ncpus))
        queue.ncpu = args.ncpus
    if args.memory:
        logger.info('Overriding memory to {}'.format(args.memory))
        queue.memory = args.memory

    # Create a QM object
    from htmd.qm import Psi4, Gaussian, FakeQM2

    if args.qmml:
        import importlib
        from htmd.qm.custom import CustomQM
        qm = CustomQM()
        qmml_module = importlib.import_module(args.qmml)
        logger.info('QMML module: {}'.format(qmml_module))
        qmml_calculator = qmml_module.get_calculator()
        logger.info('QMML calculator: {}'.format(qmml_calculator))
        qm.calculator = qmml_calculator
    else:
        if args.code == 'Psi4':
            qm = Psi4()
        elif args.code == 'Gaussian':
            qm = Gaussian()
        else:
            raise NotImplementedError

    # This is for debugging only!
    if args.fake_qm:
        qm = FakeQM2()
        logger.warning('Using FakeQM')

    # Start processing
    from htmd.parameterization.fftype import fftype
    from htmd.parameterization.util import getEquivalentsAndDihedrals, minimize, fitDihedrals, _qm_method_name
    from htmd.parameterization.parameterset import recreateParameters, createMultitermDihedralTypes, inventAtomTypes
    from htmd.parameterization.writers import writeParameters

    # Get rotatable dihedral angles
    mol, equivalents, all_dihedrals = getEquivalentsAndDihedrals(mol)

    if args.list:
        print('\n === Parameterizable dihedral angles of {} ===\n'.format(args.filename))
        with open('torsions.txt', 'w') as fh:
            for dihedral in all_dihedrals:
                dihedral_name = '-'.join(mol.name[list(dihedral[0])])
                print('  {}'.format(dihedral_name))
                fh.write(dihedral_name+'\n')
        print()
        sys.exit(0)

    # Print arguments
    print('\n === Arguments ===\n')
    for key, value in vars(args).items():
        print('{:>12s}: {:s}'.format(key, str(value)))

    # Get the molecular charge
    charge = int(round(np.sum(mol.charge)))
    if args.charge is None:
        args.charge = charge
        logger.info('Molecular charge is set to {} by adding up the atomic charges in {}'.format(args.charge,
                                                                                              args.filename))
    else:
        logger.info('Molecular charge is set to {}'.format(args.charge))
        if args.charge_type == 'None' and args.charge != charge:
            raise ValueError(
                'The molecular charge is set to {}, but the partial atomic charges in {} add up to {}'.format(
                    args.charge, args.filename, charge))

    # Select which dihedrals to fit
    parameterizable_dihedrals = [list(dih[0]) for dih in all_dihedrals]
    if len(args.dihedral) > 0:
        all_dihedral_names = ['-'.join(mol.name[list(dihedral[0])]) for dihedral in all_dihedrals]
        parameterizable_dihedrals = []
        for dihedral_name in args.dihedral:
            if dihedral_name not in all_dihedral_names:
                raise ValueError('%s is not recognized as a rotatable dihedral angle' % dihedral_name)
            parameterizable_dihedrals.append(list(all_dihedrals[all_dihedral_names.index(dihedral_name)][0]))

    # Set up the QM object
    qm.theory = args.theory
    qm.basis = args.basis
    qm.solvent = args.environment
    qm.queue = queue
    qm.charge = args.charge

    print('\n === Parameterizing %s ===\n' % args.filename)
    for method in args.forcefield:

        print(" === Fitting for %s ===\n" % method)
        printReport(mol, args.charge, equivalents, all_dihedrals)

        _charge = mol.charge.copy()
        parameters, mol = fftype(mol, method=method, rtfFile=rtfFile, prmFile=prmFile, netcharge=args.charge)
        assert np.all(mol.charge == _charge), 'fftype is meddling with charges!'

        if isinstance(qm, FakeQM2):
            qm._parameters = parameters

        # Copy the molecule to preserve initial coordinates
        orig_coor = mol.coords.copy()

        # Update B3LYP to B3LYP-D3
        # TODO: this is silent and not documented stuff
        if qm.theory == 'B3LYP':
            qm.correction = 'D3'

        # Update basis sets
        # TODO: this is silent and not documented stuff
        if args.charge < 0 and qm.solvent == 'vacuum':
            if qm.basis == '6-31G*':
                qm.basis = '6-31+G*'
            if qm.basis == '6-311G**':
                qm.basis = '6-311++G**'
            if qm.basis == 'cc-pVDZ':
                qm.basis = 'aug-cc-pVDZ'
            logger.info('Changing basis sets to %s' % qm.basis)

        # Minimize molecule
        if args.minimize:
            print('\n == Minimizing ==\n')
            mol = minimize(mol, qm, args.outdir)

        # Fit charges
        mol = _fit_charges(mol, args, qm)

        # Fit dihedral angle parameters
        if args.fit_dihedral:
            print('\n == Fitting dihedral angle parameters ==\n')

            # Set random number generator seed
            if args.seed:
                np.random.seed(args.seed)

            # Invent new atom types for dihedral atoms
            mol, originaltypes = inventAtomTypes(mol, parameterizable_dihedrals, equivalents)
            parameters = recreateParameters(mol, originaltypes, parameters)
            parameters = createMultitermDihedralTypes(parameters)
            if isinstance(qm, FakeQM2):
                qm._parameters = parameters

            # Fit the parameters
            fitDihedrals(mol, qm, method, parameters, all_dihedrals, parameterizable_dihedrals, args.outdir,
                         geomopt=args.optimize_dihedral)

        # Output the FF parameters
        print('\n == Writing results ==\n')
        writeParameters(mol, parameters, qm, method, args.charge, args.outdir, original_coords=orig_coor)

        # Write energy file
        energyFile = os.path.join(args.outdir, 'parameters', method, _qm_method_name(qm), 'energies.txt')
        printEnergies(mol, parameters, energyFile)
        logger.info('Write energy file: %s' % energyFile)


if __name__ == "__main__":

    args = sys.argv[1:] if len(sys.argv) > 1 else ['-h']
    main_parameterize(arguments=args)
