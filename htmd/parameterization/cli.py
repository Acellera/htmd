# (c) 2015-2017 Acellera Ltd http://www.acellera.com
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
from htmd.queues.localqueue import LocalCPUQueue
from htmd.queues.slurmqueue import SlurmQueue
from htmd.queues.lsfqueue import LsfQueue
from htmd.queues.pbsqueue import PBSQueue
from htmd.queues.acecloudqueue import AceCloudQueue
from htmd.qm import Psi4, Gaussian, FakeQM2
from htmd.parameterization.ffmolecule import FFMolecule
from htmd.parameterization.fftype import FFTypeMethod
from htmd.ffevaluation.ffevaluate import FFEvaluate

logger = logging.getLogger(__name__)


def getArgumentParser():

    parser = argparse.ArgumentParser(description='Acellera small molecule parameterization tool')

    parser.add_argument('filename', help='Molecule file in MOL2 format')
    parser.add_argument('-c', '--charge', type=int,
                        help='Total charge of the molecule (default: sum of partial charges)')
    parser.add_argument('-l', '--list', action='store_true', help='List parameterizable dihedral angles')
    parser.add_argument('--rtf-prm', nargs=2, metavar='<filename>', help='CHARMM RTF and PRM files')
    parser.add_argument('-ff', '--forcefield', nargs='+', default=['GAFF2'], choices=['GAFF', 'GAFF2', 'CGENFF'],
                        help='Inital force field guess (default: %(default)s)')
    parser.add_argument('--fix-charge', nargs='+', default=[], metavar='<atom name>',
                        help='Fix atomic charge during charge fitting (default: none)')
    parser.add_argument('-d', '--dihedral', nargs='+', default=[], metavar='A1-A2-A3-A4',
                        help='Select dihedral angle to parameterize (default: all parameterizable dihedral angles)')
    parser.add_argument('--code', default='Psi4', choices=['Psi4', 'Gaussian'], help='QM code (default: %(default)s)')
    parser.add_argument('--theory', default='B3LYP', choices=['HF', 'B3LYP'],
                        help='QM level of theory (default: %(default)s)')
    parser.add_argument('--basis', default='cc-pVDZ', choices=['6-31G*', '6-31+G*', 'cc-pVDZ', 'aug-cc-pVDZ'],
                        help='QM basis set (default: %(default)s)')
    parser.add_argument('--environment', default='vacuum', choices=['vacuum', 'PCM'],
                        help='QM environment (default: %(default)s)')
    parser.add_argument('--no-min', action='store_false', dest='minimize',
                        help='Do not perform QM structure minimization')
    parser.add_argument('--no-esp', action='store_false', dest='fit_charges', help='Do not perform QM charge fitting')
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

    return parser


def printEnergies(molecule, filename):

    assert molecule.numFrames == 1
    energies = FFEvaluate(molecule).run(molecule.coords[:, :, 0])

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


def main_parameterize(arguments=None):

    args = getArgumentParser().parse_args(args=arguments)

    if not os.path.exists(args.filename):
        raise ValueError('File %s cannot be found' % args.filename)

    method_map = {'GAFF': FFTypeMethod.GAFF, 'GAFF2': FFTypeMethod.GAFF2, 'CGENFF': FFTypeMethod.CGenFF_2b6}
    methods = [method_map[method] for method in args.forcefield]  # TODO: move into FFMolecule

    # Get RTF and PRM file names
    rtf, prm = None, None
    if args.rtf_prm:
        rtf, prm = args.rtf_prm

    # Create a queue for QM
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

    if hasattr(queue, 'environment'): # TODO: LocalCPUQueue does not have it
        queue.environment = 'PATH,LD_LIBRARY_PATH,PSI_SCRATCH' # Use Psi4 from an active conda environment

    # Override default ncpus
    if args.ncpus:
        logger.info('Overriding ncpus to {}'.format(args.ncpus))
        queue.ncpu = args.ncpus
    if args.memory:
        logger.info('Overriding memory to {}'.format(args.memory))
        queue.memory = args.memory

    # Create a QM object
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

    # Set up the QM object
    qm.theory = args.theory
    qm.basis = args.basis
    qm.solvent = args.environment
    qm.queue = queue

    # List rotatable dihedral angles
    if args.list:

        mol = FFMolecule(args.filename, method=methods[0], netcharge=args.charge, rtf=rtf, prm=prm, qm=qm,
                         outdir=args.outdir)
        print('\n === Parameterizable dihedral angles of %s ===\n' % args.filename)
        with open('torsions.txt', 'w') as fh:
            for dihedral in mol.getRotatableDihedrals():
                dihedral_name = '-'.join(mol.name[dihedral])
                print('  '+dihedral_name)
                fh.write(dihedral_name+'\n')
        print()
        sys.exit(0)

    # Print arguments
    print('\n === Arguments ===\n')
    for key, value in vars(args).items():
        print('{:>12s}: {:s}'.format(key, str(value)))

    print('\n === Parameterizing %s ===\n' % args.filename)
    for method in methods:
        print(" === Fitting for %s ===\n" % method.name)

        # Create the molecule
        mol = FFMolecule(args.filename, method=method, netcharge=args.charge, rtf=rtf, prm=prm, qm=qm,
                         outdir=args.outdir)
        mol.printReport()

        # Copy the molecule to preserve initial coordinates
        mol_orig = mol.copy()

        # Update B3LYP to B3LYP-D3
        # TODO: this is silent and not documented stuff
        if qm.theory == 'B3LYP':
            qm.correction = 'D3'

        # Update basis sets
        # TODO: this is silent and not documented stuff
        if mol.netcharge < 0 and qm.solvent == 'vacuum':
            if qm.basis == '6-31G*':
                qm.basis = '6-31+G*'
            if qm.basis == 'cc-pVDZ':
                qm.basis = 'aug-cc-pVDZ'
            logger.info('Changing basis sets to %s' % qm.basis)

        # Minimize molecule
        if args.minimize:
            print('\n == Minimizing ==\n')
            mol.minimize()

        # Fit ESP charges
        if args.fit_charges:
            print('\n == Fitting ESP charges ==\n')

            # Set random number generator seed
            if args.seed:
                np.random.seed(args.seed)

            # Select the atoms with fixed charges
            fixed_atom_indices = []
            for fixed_atom_name in args.fix_charge:

                if fixed_atom_name not in mol.name:
                    raise ValueError('Atom %s is not found. Check --fix-charge arguments' % fixed_atom_name)

                for aton_index in range(mol.numAtoms):
                    if mol.name[aton_index] == fixed_atom_name:
                        fixed_atom_indices.append(aton_index)
                        logger.info('Charge of atom %s is fixed to %f' % (fixed_atom_name, mol.charge[aton_index]))

            # Fit ESP charges
            _, qm_dipole = mol.fitCharges(fixed=fixed_atom_indices)

            # Copy the new charges to the original molecule
            mol_orig.charge[:] = mol.charge

            # Print dipoles
            logger.info('QM dipole: %f %f %f; %f' % tuple(qm_dipole))
            mm_dipole = mol.getDipole()
            if np.all(np.isfinite(mm_dipole)):
                logger.info('MM dipole: %f %f %f; %f' % tuple(mm_dipole))
            else:
                logger.warning('MM dipole cannot be computed. Check if elements are detected correctly.')

        # Fit dihedral angle parameters
        if args.fit_dihedral:
            print('\n == Fitting dihedral angle parameters ==\n')

            # Set random number generator seed
            if args.seed:
                np.random.seed(args.seed)

            # Get all rotatable dihedrals
            all_dihedrals = mol.getRotatableDihedrals()

            # Choose which dihedrals to fit
            dihedrals = []
            all_dihedral_names = ['-'.join(mol.name[dihedral]) for dihedral in all_dihedrals]
            for dihedral_name in args.dihedral:
                if dihedral_name not in all_dihedral_names:
                    raise ValueError('%s is not recognized as a rotatable dihedral angle' % dihedral_name)
                dihedrals.append(all_dihedrals[all_dihedral_names.index(dihedral_name)])
            dihedrals = dihedrals if len(dihedrals) > 0 else all_dihedrals  # Set default to all dihedral angles

            # Fit the parameters
            mol.fitDihedrals(dihedrals, args.optimize_dihedral)

        # Output the FF parameters
        print('\n == Writing results ==\n')
        mol.writeParameters(mol_orig)

        # Write energy file
        energyFile = os.path.join(mol.outdir, 'parameters', method.name, mol.qm_method_name(), 'energies.txt')
        printEnergies(mol, energyFile)
        logger.info('Write energy file: %s' % energyFile)

        
if __name__ == "__main__":

    args = sys.argv[1:] if len(sys.argv) > 1 else ['-h']
    main_parameterize(arguments=args)
