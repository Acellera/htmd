# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
import argparse


def cli_parser():
    ncpus = os.cpu_count()
    try:
        ncpus = int(os.getenv("NCPUS"))
    except:
        pass

    parser = argparse.ArgumentParser(description="Acellera Small Molecule Parameterization Tool")
    parser.add_argument("-m", "--mol2", help="Molecule to parameterise, in mol2 format", required=True, type=str,
                        default=None, metavar="<input.mol2>", action="store", dest="mol")
    parser.add_argument("-l", "--list", "--list-torsions", help="List parameterisable torsions", action="store_true",
                        default=False,
                        dest="list")
    parser.add_argument("-c", "--charge", help="Net charge on molecule (default: sum of the partial charges on the "
                                               ".mol2 file)", type=int, default=None, action="store", dest="charge")
    parser.add_argument("--rtf", help="Inital RTF parameters (req --prm)", type=str, default=None, dest="rtf")
    parser.add_argument("--prm", help="Inital PRM parameters (req --rtf)", type=str, default=None, dest="prm")
    parser.add_argument("-o", "--outdir", help="Output directory (default: %(default)s)", type=str, default="./",
                        dest="outdir")
    parser.add_argument("-t", "--torsion", metavar="A1-A2-A3-A4", help="Torsion to parameterise (default: %(default)s)",
                        default="all", dest="torsion")
    parser.add_argument("-n", "--ncpus", help="Number of CPUs to use (default: %(default)s)", default=ncpus,
                        dest="ncpus")
    parser.add_argument("-f", "--forcefield", help="Inital FF guess to use (default: %(default)s)",
                        choices=["GAFF", "GAFF2", "CGENFF", "all"], default="all")
    parser.add_argument("-b", "--basis", help="QM Basis Set (default: %(default)s)", choices=["6-31G*", "cc-pVDZ"],
                        default="cc-pVDZ", dest="basis")
    parser.add_argument("--theory", help="QM Theory (default: %(default)s)", choices=["HF", "B3LYP"],
                        default="B3LYP", dest="theory")
    parser.add_argument("--vacuum", help="Perform QM calculations in vacuum (default: %(default)s)",
                        action="store_false", dest="vacuum",
                        default=True)
    parser.add_argument("--no-min", help="Do not perform QM minimisation (default: %(default)s)", action="store_true",
                        dest="nomin", default=False)
    parser.add_argument("--no-esp", help="Do not perform QM charge fitting (default: %(default)s)", action="store_true",
                        dest="noesp", default=False)
    parser.add_argument("--no-torsions", help="Do not perform torsion fitting (default: %(default)s)",
                        action="store_true", dest="notorsion", default=False)
    parser.add_argument("-e", "--exec", help="Mode of execution for the QM calculations (default: %(default)s)",
                        choices=["inline", "LSF", "Slurm"], default="inline", dest="exec")
    parser.add_argument("--qmcode", help="QM code (default: %(default)s)", choices=["Psi4", "Gaussian"],
                        default="Psi4", dest="qmcode")
    parser.add_argument("--freeze-charge", metavar="A1",
                        help="Freeze the charge of the named atom (default: %(default)s)", action="append",
                        default=None, dest="freezeq")
    parser.add_argument("--no-geomopt", help="Do not perform QM geometry optimisation when fitting torsions (default: "
                                             "%(default)s)", action="store_false", dest="geomopt", default=True)
    parser.add_argument("--seed", help="Random number generator seed (default: %(default)s)", type=int,
                        default=20170920, dest="seed")
    return parser


def main_parameterize(arguments=None):

    import numpy as np

    from htmd.parameterization.ffmolecule import FFMolecule, FFEvaluate
    from htmd.parameterization.fftype import FFTypeMethod
    from htmd.queues.localqueue import LocalCPUQueue
    from htmd.queues.lsfqueue import LsfQueue
    from htmd.queues.slurmqueue import SlurmQueue
    from htmd.qm import Psi4, Gaussian

    args = cli_parser().parse_args(args=arguments)

    def printEnergies(m, filename):
        fener = open(filename, "w")
        ffe = FFEvaluate(m)
        energies = ffe.run(m.coords[:, :, 0])
        for out in sys.stdout, fener:
            print('''
== Diagnostic Energies ==

Bond     : {BOND_ENERGY}
Angle    : {ANGLE_ENERGY}
Dihedral : {DIHEDRAL_ENERGY}
Improper : {IMPROPER_ENERGY}
Electro  : {ELEC_ENERGY}
VdW      : {VDW_ENERGY}

'''.format(BOND_ENERGY=energies['bond'], ANGLE_ENERGY=energies['angle'], DIHEDRAL_ENERGY=energies['dihedral'],
           IMPROPER_ENERGY=energies['improper'], ELEC_ENERGY=energies['elec'], VDW_ENERGY=energies['vdw']),
                  file=out)
        fener.close()

    filename = args.mol
    if not os.path.exists(filename):
        print("File {} not found. Please check that the file exists and that the path is correct.".format(filename))
        sys.exit(0)

    # Create a queue for QM
    if args.exec == 'inline':
        queue = LocalCPUQueue()
    elif args.exec == 'Slurm':
        queue = SlurmQueue()
        queue.partition = SlurmQueue._defaults['cpu_partition']
    elif args.exec == 'LSF':
        queue = LsfQueue()
        queue.queue = LsfQueue._defaults['cpu_queue']
    else:
        raise NotImplementedError

    # Communicate the # of CPUs to use to the QM engine via environment variable
    # TODO this should be set up via a queue
    os.environ['NCPUS'] = str(args.ncpus)

    # Create a QM object
    if args.qmcode == 'Psi4':
        qm = Psi4()
    elif args.qmcode == 'Gaussian':
        qm = Gaussian()
    else:
        raise NotImplementedError

    # Set up the QM object
    qm.theory = args.theory
    qm.basis = args.basis
    qm.solvent = 'vacuum' if args.vacuum else 'PCM'
    qm.queue = queue

    if args.forcefield == "CGENFF":
        methods = [FFTypeMethod.CGenFF_2b6]
    elif args.forcefield == "GAFF":
        methods = [FFTypeMethod.GAFF]
    elif args.forcefield == "GAFF2":
        methods = [FFTypeMethod.GAFF2]
    elif args.forcefield == "all":
        methods = [FFTypeMethod.CGenFF_2b6, FFTypeMethod.GAFF2]
    else:
        print("Unknown initial guess force-field: {}".format(args.forcefield))
        sys.exit(1)

    # Just list torsions?
    if args.list:
        print(" === Listing soft torsions of {} ===\n".format(filename))
        mol = FFMolecule(filename=filename, method=methods[0], netcharge=args.charge, rtf=args.rtf, prm=args.prm,
                         qm=qm, outdir=args.outdir)
        print('Detected soft torsions:')
        with open('torsions.txt', 'w') as fh:
            for dihedral in mol.getSoftTorsions():
                name = "%s-%s-%s-%s" % tuple(mol.name[dihedral])
                print('\t'+name)
                fh.write(name+'\n')
        sys.exit(0)

    # Small report
    print(" === List of arguments used ===\n")
    for i in vars(args):
        print('{:>10s}:  {:<10s}'.format(i, str(vars(args)[i])))

    print("\n === Parameterizing {} ===\n".format(filename))
    for method in methods:
        sys.stdout.flush()
        print(" === Fitting for FF %s ===\n" % method.name)

        mol = FFMolecule(filename=filename, method=method, netcharge=args.charge, rtf=args.rtf, prm=args.prm,
                         qm=qm, outdir=args.outdir)

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

        mol_orig = mol.copy()

        if not args.nomin:
            print("\n == Minimizing ==\n")
            mol.minimize()

        sys.stdout.flush()
        if not args.noesp:
            print("\n == Charge fitting ==\n")

            # Set random number generator seed
            if args.seed:
                np.random.seed(args.seed)

            # Select the atoms that are to have frozen charges in the fit
            fixed_charges = []
            if args.freezeq:
                for i in args.freezeq:
                    found = False
                    for d in range(len(mol.name)):
                        if mol.name[d] == i:
                            ni = d
                            print("Fixing charge for atom %s to %f" % (i, mol.charge[ni]))
                            fixed_charges.append(ni)
                            found = True
                    if not found:
                        raise ValueError(" No atom named %s (--freeze-charge)" % i)

            # Fit ESP charges
            score, qm_dipole = mol.fitCharges(fixed=fixed_charges)

            rating = "GOOD"
            if score > 1:
                rating = "CHECK"
            if score > 10:
                rating = "BAD"
            print("Charge Chi^2 score : %f : %s\n" % (score, rating))
            print("QM Dipole   : %f %f %f ; %f" % tuple(qm_dipole))

            # Compare dipoles
            mm_dipole = mol.getDipole()
            score = np.sum((qm_dipole[:3] - mm_dipole[:3])**2)

            rating = "GOOD"
            if score > 1:
                rating = "CHECK"
            print("MM Dipole   : %f %f %f ; %f" % tuple(mm_dipole))
            print("Dipole Chi^2 score : %f : %s" % (score, rating))
            print("")

        sys.stdout.flush()
        # Iterative dihedral fitting
        if not args.notorsion:
            print("\n == Torsion fitting ==\n")

            # Set random number generator seed
            if args.seed:
                np.random.seed(args.seed)

            all_dihedrals = mol.getSoftTorsions()

            # Choose which dihedrals to fit
            dihedrals = []
            if args.torsion == 'all':
                dihedrals = all_dihedrals
            else:
                all_names = ['%s-%s-%s-%s' % tuple(mol.name[dihedral]) for dihedral in all_dihedrals]
                for name in args.torsion.split(','):
                    if name in all_names:
                        dihedrals.append(all_dihedrals[all_names.index(name)])
                    else:
                        raise ValueError('%s is not recognized as a soft torsion\n' % name)

            mol.fitSoftTorsions(dihedrals, args.geomopt)

        printEnergies(mol, 'energies.txt')

        # Output the FF parameters
        paramdir = os.path.join(args.outdir, 'parameters', method.name, mol.output_directory_name())
        os.makedirs(paramdir, exist_ok=True)
        print("\n == Output to {} ==\n".format(paramdir))

        if method.name == "CGenFF_2b6":
            try:
                mol._rtf.write(os.path.join(paramdir, "mol.rtf"))
                mol._prm.write(os.path.join(paramdir, "mol.prm"))
                for ext in ['psf', 'xyz', 'coor', 'mol2', 'pdb']:
                    mol.write(os.path.join(paramdir, "mol." + ext))
                mol_orig.write(os.path.join(paramdir, "mol-orig.mol2"))
                f = open(os.path.join(paramdir, "input.namd"), "w")
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
            except ValueError as e:
                print("Not writing CHARMM PRM: {}".format(str(e)))

        elif method.name == "GAFF" or method.name == "GAFF2":
            try:
                # types need to be remapped because Amber FRCMOD format limits the type to characters
                # writeFrcmod does this on the fly and returns a mapping that needs to be applied to the mol
                typemap = mol._prm.writeFrcmod(mol._rtf, os.path.join(paramdir, "mol.frcmod"))
                for ext in ['coor', 'mol2', 'pdb']:
                    mol.write(os.path.join(paramdir, "mol." + ext), typemap=typemap)
                mol_orig.write(os.path.join(paramdir, "mol-orig.mol2"), typemap=typemap)
                f = open(os.path.join(paramdir, "tleap.in"), "w")
                tmp = '''loadAmberParams mol.frcmod
A = loadMol2 mol.mol2
saveAmberParm A structure.prmtop mol.crd
quit'''
                print(tmp, file=f)
                f.close()
                f = open(os.path.join(paramdir, "input.namd"), "w")
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
            except ValueError as e:
                print("Not writing Amber FRCMOD: {}".format(str(e)))


if __name__ == "__main__":

    main_parameterize(arguments=['-h'])
    sys.exit(0)
