# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import argparse
import sys


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
    parser.add_argument("-b", "--basis", help="QM Basis Set (default: %(default)s)", choices=["6-31g-star", "cc-pVDZ"],
                        default="cc-pVDZ", dest="basis")
    parser.add_argument("--theory", help="QM Theory (default: %(default)s)", choices=["RHF", "B3LYP"],
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
    parser.add_argument("--qmcode", help="QM code (default: %(default)s)", choices=["Gaussian", "PSI4", "TeraChem"],
                        default="PSI4", dest="qmcode")
    parser.add_argument("--freeze-charge", metavar="A1",
                        help="Freeze the charge of the named atom (default: %(default)s)", action="append",
                        default=None, dest="freezeq")
    parser.add_argument("--no-geomopt", help="Do not perform QM geometry optimisation when fitting torsions (default: "
                                             "%(default)s)", action="store_false", dest="geomopt", default=True)
    return parser


def main_parameterize(arguments=None):
    args = cli_parser().parse_args(args=arguments)

    from htmd.parameterization.ffmolecule import FFMolecule, FFEvaluate
    from htmd.parameterization.fftype import FFTypeMethod
    from htmd.qm.qmcalculation import Theory, BasisSet, Execution, Code
    import numpy as np
    import math

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

    # Communicate the # of CPUs to use to the QM engine via environment variable
    os.environ['NCPUS'] = str(args.ncpus)

    filename = args.mol
    if not os.path.exists(filename):
        print("File {} not found. Please check that the file exists and that the path is correct.".format(filename))
        sys.exit(0)

    if args.qmcode == "Gaussian":
        code = Code.Gaussian
    elif args.qmcode == "PSI4":
        code = Code.PSI4
    elif args.qmcode == "TeraChem":
        code = Code.TeraChem
    else:
        print("Unknown QM code: {}".format(args.qmcode))
        sys.exit(1)

    if args.exec == "inline":
        execution = Execution.Inline
    elif args.exec == "LSF":
        execution = Execution.LSF
    elif args.exec == "Slurm":
        execution = Execution.Slurm
    else:
        print("Unknown execution mode: {}".format(args.exec))
        sys.exit(1)

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

    if args.basis == "6-31g-star":
        basis = BasisSet._6_31G_star
    elif args.basis == "cc-pVDZ":
        basis = BasisSet._cc_pVDZ
    else:
        print("Unknown basis {}".format(args.basis))
        sys.exit(1)

    if args.theory == "RHF":
        theory = Theory.RHF
    elif args.theory == "B3LYP":
        theory = Theory.B3LYP
    else:
        print("Unknown theory %s".format(args.theory))
        sys.exit(1)

    if args.vacuum:
        solvent = False
    else:
        solvent = True

    # Just list torsions?
    if args.list:
        print(" === Listing soft torsions of {} ===\n".format(filename))
        mol = FFMolecule(filename=filename, method=methods[0], netcharge=args.charge, rtf=args.rtf, prm=args.prm,
                         basis=basis, theory=theory, solvent=solvent, execution=execution, qmcode=code,
                         outdir=args.outdir)
        dihedrals = mol.getSoftTorsions()
        print("Detected soft torsions:")
        fh = open("torsions.txt", "w")
        for d in dihedrals:
            print("\t{}-{}-{}-{}".format(mol.name[d[0]], mol.name[d[1]], mol.name[d[2]], mol.name[d[3]]))
            print("{}-{}-{}-{}".format(mol.name[d[0]], mol.name[d[1]], mol.name[d[2]], mol.name[d[3]]), file=fh)
        fh.close()
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
                         basis=basis, theory=theory, solvent=solvent, execution=execution, qmcode=code,
                         outdir=args.outdir)
        dihedrals = mol.getSoftTorsions()
        mol_orig = mol.copy()

        if not args.nomin:
            print("\n == Minimizing ==\n")
            mol.minimize()

        sys.stdout.flush()
        if not args.noesp:
            print("\n == Charge fitting ==\n")

            # Select the atoms that are to have frozen charges in the fit
            fixq = []
            if args.freezeq:
                for i in args.freezeq:
                    found = False
                    for d in range(len(mol.name)):
                        if mol.name[d] == i:
                            ni = d
                            print("Fixing charge for atom %s to %f" % (i, mol.charge[ni]))
                            fixq.append(ni)
                            found = True
                    if not found:
                        raise ValueError(" No atom named %s (--freeze-charge)" % i)

            (score, qm_dipole, mm_dipole) = mol.fitCharges(fixed=fixq)

            rating = "GOOD"
            if score > 1:
                rating = "CHECK"
            if score > 10:
                rating = "BAD"

            print("Charge Chi^2 score : %f : %s" % (score, rating))
            print("QM Dipole   : %f %f %f ; %f" % (qm_dipole[0], qm_dipole[1], qm_dipole[2], qm_dipole[3]))
            print("MM Dipole   : %f %f %f ; %f" % (mm_dipole[0], mm_dipole[1], mm_dipole[2], mm_dipole[3]))
            d = 0.
            for i in range(3):
                x = qm_dipole[i] - mm_dipole[i]
                d = d + x * x

            rating = "GOOD"
            if score > 1:
                rating = "CHECK"
            print("Dipole Chi^2 score : %f : %s" % (d, rating))
            print("")

        sys.stdout.flush()
        # Iterative dihedral fitting
        if not args.notorsion:
            print("\n == Torsion fitting ==\n")

            scores = np.ones(len(dihedrals))
            converged = False
            iteration = 1
            ref_mm = dict()
            while not converged:
                rets = []

                print("\nIteration %d" % iteration)

                last_scores = scores
                scores = np.zeros(len(dihedrals))
                idx = 0
                for d in dihedrals:
                    name = "%s-%s-%s-%s" % (mol.name[d[0]], mol.name[d[1]], mol.name[d[2]], mol.name[d[3]])
                    if args.torsion == 'all' or name in args.torsion.split(','):
                        print("\n == Fitting torsion {} ==\n".format(name))
                        try:
                            ret = mol.fitSoftTorsion(d, geomopt=args.geomopt)
                            rets.append(ret)
                            if iteration == 1:
                                ref_mm[name] = ret
                            rating = "GOOD"
                            if ret.chisq > 10:
                                rating = "CHECK"
                            if ret.chisq > 100:
                                rating = "BAD"
                            print("Torsion %s Chi^2 score : %f : %s" % (name, ret.chisq, rating))
                            sys.stdout.flush()
                            scores[idx] = ret.chisq
                            # Always use the mm_orig from first iteration (unmodified)
                            ret.mm_original = ref_mm[name].mm_original
                            phi_original = ref_mm[name].phi
                            fn = mol.plotTorsionFit(ret, phi_original, show=False)
                        except Exception as e:
                            print("Error in fitting")
                            print(str(e))
                            raise
                            scores[idx] = 0.
                            pass
                            # print(fn)
                    idx += 1
                # print(scores)
                if iteration > 1:
                    converged = True
                    for j in range(len(scores)):
                        # Check convergence
                        try:
                            relerr = (scores[j] - last_scores[j]) / last_scores[j]
                        except:
                            relerr = 0.
                        if math.isnan(relerr):
                            relerr = 0.
                        convstr = "- converged"
                        if math.fabs(relerr) > 1.e-2:
                            convstr = ""
                            converged = False
                        print(" Dihedral %d relative error : %f %s" % (j, relerr, convstr))

                iteration += 1

            print(" Fitting converged at iteration %d" % (iteration - 1))
            if len(rets):
                fit = mol.plotConformerEnergies(rets, show=False)
                print("\n Fit of conformer energies: RMS %f Variance %f" % (fit[0], fit[1]))

        printEnergies(mol, 'energies.txt')

        # Output the ff parameters
        paramdir = os.path.join(args.outdir, "parameters", method.name, mol.output_directory_name())
        print("\n == Output to {} ==\n".format(paramdir))
        try:
            os.makedirs(paramdir, exist_ok=True)
        except:
            raise OSError('Directory {} could not be created. Check if you have permissions.'.format(paramdir))

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
    sys.exit(0)


if __name__ == "__main__":

    main_parameterize(arguments=['-h'])
    sys.exit(0)
