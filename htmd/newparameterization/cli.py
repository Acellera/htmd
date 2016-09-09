# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import matplotlib

matplotlib.use('Agg')

import argparse

from htmd.newparameterization.ffmolecule import FFMolecule, FFEvaluate
from htmd.newparameterization.fftype import FFTypeMethod
from htmd.qm.qmcalculation import BasisSet
import sys
import os

def printEnergies( mol ):
  print(" == Diagnostic Energies == ")
  ffe = FFEvaluate( mol )
  energies = ffe.evaluate( mol.coords[:,:,0] )
  print( "" )
  print( " Bond     : %f" % ( energies['bond'] ) )
  print( " Angle    : %f" % ( energies['angle'] ) )
  print( " Dihedral : %f" % ( energies['dihedral'] ) )
  print( " Improper : %f" % ( energies['improper'] ) )
  print( " Electro  : %f" % ( energies['elec'] ) )
  print( " VdW      : %f" % ( energies['vdw'] ) )
  print( "" )
 

def main_parameterize():
    ncpus = os.cpu_count()
    try:
        ncpus = int(os.getenv("NCPUS"))
    except:
        pass

    parser = argparse.ArgumentParser(description="Acellera Small Molecule Parameterisation Version 2.0")
    parser.add_argument("-m", "--mol2", help="Molecule to parameterise, in mol2 format", type=str, default="input.mol2",
                        action="store", dest="mol")
    parser.add_argument("-c", "--charge", help="Net charge on molecule", type=int, default=None, action="store",
                        dest="charge")
    parser.add_argument("--rtf", help="Inital RTF parameters (req --prm)", type=str, default=None, dest="rtf")
    parser.add_argument("--prm", help="Inital PRM parameters (req --rtf)", type=str, default=None, dest="prm")
    parser.add_argument("-o", "--output", metavar="DIRECTORY", help="Output directory", type=str, default="parameters")
    parser.add_argument("-l", "--list", help="List parameterisable torsions", action="store_true", default=False,
                        dest="list")
    parser.add_argument("-t", "--torsion", metavar="A1-A2-A3-A4", help="Torsion to parameterise (default all)",
                        action="append", default=None, dest="torsion")
    parser.add_argument("-n", "--ncpus", help="Number of CPUs to use (default %d)" % (ncpus), default=ncpus, dest="ncpus")
    parser.add_argument( "-f", "--forcefield", help="Inital FF guess to use", choices=[ "GAFF", "GAFF2", "CGENFF"], default="CGENFF" )
    parser.add_argument ( "-b", "--basis", help="QM Basis Set", choices=[ "6-31g-star", "cc-pVTZ" ], default="6-31g-star", dest="basis") 

    args = parser.parse_args()

    # Communicate the # of CPUs to use to the QM engine via environment variable
    os.putenv("NCPUS", str(args.ncpus))

    filename = args.mol

    if not os.path.exists(filename):
        print("file %s not found" % filename)
        sys.exit(0)

    # Make outdir directory [outdir]/[ff-name]
    outdir = os.path.join( args.output, args.forcefield, args.basis )
    try:
        os.mkdir(args.output)
    except: 
        pass
    try:
        os.mkdir( os.path.join( args.output, args.forcefield ) )
    except:
        pass
    try:
        os.mkdir( os.path.join( args.output, args.forcefield, args.basis ) )
    except:
        pass

    basis = BasisSet._6_31G_star

    if args.basis == "6-31g-star"  : basis = BasisSet._6_31G_star
    if args.basis == "cc_pVTZ" : basis = BasisSet._cc_pVTZ

    print(" === Parameterizing %s ===\n" % filename)

    method = FFTypeMethod.CGenFF_2b6
    if  args.forcefield == "GAFF" : method = FFTypeMethod.GAFF
    if  args.forcefield == "GAFF2": method = FFTypeMethod.GAFF2

    mol = FFMolecule(filename=filename, method=method, netcharge=args.charge, rtf=args.rtf, prm=args.prm, basis = basis )

    dihedrals = mol.getSoftDihedrals()
    
    if args.list:
        print("Detected soft torsions:")
        for d in dihedrals:
            print("\t%s-%s-%s-%s" % (mol.name[d[0]], mol.name[d[1]], mol.name[d[2]], mol.name[d[3]]))
        sys.exit(0)

    print("\n == Minimizing ==\n")
    mol.minimize()

    print("\n == Charge fitting ==\n")
    if 1:
      (score, qm_dipole, mm_dipole) = mol.fitCharges()

      rating="GOOD"
      if score > 1:  rating="CHECK"
      if score > 10: rating="BAD"

      print("Charge Chi^2 score : %f : %s" % ( score, rating ))
      print("QM Dipole   : %f %f %f ; %f" % (qm_dipole[0], qm_dipole[1], qm_dipole[2], qm_dipole[3]))
      print("MM Dipole   : %f %f %f ; %f" % (mm_dipole[0], mm_dipole[1], mm_dipole[2], mm_dipole[3]))
      d = qm_dipole[3] - mm_dipole[3]
      d = d * d
      rating="GOOD"
      if score > 1:  rating="CHECK"
      print("Dipole Chi^2 score : %f : %s" % ( d, rating ) )
      print("")

    if 1:
      for d in dihedrals:
        name = "%s-%s-%s-%s" % (mol.name[d[0]], mol.name[d[1]], mol.name[d[2]], mol.name[d[3]])
        if not args.torsion or name in args.torsion:
            print("\n == Fitting torsion %s ==\n" % (name))
            try:
                ret = mol.fitSoftDihedral(d)

                rating="GOOD"
                if ret.chisq > 10:  rating="CHECK"
                if ret.chisq > 100: rating="BAD"
                print("Torsion %s Chi^2 score : %f : %s" % (name, ret.chisq, rating))

                fn = mol.plotDihedralFit(ret, show=False, directory= os.path.join( outdir, "plots") )
            except:
                pass
                # print(fn)

    printEnergies( mol )



    print("\n == Output to %s ==\n" % (outdir) )

    try:
      if args.forcefield == "CGENFF":
        mol._rtf.write( os.path.join( outdir, "mol.rtf" ) )
        mol._prm.write( os.path.join( outdir, "mol.prm" ) ) 
        mol.write( os.path.join( outdir, "mol.psf" ) )
        mol.write( os.path.join( outdir, "mol.xyz" ) )
        mol.write( os.path.join( outdir, "mol.coor" ) )
        mol.write( os.path.join( outdir, "mol.mol2" ) )
        mol.write( os.path.join( outdir, "mol.pdb" ) )
        f = open( os.path.join(outdir, "input.namd" ), "w"  )
        print("parameters mol.prm", file=f )
        print("paraTypeCharmm on", file=f )
        print("coordinates mol.pdb", file=f )
        print("bincoordinates mol.coor", file=f )
        print("temperature 0", file=f )
        print("timestep 0", file=f )
        print("1-4scaling 1.0", file=f )
        print("exclude scaled1-4", file=f )
 
        print("outputname .out", file=f )
        print("outputenergies 1", file=f )
        print("structure mol.psf", file=f )
        print("cutoff 20.", file=f )
        print("switching off", file=f )
        print("stepsPerCycle 1", file=f )
        print("rigidbonds none", file=f )
 
        print("cellBasisVector1 50. 0. 0.", file=f )
        print("cellBasisVector2 0. 50. 0.", file=f )
        print("cellBasisVector3 0. 0. 50.", file=f )
        print("run 0", file=f )
        f.close()
    except ValueError as e:
      print("Not writing CHARMM PRM: %s" % ( str(e) ) )

    try:
     if args.forcefield == "GAFF" or args.forcefield == "GAFF2":
       # types need to be remapped because Amber FRCMOD format limits the type to characters
       # writeFrcmod does this on the fly and returns a mapping that needs to be applied to the mol
       typemap = mol._prm.writeFrcmod( mol._rtf, os.path.join( outdir, "mol.frcmod") )
       mol.write( os.path.join( outdir, "mol.mol2" ), typemap = typemap )
       mol.write( os.path.join( outdir, "mol.pdb" ), typemap = typemap )
       mol.write( os.path.join( outdir, "mol.coor" ), typemap = typemap )
       f = open( os.path.join( outdir, "tleap.in" ) , "w" )
       print("loadAmberParams mol.frcmod", file=f )
       print("A = loadMol2 mol.mol2", file=f )
       print("saveAmberParm A structure.prmtop mol.crd", file=f )
       print("quit", file=f )
       f.close()
       f = open( os.path.join(outdir, "input.namd" ), "w"  )
       print("parmfile structure.prmtop", file=f )
       print("amber on", file=f )
       print("coordinates mol.pdb", file=f )
       print("bincoordinates mol.coor", file=f )
       print("temperature 0", file=f )
       print("timestep 0", file=f )
       print("1-4scaling 0.83333333", file=f )
       print("exclude scaled1-4", file=f )

       print("outputname .out", file=f )
       print("outputenergies 1", file=f )
       print("cutoff 20.", file=f )
       print("switching off", file=f )
       print("stepsPerCycle 1", file=f )
       print("rigidbonds none", file=f )

       print("cellBasisVector1 50. 0. 0.", file=f )
       print("cellBasisVector2 0. 50. 0.", file=f )
       print("cellBasisVector3 0. 0. 50.", file=f )
       print("run 0", file=f )
       f.close()
    except ValueError as e: 
      print("Not writing Amber FRCMOD: %s" % (str(e)) )
    sys.exit(0)


if __name__ == "__main__":
    main_parameterize()
