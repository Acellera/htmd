# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import matplotlib
matplotlib.use('Agg')

import argparse

from htmd.newparameterization.ffmolecule import FFMolecule
from htmd.newparameterization.fftype import FFTypeMethod
import sys
import os


def main_parameterize():
  ncpus = os.cpu_count()
  try:
    ncpus = int( os.getenv("NCPUS") )
  except:
    pass

  parser = argparse.ArgumentParser( description="Acellera Small Molecule Parameterisation Version 2.0" )
  parser.add_argument( "-m", "--mol2", help="Molecule to parameterise, in mol2 format", type=str, default="input.mol2", action="store", dest="mol" )
  parser.add_argument( "-c", "--charge", help="Net charge on molecule", type=int, default=None, action="store", dest="charge" )
  parser.add_argument( "--rtf", help="Inital RTF parameters (req --prm)", type=str, default=None, dest="rtf" )
  parser.add_argument( "--prm", help="Inital PRM parameters (req --rtf)", type=str, default=None, dest="prm" )
  parser.add_argument( "-o", "--output", metavar="DIRECTORY", help="Output directory", type=str, default="parameters" )
  parser.add_argument( "-l", "--list",help="List parameterisable torsions", action="store_true", default=False, dest="list" )
  parser.add_argument( "-t", "--torsion", metavar="A1-A2-A3-A4",  help="Torsion to parameterise (default all)", action="append", default=None, dest="torsion" ) 
  parser.add_argument( "-n", "--ncpus",  help="Number of CPUs to use (default %d)" % (ncpus),  default=ncpus, dest="ncpus" ) 

  args =parser.parse_args()

  # Communicate the # of CPUs to use to the QM engine via environment variable
  os.putenv("NCPUS", args.ncpus )

  filename = args.mol

  if( not os.path.exists( filename ) ):
    print("file %s not found" % ( filename ) )
    sys.exit(1)
 

  print(" === Parameterizing %s ===\n" % ( filename )) 

  mol = FFMolecule( filename=filename, method=FFTypeMethod.CGenFF_2b6, netcharge=args.charge, rtf=args.rtf, prm=args.prm )
  

  dihedrals = mol.getSoftDihedrals()

  if args.list:
    print("Detected soft torsions:" )
    for d in dihedrals:
      print("\t%s-%s-%s-%s" % ( mol.name[ d[0] ], mol.name[ d[1] ], mol.name[ d[2] ], mol.name[ d[3] ]  ))
    sys.exit(0)
 
  print(" == Minimizing ==\n")
  mol.minimize()

  print(" == Charge fitting ==\n")
#  (score, qm_dipole, mm_dipole) = mol.fitCharges()


  print( "\tCharge Chi^2 score : %f" %( score ) )
  print( "\tQM Dipole   : %f %f %f ; %f" % ( qm_dipole[0], qm_dipole[1], qm_dipole[2], qm_dipole[3] ) )
  print( "\tMM Dipole   : %f %f %f ; %f" % ( mm_dipole[0], mm_dipole[1], mm_dipole[2], mm_dipole[3] ) )
  print( "" )	


  for d in dihedrals:
    name="%s-%s-%s-%s" % ( mol.name[ d[0] ], mol.name[ d[1] ], mol.name[ d[2] ], mol.name[ d[3] ]  )
    if not args.torsion or name in args.torsion:
      print(" == Fitting torsion %s ==\n" % (name ) )
      try:
        ret = mol.fitSoftDihedral( d )

        print("\tTorsion %s Chi^2 score : %f" % ( name, ret.chisq ) )

        fn  = mol.plotDihedralFit( ret, show=False, directory="plots" )
      except:
        pass
    #print(fn)

  print(" == Output to %s ==\n", args.output );

  try:
    os.mkdir( args.output )
  except:
    pass

  mol._rtf.write( "parameters/mol.rtf" )
  mol._prm.write( "parameters/mol.prm" )
  mol.write( "parameters/mol.psf" )
  mol.write( "parameters/mol.xyz" )
  mol.write( "parameters/mol.pdb" )
  mol.write( "parameters/mol.mol2" )
  
  sys.exit(0)	


if __name__ == "__main__":
  main_parameterize()
