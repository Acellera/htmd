# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.newparameterization import ffmolecule
from htmd.newparameterization.fftype import FFTypeMethod
import sys
import os


def main_parameterize():
  if(len(sys.argv)!=2):
    print("Syntax: %s input.mol" % ( sys.argv[0] ) )
    sys.exit(0)
 
  filename = sys.argv[1]

  print("Parameterizing %s" % ( filename )) 

  mol = FFMolecule( filename=filename, method=FFTypeMethod.CGenFF_2b6 )
  

  print("Minimizing")
  mol.minimize()

  print("Charge fitting")
  (score, qm_dipole, mm_dipole) = mol.fitCharges()

  print( "Chi^2 score : %f" %( score ) )
  print( "QM Dipole   : %f %f %f ; %f" % ( qm_dipole[0], qm_dipole[1], qm_dipole[2], qm_dipole[3] ) )
  print( "MM Dipole   : %f %f %f ; %f" % ( mm_dipole[0], mm_dipole[1], mm_dipole[2], mm_dipole[3] ) )

  dihedrals = mol.getSoftDihedrals()
  for d in dihedrals:
    print("\nFitting dihedral %s-%s-%s-%s" % ( mol.name[ d[0] ], mol.name[ d[1] ], mol.name[ d[2] ], mol.name[ d[3] ]  ))
    ret = mol.fitSoftDihedral( d )

    print("Chi^2 score : %f" % ( ret.chisq ) )

    fn  = mol.plotDihedralFit( ret, show=False, directory="plots" )
    #print(fn)

  try:
    os.mkdir("parameters")
  except:
    pass
  mol._rtf.write( "parameters/mol.rtf" )
  mol._prm.write( "parameters/mol.prm" )
  mol.write( "parameters/mol.psf" )
  mol.write( "parameters/mol.xyz" )
  mol.write( "parameters/mol.pdb" )
  



if __name__ == "__main__":
  main_parameterize()
