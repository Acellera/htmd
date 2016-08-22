# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.newparameterization.ffmolecule import FFMolecule
from htmd.newparameterization.fftype import FFTypeMethod
import sys
import os


def main_scan():
  if(len(sys.argv)!=2):
    print("Syntax: %s input.mol" % ( sys.argv[0] ) )
    sys.exit(0)
 
  filename = sys.argv[1]

  print("Scan %s" % ( filename )) 

  mol = FFMolecule( filename=filename, method=FFTypeMethod.CGenFF_2b6 )

  dihedrals = mol.getSoftDihedrals()
  for d in dihedrals:
    print("\nScanning dihedral %s-%s-%s-%s" % ( mol.name[ d[0] ], mol.name[ d[1] ], mol.name[ d[2] ], mol.name[ d[3] ]  ))
    qmset =mol.scanSoftDihedral( d, directory="scan", step=20 );
    print(qmset)


if __name__ == "__main__":
  main_scan()
