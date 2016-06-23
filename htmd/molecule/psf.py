# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import htmd
import numpy;
import os
from htmd.molecule.support import *


class PSF:
    charges = None
    bonds = None
    masses = None
    resid   = None
    resname = None   
    atomname = None
    atomtype = None
    serial = None
    segid = None
    insertion = None
    angles = None
    dihedrals = None
    impropers = None

 
def PSFread(filename):
    import re
    residinsertion = re.compile('(\d+)([a-zA-Z])')
    psf = PSF
    f = open(filename, 'r')
    mode = None
    c = 0

    for line in f:
        if line.strip() == "":
            mode = None

        if mode == 'atom':
            l = line.split()
            psf.serial[c]   = l[0]
            psf.segid[c]    = l[1]
            match = residinsertion.findall(l[2])
            if match:
                resid = int(match[0][0])
                insertion = match[0][1]
            else:
                resid = int(l[2])
                insertion = ''
            psf.resid[c]    = resid
            psf.insertion[c] = insertion
            psf.resname[c]  = l[3]
            psf.atomname[c] = l[4]
            psf.atomtype[c] = l[5]
            psf.charges[c] = float(l[6])
            psf.masses[c] = float(l[7])
            c += 1
        # psf.charges = numpy.append( psf.charges, l[6] )
        #	psf.masses  = numpy.append( psf.masses,  l[7] )
        elif mode == 'bond':
            l = line.split()
            for x in range(0, len(l), 2):
                psf.bonds[c, 0] = int(l[x]) - 1
                psf.bonds[c, 1] = int(l[x + 1]) - 1
                c += 1
        elif mode == 'angle':
            l = line.split()
            for x in range(0, len(l), 3):
                psf.angles[c, 0] = int(l[x]) - 1
                psf.angles[c, 1] = int(l[x + 1]) - 1
                psf.angles[c, 2] = int(l[x + 2]) - 1
                c += 1
        elif mode == 'dihedral':
            l = line.split()
            for x in range(0, len(l), 4):
                psf.dihedrals[c, 0] = int(l[x]) - 1
                psf.dihedrals[c, 1] = int(l[x + 1]) - 1
                psf.dihedrals[c, 2] = int(l[x + 2]) - 1
                psf.dihedrals[c, 3] = int(l[x + 3]) - 1
                c += 1
        elif mode == 'improper':
            l = line.split()
            for x in range(0, len(l), 4):
                psf.impropers[c, 0] = int(l[x]) - 1
                psf.impropers[c, 1] = int(l[x + 1]) - 1
                psf.impropers[c, 2] = int(l[x + 2]) - 1
                psf.impropers[c, 3] = int(l[x + 3]) - 1
                c += 1







        if '!NATOM' in line:
            mode = 'atom'
            l = line.split()
            psf.segid    = numpy.empty([(int(l[0]))], dtype=object )
            psf.serial   = numpy.zeros([int(l[0])], dtype=numpy.uint32)
            psf.resid    = numpy.zeros([int(l[0])], dtype=numpy.uint32)
            psf.resname  = numpy.empty([(int(l[0]))], dtype=object )
            psf.atomname = numpy.empty([(int(l[0]))], dtype=object )
            psf.atomtype = numpy.empty([(int(l[0]))], dtype=object )
            psf.charges = numpy.zeros([int(l[0])], dtype=numpy.float32)
            psf.masses = numpy.zeros([int(l[0])], dtype=numpy.float32)
            psf.insertion = numpy.empty([(int(l[0]))], dtype=object)
            
            c = 0
        elif '!NBOND' in line:
            mode = 'bond'
            l = line.split()
            psf.bonds = numpy.zeros([int(l[0]), 2], dtype=numpy.uint32)
            c = 0
        elif '!NTHETA' in line:
            mode = 'angle'
            l = line.split()
            psf.angles = numpy.zeros([int(l[0]), 3], dtype=numpy.uint32)
            c = 0
        elif '!NPHI' in line:
            mode = 'dihedral'
            l = line.split()
            psf.dihedrals = numpy.zeros([int(l[0]), 4], dtype=numpy.uint32)
            c = 0
        elif '!NIMPHI' in line:
            mode = 'improper'
            l = line.split()
            psf.impropers = numpy.zeros([int(l[0]), 4], dtype=numpy.uint32)
            c = 0





    f.close()
    return psf
    pass


def PSFwrite(molecule, filename):
    m = molecule

    f = open(filename, 'w')
    print("PSF\n", file=f)
    print("%8d !TITLE\n" % (0), file=f)
    print("%8d !NATOM" % (len(m.serial)), file=f)
    for i in range(len(m.serial)):
        charge = 0
        mass = 1
        atomtype = ""
        if (m.masses is not None) and (i < len(m.masses)):
            mass = m.masses[i]
        if (m.charge is not None) and (i < len(m.charge)):
            charge = m.charge[i]
        if( m.atomtype is not None ) and ( i<len(m.atomtype)):
            atomtype   = m.atomtype[i]
        print("%8d %-4s %-5s%-4s %-4s %-6s %-2s %10.6f  %8.6f  %10d" %
              (int(m.serial[i]),
               m.segid[i],
               str(m.resid[i])+m.insertion[i],
               (m.resname[i]),
               m.name[i],
               atomtype,
               m.element[i],
               charge,
               mass,
               0
               ),
              file=f)
    print("\n\n", file=f)
    print(" %8d !NBOND: bonds" % (m.bonds.shape[0]), file=f)
    for i in range(m.bonds.shape[0]):
        if (i and not (i % 4)):
            print("", file=f)
        print("%10d%10d" % (m.bonds[i, 0] + 1, m.bonds[i, 1] + 1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NTHETA: angles\n" % (m.angles.shape[0]), file=f)
    for i in range(m.angles.shape[0]):
        if (i and not (i % 3)):
            print("", file=f)
        print("%10d%10d%10d" % (m.angles[i, 0] + 1, m.angles[i, 1] + 1, m.angles[i,2]+1), file=f, end="")


    print("%10d !NPHI: dihedrals\n" % (m.dihedrals.shape[0]), file=f)
    for i in range(m.dihedrals.shape[0]):
        if (i and not (i % 3)):
            print("", file=f)
        print("%10d%10d%10d%10d" % (m.dihedrals[i, 0] + 1, m.dihedrals[i, 1] + 1, m.dihedrals[i,2]+1, m.dihedrals[i,3]+1), file=f, end="")



    print("%10d !NIMPHI: impropers\n" % (m.impropers.shape[0]), file=f)
    for i in range(m.impropers.shape[0]):
        if (i and not (i % 3)):
            print("", file=f)
        print("%10d%10d%10d%10d" % (m.impropers[i, 0] + 1, m.impropers[i, 1] + 1, m.impropers[i,2]+1, m.impropers[i,3]+1), file=f, end="")



    print("%10d !NDON: donors\n" % (0), file=f)
    print("%10d !NACC: acceptors\n" % (0), file=f)
    print("%10d !NNB: acceptors\n" % (0), file=f)
    print("%10d %10d !NGRP \n" % (0, 0), file=f)
    f.close()
    pass


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    import sys

    m =  Molecule(path.join(home(), 'data', 'dhfr', 'dhfr.pdb')) 
    m.write( 'test.psf' )
    m =  Molecule(path.join(home(), 'data', 'dhfr', 'dhfr.psf'))
    m.write( 'test.psf' )

    from htmd.molecule.vmdparser import guessAnglesAndDihedrals
    (angles, dihedrals ) = guessAnglesAndDihedrals( m.bonds )

    if( angles.shape != m.angles.shape ): 
      print( angles.shape )
      print( m.angles.shape )
      print("Mismatch in Guessed Angles" ) 
      sys.exit(1)

    if( dihedrals.shape != m.dihedrals.shape ): 
      print( dihedrals.shape )
      print( m.dihedrals.shape )
      print("Mismatch in Guessed dihedrals" ) 
      sys.exit(1)

    m.write('test.psf')

#(2496, 3)
#(11584, 3)
#(1872, 4)
#(6701, 4)

