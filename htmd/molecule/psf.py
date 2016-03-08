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


def PSFread(filename):
    psf = PSF
    f = open(filename, 'r')
    mode = None
    c = 0

    for line in f:
        if line.strip() == "":
            mode = None

        if mode == 'atom':
            l = line.split()
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

        if '!NATOM' in line:
            mode = 'atom'
            l = line.split()
            psf.charges = numpy.zeros([int(l[0])])
            psf.masses = numpy.zeros([int(l[0])])
            c = 0
        elif '!NBOND' in line:
            mode = 'bond'
            l = line.split()
            psf.bonds = numpy.zeros([int(l[0]), 2], dtype=int)
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
        if (m.masses is not None) and (i < len(m.masses)):
            mass = m.masses[i]
        if (m.charge is not None) and (i < len(m.charge)):
            charge = m.charge[i]
        print("%8d P1   %-5d%-4s %-4s %-5s %.6f  %8.5f  %10d" %
              (int(m.serial[i]),
               int(m.resid[i]),
               (m.resname[i]),
               m.name[i],
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
    print("%10d !NTHETA: angles\n" % (0), file=f)
    print("%10d !NPHI: dihedrals\n" % (0), file=f)
    print("%10d !NIMPHI: impropers\n" % (0), file=f)
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

    m = Molecule(path.join(home(), 'data', 'dhfr', 'dhfr.pdb'))
    m.read(path.join(home(), 'data', 'dhfr', 'dhfr.psf'))
    m.write(path.join(home(), 'data', 'dhfr', 'test.psf'))

