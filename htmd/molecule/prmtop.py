# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import numpy;
import os
from htmd.molecule.support import *


class PRMTOP:
    charges = None
    bonds = None


def PRMTOPread(filename):
    psf = PRMTOP
    f = open(filename, 'r')
    mode = None
    c = 0

    lines = []
    for line in f:
        lines.append(line)

    p_idx = -1
    c_idx = -1
    b1_idx = -1
    b2_idx = -1
    for i in range(0, len(lines)):
        #		print(lines[i])
        if lines[i].startswith("%FLAG POINTERS"):
            p_idx = i
        if lines[i].startswith("%FLAG CHARGE"):
            c_idx = i
        if lines[i].startswith("%FLAG BONDS_INC_HYDROGEN"):
            b1_idx = i
        if lines[i].startswith("%FLAG BONDS_WITHOUT_HYDROGEN"):
            b2_idx = i

    if (p_idx == -1 or c_idx == -1 or b1_idx == -1 or b2_idx == -1):
        raise NameError("Syntax error in PRMTOP file")

    l = lines[p_idx + 2].split()
    natoms = int(l[0])
    nb1 = int(l[2])
    nb2 = int(l[3])
    nc = natoms

    prmtop = PRMTOP()

    prmtop.charges = numpy.zeros((natoms))
    prmtop.bonds = numpy.zeros((nb1 + nb2, 2), dtype=int)

    i = c_idx + 2
    idx = 0
    while not lines[i].startswith('%'):
        v = lines[i].split()
        for c in v:
            prmtop.charges[idx] = float(c) / 18.2223  # Scaling factor for charges
            idx += 1
        i += 1

    idx = 0
    v = []

    i = b1_idx + 2
    while not lines[i].startswith('%'):
        v.extend(lines[i].split())
        i += 1
    i = b2_idx + 2
    while not lines[i].startswith('%'):
        v.extend(lines[i].split())
        i += 1

    for c in range(0, len(v), 3):
        prmtop.bonds[idx, 0] = int(v[c + 0]) / 3
        prmtop.bonds[idx, 1] = int(v[c + 1]) / 3
        idx += 1

    return prmtop
    pass


def PRMTOPwrite(filename):
    pass


if __name__ == "__main__":
    from htmd.home import home
    from os import path
    from htmd.molecule.molecule import Molecule

    m = Molecule(path.join(home(), 'data', 'amber', 'amber.pdb'))
    m.read(path.join(home(), 'data', 'amber', 'amber.prmtop') )
    print( m.charge )
    print( m.bonds )

