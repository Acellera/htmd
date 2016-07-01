# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import numpy;
import os
from htmd.molecule.support import *


def PRMTOPread(filename):
    f = open(filename, 'r')
    names = []
    charges = []
    masses = []
    uqresnames = []
    residx = []
    bondsidx = []

    section = None
    for line in f:
        if line.startswith('%FLAG POINTERS'):
            section = 'pointers'
        elif line.startswith('%FLAG ATOM_NAME'):
            section = 'names'
        elif line.startswith('%FLAG CHARGE'):
            section = 'charges'
        elif line.startswith('%FLAG MASS'):
            section = 'masses'
        elif line.startswith('%FLAG ATOM_TYPE_INDEX'):
            section = 'type'
        elif line.startswith('%FLAG RESIDUE_LABEL'):
            section = 'resname'
        elif line.startswith('%FLAG RESIDUE_POINTER'):
            section = 'resstart'
        elif line.startswith('%FLAG BONDS_INC_HYDROGEN') or line.startswith('%FLAG BONDS_WITHOUT_HYDROGEN'):
            section = 'bonds'
        elif line.startswith('%FLAG BOX_DIMENSIONS'):
            section = 'box'
        elif line.startswith('%FLAG'):
            section = None

        if line.startswith('%'):
            continue

        if section == 'pointers':
            pass
        elif section == 'names':
            fieldlen = 4
            names += [line[i:i + fieldlen].strip() for i in range(0, len(line), fieldlen)
                      if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'charges':
            fieldlen = 16
            charges += [float(line[i:i + fieldlen].strip()) / 18.2223 for i in range(0, len(line), fieldlen)
                        if len(line[i:i + fieldlen].strip()) != 0]  # 18.2223 = Scaling factor for charges
        elif section == 'masses':
            fieldlen = 16
            masses += [float(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                       if len(line[i:i + fieldlen].strip()) != 0]  # 18.2223 = Scaling factor for charges
        elif section == 'resname':
            fieldlen = 4
            uqresnames += [line[i:i + fieldlen].strip() for i in range(0, len(line), fieldlen)
                           if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'resstart':
            fieldlen = 8
            residx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                       if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'bonds':
            fieldlen = 8
            bondsidx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                         if len(line[i:i + fieldlen].strip()) != 0]

    # Replicating unique resnames according to their start and end indeces
    residx.append(len(names)+1)
    resnames = []
    resid = []
    for i in range(len(residx) - 1):
        numresatoms = residx[i+1] - residx[i]
        resnames += [uqresnames[i]] * numresatoms
        resid += [i+1] * numresatoms

    # Processing bond triplets
    bonds = []
    for i in range(0, len(bondsidx), 3):
        bonds.append([int(bondsidx[i] / 3), int(bondsidx[i+1] / 3)])

    return names, charges, masses, resnames, resid, bonds


def PRMTOPwrite(filename):
    pass


if __name__ == "__main__":
    from htmd.home import home
    from os import path
    from htmd.molecule.molecule import Molecule

    m = Molecule(path.join(home(), 'data', 'amber', 'test-prmtop.pdb'))
    m.read(path.join(home(), 'data', 'amber', 'test-prmtop.prmtop') )
    print( m.charge )
    print( m.bonds )

