# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import numpy
import os
from htmd.molecule.support import *
import struct


class BINCOOR:
    charges = None
    bonds = None
    masses = None


def BINCOORread(filename):
    f = open(filename, 'rb')
    dat = f.read(4)
    fmt = 'i'
    natoms = struct.unpack(fmt, dat)[0]
    dat = f.read(natoms * 3 * 8)
    fmt = 'd' * (natoms * 3)
    coords = struct.unpack(fmt, dat)
    coords = numpy.array(coords).reshape((natoms, 3, 1))
    f.close() 
    return coords


def BINCOORwrite(coords, filename):
    natoms = numpy.array([coords.shape[0]])
    f = open(filename, 'wb')

    dat = coords[:, :, 0]
    dat = dat.reshape(dat.shape[0] * 3).astype(numpy.float64)

    fmt1 = 'i' * natoms.shape[0]
    bin1 = struct.pack(fmt1, *natoms)
    fmt2 = 'd' * dat.shape[0]
    bin2 = struct.pack(fmt2, *dat)
    f.write(bin1)
    f.write(bin2)
    f.close()


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path

    m = Molecule(path.join(home(), 'data', 'dhfr', 'dhfr.pdb'))
    m.read(path.join(home(), 'data', 'dhfr', 'dhfr.psf'))
    m.write(path.join(home(), 'data', 'dhfr', 'test.psf'))
    m.write('/tmp/test.coor')
    m.read('/tmp/test.coor')
    m.write('/tmp/test2.coor')


