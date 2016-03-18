# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os

import htmd.lib
from htmd.molecule.support import *


def vmdselection(selection, coordinates, atomname, atomtype, resname, resid, chain=None, segname=None, insert=None,
                 altloc=None, beta=None, occupancy=None, bonds=None):
    import platform
    libdir = htmd.lib.path()

    if coordinates.ndim == 2:
        coordinates = numpy.atleast_3d(coordinates)

    if coordinates.shape[1] != 3:
        print(coordinates.shape)
        raise NameError("Coordinates needs to be natoms x 3 x nframes")

    natoms = coordinates.shape[0]
    nframes = coordinates.shape[2]
    # Sanity check the inputs

    #	print(natoms)
    #	print (len(atomname))

    if bonds and bonds.shape[1] != 2:
        raise NameError("'bonds' not nbonds x 2 in length")
    if len(atomname) != natoms:
        #        print(natoms)
        #        print(len(atomname))
        raise NameError("'atomname' not natoms in length")
    if len(atomtype) != natoms:
        raise NameError("'atomtype' not natoms in length")
    if len(resname) != natoms:
        raise NameError("'resname' not natoms in length")
    if len(resid) != natoms:
        raise NameError("'resid' not natoms in length")
    if chain is not None and len(chain) != natoms:
        raise NameError("'chain' not natoms in length")
    if segname is not None and len(segname) != natoms:
        raise NameError("'segname' not natoms in length")
    if insert is not None and len(insert) != natoms:
        raise NameError("'insert' not natoms in length")
    if altloc is not None and len(altloc) != natoms:
        raise NameError("'altloc' not natoms in length")
    if beta is not None and len(beta) != natoms:
        raise NameError("'beta' not natoms in length")
    if occupancy is not None and len(occupancy) != natoms:
        raise NameError("'occupancy' not natoms in length")

    if platform.system() == "Windows":
        cdll.LoadLibrary(os.path.join(libdir, "libgcc_s_seh-1.dll"))
        if (os.path.exists(os.path.join(libdir, "psprolib.dll"))):
            cdll.LoadLibrary(os.path.join(libdir, "psprolib.dll"))

    parser = cdll.LoadLibrary(os.path.join(libdir, "libvmdparser.so"))

    c_selection = create_string_buffer(selection.encode('ascii'), len(selection) + 1)
    c_natoms = c_int(natoms)
    c_nframes = c_int(nframes)
    c_atomname = pack_string_buffer(atomname)
    c_atomtype = pack_string_buffer(atomtype)
    c_resname = pack_string_buffer(resname)
    c_resid = pack_int_buffer(resid)
    c_chain = None
    c_segname = None
    c_insert = None
    c_altloc = None
    c_beta = None
    c_occupancy = None

    c_coords = None

    c_nbonds = None
    c_bonds = None

    if chain is not None:
        c_chain = pack_string_buffer(chain)
    if segname is not None:
        c_segname = pack_string_buffer(segname)
    if insert is not None:
        c_insert = pack_string_buffer(insert)
    if altloc is not None:
        c_altloc = pack_string_buffer(altloc)
    if beta is not None:
        c_beta = pack_double_buffer(beta)
    if occupancy is not None:
        c_occupancy = pack_double_buffer(occupancy)

    c_bonds = None
    nbonds = 0

    if bonds:  # TODO: Replace the loops for bonds with ravel
        nbonds = bonds.shape[0]
        if nbonds > 0:
            ll = nbonds * 2
            c_bonds = (c_int * ll)()

    for z in range(0, nbonds):
        for y in [0, 1]:
            c_bonds[z * 2 + y] = bonds[z, y]

    c_nbonds = c_int(nbonds)

    ll = natoms * nframes
    c_output_buffer = (c_int * ll)()

    lenv = natoms * 3 * nframes
    c_coords = coordinates.ctypes.data_as(POINTER(c_float))

    retval = parser.atomselect(
        c_selection,
        c_natoms,
        c_beta,
        c_occupancy,
        c_atomtype,
        c_atomname,
        c_resname,
        c_resid,
        c_chain,
        c_segname,
        c_insert,
        c_altloc,
        c_coords,
        c_nframes,
        c_nbonds,
        c_bonds,
        c_output_buffer)

    if retval != 0:
        raise NameError('Could not parse selection "' + selection + '". Is the selection a valid VMD atom selection?')

    retval = numpy.empty((natoms, nframes), dtype=numpy.bool_)

    for frame in range(nframes):
        for atom in range(natoms):
            retval[atom, frame] = c_output_buffer[frame * natoms + atom]

    return numpy.squeeze(retval)


#    return (retval.reshape(natoms, nframes))


def guessbonds(coordinates, atomname, atomtype, resname, resid, chain, segname, insertion, altloc):
    # if it's a single frame, resize to be a 3d array
    if coordinates.ndim == 2:
        c = coordinates.shape
        coordinates = coordinates.reshape((c[0], c[1], 1))

    #    print(coordinates.shape)
    natoms = coordinates.shape[0]
    nframes = coordinates.shape[2]
    # Sanity check the inputs

    #    print(natoms)
    #    print(len(atomname))

    if len(atomname) != natoms:
        raise NameError("'atomname' not natoms in length")
    if len(atomtype) != natoms:
        raise NameError("'atomtype' not natoms in length")
    if len(resname) != natoms:
        raise NameError("'resname' not natoms in length")
    if len(resid) != natoms:
        raise NameError("'resid' not natoms in length")

    libdir = htmd.lib.path()

    parser = cdll.LoadLibrary(os.path.join(libdir, "libvmdparser.so"))

    c_natoms = c_int(natoms)
    c_atomname = pack_string_buffer(atomname)
    c_atomtype = pack_string_buffer(atomtype)
    c_resname = pack_string_buffer(resname)

    c_chain = pack_string_buffer(chain)
    c_segname = pack_string_buffer(segname)
    c_insert = pack_string_buffer(insertion)
    c_altLoc = pack_string_buffer(altloc)
    c_resid = pack_int_buffer(resid)
    c_nframes = c_int(nframes)
    c_coords = None

    c_nbonds = (c_int * 1)()
    lenv = natoms * 8  # some dumb guess about the max # of bonds likely to be created -- natoms*4
    c_bonds = (c_int * lenv)()

    z = 0

    c_coords = coordinates.ctypes.data_as(POINTER(c_float))

    retval = fn = parser.guessbonds(
        c_natoms,
        c_nframes,
        c_atomtype,
        c_atomname,
        c_resname,
        c_resid,
        c_chain,
        c_segname,
        c_insert,
        c_altLoc,
        c_coords,
        c_nbonds,
        c_bonds
    )
    nbonds = c_nbonds[0]
    bonds = numpy.empty((nbonds, 2), dtype=numpy.uint32)
    for y in range(0, nbonds):
        for x in range(0, 2):
            bonds[y, x] = int(c_bonds[y * 2 + x])

    retval = bonds
    return retval.reshape(nbonds, 2)
