# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import ctypes as ct

import numpy
import numpy as np

from htmd.home import home


def wrap(coordinates, bonds, box, centersel=None):
    """Wrap the coordinates back into the unit cell. Molecules will remain continuous, so may escape the bounds of the prinary unit cell.

    Parameters
    ----------

    coordinates :
    bonds       :
    box         :

    Return
    ------
 
    coordinates  
    """

    import platform
    libdir = home(libDir=True)

    if coordinates.dtype != np.float32:
        raise ValueError("Coordinates is not float32")

    if coordinates.ndim == 2:
        c = coordinates.shape
        coordinates = coordinates.reshape((c[0], c[1], 1))

    if coordinates.shape[1] != 3:
        #    print(coordinates.shape)
        raise NameError("Coordinates needs to be natoms x 3 x nframes")

    z = coordinates.shape[2]
    if coordinates.strides[0] != 12 * z or coordinates.strides[1] != 4 * z:
        # It's a view -- need to make a copy to ensure contiguity of memory
        coordinates = numpy.array(coordinates, dtype=numpy.float32)
    if coordinates.strides[0] != 12 * z or coordinates.strides[1] != 4 * z:
        raise ValueError("Coordinates is a view with unsupported strides")

    natoms = coordinates.shape[0]
    nframes = coordinates.shape[2]
    # print(box.shape)
    if numpy.size(bonds, 1) != 2:
        raise NameError("'bonds' not nbonds x 2 in length")
    if numpy.size(box, 0) != 3:
        raise NameError("'box' not nframes x 3 in length")
    if numpy.size(box, 1) != nframes:
        raise NameError("'box' not nframes x 3 in length")

    if platform.system() == "Windows":
        ct.cdll.LoadLibrary(os.path.join(libdir, "libgcc_s_seh-1.dll"))
        if os.path.exists(os.path.join(libdir, "psprolib.dll")):
            ct.cdll.LoadLibrary(os.path.join(libdir, "psprolib.dll"))

    lib = ct.cdll.LoadLibrary(os.path.join(libdir, "libvmdparser.so"))

    nbonds = bonds.shape[0]
    ll = 3 * nbonds
    c_bonds = (ct.c_int * ll)()

    for z in range(0, nbonds):
        for y in [0, 1]:
            c_bonds[z * 2 + y] = bonds[z, y]

    ll = 3 * nframes
    c_box = (ct.c_double * ll)()

    z = 0
    for i in range(0, nframes):
        for j in range(0, 3):
            c_box[z] = box[j][i]
            z = z + 1

    c_nbonds = ct.c_int(nbonds)
    c_natoms = ct.c_int(natoms)
    c_nframes = ct.c_int(nframes)
    lenv = natoms * 3 * nframes

    if centersel is None:
        centersel = numpy.array([-1], dtype=numpy.int32)
    centersel = numpy.append(centersel, numpy.array([-1], dtype=numpy.int32))
    c_centersel = centersel.ctypes.data_as(ct.POINTER(ct.c_int))
    c_coords = coordinates.ctypes.data_as(ct.POINTER(ct.c_float))
    lib.wrap(
        c_bonds,
        c_coords,
        c_box, c_nbonds, c_natoms, c_nframes, c_centersel
    )

    return coordinates
