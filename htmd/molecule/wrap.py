# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from ctypes import *

import numpy

import htmd.home


def wrap( coordinates, bonds, box, centersel=None ):
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
    libdir = htmd.home(libDir=True)

    if coordinates.ndim == 2:
        c = coordinates.shape
        coordinates = coordinates.reshape((c[0], c[1], 1))

    if coordinates.shape[1] != 3:
    #    print(coordinates.shape)
        raise NameError("Coordinates needs to be natoms x 3 x nframes")

    natoms = coordinates.shape[0]
    nframes = coordinates.shape[2]
    #print(box.shape)
    if numpy.size(bonds, 1) != 2:
        raise NameError("'bonds' not nbonds x 2 in length")
    if numpy.size(box, 0) != 3:
        raise NameError("'box' not nframes x 3 in length")
    if numpy.size(box, 1) != nframes:
        raise NameError("'box' not nframes x 3 in length")


    if platform.system() == "Windows":
      cdll.LoadLibrary( os.path.join( libdir, "libgcc_s_seh-1.dll" ) )
      if( os.path.exists( os.path.join( libdir, "psprolib.dll" ) ) ):
        cdll.LoadLibrary( os.path.join( libdir, "psprolib.dll" ) )

    lib = cdll.LoadLibrary( os.path.join( libdir , "libvmdparser.so") )

    nbonds = bonds.shape[0]
    ll = 3 * nbonds
    c_bonds = (c_int * ll)()

    for z in range(0, nbonds):
        for y in [0, 1]:
            c_bonds[z * 2 + y] = bonds[z, y]


    ll = 3 * nframes
    c_box = (c_double * ll)()

    z=0
    for i in range(0, nframes ):
       for j in range(0,3):
          c_box[z] = box[j][i]
          z=z+1

    c_nbonds = c_int(nbonds)
    c_natoms = c_int(natoms)
    c_nframes= c_int(nframes)
    lenv = natoms * 3 * nframes

    c_coords = coordinates.ctypes.data_as(POINTER(c_float))
    lib.wrap(
        c_bonds,
        c_coords,
        c_box, c_nbonds, c_natoms, c_nframes
    )

    if centersel is not None:
      for frame in range(nframes):
        center=numpy.array([0,0,0])
        center[0]  = numpy.mean(coordinates[centersel,0,frame])
        center[1]  = numpy.mean(coordinates[centersel,1,frame])
        center[2]  = numpy.mean(coordinates[centersel,2,frame])
        coordinates[:,:,frame] = coordinates[:,:,frame] - numpy.array([center[0], center[1], center[2]])   
    return coordinates

