import numpy as np
import ctypes as ct
import os
from htmd.molecule.support import xtc_lib


def XTCwrite(coords, box, filename, time=None, step=None):
    nframes = np.size(coords, 2)
    if np.size(time) != nframes:
        time = np.zeros(nframes)
    if np.size(step) != nframes:
        step = np.zeros(nframes, dtype=int)

    if os.path.isfile(filename):
        os.unlink(filename)

    lib = xtc_lib()
    bbox = (ct.c_float * 3)()
    natoms = ct.c_int(coords.shape[0])
    cstep = ct.c_int()
    # print(coords.shape)
    for f in range(coords.shape[2]):
        cstep = ct.c_int(step[f])
        ctime = ct.c_float(time[f])  # TODO FIXME
        # print ( step )
        # print ( time )
        bbox[0] = box[0, f] * 0.1
        bbox[1] = box[1, f] * 0.1
        bbox[2] = box[2, f] * 0.1

        data = coords[:, :, f].astype(np.float32) * 0.1  # Convert from A to nm
        pos = data.ctypes.data_as(ct.POINTER(ct.c_float))
        lib['libxtc'].xtc_write(
            ct.c_char_p(filename.encode("ascii")),
            natoms,
            cstep,
            ctime,
            pos,
            bbox)