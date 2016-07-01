import ctypes as ct
import numpy as np
from htmd.molecule.support import pack_double_buffer, pack_int_buffer, pack_string_buffer, pack_ulong_buffer, xtc_lib


class Trajectory:  # TODO: Remove this class
    box = np.array((0, 0))
    natoms = 0
    nframes = 0
    time = np.array(0)
    step = np.array(0)
    coords = np.array((0, 3, 0))

    def __str__(self):
        return "Trajectory with " + str(self.nframes) + " frames, each with " + str(
            self.natoms) + " atoms and timestep of " + str(self.time)


def XTCread(filename, frames=None):
    class __xtc(ct.Structure):
        _fields_ = [("box", (ct.c_float * 3)),
                    ("natoms", ct.c_int),
                    ("step", ct.c_ulong),
                    ("time", ct.c_double),
                    ("pos", ct.POINTER(ct.c_float))]

    lib = xtc_lib()
    nframes = pack_ulong_buffer([0])
    natoms = pack_int_buffer([0])
    deltastep = pack_int_buffer([0])
    deltat = pack_double_buffer([0])

    lib['libxtc'].xtc_read.restype = ct.POINTER(__xtc)
    lib['libxtc'].xtc_read_frame.restype = ct.POINTER(__xtc)

    if frames is None:
        retval = lib['libxtc'].xtc_read(
            ct.c_char_p(filename.encode("ascii")),
            natoms,
            nframes, deltat, deltastep)

        if not retval:
            raise RuntimeError('XTC file {} possibly corrupt.'.format(filename))

        frames = range(nframes[0])
        t = Trajectory()
        t.natoms = natoms[0]
        t.nframes = len(frames)
        t.coords = np.zeros((natoms[0], 3, t.nframes), dtype=np.float32)
        t.step = np.zeros(t.nframes, dtype=np.uint64)
        t.time = np.zeros(t.nframes, dtype=np.float32)
        t.box = np.zeros((3, t.nframes), dtype=np.float32)

        for i, f in enumerate(frames):
            if f >= nframes[0] or f < 0:
                raise NameError('Frame index out of range in XTCread with given frames')
            t.step[i] = retval[f].step
            t.time[i] = retval[f].time
            t.box[0, i] = retval[f].box[0]
            t.box[1, i] = retval[f].box[1]
            t.box[2, i] = retval[f].box[2]
            #		print( t.coords[:,:,f].shape)
            #		print ( t.box[:,f] )
            #   t.step[i] = deltastep[0] * i
            t.coords[:, :, i] = np.ctypeslib.as_array(retval[f].pos, shape=(natoms[0], 3))

        for f in range(len(frames)):
            lib['libc'].free(retval[f].pos)
        lib['libc'].free(retval)

    else:
        if not isinstance(frames, list) and not isinstance(frames, np.ndarray):
            frames = [frames]
        t = Trajectory()
        t.natoms = 0
        t.nframes = len(frames)
        t.coords = None
        t.step = None
        t.time = None
        t.box = None

        nframes = len(frames)
        i = 0
        for f in frames:
            retval = lib['libxtc'].xtc_read_frame(
                ct.c_char_p(filename.encode("ascii")),
                natoms,
                ct.c_int(f))
            if t.coords is None:
                t.natoms = natoms[0]
                t.coords = np.zeros((natoms[0], 3, nframes), dtype=np.float32)
                t.step = np.zeros(nframes, dtype=np.uint64)
                t.time = np.zeros(nframes, dtype=np.float32)
                t.box = np.zeros((3, nframes), dtype=np.float32)

            t.step[i] = retval[0].step
            t.time[i] = retval[0].time
            t.box[0, i] = retval[0].box[0]
            t.box[1, i] = retval[0].box[1]
            t.box[2, i] = retval[0].box[2]
            t.coords[:, :, i] = np.ctypeslib.as_array(retval[0].pos, shape=(natoms[0], 3))
            i += 1

            lib['libc'].free(retval[0].pos)
            lib['libc'].free(retval)

    if np.size(t.coords, 2) == 0:
        raise NameError('Malformed XTC file. No frames read from: {}'.format(filename))
    if np.size(t.coords, 0) == 0:
        raise NameError('Malformed XTC file. No atoms read from: {}'.format(filename))

    # print( t.step )
    # print( t.time )
    #	print( t.coords[:,:,0] )
    # print(t.coords.shape)
    t.coords *= 10.  # Convert from nm to Angstrom
    t.box *= 10.  # Convert from nm to Angstrom
    return t


def CRDread(filename):
    f = open(filename, 'r')
    coords = []

    fieldlen = 12
    k = 0
    for line in f:
        k += 1
        if k < 3:
            continue

        coords += [float(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                   if len(line[i:i + fieldlen].strip()) != 0]

    return [coords[i:i+3] for i in range(0, len(coords), 3)]