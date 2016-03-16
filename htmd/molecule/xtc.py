# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from ctypes import *
import os
import inspect
from htmd.molecule.support import *


class Trajectory:
    box = numpy.array((0, 0))
    natoms = 0
    nframes = 0
    time = numpy.array(0)
    step = numpy.array(0)
    coords = numpy.array((0, 3, 0))

    def __str__(self):
        return "Trajectory with " + str(self.nframes) + " frames, each with " + str(
            self.natoms) + " atoms and timestep of " + str(self.time)


class __xtc(Structure):
    _fields_ = [("box", (c_float * 3)),
                ("natoms", c_int),
                ("step", c_ulong),
                ("time", c_double),
                ("pos", POINTER(c_float))]


def xtc_lib():
    import platform
    libdir = os.path.join(os.path.dirname(inspect.getfile(XTCread)), "..", "lib", platform.system())

    lib = {}


    if platform.system() == "Windows":
      lib['libc'] = cdll.msvcrt
      cdll.LoadLibrary( os.path.join( libdir, "libgcc_s_seh-1.dll" ) )
    else:
      lib['libc'] = cdll.LoadLibrary("libc.so.6")

    lib['libxtc'] = cdll.LoadLibrary(os.path.join(libdir, "libxtc.so"))
    return lib


def XTCread(filename, frames=None):
    lib = xtc_lib()
    nframes = pack_ulong_buffer([0])
    natoms = pack_int_buffer([0])
    deltastep = pack_int_buffer([0])
    deltat    = pack_double_buffer([0])

    lib['libxtc'].xtc_read.restype = POINTER(__xtc)
    lib['libxtc'].xtc_read_frame.restype = POINTER(__xtc)

    if frames is None:
        retval = lib['libxtc'].xtc_read(
            c_char_p(filename.encode("ascii")),
            natoms,
            nframes, deltat, deltastep)
        frames = range( nframes[0] )
        t = Trajectory()
        t.natoms = natoms[0]
        t.nframes = len(frames)
        t.coords = numpy.zeros((natoms[0], 3, t.nframes), dtype=numpy.float32)
        t.step = numpy.zeros(t.nframes, dtype=numpy.uint64)
        t.time = numpy.zeros(t.nframes, dtype=numpy.float32)
        t.box = numpy.zeros((3, t.nframes), dtype=numpy.float32)

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
            t.coords[:, :, i] = numpy.ctypeslib.as_array(retval[f].pos, shape=(natoms[0], 3))

        for f in range(len(frames)):
            lib['libc'].free(retval[f].pos)
        lib['libc'].free(retval)

    else:
        if not isinstance(frames, list) and not isinstance(frames, numpy.ndarray):
            frames = [frames]
        t = Trajectory()
        t.natoms  = 0
        t.nframes = len(frames)
        t.coords = None
        t.step   = None  
        t.time   = None
        t.box    = None

        nframes = len(frames)
        i=0
        for f in frames:
            retval = lib['libxtc'].xtc_read_frame(
                c_char_p(filename.encode("ascii")),
                natoms,
                c_int(f))
            if t.coords is None:
                t.natoms = natoms[0]
                t.coords = numpy.zeros((natoms[0], 3, nframes), dtype=numpy.float32)
                t.step = numpy.zeros(nframes, dtype=numpy.uint64)
                t.time = numpy.zeros(nframes, dtype=numpy.float32)
                t.box = numpy.zeros((3, nframes), dtype=numpy.float32)

            t.step[i] = retval[0].step
            t.time[i] = retval[0].time
            t.box[0, i] = retval[0].box[0]
            t.box[1, i] = retval[0].box[1]
            t.box[2, i] = retval[0].box[2]
            t.coords[:, :, i] = numpy.ctypeslib.as_array(retval[0].pos, shape=(natoms[0], 3))
            i += 1

            lib['libc'].free(retval[0].pos)
            lib['libc'].free(retval)

    if np.size(t.coords, 2) == 0:
        raise NameError('Malformed XTC file. No frames read from: {}'.format(filename))
    if np.size(t.coords, 0) == 0:
        raise NameError('Malformed XTC file. No atoms read from: {}'.format(filename))

    #print( t.step )
    #print( t.time )
    #	print( t.coords[:,:,0] )
    #print(t.coords.shape)
    t.coords *= 10.  # Convert from nm to Angstrom
    t.box *= 10.  # Convert from nm to Angstrom
    return t


def XTCwrite(coords, box, filename, time=None, step=None):
    nframes = np.size(coords, 2)
    if np.size(time) != nframes:
        time = np.zeros(nframes)
    if np.size(step) != nframes:
        step = np.zeros(nframes, dtype=int)

    if os.path.isfile(filename):
        os.unlink(filename)

    lib = xtc_lib()
    bbox = (c_float * 3)()
    natoms = c_int(coords.shape[0])
    cstep = c_int()
    #print(coords.shape)
    for f in range(coords.shape[2]):
        cstep = c_int(step[f])
        ctime = c_float(time[f])  # TODO FIXME
        #print ( step )
        #print ( time )
        bbox[0] = box[0, f] * 0.1
        bbox[1] = box[1, f] * 0.1
        bbox[2] = box[2, f] * 0.1

        data = coords[:, :, f].astype(numpy.float32) * 0.1  # Convert from A to nm
        pos = data.ctypes.data_as(POINTER(c_float))
        lib['libxtc'].xtc_write(
            c_char_p(filename.encode("ascii")),
            natoms,
            cstep,
            ctime,
            pos,
            bbox)
    pass


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path

    m = Molecule(path.join(home(), 'data', 'dhfr', 'dhfr.pdb'))
    m.read(path.join(home(), 'data', 'dhfr', 'dhfr.xtc'))
    m.write('/tmp/test.xtc')
