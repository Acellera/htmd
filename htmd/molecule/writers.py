import numpy as np
import ctypes as ct
import os
from htmd.molecule.support import xtc_lib


def XTCwrite(coords, box, filename, time=None, step=None):
    numFrames = coords.shape[2]
    if np.size(box, 1) != numFrames:  # Box should have as many frames as trajectory
        box = np.tile(box, (1, numFrames))

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


def BINCOORwrite(coords, filename):
    import struct
    natoms = np.array([coords.shape[0]])
    f = open(filename, 'wb')

    dat = coords[:, :, 0]
    dat = dat.reshape(dat.shape[0] * 3).astype(np.float64)

    fmt1 = 'i' * natoms.shape[0]
    bin1 = struct.pack(fmt1, *natoms)
    fmt2 = 'd' * dat.shape[0]
    bin2 = struct.pack(fmt2, *dat)
    f.write(bin1)
    f.write(bin2)
    f.close()


def PSFwrite(molecule, filename):
    m = molecule

    f = open(filename, 'w')
    print("PSF\n", file=f)
    print("%8d !TITLE\n" % (0), file=f)
    print("%8d !NATOM" % (len(m.serial)), file=f)
    for i in range(len(m.serial)):
        charge = 0
        mass = 1
        atomtype = ""
        if (m.masses is not None) and (i < len(m.masses)):
            mass = m.masses[i]
        if (m.charge is not None) and (i < len(m.charge)):
            charge = m.charge[i]
        if( m.atomtype is not None ) and ( i<len(m.atomtype)):
            atomtype   = m.atomtype[i]
        print("%8d %-4s %-5s%-4s %-4s %-6s %-2s %10.6f  %8.6f  %10d" %
              (int(m.serial[i]),
               m.segid[i],
               str(m.resid[i])+str(m.insertion[i]),
               (m.resname[i]),
               m.name[i],
               atomtype,
               m.element[i],
               charge,
               mass,
               0
               ),
              file=f)
    print("\n\n", file=f)
    print(" %8d !NBOND: bonds" % (m.bonds.shape[0]), file=f)
    for i in range(m.bonds.shape[0]):
        if i and not (i % 4):
            print("", file=f)
        print("%10d%10d" % (m.bonds[i, 0] + 1, m.bonds[i, 1] + 1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NTHETA: angles" % (m.angles.shape[0]), file=f)
    for i in range(m.angles.shape[0]):
        if i and not (i % 3):
            print("", file=f)
        print("%10d%10d%10d" % (m.angles[i, 0] + 1, m.angles[i, 1] + 1, m.angles[i,2]+1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NPHI: dihedrals" % (m.dihedrals.shape[0]), file=f)
    for i in range(m.dihedrals.shape[0]):
        if i and not (i % 3):
            print("", file=f)
        print("%10d%10d%10d%10d" % (m.dihedrals[i, 0] + 1, m.dihedrals[i, 1] + 1, m.dihedrals[i,2]+1, m.dihedrals[i,3]+1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NIMPHI: impropers" % (m.impropers.shape[0]), file=f)
    for i in range(m.impropers.shape[0]):
        if i and not (i % 3):
            print("", file=f)
        print("%10d%10d%10d%10d" % (m.impropers[i, 0] + 1, m.impropers[i, 1] + 1, m.impropers[i,2]+1, m.impropers[i,3]+1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NDON: donors\n" % (0), file=f)
    print("%10d !NACC: acceptors\n" % (0), file=f)
    print("%10d !NNB: acceptors\n" % (0), file=f)
    print("%10d %10d !NGRP \n" % (0, 0), file=f)
    f.close()


def XYZwrite(src, filename):
    import re
    fh = open(filename, "w")
    natoms = len(src.record)
    print("%d\n" % (natoms), file=fh)
    for i in range(natoms):
        e = src.element[i].strip()
        if not len(e):
            e = re.sub("[1234567890]*", "", src.name[i])
        print("%s   %f   %f    %f" % (e, src.coords[i, 0, src.frame], src.coords[i, 1, src.frame], src.coords[i, 2, src.frame]), file=fh)
    fh.close()
