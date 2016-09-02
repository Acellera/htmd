import numpy as np
import ctypes as ct
import os
from htmd.molecule.support import xtc_lib
import collections
import logging
logger = logging.getLogger(__name__)


_Pair = collections.namedtuple('Atom', 'resname name')
# The following is taken from MDAnalysis to align atom names based on the strict PDB formatting.
# Some programs like AutoDock require this strict formatting to guess atom types.

# These attributes are used to deduce how to format the atom name.
_ions = ('FE', 'AS', 'ZN', 'MG', 'MN', 'CO', 'BR',
        'CU', 'TA', 'MO', 'AL', 'BE', 'SE', 'PT',
        'EU', 'NI', 'IR', 'RH', 'AU', 'GD', 'RU')
# Mercurial can be confused for hydrogen gamma. Yet, mercurial is
# rather rare in the PDB. Here are all the residues that contain
# mercurial.
_special_hg = ('CMH', 'EMC', 'MBO', 'MMC', 'HGB', 'BE7', 'PMB')
# Chloride can be confused for a carbon. Here are the residues that
# contain chloride.
_special_cl = ('0QE', 'CPT', 'DCE', 'EAA', 'IMN', 'OCZ', 'OMY', 'OMZ',
              'UN9', '1N1', '2T8', '393', '3MY', 'BMU', 'CLM', 'CP6',
              'DB8', 'DIF', 'EFZ', 'LUR', 'RDC', 'UCL', 'XMM', 'HLT',
              'IRE', 'LCP', 'PCI', 'VGH')
# In these pairs, the atom name is aligned on the first column
# (column 13).
_include_pairs = (_Pair('OEC', 'CA1'),
                  _Pair('PLL', 'PD'),
                  _Pair('OEX', 'CA1'))
# In these pairs, the atom name is aligned on the second column
# (column 14), but other rules would align them on the first column.
_exclude_pairs = (_Pair('C14', 'C14'), _Pair('C15', 'C15'),
                  _Pair('F9F', 'F9F'), _Pair('OAN', 'OAN'),
                  _Pair('BLM', 'NI'), _Pair('BZG', 'CO'),
                  _Pair('BZG', 'NI'), _Pair('VNL', 'CO1'),
                  _Pair('VNL', 'CO2'), _Pair('PF5', 'FE1'),
                  _Pair('PF5', 'FE2'), _Pair('UNL', 'UNL'))


def _deduce_PDB_atom_name(name, resname):
    """Deduce how the atom name should be aligned.
    Atom name format can be deduced from the atom type, yet atom type is
    not always available. This function uses the atom name and residue name
    to deduce how the atom name should be formatted. The rules in use got
    inferred from an analysis of the PDB. See gist at
    <https://gist.github.com/jbarnoud/37a524330f29b5b7b096> for more
    details.
    """
    if len(name) >= 4:
        return name[:4]
    elif len(name) == 1:
        return ' {}  '.format(name)
    elif ((resname == name
           or name[:2] in _ions
           or name == 'UNK'
           or (resname in _special_hg and name[:2] == 'HG')
           or (resname in _special_cl and name[:2] == 'CL')
           or _Pair(resname, name) in _include_pairs)
          and _Pair(resname, name) not in _exclude_pairs):
        return '{:<4}'.format(name)
    return ' {:<3}'.format(name)


def checkTruncations(mol):
    fieldsizes = {'record': 6, 'serial': 5, 'name': 4, 'altloc': 1, 'resname': 4, 'chain': 1, 'resid': 4,
                  'insertion': 1, 'segid': 4, 'element': 2}
    for f in fieldsizes:
        if np.any([True if len(x) > fieldsizes[f] else False for x in mol.__dict__[f].astype('str')]):
            if fieldsizes[f] == 1:
                logger.warning('Field "{}" of PDB overflows. Your data will be truncated to 1 character.'.format(f))
            else:
                logger.warning('Field "{}" of PDB overflows. Your data will be truncated to {} characters.'.format(f, fieldsizes[f]))


def PDBwrite(mol, filename):
    def format83(f):
        """Format a single float into a string of width 8, with ideally 3 decimal places of precision. If the number is
        a little too large, we can gracefully degrade the precision by lopping off some of the decimal places. If it's
        much too large, we throw a NameError"""
        if -999.999 < f < 9999.999:
            return '%8.3f' % f
        if -9999999 < f < 99999999:
            return ('%8.3f' % f)[:8]
        raise NameError('coordinate "%s" could not be represented '
                        'in a width-8 field' % f)

    def format62(f):
        if -9.999 < f < 99.999:
            return '%6.2f' % f
        if -99999 < f < 999999:
            return ('%6.2f' % f)[:6]
        raise NameError('coordinate "%s" could not be represented '
                        'in a width-6 field' % f)

    checkTruncations(mol)
    coords = np.atleast_3d(mol.coords[:, :, mol.frame])
    numFrames = coords.shape[2]
    serial = np.arange(1, np.size(coords, 0) + 1)

    fh = open(filename, 'w')
    # TODO FIXME  -- should take box from traj frame
    box = mol.box
    if box is not None:
        box = np.atleast_2d(np.atleast_2d(box)[:, mol.frame])
        print("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % (box[0, 0], box[0, 1], box[0, 2], 90, 90, 90),
            file=fh)

    for frame in range(numFrames):
        print("MODEL    %5d" % (frame + 1), file=fh)
        for i in range(0, len(mol.record)):
            name = _deduce_PDB_atom_name(mol.name[i], mol.resname[i])

            if mol.serial[i] < 100000:
                ser = str(int(mol.serial[i]))
            else:
                ser = '*****'

            print(
                "{!s:6.6}{!s:>5.5} {}{!s:>1.1}{!s:4.4}{!s:>1.1}{!s:>4.4}{!s:>1.1}   {}{}{}{}{}      {!s:4.4}{!s:>2.2}  ".format(
                    mol.record[i],
                    ser, name, mol.altloc[i],
                    mol.resname[i], mol.chain[i],
                    mol.resid[i],
                    mol.insertion[i],
                    format83(mol.coords[i, 0, frame]),
                    format83(mol.coords[i, 1, frame]),
                    format83(mol.coords[i, 2, frame]),
                    format62(mol.occupancy[i]),
                    format62(mol.beta[i]),
                    mol.segid[i],
                    mol.element[i]
                ), file=fh
            )
            # TODO : convert charges to ints if we ever write them
            if i < len(mol.record) - 1 and mol.segid[i] != mol.segid[i + 1]:
                print("TER", file=fh)

        if mol.bonds is not None and len(mol.bonds) != 0:
            bondedatoms = np.unique(mol.bonds)

            for a in bondedatoms:
                idx = mol.bonds[mol.bonds[:, 0] == a, 1]
                idx = np.unique(np.append(idx, mol.bonds[mol.bonds[:, 1] == a, 0]))
                # I need to support multi-line printing of atoms with more than 4 bonds
                for k in range(int(np.ceil(len(idx) / 4))):
                    print("CONECT%5d" % (a + 1), file=fh, end="")
                    for j in range((k * 4), np.min((len(idx), (k + 1) * 4))):
                        print("%5d" % (idx[j] + 1), file=fh, end="")
                    print("", file=fh)

        print("ENDMDL", file=fh)
    print("END", file=fh)

    fh.close()


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
