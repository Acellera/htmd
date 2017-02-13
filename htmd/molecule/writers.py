# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import ctypes as ct
import os
from htmd.molecule.support import xtc_lib
from htmd.util import ensurelist
import collections
import logging
import numbers

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
                logger.warning('Field "{}" of PDB overflows. Your data will be truncated to {} characters.'.format(f,
                                                                                                                   fieldsizes[
                                                                                                                       f]))


def PDBwrite(mol, filename, frame):
    frame = ensurelist(frame)

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
    coords = np.atleast_3d(mol.coords[:, :, frame])
    numFrames = coords.shape[2]
    serial = np.arange(1, np.size(coords, 0) + 1)

    fh = open(filename, 'w')
    # TODO FIXME  -- should take box from traj frame
    box = mol.box
    if box is not None and not np.all(mol.box == 0):
        box = np.atleast_2d(np.atleast_2d(box)[:, mol.frame])
        print("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % (box[0, 0], box[0, 1], box[0, 2], 90, 90, 90),
              file=fh)

    for f in frame:
        print("MODEL    %5d" % (f + 1), file=fh)
        for i in range(0, len(mol.record)):
            name = _deduce_PDB_atom_name(mol.name[i], mol.resname[i])

            if serial[i] < 100000:
                ser = str(int(serial[i]))
            else:
                ser = '*****'

            print(
                "{!s:6.6}{!s:>5.5} {}{!s:>1.1}{!s:4.4}{!s:>1.1}{!s:>4.4}{!s:>1.1}   {}{}{}{}{}      {!s:4.4}{!s:>2.2}  ".format(
                    mol.record[i],
                    ser, name, mol.altloc[i],
                    mol.resname[i], mol.chain[i],
                    mol.resid[i],
                    mol.insertion[i],
                    format83(mol.coords[i, 0, f]),
                    format83(mol.coords[i, 1, f]),
                    format83(mol.coords[i, 2, f]),
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
            bondedatoms = bondedatoms[bondedatoms < 99998]  # Don't print bonds over 99999 as it overflows the field

            for a in bondedatoms:
                partners = mol.bonds[mol.bonds[:, 0] == a, 1]
                partners = np.unique(np.append(partners, mol.bonds[mol.bonds[:, 1] == a, 0]))
                partners = partners[partners < 99998] + 1  # Don't print bonds over 99999 as it overflows the field
                # I need to support multi-line printing of atoms with more than 4 bonds
                while len(partners) >= 3:  # Write bonds as long as they are more than 3 in fast more
                    print("CONECT%5d%5d%5d%5d" % (a + 1, partners[0], partners[1], partners[2]), file=fh)
                    partners = partners[3:]
                if len(partners) > 0:  # Write the rest of the bonds
                    line = "CONECT%5d" % (a + 1)
                    for p in partners:
                        line = "%s%5d" % (line, p)
                    print(line, file=fh)

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

    import string
    # Letters to be used for default segids, if free: 0123456789abcd...ABCD..., minus chain symbols already used
    used_segids = set(m.segid)
    segid_alphabet = list(string.digits + string.ascii_letters)
    available_segids = [x for x in segid_alphabet if x not in used_segids]

    f = open(filename, 'w')
    print("PSF\n", file=f)
    print("%8d !NTITLE\n" % (0), file=f)
    print("%8d !NATOM" % (len(m.serial)), file=f)
    for i in range(len(m.serial)):
        charge = 0
        mass = 1
        segid = available_segids[0]
        atomtype = ""
        if (m.masses is not None) and (i < len(m.masses)):
            mass = m.masses[i]
        if (m.charge is not None) and (i < len(m.charge)):
            charge = m.charge[i]
        if (m.atomtype is not None) and (i < len(m.atomtype)):
            atomtype = m.atomtype[i]
        if (m.segid is not None) and (i < len(m.segid)) and m.segid[i] != '':
            segid = m.segid[i]
        elif m.segid[i] == '' and m.chain[i] != '':
            segid = m.chain[i]
        print("%8d %-4s %-5s%-4s %-4s %-6s %-2s %10.6f  %8.6f  %10d" %
              (int(m.serial[i]),
               segid,
               str(m.resid[i]) + str(m.insertion[i]),
               (m.resname[i]),
               m.name[i],
               atomtype,
               "",  # m.element[i], # NAMD barfs if this is set
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
        print("%10d%10d%10d" % (m.angles[i, 0] + 1, m.angles[i, 1] + 1, m.angles[i, 2] + 1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NPHI: dihedrals" % (m.dihedrals.shape[0]), file=f)
    for i in range(m.dihedrals.shape[0]):
        if i and not (i % 2):
            print("", file=f)
        print("%10d%10d%10d%10d" % (
        m.dihedrals[i, 0] + 1, m.dihedrals[i, 1] + 1, m.dihedrals[i, 2] + 1, m.dihedrals[i, 3] + 1), file=f, end="")

    print("\n\n", file=f)
    print("%10d !NIMPHI: impropers" % (m.impropers.shape[0]), file=f)
    for i in range(m.impropers.shape[0]):
        if i and not (i % 2):
            print("", file=f)
        print("%10d%10d%10d%10d" % (
        m.impropers[i, 0] + 1, m.impropers[i, 1] + 1, m.impropers[i, 2] + 1, m.impropers[i, 3] + 1), file=f, end="")

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
        print("%s   %f   %f    %f" % (
        e, src.coords[i, 0, src.frame], src.coords[i, 1, src.frame], src.coords[i, 2, src.frame]), file=fh)
    fh.close()


def MOL2write(mol, filename):
    with open(filename, "w") as f:
        print("@<TRIPOS>MOLECULE", file=f)
        print("    MOL", file=f)
        unique_bonds = np.array([list(t) for t in set(map(tuple, [sorted(x) for x in mol.bonds]))])
        print("%5d %5d %5d %5d %5d" % (mol.numAtoms, unique_bonds.shape[0], 0, 0, 0), file=f)
        print("SMALL\nUSER_CHARGES\n\n", file=f)
        '''
        @<TRIPOS>ATOM
        Each data record associated with this RTI consists of a single data line. This
        data line contains all the information necessary to reconstruct one atom
        contained within the molecule. The atom ID numbers associated with the atoms
        in the molecule will be assigned sequentially when the .mol2 file is read into
        SYBYL.
        Format:
        atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]

        • atom_id (integer) = the ID number of the atom at the time the file was
        created. This is provided for reference only and is not used when the
        .mol2 file is read into SYBYL.
        • atom_name (string) = the name of the atom.
        • x (real) = the x coordinate of the atom.
        • y (real) = the y coordinate of the atom.
        • z (real) = the z coordinate of the atom.
        • atom_type (string) = the SYBYL atom type for the atom.
        • subst_id (integer) = the ID number of the substructure containing the
        atom.
        • subst_name (string) = the name of the substructure containing the atom.
        • charge (real) = the charge associated with the atom.
        • status_bit (string) = the internal SYBYL status bits associated with the
        atom. These should never be set by the user. Valid status bits are
        DSPMOD, TYPECOL, CAP, BACKBONE, DICT, ESSENTIAL, WATER and
        DIRECT.
        '''
        print("@<TRIPOS>ATOM", file=f)
        for i in range(mol.coords.shape[0]):
            print('{:7d} {:8s} {:9.4f} {:9.4f} {:9.4f} {:8s} '.format(i + 1, mol.name[i], mol.coords[i, 0, mol.frame],
                                                                      mol.coords[i, 1, mol.frame],
                                                                      mol.coords[i, 2, mol.frame],
                                                                      mol.element[i]),  # TODO: implement SYBYL atom types
                  end='',
                  file=f)
            if isinstance(mol.resid[i], numbers.Integral):
                print('{:3d} '.format(mol.resid[i]), end='', file=f)
                if mol.resname[i] != '':
                    print('{:4s} '.format(mol.resname[i]), end='', file=f)
                    if isinstance(mol.charge[i], numbers.Real):
                        print('{:12.6f}'.format(mol.charge[i]), end='', file=f)
            print('', file=f)
        print("@<TRIPOS>BOND", file=f)
        for i in range(unique_bonds.shape[0]):
            print("%6d %4d %4d un" % (i + 1, unique_bonds[i, 0] + 1, unique_bonds[i, 1] + 1), file=f) # TODO: implement SYBYL bond types
        print("", file=f)


def GROwrite(mol, filename):
    import pandas as pd
    from collections import OrderedDict
    coor = mol.coords[:, :, mol.frame] / 10  # Convert to nm
    box = mol.box[:, mol.frame] / 10  # Convert to nm
    datadict = OrderedDict([('resid', mol.resid), ('resname', mol.resname), ('name', mol.name), ('serial', mol.serial),
                            ('posx', coor[:, 0]), ('posy', coor[:, 1]), ('posz', coor[:, 2])])
    a = pd.DataFrame(data=datadict)
    with open(filename, 'wb') as fh:
        fh.write(b'Generated with HTMD, t= %f\n' % (mol.fstep * 1000))
        fh.write(b'%5d\n' % mol.numAtoms)
        np.savetxt(fh, a.values, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f')
        fh.write(b'%f %f %f 0 0 0 0 0 0' % (box[0], box[1], box[2]))
