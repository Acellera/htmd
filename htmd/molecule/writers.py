# (c) 2015-2018 Acellera Ltd http://www.acellera.com
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


def _getPDBElement(name, element, lowersecond=True):
    """
    Given a PDB atom name of 4 characters (including spaces), get the element
    """
    import re
    regH_old = re.compile('H.[123][123]') # Matches i.e. HE13
    regH_inv = re.compile('[123]H')  # Matches i.e. 2H
    element_backup = element.strip()
    if not element.isalpha():
        element = name[0:2].strip()
        if element and element[0].isdigit():
            if element_backup:
                element = element_backup
            else:
                element = name[1]
        if element and len(element) > 1 and element[1].isdigit():
            if element_backup:
                element = element_backup
            else:
                element = name[0]
    if element:
        element = element.strip()
    if regH_old.match(name.strip()) or regH_inv.match(name.strip()):
        element = 'H'
    if len(element) == 2:
        if lowersecond:
            element = element[0] + element[1].lower()
        else:
            element = element[0] + element[1]
    return element


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


def PDBwrite(mol, filename, frames=None, writebonds=True):
    if frames is None:
        frames = mol.frame
    frames = ensurelist(frames)

    checkTruncations(mol)
    coords = np.atleast_3d(mol.coords[:, :, frames])
    numFrames = coords.shape[2]
    nAtoms = coords.shape[0]

    serial = np.arange(1, np.size(coords, 0) + 1).astype(object)
    serial[serial > 99999] = '*****'
    serial = serial.astype('U5')

    if nAtoms > 0:
        if coords.max() >= 1E8 or coords.min() <= -1E7:
            raise RuntimeError('Cannot write PDB coordinates with values smaller than -1E7 or larger than 1E8')
        if mol.occupancy.max() >= 1E6 or mol.occupancy.min() <= -1E5:
            raise RuntimeError('Cannot write PDB occupancy with values smaller than -1E5 or larger than 1E6')
        if mol.beta.max() >= 1E6 or mol.beta.min() <= -1E5:
            raise RuntimeError('Cannot write PDB beta/temperature with values smaller than -1E5 or larger than 1E6')

    fh = open(filename, 'w')
    box = mol.box[:, frames[0]]
    if box is not None and not np.all(mol.box == 0):
        fh.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n" % (box[0], box[1], box[2], 90, 90, 90))

    for f in range(numFrames):
        fh.write("MODEL    %5d\n" % (frames[f] + 1))
        for i in range(0, len(mol.record)):
            name = _deduce_PDB_atom_name(mol.name[i], mol.resname[i])

            fh.write(
                "{!s:6.6}{!s:>5.5} {}{!s:>1.1}{!s:4.4}{!s:>1.1}{!s:>4.4}{!s:>1.1}   {}{}{}{}{}      {!s:4.4}{!s:>2.2}  \n".format(
                    mol.record[i],
                    serial[i], name, mol.altloc[i],
                    mol.resname[i], mol.chain[i],
                    mol.resid[i],
                    mol.insertion[i],
                    '{:8.3f}'.format(coords[i, 0, f])[:8],
                    '{:8.3f}'.format(coords[i, 1, f])[:8],
                    '{:8.3f}'.format(coords[i, 2, f])[:8],
                    '{:6.2f}'.format(mol.occupancy[i])[:6],
                    '{:6.2f}'.format(mol.beta[i])[:6],
                    mol.segid[i],
                    mol.element[i]
                )
            )
            # TODO : convert charges to ints if we ever write them
            if i < len(mol.record) - 1 and mol.segid[i] != mol.segid[i + 1]:
                fh.write("TER\n")

        if writebonds and mol.bonds is not None and len(mol.bonds) != 0:
            bondedatoms = np.unique(mol.bonds)
            bondedatoms = bondedatoms[bondedatoms < 99998]  # Don't print bonds over 99999 as it overflows the field

            for a in bondedatoms:
                partners = mol.bonds[mol.bonds[:, 0] == a, 1]
                partners = np.unique(np.append(partners, mol.bonds[mol.bonds[:, 1] == a, 0]))
                partners = partners[partners < 99998] + 1  # Don't print bonds over 99999 as it overflows the field
                # I need to support multi-line printing of atoms with more than 4 bonds
                while len(partners) >= 3:  # Write bonds as long as they are more than 3 in fast more
                    fh.write("CONECT%5d%5d%5d%5d\n" % (a + 1, partners[0], partners[1], partners[2]))
                    partners = partners[3:]
                if len(partners) > 0:  # Write the rest of the bonds
                    line = "CONECT%5d" % (a + 1)
                    for p in partners:
                        line = "%s%5d" % (line, p)
                    fh.write(line)
                    fh.write('\n')

        fh.write("ENDMDL\n")
    fh.write("END\n")

    fh.close()


def XTCwrite(mol, filename):
    coords = mol.coords
    box = mol.box
    time = mol.time
    step = mol.step
    numFrames = mol.numFrames

    if np.size(box, 1) != numFrames:  # Box should have as many frames as trajectory
        box = np.tile(box, (1, numFrames))

    nframes = np.size(coords, 2)
    if np.size(time) != nframes:
        time = np.zeros(nframes)
    if np.size(step) != nframes:
        step = np.zeros(nframes, dtype=int)

    if os.path.isfile(filename):
        os.unlink(filename)

    box = box.astype(np.float32) * 0.1
    step = step.astype(np.int32)
    time = time.astype(np.float32)
    coords = coords.astype(np.float32) * 0.1  # Convert from A to nm
    if not box.flags['C_CONTIGUOUS']:
        box = np.ascontiguousarray(box)
    if not step.flags['C_CONTIGUOUS']:
        step = np.ascontiguousarray(step)
    if not time.flags['C_CONTIGUOUS']:
        time = np.ascontiguousarray(time)
    if not coords.flags['C_CONTIGUOUS']:
        coords = np.ascontiguousarray(coords)

    lib = xtc_lib()
    natoms = ct.c_int(coords.shape[0])
    nframes = ct.c_int(coords.shape[2])

    cstep = step.ctypes.data_as(ct.POINTER(ct.c_int))
    ctime = time.ctypes.data_as(ct.POINTER(ct.c_float))
    cbox = box.ctypes.data_as(ct.POINTER(ct.c_float))
    ccoords = coords.ctypes.data_as(ct.POINTER(ct.c_float))
    lib['libxtc'].xtc_write(
        ct.c_char_p(filename.encode("ascii")),
        natoms,
        nframes,
        cstep,
        ctime,
        ccoords,
        cbox)


def BINCOORwrite(mol, filename):
    import struct
    natoms = np.array([mol.numAtoms])
    f = open(filename, 'wb')

    dat = mol.coords[:, :, mol.frame].copy()
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
        print("{!s:>8.8} {!s:4.4} {!s:5.5}{!s:4.4} {!s:4.4} {!s:6.6} {!s:2.2} {:10.6}  {:8.6}  {!s:>10.10}".format(
            int(m.serial[i]),
            segid,
            str(m.resid[i]) + str(m.insertion[i]),
            (m.resname[i]),
            m.name[i],
            atomtype,
            "",  # m.element[i], # NAMD barfs if this is set
            charge,
            mass,
            0
        ), file=f)
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
    # According ParmEd, CHARMM PSF has to have an extra blank line after NNB
    # https: // github.com / ParmEd / ParmEd / blob / master / parmed / charmm / psf.py#L151
    print("%10d !NNB\n\n" % (0), file=f)
    print("%10d %10d !NGRP\n" % (0, 0), file=f)
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
    uqresname = np.unique(mol.resname)
    if len(uqresname) > 1:
        raise RuntimeError('MOL2 file can only be written for a single residue. We detected {} resnames in the '
                           'Molecule.'.format(len(uqresname)))
    if len(uqresname[0]) == 0:
        raise RuntimeError('MOL2 file can only be written if a resname is defined for the Molecule. Currently the '
                           'resname is empty.')

    with open(filename, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("    {}\n".format(uqresname[0]))
        unique_bonds = [list(t) for t in set(map(tuple, [sorted(x) for x in mol.bonds]))]
        unique_bonds = np.array(sorted(unique_bonds, key=lambda x: (x[0], x[1])))
        f.write("%5d %5d %5d %5d %5d\n" % (mol.numAtoms, unique_bonds.shape[0], 0, 0, 0))
        f.write("SMALL\nUSER_CHARGES\n\n\n")
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
        f.write("@<TRIPOS>ATOM\n")
        for i in range(mol.coords.shape[0]):
            f.write('{:7d} {:8s} {:9.4f} {:9.4f} {:9.4f} {:8s} '.format(i + 1, mol.name[i], mol.coords[i, 0, mol.frame],
                                                                      mol.coords[i, 1, mol.frame],
                                                                      mol.coords[i, 2, mol.frame],
                                                                      mol.atomtype[i] if mol.atomtype[i] != ''
                                                                      else mol.element[i])
                    )
            if isinstance(mol.resid[i], numbers.Integral):
                f.write('{:3d} '.format(mol.resid[i]))
                if mol.resname[i] != '':
                    f.write('{:4s} '.format(mol.resname[i]))
                    if isinstance(mol.charge[i], numbers.Real):
                        f.write('{:12.6f}'.format(mol.charge[i]))
            f.write('\n')
        f.write("@<TRIPOS>BOND\n")
        for i in range(unique_bonds.shape[0]):
            bt = 'un'
            if len(mol.bondtype) > 0:
                idx = (mol.bonds[:, 0] == unique_bonds[i, 0]) & (mol.bonds[:, 1] == unique_bonds[i, 1])
                idx |= (mol.bonds[:, 0] == unique_bonds[i, 1]) & (mol.bonds[:, 1] == unique_bonds[i, 0])
                tmp = np.unique(mol.bondtype[idx])
                assert len(tmp) == 1, 'There should only exist one bond type for atoms {} {}'.format(unique_bonds[i, 0], unique_bonds[i, 1])
                bt = tmp[0]
            f.write("{:6d} {:4d} {:4d} {}\n".format(i + 1, unique_bonds[i, 0] + 1, unique_bonds[i, 1] + 1, bt))


def SDFwrite(mol, filename):
    import datetime
    mol2bonds = {'1': 1, '2': 2, '3': 3, 'ar': 4}
    with open(filename, 'w') as fh:
        fh.write('{}\n'.format(mol.viewname))
        currtime = datetime.datetime.now().strftime('%m%d%y%H%M')
        fh.write('  HTMD    {}3D  {scaling:10s}{energy:12s}{registry:6s}\n'.format(currtime, scaling='', energy='', registry=''))
        fh.write('{comments}\n'.format(comments=''))
        fh.write('  0  0  0     0  0{:12s}999 V3000\n'.format(''))
        fh.write('M  V30 BEGIN CTAB\n')
        fh.write('M  V30 COUNTS {na} {nb} {nsg} {n3d} {chiral}\n'.format(na=mol.numAtoms, nb=mol.bonds.shape[0], nsg=0, n3d=0, chiral=0))
        coor = mol.coords[:, :, mol.frame]
        fh.write('M  V30 BEGIN ATOM\n')
        for i in range(mol.numAtoms):
            atype = mol.atomtype[i]
            if atype == '':
                atype = mol.element[i]
            if atype == '':
                atype = mol.name[i]
            fh.write('M  V30 {index} {type} {x:9.4f} {y:9.4f} {z:9.4f} {aamap} {charge}\n'.format(index=i, type=atype, x=coor[i, 0], y=coor[i, 1], z=coor[i, 2], aamap=0, charge=mol.charge[i]))
        fh.write('M  V30 END ATOM\n')

        fh.write('M  V30 BEGIN BOND\n')
        for i in range(mol.bonds.shape[0]):

            if mol.bondtype[i] != '':
                btype = mol2bonds[mol.bondtype[i].strip()]
            else:
                btype = 1
            fh.write('M  V30 {index} {type} {atom1} {atom2}\n'.format(index=i, type=btype, atom1=mol.bonds[i, 0], atom2=mol.bonds[i, 1]))
        fh.write('M  V30 END BOND\n')
        fh.write('M  V30 END CTAB\n')
        fh.write('M  END\n')


def GROwrite(mol, filename):
    import pandas as pd
    from collections import OrderedDict
    coor = mol.coords[:, :, mol.frame] / 10  # Convert to nm
    box = mol.box[:, mol.frame] / 10  # Convert to nm
    datadict = OrderedDict([('resid', mol.resid), ('resname', mol.resname), ('name', mol.name), ('serial', mol.serial),
                            ('posx', coor[:, 0]), ('posy', coor[:, 1]), ('posz', coor[:, 2])])
    a = pd.DataFrame(data=datadict)
    with open(filename, 'wb') as fh:
        if mol.fstep is not None:
            fh.write(b'Generated with HTMD, t= %f\n' % (mol.fstep * 1000))
        else:
            fh.write(b'Generated with HTMD\n')
        fh.write(b'%5d\n' % mol.numAtoms)
        np.savetxt(fh, a.values, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f')
        fh.write(b'%f %f %f 0 0 0 0 0 0' % (box[0], box[1], box[2]))


# Taken from trajectory.py Trajectory()._savers() method of MDtraj
_MDTRAJ_TOPOLOGY_SAVERS = ('pdb', 'pdb.gz', 'xyz', 'xyz.gz')

_MDTRAJ_TRAJECTORY_SAVERS = ('xtc', 'trr', 'dcd', 'h5', 'binpos', 'nc', 'netcdf', 'ncrst', 'crd', 'mdcrd', 'ncdf',
                             'lh5', 'lammpstrj', 'gro', 'rst7', 'tng')

_MDTRAJ_SAVERS = _MDTRAJ_TRAJECTORY_SAVERS + _MDTRAJ_TOPOLOGY_SAVERS


def MDTRAJwrite(mol, filename):
    try:
        import mdtraj as md
        from htmd.util import tempname
        ext = os.path.splitext(filename)[1][1:]
        if ext == 'gz':
            pieces = filename.split('.')
            ext = '{}.{}'.format(pieces[-2], pieces[-1])

        if ext in _MDTRAJ_TOPOLOGY_SAVERS:
            tmppdb = tempname(suffix='.pdb')
            mol.write(tmppdb)
            traj = md.load(tmppdb)
            os.remove(tmppdb)
        elif ext in _MDTRAJ_TRAJECTORY_SAVERS:
            tmppdb = tempname(suffix='.pdb')
            tmpxtc = tempname(suffix='.xtc')
            mol.write(tmppdb)
            mol.write(tmpxtc)
            traj = md.load(tmpxtc, top=tmppdb)
            os.remove(tmppdb)
            os.remove(tmpxtc)
        else:
            raise ValueError('Unknown file type for file {}'.format(filename))
        # traj.xyz = np.swapaxes(np.swapaxes(self.coords, 1, 2), 0, 1) / 10
        # traj.time = self.time
        # traj.unitcell_lengths = self.box.T / 10
        traj.save(filename)
    except Exception as e:
        raise ValueError('MDtraj reader failed for file {} with error "{}"'.format(filename, e))


_WRITERS = {'psf': PSFwrite,
            'pdb': PDBwrite,
            'mol2': MOL2write,
            'sdf': SDFwrite,
            'xyz': XYZwrite,
            'gro': GROwrite,
            'coor': BINCOORwrite,
            'xtc': XTCwrite}


for ext in _MDTRAJ_SAVERS:
    if ext not in _WRITERS:
        _WRITERS[ext] = MDTRAJwrite


if __name__ == '__main__':
    from htmd.home import home
    from htmd.molecule.molecule import Molecule, mol_equal
    from htmd.util import tempname
    import numpy as np
    import os
    testfolder = home(dataDir='metricdistance')
    mol = Molecule(os.path.join(testfolder, 'filtered.pdb'))
    mol.coords = np.tile(mol.coords, (1, 1, 2))
    mol.filter('protein and resid 1 to 20')
    mol.boxangles = np.ones((3, 2), dtype=np.float32) * 90
    mol.box = np.ones((3, 2), dtype=np.float32) * 15
    mol.step = np.arange(2)
    mol.time = np.arange(2) * 1E5
    mol.fileloc = [mol.fileloc[0], mol.fileloc[0]]

    for ext in _WRITERS:
        tmp = tempname(suffix='.'+ext)
        if ext == 'mol2':
            mol.write(tmp, sel='resid 1')
        else:
            mol.write(tmp)
        print('Can write {} files'.format(ext))


    # from difflib import Differ
    # d = Differ()
    #
    # with open(tmp, 'rb') as f1, open(a, 'rb') as f2:
    #     res = d.compare([x.decode('utf8') for x in f1.readlines()], [x.decode('utf8') for x in f2.readlines()])
    #
    #     try:
    #         next(res)
    #     except StopIteration:
    #         print('Same files')
    #     else:
    #         assert('Error in MDtraj writer.')
    #
    # assert filecmp.cmp(tmp, os.path.join(testfolder, '3PTB.h5'))
