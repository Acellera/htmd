# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import numpy
import numpy as np
import logging
import collections
import re

logger = logging.getLogger(__name__)

# If PDB spec says "COLUMNS 18-20" this means line[17:20]
Pair = collections.namedtuple('Atom', 'resname name')


class PDBParser:
    # The following is taken from MDAnalysis to align atom names based on the strict PDB formatting.
    # Some programs like AutoDock require this strict formatting to guess atom types.

    # These attributes are used to deduce how to format the atom name.
    ions = ('FE', 'AS', 'ZN', 'MG', 'MN', 'CO', 'BR',
            'CU', 'TA', 'MO', 'AL', 'BE', 'SE', 'PT',
            'EU', 'NI', 'IR', 'RH', 'AU', 'GD', 'RU')
    # Mercurial can be confused for hydrogen gamma. Yet, mercurial is
    # rather rare in the PDB. Here are all the residues that contain
    # mercurial.
    special_hg = ('CMH', 'EMC', 'MBO', 'MMC', 'HGB', 'BE7', 'PMB')
    # Chloride can be confused for a carbon. Here are the residues that
    # contain chloride.
    special_cl = ('0QE', 'CPT', 'DCE', 'EAA', 'IMN', 'OCZ', 'OMY', 'OMZ',
                  'UN9', '1N1', '2T8', '393', '3MY', 'BMU', 'CLM', 'CP6',
                  'DB8', 'DIF', 'EFZ', 'LUR', 'RDC', 'UCL', 'XMM', 'HLT',
                  'IRE', 'LCP', 'PCI', 'VGH')
    # In these pairs, the atom name is aligned on the first column
    # (column 13).
    include_pairs = (Pair('OEC', 'CA1'),
                     Pair('PLL', 'PD'),
                     Pair('OEX', 'CA1'))
    # In these pairs, the atom name is aligned on the second column
    # (column 14), but other rules would align them on the first column.
    exclude_pairs = (Pair('C14', 'C14'), Pair('C15', 'C15'),
                     Pair('F9F', 'F9F'), Pair('OAN', 'OAN'),
                     Pair('BLM', 'NI'), Pair('BZG', 'CO'),
                     Pair('BZG', 'NI'), Pair('VNL', 'CO1'),
                     Pair('VNL', 'CO2'), Pair('PF5', 'FE1'),
                     Pair('PF5', 'FE2'), Pair('UNL', 'UNL'))

    def __init__(self, filename=None, mode='pdb'):
        self.record = []
        self.serial = []
        self.name = []
        self.altloc = []
        self.resname = []
        self.chain = []
        self.resid = []
        self.insertion = []
        self.coords = None
        self.occupancy = []
        self.beta = []
        self.segid = []
        self.element = []
        self.charge = []
        self.bonds = []
        self.ssbonds = []

        self.box = None

        self.fieldsizes = {'record': 6, 'serial': 5, 'name': 4, 'altloc': 1, 'resname': 4, 'chain': 1, 'resid': 4,
                           'insertion': 1, 'segid': 4, 'element': 2}

        if filename:
            self.readPDB(filename, mode)
            self._readPDBcoords(filename)

    def checkTruncations(self):
        for f in self.fieldsizes:
            if np.any([True if len(x) > self.fieldsizes[f] else False for x in self.__dict__[f].astype('str')]):
                if self.fieldsizes[f] == 1:
                    logger.warning('Field "{}" of PDB overflows. Your data will be truncated to 1 character.'.format(f))
                else:
                    logger.warning(
                        'Field "{}" of PDB overflows. Your data will be truncated to {} characters.'.format(f,
                                                                                                            self.fieldsizes[
                                                                                                                f]))

    def writePDB(self, filename):
        self.checkTruncations()
        fh = open(filename, 'w')

        # TODO FIXME  -- should take box from traj frame
        box = self.box
        if box is not None:
            print(
                "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % (box[0, 0], box[1, 0], box[2, 0], 90, 90, 90),
                file=fh)

        for frame in range(0, self.coords.shape[2]):
            print("MODEL    %5d" % (frame + 1), file=fh)
            # print ( len(self.record) )
            # print ( self.coords[0].shape )
            for i in range(0, len(self.record)):
                name = self._deduce_PDB_atom_name(self.name[i], self.resname[i])

                if self.serial[i] < 100000:
                    ser = str(int(self.serial[i]))
                else:
                    ser = '*****'

                print(
                    "{!s:6.6}{!s:>5.5} {}{!s:>1.1}{!s:4.4}{!s:>1.1}{!s:>4.4}{!s:>1.1}   {}{}{}{}{}      {!s:4.4}{!s:>2.2}  ".format(
                        self.record[i],
                        ser, name, self.altloc[i],
                        self.resname[i], self.chain[i],
                        self.resid[i],
                        self.insertion[i],
                        self.format83(self.coords[i, 0, frame]),
                        self.format83(self.coords[i, 1, frame]),
                        self.format83(self.coords[i, 2, frame]),
                        self.format62(self.occupancy[i]),
                        self.format62(self.beta[i]),
                        self.segid[i],
                        self.element[i]
                    ), file=fh
                )
                # TODO : convert charges to ints if we ever write them
                if i < len(self.record) - 1 and self.segid[i] != self.segid[i + 1]:
                    print("TER", file=fh)

            if self.bonds is not None and len(self.bonds) != 0:
                bondedatoms = np.unique(self.bonds)

                for a in bondedatoms:
                    idx = self.bonds[self.bonds[:, 0] == a, 1]
                    idx = np.unique(np.append(idx, self.bonds[self.bonds[:, 1] == a, 0]))
                    # I need to support multi-line printing of atoms with more than 4 bonds
                    for k in range(int(np.ceil(len(idx) / 4))):
                        print("CONECT%5d" % (a + 1), file=fh, end="")
                        for j in range((k * 4), np.min((len(idx), (k + 1) * 4))):
                            print("%5d" % (idx[j] + 1), file=fh, end="")
                        print("", file=fh)

            print("ENDMDL", file=fh)
        print("END", file=fh)

        fh.close()

    def _readPDBcoords(self, filename):
        if self.coords is None:
            raise NameError("Can't read coords until a full PDB had been read")
        fh = open(filename, 'r')

        natoms = len(self.coords)

        coords = np.zeros((natoms, 3, 1))
        first = 1
        n = 0
        global_line_counter = 0
        for line in fh:
            global_line_counter += 1
            record_type = line[0:6]

            if record_type == "MODEL ":
                if first == 0:
                    coords = numpy.append(coords, np.zeros((natoms, 3, 1)), axis=2)
                first = 0
            if record_type == "ENDMDL":
                if n != natoms: raise ValueError("Mismatch in number of atoms in PDB frame " + str(f))
                n = 0
            if record_type == "ATOM  " or record_type == "HETATM":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except:
                    # Should we allow parsing to continue in permissive mode?
                    # If so, what coordinates should we default to?  Easier to abort!
                    raise NameError("Invalid or missing coordinate(s) at line %i." % global_line_counter)
                f = coords.shape[2] - 1
                coords[n, 0, f] = x
                coords[n, 1, f] = y
                coords[n, 2, f] = z
                n = n + 1

        self.coords = coords
        fh.close()

    def readPDB(self, filename, mode='pdb'):
        chargeregex = re.compile('\d[\+\-]')

        self.box = numpy.zeros((3, 1))
        fh = open(filename, 'r')
        global_line_counter = 0
        currter = 0
        tergroups = []
        for line in fh:
            global_line_counter += 1
            record_type = line[0:6]
            if record_type == "ATOM  " or record_type == "HETATM":
                self.record.append(record_type.strip())

                fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                # if len(split_list) != 1:
                #    name = fullname
                # else:
                name = split_list[0].strip()
                altloc = line[16].strip()
                resname = line[17:21].strip()
                chainid = line[21].strip()

                try:
                    serial_number = int(line[6:11])
                except:
                    serial_number = 0

                resseq = int(line[22:26].split()[0])  # sequence identifier
                icode = line[26].strip()  # insertion code

                if record_type == "HETATM":  # hetero atom flag
                    if resname == "HOH" or resname == "WAT":
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H"
                else:
                    hetero_flag = " "

                residue_id = resseq  # (hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except:
                    # Should we allow parsing to continue in permissive mode?
                    # If so, what coordinates should we default to?  Easier to abort!
                    raise NameError("Invalid or missing coordinate(s) at line %i." % global_line_counter)

                coord = [x, y, z]
                # occupancy & B factor
                try:
                    occupancy = float(line[54:60])
                except:
                    occupancy = 0

                try:
                    bfactor = float(line[60:66])
                except:
                    bfactor = 0.0  # The PDB use a default of zero if the data is missing

                charge = ''
                if mode == 'pdb':
                    segid = line[72:76].strip()
                    element = line[76:78].strip()
                    charge = line[78:80].strip()
                elif mode == 'pdbqt':
                    segid = ''
                    charge = line[70:76].strip()
                    element = line[77:79].strip()

                if len(charge) != 0:
                    # PDB supports charges of style: 1+ 2- Need to swap the order to convert to float
                    if len(charge) == 2 and chargeregex.match(charge):
                        charge = charge[1] + charge[0]
                    try:
                        charge = float(charge)
                    except:
                        logger.warning('Could not convert charge "{}" to float'.format(charge))
                        charge = 0
                else:
                    charge = 0

                self.name.append(name)
                self.serial.append(serial_number)
                self.altloc.append(altloc)
                self.resname.append(resname)
                self.insertion.append(icode)
                self.chain.append(chainid)
                self.resid.append(residue_id)
                self.occupancy.append(occupancy)
                self.beta.append(bfactor)
                self.segid.append(segid)
                self.element.append(element)
                self.charge.append(charge)
                if self.coords is None:
                    self.coords = [coord]
                else:
                    self.coords.append(coord)
                tergroups.append(str(currter))

            elif record_type == "CONECT":
                # If PDB spec says "COLUMNS 18-20" this means line[17:20]
                a1 = a2 = a3 = a4 = a5 = None
                try:
                    a1 = self.serial.index(int(line[6:11]))
                    a2 = self.serial.index(int(line[11:16]))
                    a3 = self.serial.index(int(line[16:21]))
                    a4 = self.serial.index(int(line[21:26]))
                    a5 = self.serial.index(int(line[26:31]))
                except:
                    pass
                if a1:
                    if a2:
                        self.bonds.append([a1, a2])
                    if a3:
                        self.bonds.append([a1, a3])
                    if a4:
                        self.bonds.append([a1, a4])
                    if a5:
                        self.bonds.append([a1, a5])
            elif record_type == "SSBOND":
                # If PDB spec says "COLUMNS 18-20" this means line[17:20]
                a1 = a2 = a3 = a4 = a5 = None
                try:
                    a1 = self.serial.index(int(line[6:11]))
                    a2 = self.serial.index(int(line[11:16]))
                    a3 = self.serial.index(int(line[16:21]))
                    a4 = self.serial.index(int(line[21:26]))
                    a5 = self.serial.index(int(line[26:31]))
                except:
                    pass
                if a1:
                    if a2:
                        self.ssbonds.append([a1, a2])
                    if a3:
                        self.ssbonds.append([a1, a3])
                    if a4:
                        self.ssbonds.append([a1, a4])
                    if a5:
                        self.ssbonds.append([a1, a5])

            elif record_type == "CRYST1":
                a = b = c = 0.
                try:
                    a = float(line[6:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                except:
                    pass
                    alpha = beta = gamma = 90.
                try:
                    alpha = float(line[33:40])
                    beta = float(line[40:47])
                    gamma = float(line[47:54])
                except:
                    pass

                self.box[0, 0] = a
                self.box[1, 0] = b
                self.box[2, 0] = c

            elif record_type[0:3] == "TER":
                currter += 1

            elif record_type == "END" or record_type == "ENDMDL":
                break

        fh.close()

        # If no segids were read I can use TER groups to set segids
        if all(x == '' for x in self.segid) and currter != 0:
            self.segid = tergroups

    def format83(self, f):
        """Format a single float into a string of width 8, with ideally 3 decimal
        places of precision. If the number is a little too large, we can
        gracefully degrade the precision by lopping off some of the decimal
        places. If it's much too large, we throw a NameError"""
        if -999.999 < f < 9999.999:
            return '%8.3f' % f
        if -9999999 < f < 99999999:
            return ('%8.3f' % f)[:8]
        raise NameError('coordinate "%s" could not be represented '
                        'in a width-8 field' % f)

    def format62(self, f):
        """Format a single float into a string of width 8, with ideally 3 decimal
        places of precision. If the number is a little too large, we can
        gracefully degrade the precision by lopping off some of the decimal
        places. If it's much too large, we throw a NameError"""
        if -9.999 < f < 99.999:
            return '%6.2f' % f
        if -99999 < f < 999999:
            return ('%6.2f' % f)[:6]
        raise NameError('coordinate "%s" could not be represented '
                        'in a width-6 field' % f)

    def _deduce_PDB_atom_name(self, name, resname):
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
               or name[:2] in self.ions
               or name == 'UNK'
               or (resname in self.special_hg and name[:2] == 'HG')
               or (resname in self.special_cl and name[:2] == 'CL')
               or Pair(resname, name) in self.include_pairs)
              and Pair(resname, name) not in self.exclude_pairs):
            return '{:<4}'.format(name)
        return ' {:<3}'.format(name)

# readPDB('test.pdb')
