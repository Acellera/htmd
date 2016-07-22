import ctypes as ct
import numpy as np
from htmd.molecule.support import pack_double_buffer, pack_int_buffer, pack_string_buffer, pack_ulong_buffer, xtc_lib
import logging
logger = logging.getLogger(__name__)


class Topology:
    def __init__(self, pandasdata=None):
        self.record = []
        self.serial = []
        self.name = []
        self.altloc = []
        self.element = []
        self.resname = []
        self.chain = []
        self.resid = []
        self.insertion = []
        self.occupancy = []
        self.beta = []
        self.segid = []
        self.bonds = []
        self.charge = []
        self.masses = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.atomtype = []

        if pandasdata is not None:
            for field in self.__dict__:
                fielddat = pandasdata.get(field)
                if fielddat is not None and not np.all(fielddat.isnull()):
                    if fielddat.dtype == object:  # If not all are NaN replace NaNs with default values
                        pandasdata.loc[fielddat.isnull(), field] = ''
                    else:
                        pandasdata.loc[fielddat.isnull(), field] = 0
                    self.__dict__[field] = fielddat.tolist()

    @property
    def atominfo(self):
        return ['record', 'serial', 'name', 'altloc', 'element', 'resname', 'chain', 'resid', 'insertion',
                     'occupancy', 'beta', 'segid', 'charge', 'masses', 'atomtype']


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
    #default_name
    #  7196
    #  -7.0046035  10.4479194  20.8320000  -7.3970000   9.4310000  20.8320000
    #  -7.0486898   8.9066002  21.7218220  -7.0486899   8.9065995  19.9421780

    f = open(filename, 'r')
    coords = []

    fieldlen = 12
    k = 0
    for line in f:
        k += 1
        if k < 3:  # skip first 2 lines
            continue

        coords += [float(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                   if len(line[i:i + fieldlen].strip()) != 0]

    return [coords[i:i+3] for i in range(0, len(coords), 3)]


def XYZread(filename):
    topo = Topology()
    coords = []

    f = open(filename, 'r')
    natoms = int(f.readline().split()[0])
    f.readline()
    for i in range(natoms):
        s = f.readline().split()
        topo.record.append('HETATM')
        topo.serial.append(i + 1)
        topo.element.append(s[0])
        topo.name.append(s[0])
        topo.resname.append('MOL')
        coords.append(s[1:4])

    return topo, coords


def GJFread(filename):
    # $rungauss
    # %chk=ts_rhf
    # %mem=2000000
    # #T RHF/6-31G(d) TEST
    #
    # C9H8O4
    #
    # 0,1
    # C1,2.23927,-0.379063,0.262961
    # C2,0.842418,1.92307,-0.424949
    # C3,2.87093,0.845574,0.272238

    #import re
    #regex = re.compile('(\w+),([-+]?[0-9]*\.?[0-9]+),([-+]?[0-9]*\.?[0-9]+),([-+]?[0-9]*\.?[0-9]+)')

    topo = Topology()
    coords = []

    f = open(filename, "r")
    for line in f:
        pieces = line.split(sep=',')
        if len(pieces) == 4 and not line.startswith('$') and not line.startswith('%') and not line.startswith('#'):
            topo.record.append('HETATM')
            topo.element.append(pieces[0])
            topo.name.append(pieces[0])
            topo.resname = 'MOL'
            coords.append([float(s) for s in pieces[1:4]])
    topo.serial = range(len(topo.record))

    return topo, coords


def MOL2read(filename):
    import re

    topo = Topology()
    coords = []

    f = open(filename, "r")
    l = f.readlines()
    f.close()
    start = None
    end = None
    bond = None
    for i in range(len(l)):
        if l[i].startswith("@<TRIPOS>ATOM"): start = i + 1
        if l[i].startswith("@<TRIPOS>BOND"):
            end = i - 1
            bond = i + 1

    if not start or not end:
        raise ValueError("File cannot be read")

    natoms = end - start + 1
    for i in range(natoms):
        s = l[i + start].strip().split()
        topo.record.append("HETATM")
        topo.serial.append(int(s[0]))
        topo.element.append(re.sub("[0123456789]*", "", s[1]))
        topo.name.append(s[1])
        coords.append([float(x) for x in s[2:5]])
        topo.charge.append(float(s[8]))
        topo.atomtype.append(s[5])
        topo.resname.append(s[6])
    if bond:
        for i in range(bond, len(l)):
            b = l[i].split()
            if len(b) != 4:
                break
            topo.bonds.append([int(b[1]) - 1, int(b[2]) - 1])
    return topo, coords


def MAEread(fname):
    """ Reads maestro files.

    Parameters
    ----------
    fname : str
        .mae file

    Returns
    -------
    topo : Topology
    coords : list of lists
    """
    section = None
    section_desc = False
    section_data = False

    topo = Topology()
    coords = []
    heteros = []

    import csv
    with open(fname, newline='') as fp:
        reader = csv.reader(fp, delimiter=' ', quotechar='"', skipinitialspace=True)
        for row in reader:
            if len(row) == 0:
                continue

            if row[0].startswith('m_atom'):
                section = 'atoms'
                section_desc = True
                section_cols = []
            elif row[0].startswith('m_bond'):
                section = 'bonds'
                section_desc = True
                section_cols = []
            elif row[0].startswith('m_PDB_het_residues'):
                section = 'hetresidues'
                section_desc = True
                section_cols = []
            elif section_desc and row[0] == ':::':  # Once the section description has finished create a map from names to columns
                section_dict = dict(zip(section_cols, range(len(section_cols))))
                section_desc = False
                section_data = True
            elif section_data and (row[0] == ':::' or row[0] == '}'):
                section_data = False
            else:  # It's actual data
                if section_desc:
                    section_cols.append(row[0])

                # Reading the data of the atoms section
                if section == 'atoms' and section_data:
                    topo.record.append('ATOM')
                    row = np.array(row)
                    if len(row) != len(section_dict) +1:  # TODO: fix the reader
                        raise RuntimeError('{} has {} fields in the m_atom section description, but {} fields in the '
                                           'section data. Please check for missing fields in the mae file.'
                                           .format(fname, len(section_dict), len(row)))
                    row[row == '<>'] = 0
                    if 'i_pdb_PDB_serial' in section_dict:
                        topo.serial.append(row[section_dict['i_pdb_PDB_serial']])
                    if 's_m_pdb_atom_name' in section_dict:
                        topo.name.append(row[section_dict['s_m_pdb_atom_name']].strip())
                    if 's_m_pdb_residue_name' in section_dict:
                        topo.resname.append(row[section_dict['s_m_pdb_residue_name']].strip())
                    if 'i_m_residue_number' in section_dict:
                        topo.resid.append(int(row[section_dict['i_m_residue_number']]))
                    if 's_m_chain_name' in section_dict:
                        topo.chain.append(row[section_dict['s_m_chain_name']])
                    if 's_pdb_segment_id' in section_dict:
                        topo.segid.append(row[section_dict['s_pdb_segment_id']])
                    if 'r_m_pdb_occupancy' in section_dict:
                        topo.occupancy.append(float(row[section_dict['r_m_pdb_occupancy']]))
                    if 'r_m_pdb_tfactor' in section_dict:
                        topo.beta.append(float(row[section_dict['r_m_pdb_tfactor']]))
                    if 's_m_insertion_code' in section_dict:
                        topo.insertion.append(row[section_dict['s_m_insertion_code']].strip())
                    if '' in section_dict:
                        topo.element.append('')  # TODO: Read element
                    if '' in section_dict:
                        topo.altloc.append('')  # TODO: Read altloc. Quite complex actually. Won't bother.
                    if 'r_m_x_coord' in section_dict:
                        coords.append(
                            [float(row[section_dict['r_m_x_coord']]), float(row[section_dict['r_m_y_coord']]),
                             float(row[section_dict['r_m_z_coord']])])
                    topo.masses.append(0)

                # Reading the data of the bonds section
                if section == 'bonds' and section_data:
                    topo.bonds.append([int(row[section_dict['i_m_from']]) - 1, int(row[section_dict['i_m_to']]) - 1])  # -1 to conver to 0 indexing

                # Reading the data of the hetero residue section
                if section == 'hetresidues' and section_data:
                    heteros.append(row[section_dict['s_pdb_het_name']].strip())

    for h in heteros:
        topo.record[topo.resname == h] = 'HETATM'
    return topo, coords


# def MAEreadPandas(fname):
#     """ Reads maestro files.
#
#     Parameters
#     ----------
#     fname : str
#         .mae file
#
#     Returns
#     -------
#     topo : Topology
#     coords : list of lists
#     """
#     import io
#     from pandas import read_csv
#
#     section = None
#     section_desc = False
#     section_data = False
#
#     topo = Topology()
#     coords = []
#     data = None
#
#     f = open(fname, 'r')
#     for line in f:
#         if len(line) == 0:
#             continue
#
#         if line.strip().startswith('m_atom'):
#             section = 'atoms'
#             section_desc = True
#             section_cols = []
#         elif line.strip().startswith('m_bond'):
#             section = 'bonds'
#             section_desc = True
#             section_cols = []
#         elif line.strip().startswith('m_PDB_het_residues'):
#             section = 'hetresidues'
#             section_desc = True
#             section_cols = []
#         elif section_desc and line.strip() == ':::':  # Once the section description has finished create a map from names to columns
#             section_desc = False
#             section_data = True
#         elif section_data and (line.strip() == ':::' or line.strip() == '}'):
#             data.seek(0)
#             if section == 'atoms':
#                 atomdata = read_csv(data, delim_whitespace=True, skipinitialspace=True)
#             elif section == 'bonds':
#                 bonddata = read_csv(data, delim_whitespace=True, skipinitialspace=True)
#             elif section == 'hetresidues':
#                 hetdata = read_csv(data, delim_whitespace=True, skipinitialspace=True)
#             data = None
#             section = None
#             section_data = False
#         else:  # Descriptions or data
#             # Reading the section description column names
#             if section_desc:
#                 section_cols.append(line.strip())
#             # Reading the data
#             if section_data and (section == 'atoms' or section == 'bonds' or section =='hetresidues'):
#                 if data is None:
#                     data = io.StringIO()
#                     data.write(' '.join(section_cols) + '\n')
#                 data.write(line)
#
#     from IPython.core.debugger import Tracer
#     Tracer()()
#
#     for h in heteros:
#         topo.record[topo.resname == h] = 'HETATM'
#     return topo, coords

def PDBread(filename, mode='pdb'):
    from pandas import read_fwf
    import io

    """
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    """
    if mode == 'pdb':
        topocolspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 21), (21, 22), (22, 26), (26, 27),
                        (54, 60), (60, 66), (72, 76), (76, 78), (78, 79), (79, 80)]
        toponames = ('record', 'serial', 'name', 'altloc', 'resname', 'chain', 'resid', 'insertion',
                     'occupancy', 'beta', 'segid', 'element', 'charge', 'chargesign')
    elif mode == 'pdbqt':
        topocolspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 21), (21, 22), (22, 26), (26, 27),
                        (54, 60), (60, 66), (70, 76), (77, 79)]
        toponames = ('record', 'serial', 'name', 'altloc', 'resname', 'chain', 'resid', 'insertion',
                     'occupancy', 'beta', 'charge', 'element')
    topodtypes = {
        'record': str,
        'serial': np.int,
        'name': str,
        'altloc': str,
        'resname': str,
        'chain': str,
        'resid': np.int,
        'insertion': str,
        'occupancy': np.float32,
        'beta': np.float32,
        'segid': str,
        'element': str,
        'charge': np.float32,
        'chargesign': str,
    }
    coordcolspecs = [(30, 38), (38, 46), (46, 54)]
    coordnames = ('x', 'y', 'z')

    """
    COLUMNS       DATA  TYPE      FIELD        DEFINITION
    -------------------------------------------------------------------------
     1 -  6        Record name    "CONECT"
     7 - 11        Integer        serial       Atom  serial number
    12 - 16        Integer        serial       Serial number of bonded atom
    17 - 21        Integer        serial       Serial  number of bonded atom
    22 - 26        Integer        serial       Serial number of bonded atom
    27 - 31        Integer        serial       Serial number of bonded atom
    """
    bondcolspecs = [(6, 11), (11, 16), (16, 21), (21, 26), (26, 31)]
    bondnames = ('serial1', 'serial2', 'serial3', 'serial4', 'serial5')

    """
    COLUMNS       DATA  TYPE    FIELD          DEFINITION
    -------------------------------------------------------------
     1 -  6       Record name   "CRYST1"
     7 - 15       Real(9.3)     a              a (Angstroms).
    16 - 24       Real(9.3)     b              b (Angstroms).
    25 - 33       Real(9.3)     c              c (Angstroms).
    34 - 40       Real(7.2)     alpha          alpha (degrees).
    41 - 47       Real(7.2)     beta           beta (degrees).
    48 - 54       Real(7.2)     gamma          gamma (degrees).
    56 - 66       LString       sGroup         Space  group.
    67 - 70       Integer       z              Z value.
    """
    boxcolspecs = [(6, 15), (15, 24), (24, 33)]
    boxnames = ('a', 'b', 'c')

    def concatCoords(coords, coorddata):
        if coorddata.tell() != 0:  # Not empty
            coorddata.seek(0)
            parsedcoor = read_fwf(coorddata, colspecs=coordcolspecs, names=coordnames)
            if coords is None:
                coords = np.zeros((len(parsedcoor), 3, 0), dtype=np.float32)
            currcoords = np.vstack((parsedcoor.x, parsedcoor.y, parsedcoor.z)).T
            coords = np.append(coords, currcoords[:, :, np.newaxis], axis=2)
        return coords

    teridx = []
    currter = 0
    topoend = False

    crystdata = io.StringIO()
    topodata = io.StringIO()
    conectdata = io.StringIO()
    coorddata = io.StringIO()

    coords = None

    f = open(filename, 'r')
    for line in f:
        if line.startswith('CRYST1'):
            crystdata.write(line)
        if line.startswith('ATOM') or line.startswith('HETATM'):
            coorddata.write(line)
        if (line.startswith('ATOM') or line.startswith('HETATM')) and not topoend:
            topodata.write(line)
            teridx.append(str(currter))
        if line.startswith('TER'):
            currter += 1
        if line.startswith('END'):
            topoend = True
        if line.startswith('CONECT'):
            conectdata.write(line)
        if line.startswith('MODEL'):
            coords = concatCoords(coords, coorddata)
            coorddata = io.StringIO()
    crystdata.seek(0)
    topodata.seek(0)
    conectdata.seek(0)

    coords = concatCoords(coords, coorddata)

    parsedbonds = read_fwf(conectdata, colspecs=bondcolspecs, names=bondnames)
    parsedbox = read_fwf(crystdata, colspecs=boxcolspecs, names=boxnames)
    parsedtopo = read_fwf(topodata, colspecs=topocolspecs, names=toponames) #, dtype=topodtypes)
    if not np.all(parsedtopo.chargesign.isnull()):
        parsedtopo.loc[parsedtopo.chargesign == '-', 'charge'] *= -1

    if len(parsedtopo) > 99999:
        logger.warning('Reading PDB file with more than 99999 atoms. Bond information can be wrong.')

    topo = Topology(parsedtopo)

    # TODO: Speed this up. This is the slowest part for large PDB files. From 700ms to 7s
    serials = parsedtopo.serial.as_matrix()
    mapserials = np.ones(np.max(serials)+1) * -1
    mapserials[serials] = list(range(np.max(serials)))
    for i in range(len(parsedbonds)):
        row = parsedbonds.loc[i].tolist()
        for b in range(1, 5):
            if not np.isnan(row[b]):
                topo.bonds.append([int(row[0]), int(row[b])])
    topo.bonds = np.array(topo.bonds, dtype=np.uint32)
    topo.bonds[:] = mapserials[topo.bonds[:]]

    if len(topo.segid) == 0 and currter != 0:  # If no segid was read, use the TER rows to define segments
        topo.segid = teridx
    return topo, coords


def BINCOORread(filename):
    import struct
    f = open(filename, 'rb')
    dat = f.read(4)
    fmt = 'i'
    natoms = struct.unpack(fmt, dat)[0]
    dat = f.read(natoms * 3 * 8)
    fmt = 'd' * (natoms * 3)
    coords = struct.unpack(fmt, dat)
    coords = np.array(coords).reshape((natoms, 3, 1))
    f.close()
    return coords


def PRMTOPread(filename):
    f = open(filename, 'r')
    topo = Topology()
    uqresnames = []
    residx = []
    bondsidx = []

    section = None
    for line in f:
        if line.startswith('%FLAG POINTERS'):
            section = 'pointers'
        elif line.startswith('%FLAG ATOM_NAME'):
            section = 'names'
        elif line.startswith('%FLAG CHARGE'):
            section = 'charges'
        elif line.startswith('%FLAG MASS'):
            section = 'masses'
        elif line.startswith('%FLAG ATOM_TYPE_INDEX'):
            section = 'type'
        elif line.startswith('%FLAG RESIDUE_LABEL'):
            section = 'resname'
        elif line.startswith('%FLAG RESIDUE_POINTER'):
            section = 'resstart'
        elif line.startswith('%FLAG BONDS_INC_HYDROGEN') or line.startswith('%FLAG BONDS_WITHOUT_HYDROGEN'):
            section = 'bonds'
        elif line.startswith('%FLAG BOX_DIMENSIONS'):
            section = 'box'
        elif line.startswith('%FLAG'):
            section = None

        if line.startswith('%'):
            continue

        if section == 'pointers':
            pass
        elif section == 'names':
            fieldlen = 4
            topo.name += [line[i:i + fieldlen].strip() for i in range(0, len(line), fieldlen)
                      if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'charges':
            fieldlen = 16
            topo.charge += [float(line[i:i + fieldlen].strip()) / 18.2223 for i in range(0, len(line), fieldlen)
                        if len(line[i:i + fieldlen].strip()) != 0]  # 18.2223 = Scaling factor for charges
        elif section == 'masses':
            fieldlen = 16
            topo.masses += [float(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                       if len(line[i:i + fieldlen].strip()) != 0]  # 18.2223 = Scaling factor for charges
        elif section == 'resname':
            fieldlen = 4
            uqresnames += [line[i:i + fieldlen].strip() for i in range(0, len(line), fieldlen)
                           if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'resstart':
            fieldlen = 8
            residx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                       if len(line[i:i + fieldlen].strip()) != 0]
        elif section == 'bonds':
            fieldlen = 8
            bondsidx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                         if len(line[i:i + fieldlen].strip()) != 0]

    # Replicating unique resnames according to their start and end indeces
    residx.append(len(topo.name)+1)

    for i in range(len(residx) - 1):
        numresatoms = residx[i+1] - residx[i]
        topo.resname += [uqresnames[i]] * numresatoms
        topo.resid += [i+1] * numresatoms

    # Processing bond triplets
    for i in range(0, len(bondsidx), 3):
        topo.bonds.append([int(bondsidx[i] / 3), int(bondsidx[i+1] / 3)])

    return topo


def PSFread(filename):
    import re
    residinsertion = re.compile('(\d+)([a-zA-Z])')

    topo = Topology()

    f = open(filename, 'r')
    mode = None

    for line in f:
        if line.strip() == "":
            mode = None

        if mode == 'atom':
            l = line.split()
            topo.serial.append(l[0])
            topo.segid.append(l[1])
            match = residinsertion.findall(l[2])
            if match:
                resid = int(match[0][0])
                insertion = match[0][1]
            else:
                resid = int(l[2])
                insertion = ''
            topo.resid.append(resid)
            topo.insertion.append(insertion)
            topo.resname.append(l[3])
            topo.name.append(l[4])
            topo.atomtype.append(l[5])
            topo.charge.append(float(l[6]))
            topo.masses.append(float(l[7]))
        elif mode == 'bond':
            l = line.split()
            for x in range(0, len(l), 2):
                topo.bonds.append([int(l[x]) - 1, int(l[x + 1]) - 1])
        elif mode == 'angle':
            l = line.split()
            for x in range(0, len(l), 3):
                topo.angles.append([int(l[x]) - 1, int(l[x + 1]) - 1, int(l[x + 2]) - 1])
        elif mode == 'dihedral':
            l = line.split()
            for x in range(0, len(l), 4):
                topo.dihedrals.append([int(l[x]) - 1, int(l[x + 1]) - 1, int(l[x + 2]) - 1, int(l[x + 3]) - 1])
        elif mode == 'improper':
            l = line.split()
            for x in range(0, len(l), 4):
                topo.impropers.append([int(l[x]) - 1, int(l[x + 1]) - 1, int(l[x + 2]) - 1, int(l[x + 3]) - 1])

        if '!NATOM' in line:
            mode = 'atom'
        elif '!NBOND' in line:
            mode = 'bond'
        elif '!NTHETA' in line:
            mode = 'angle'
        elif '!NPHI' in line:
            mode = 'dihedral'
        elif '!NIMPHI' in line:
            mode = 'improper'
    f.close()
    return topo


if __name__ == '__main__':
    from htmd.home import home
    from htmd.molecule.molecule import Molecule
    import os
    testfolder = home(dataDir='molecule-readers/4RWS/')
    mol = Molecule(os.path.join(testfolder, 'structure.psf'))
    print('Can read PSF files.')
    mol.read(os.path.join(testfolder, 'traj.xtc'))
    print('Can read XTC files.')
    testfolder = home(dataDir='molecule-readers/3AM6/')
    mol = Molecule(os.path.join(testfolder, 'structure.prmtop'))
    print('Can read PRMTOP files.')
    mol.read(os.path.join(testfolder, 'structure.crd'))
    print('Can read CRD files.')
    testfolder = home(dataDir='molecule-readers/3L5E/')
    mol = Molecule(os.path.join(testfolder, 'protein.mol2'))
    mol = Molecule(os.path.join(testfolder, 'ligand.mol2'))
    print('Can read MOL2 files.')
