# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import ctypes as ct
import numpy as np
from htmd.molecule.support import pack_double_buffer, pack_int_buffer, pack_string_buffer, pack_ulong_buffer, xtc_lib
from htmd.molecule.util import sequenceID
import os
import logging
logger = logging.getLogger(__name__)

# Pandas NA values taken from https://github.com/pydata/pandas/blob/6645b2b11a82343e5f07b15a25a250f411067819/pandas/io/common.py
# Removed NA because it's natrium!
_NA_VALUES = set([
    '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A',
    'N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan', ''
])


class FormatError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


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
        self.bondtype = []
        self.crystalinfo = None

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

    def fromMolecule(self, mol):
        for field in self.__dict__:
            data = mol.__dict__[field]
            if data is None:
                continue
            if isinstance(data, np.ndarray):
                self.__dict__[field] = data.tolist()
            self.__dict__[field] = data


class Trajectory:
    def __init__(self, coords=None, box=None, boxangles=None, fileloc=None, step=None, time=None):
        self.coords = []
        self.box = []
        self.boxangles = []
        self.fileloc = []
        self.step = []
        self.time = []
        if coords is not None:
            if coords.ndim == 2:
                coords = coords[:, :, np.newaxis]
            self.coords = coords

            nframes = self.numFrames
            if box is None:
                self.box = np.zeros((3, nframes), np.float32)
            if boxangles is None:
                self.boxangles = np.zeros((3, nframes), np.float32)
            if step is None:
                self.step = np.arange(nframes, dtype=int)
            if time is None:
                self.time = np.zeros(nframes, dtype=np.float32)
        if box is not None:
            self.box = box
        if boxangles is not None:
            self.boxangles = boxangles
        if fileloc is not None:
            self.fileloc = fileloc
        if step is not None:
            self.step = step
        if time is not None:
            self.time = time

    @property
    def numFrames(self):
        return self.coords.shape[2]


class TopologyInconsistencyError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class MolFactory:
    """ This class converts Topology and Trajectory data into Molecule objects """
    @staticmethod
    def construct(topos, trajs, filename, frame):
        from htmd.molecule.molecule import Molecule

        topos = ensurelist(topos)
        trajs = ensurelist(trajs)
        if len(topos) != len(trajs):
            raise RuntimeError('Different number of topologies ({}) and trajectories ({}) were read from {}'.format(
                len(topos), len(trajs), filename))

        mols = []
        for topo, traj in zip(topos, trajs):
            natoms = MolFactory._getNumAtoms(topo, traj, filename)

            mol = Molecule()
            if topo is not None:
                mol._emptyTopo(natoms)
                MolFactory._parseTopology(mol, topo, filename)
            if traj is not None:
                mol._emptyTraj(natoms)
                MolFactory._parseTraj(mol, traj, filename, frame)

            mols.append(mol)

        if len(mols) == 1:
            return mols[0]
        else:
            return mols

    @staticmethod
    def _getNumAtoms(topo, traj, filename):
        toponatoms = None
        trajnatoms = None

        if topo is not None:
            natoms = []
            # Checking number of atoms that were read in the topology file for each field are the same
            for field in topo.atominfo:
                if len(topo.__dict__[field]) != 0:
                    natoms.append(len(topo.__dict__[field]))
            natoms = np.unique(natoms)
            if len(natoms) == 0:
                raise RuntimeError('No atoms were read from file {}.'.format(filename))
            if len(natoms) != 1:
                raise TopologyInconsistencyError('Different number of atoms read from file {} for different fields: {}.'
                                                 .format(filename, natoms))
            toponatoms = natoms[0]
        if traj is not None:
            trajnatoms = traj.coords.shape[0]

        if toponatoms is not None and trajnatoms is not None and toponatoms != trajnatoms:
            raise TopologyInconsistencyError('Different number of atoms in topology ({}) and trajectory ({}) for '
                                             'file {}'.format(toponatoms, trajnatoms, filename))

        if toponatoms is not None:
            return toponatoms
        elif trajnatoms is not None:
            return trajnatoms

    @staticmethod
    def _parseTopology(mol, topo, filename):
        from htmd.molecule.molecule import Molecule
        for field in topo.__dict__:
            if field == 'crystalinfo':
                mol.crystalinfo = topo.crystalinfo
                continue

            newfielddata = np.array(topo.__dict__[field], dtype=mol._dtypes[field])

            # Skip on empty new field data
            if newfielddata is None or len(newfielddata) == 0 or np.all([x is None for x in topo.__dict__[field]]):
                continue

            # Objects could be ints for example but we want them as str
            if mol._dtypes[field] == object and len(newfielddata) != 0:
                newfielddata = np.array([str(x) for x in newfielddata], dtype=object)

            mol.__dict__[field] = newfielddata

        if len(mol.bonds) != 0 and len(topo.bondtype) == 0:
            mol.bondtype = np.empty(mol.bonds.shape[0], dtype=Molecule._dtypes['bondtype'])
            mol.bondtype[:] = 'un'

        mol.element = mol._guessMissingElements()

        if os.path.exists(filename):
            filename = os.path.abspath(filename)
        mol.topoloc = filename
        mol.fileloc = [[filename, 0]]
        mol.viewname = os.path.basename(filename)

    @staticmethod
    def _parseTraj(mol, traj, filename, frame):
        from htmd.molecule.molecule import Molecule
        ext = os.path.splitext(filename)[1][1:]

        assert traj.coords.ndim == 3, '{} reader must return 3D coordinates array for file {}'.format(ext, filename)
        assert traj.coords.shape[1] == 3, '{} reader must return 3 values in 2nd dimension for file {}'.format(ext, filename)

        mol.coords = traj.coords.astype(Molecule._dtypes['coords'])
        if traj.box is None:
            mol.box = np.zeros((3, 1), dtype=Molecule._dtypes['box'])
        else:
            mol.box = np.array(traj.box).astype(Molecule._dtypes['box'])
            if mol.box.ndim == 1:
                mol.box = mol.box[:, np.newaxis]

        if traj.boxangles is None:
            mol.boxangles = np.zeros((3, 1), dtype=Molecule._dtypes['boxangles'])
        else:
            mol.boxangles = np.array(traj.boxangles).astype(Molecule._dtypes['boxangles'])
            if mol.boxangles.ndim == 1:
                mol.boxangles = mol.boxangles[:, np.newaxis]

        # mol.fileloc = traj.fileloc
        mol.step = np.hstack(traj.step).astype(int)
        mol.time = np.hstack(traj.time)

        if ext in _TRAJECTORY_READERS and frame is None:
            # Writing hidden index file containing number of frames in trajectory file
            if os.path.isfile(filename):
                MolFactory._writeNumFrames(filename, mol.numFrames)
            ff = range(mol.numFrames)
            # tr.step = tr.step + traj[-1].step[-1] + 1
        elif frame is None:
            ff = [0]
        elif frame is not None:
            ff = [frame]
        else:
            raise AssertionError('Should not reach here')
        mol.fileloc = [[filename, j] for j in ff]

    @staticmethod
    def _writeNumFrames(filepath, numFrames):
        """ Write the number of frames in a hidden file. Allows us to check for trajectory length issues before projecting

        Parameters
        ----------
        filepath : str
            Path to trajectory file
        numFrames : int
            Number of frames in trajectory file
        """
        filepath = os.path.abspath(filepath)
        filedir = os.path.dirname(filepath)
        basename = os.path.basename(filepath)
        numframefile = os.path.join(filedir, '.{}.numframes'.format(basename))
        if not os.path.exists(numframefile) or (os.path.exists(numframefile) and (os.path.getmtime(numframefile) < os.path.getmtime(filepath))):
            try:
                with open(numframefile, 'w') as f:
                    f.write(str(numFrames))
            except:
                pass


def XYZread(filename, frame=None, topoloc=None):
    topo = Topology()

    frames = []
    firstconf = True
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if line == '':
                break
            natoms = int(line.split()[0])
            f.readline()
            coords = []
            for i in range(natoms):
                s = f.readline().split()
                if firstconf:
                    topo.record.append('HETATM')
                    topo.serial.append(i + 1)
                    topo.element.append(s[0])
                    topo.name.append(s[0])
                    topo.resname.append('MOL')
                coords.append(s[1:4])
            frames.append(np.vstack(coords))
            firstconf = False

    coords = np.stack(frames, axis=2)
    traj = Trajectory(coords=coords)
    return MolFactory.construct(topo, traj, filename, frame)


def GJFread(filename, frame=None, topoloc=None):
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

    with open(filename, "r") as f:
        for line in f:
            pieces = line.split(sep=',')
            if len(pieces) == 4 and not line.startswith('$') and not line.startswith('%') and not line.startswith('#'):
                topo.record.append('HETATM')
                topo.element.append(pieces[0])
                topo.name.append(pieces[0])
                topo.resname = 'MOL'
                coords.append([float(s) for s in pieces[1:4]])
        topo.serial = range(len(topo.record))

    coords = np.vstack(coords)[:, :, np.newaxis]
    traj = Trajectory(coords=coords)
    return MolFactory.construct(topo, traj, filename, frame)


def MOL2read(filename, frame=None, topoloc=None, singlemol=True):
    from periodictable import elements
    element_objs = list(elements._element.values())[1:]
    element_symbols = [e.symbol for e in element_objs]
    assert len(element_symbols) == 118

    topologies = []  # Allow reading of multi-mol MOL2 files
    topologies.append(Topology())
    topo = topologies[-1]
    coordinates = [[]]
    coords = coordinates[-1]
    section = None

    molnum = 0
    unguessed = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith('@<TRIPOS>MOLECULE'):
                section = None
                molnum += 1
                if molnum > 1:  # New Molecule, create new topology
                    if singlemol:
                        break
                    topologies.append(Topology())
                    topo = topologies[-1]
                    coordinates.append([])
                    coords = coordinates[-1]
            if line.startswith('@<TRIPOS>ATOM'):
                section = 'atom'
                continue
            if line.startswith('@<TRIPOS>BOND'):
                section = 'bond'
                continue
            if line.startswith('@<TRIPOS>'):  # Skip all other sections
                section = None
                continue

            if section == 'atom':
                pieces = line.strip().split()
                topo.record.append('HETATM')
                topo.serial.append(int(pieces[0]))
                topo.name.append(pieces[1])
                coords.append([float(x) for x in pieces[2:5]])
                topo.atomtype.append(pieces[5])
                if len(pieces) > 6:
                    topo.resid.append(int(pieces[6]))
                if len(pieces) > 7:
                    topo.resname.append(pieces[7][:3])
                if len(pieces) > 8:
                    topo.charge.append(float(pieces[8]))

                element = pieces[5].split('.')[0]
                if element in element_symbols:
                    topo.element.append(element)
                else:
                    unguessed.append(pieces[5])
                    topo.element.append('')
            elif section == 'bond':
                pieces = line.strip().split()
                if len(pieces) < 4:
                    raise RuntimeError('Less than 4 values encountered in bonds definition in line {}'.format(line))
                topo.bonds.append([int(pieces[1]) - 1, int(pieces[2]) - 1])
                topo.bondtype.append(pieces[3])


    if len(unguessed) != 0:
        logger.warning('Could not guess elements for {} atoms with MOL2 atomtypes '
                       '({}).'.format(len(unguessed), ', '.join(np.unique(unguessed))))

    trajectories = []
    for cc in coordinates:
        trajectories.append(Trajectory(coords=np.vstack(cc)[:, :, np.newaxis]))

    if singlemol:
        if molnum > 1:
            logger.warning('Mol2 file {} contained multiple molecules. Only the first was read.'.format(filename))
        return MolFactory.construct(topologies[0], trajectories[0], filename, frame)
    else:
        return MolFactory.construct(topologies, trajectories, filename, frame)


def MAEread(fname, frame=None, topoloc=None):
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

    # Stripping starting and trailing whitespaces which confuse csv reader
    with open(fname, 'r') as csvfile:
        stripped = (row.strip() for row in csvfile)

        import csv
        reader = csv.reader(stripped, delimiter=' ', quotechar='"', skipinitialspace=True)
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
                    if len(row) != len(section_dict):  # TODO: fix the reader
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

    coords = np.vstack(coords)[:, :, np.newaxis]
    traj = Trajectory(coords=coords)
    return MolFactory.construct(topo, traj, fname, frame)


def _getPDB(pdbid):
    import requests
    from htmd.molecule.support import string_to_tempfile
    from htmd.home import home
    # Try loading it from the pdb data directory
    tempfile = False
    localpdb = os.path.join(home(dataDir='pdb'), pdbid.lower() + '.pdb')
    if os.path.isfile(localpdb):
        logger.info('Using local copy for {:s}: {:s}'.format(pdbid, localpdb))
        filepath = localpdb
    else:
        # or the PDB website
        logger.info('Attempting PDB query for {:s}'.format(pdbid))
        r = requests.get('https://files.rcsb.org/download/{}.pdb'.format(pdbid))
        if r.status_code == 200:
            filepath = string_to_tempfile(r.content.decode('ascii'), 'pdb')
            tempfile = True
        else:
            raise NameError('Invalid PDB code')
    return filepath, tempfile


def pdbGuessElementByName(pdtopo):
    """
    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#misalignment which states that elements
    should be right-aligned in columns 13-14 unless it's a 4 letter name when it would end up being left-aligned.
    """
    from periodictable import elements
    from collections import defaultdict
    allelements = [str(el).upper() for el in list(elements._element.values())[1:]]

    alternatives = defaultdict(list)

    if pdtopo.element.dtype == np.float64 and np.all(np.isnan(pdtopo.element)):
        noelem = np.arange(len(pdtopo.element))
    else:
        noelem = np.where((pdtopo.element.str.strip() == '') | pdtopo.element.isnull())[0]

    uqnames = np.unique(pdtopo.name[noelem])
    for name in uqnames:
        elem = None
        if name[0] == ' ':  # If there is no letter in col 13 then col 14 is the element
            elem = name[1]
        else:
            if name[-1] == ' ':  # If it's not a 4 letter name then it's a two letter element
                elem = name[0] + name[1].lower()
            else:  # With 4 letter name it could be either a 1 letter element or 2 letter element
                if name[0] in allelements:
                    elem = name[0]
                if name[:2].upper() in allelements:
                    if elem is not None:
                        altelem = name[0] + name[1].lower()
                        alternatives[(elem, altelem)].append(name)
                    else:
                        elem = name[0] + name[1].lower()
        pdtopo.loc[pdtopo.name == name, 'element'] = elem

    for elem, altelem in alternatives:
        names = np.unique(alternatives[(elem, altelem)])
        namestr = '["' + '", "'.join(names) + '"]'
        altelemname = elements.__dict__[altelem].name
        logger.warning('Atoms with names {} were guessed as element {} but could also be {} ({}). If this is a case,'
                       'you can correct them with mol.element[np.in1d(mol.name, {})] = \'{}\''.format(namestr, elem, altelem,
                                                                                            altelemname, namestr, altelem))


def PDBread(filename, mode='pdb', frame=None, topoloc=None):
    from pandas import read_fwf
    import io

    tempfile = False
    if not os.path.isfile(filename) and len(filename) == 4:  # Could be a PDB id. Try to load it from the PDB website
        filename, tempfile = _getPDB(filename)

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
                        (54, 60), (60, 66), (72, 76), (76, 78), (78, 80)]
        toponames = ('record', 'serial', 'name', 'altloc', 'resname', 'chain', 'resid', 'insertion',
                     'occupancy', 'beta', 'segid', 'element', 'charge')
    elif mode == 'pdbqt':
        # http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
        # The rigid root contains one or more PDBQT-style ATOM or HETATM records. These records resemble their
        # traditional PDB counterparts, but diverge in columns 71-79 inclusive (where the first character in the line
        # corresponds to column 1). The partial charge is stored in columns 71-76 inclusive (in %6.3f format, i.e.
        # right-justified, 6 characters wide, with 3 decimal places). The AutoDock atom-type is stored in columns 78-79
        # inclusive (in %-2.2s format, i.e. left-justified and 2 characters wide..
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
        'charge': np.float32
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
    cryst1colspecs = [(6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54), (55, 66), (66, 70)]
    cryst1names = ('a', 'b', 'c', 'alpha', 'beta', 'gamma', 'sGroup', 'z')

    """
    Guessing columns for REMARK 290 SMTRY from the example since the specs don't define them
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890
    REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
    REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
    REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
    REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       36.30027
    REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000
    REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       59.50256
    REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000
    REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       46.45545
    REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       59.50256
    REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       36.30027
    REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       46.45545
    REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000

    Guessing columns for REMARK 350   BIOMT from the example since the specs don't define them
    REMARK 350   BIOMT1   1 -0.981559  0.191159  0.000000        0.00000
    REMARK 350   BIOMT2   1 -0.191159 -0.981559  0.000000        0.00000
    REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000      -34.13878
    REMARK 350   BIOMT1   2 -0.838088  0.545535  0.000000        0.00000
    REMARK 350   BIOMT2   2 -0.545535 -0.838088  0.000000        0.00000
    REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000      -32.71633
    """
    symmetrycolspecs = [(20, 23), (23, 33), (33, 43), (43, 53), (53, 68)]
    symmetrynames = ('idx', 'rot1', 'rot2', 'rot3', 'trans')

    def concatCoords(coords, coorddata):
        if coorddata.tell() != 0:  # Not empty
            coorddata.seek(0)
            parsedcoor = read_fwf(coorddata, colspecs=coordcolspecs, names=coordnames, na_values=_NA_VALUES, keep_default_na=False)
            if coords is None:
                coords = np.zeros((len(parsedcoor), 3, 0), dtype=np.float32)
            currcoords = np.vstack((parsedcoor.x, parsedcoor.y, parsedcoor.z)).T
            if coords.shape[0] != currcoords.shape[0]:
                logger.warning('Different number of atoms read in different MODELs in the PDB file. '
                               'Keeping only the first {} model(s)'.format(coords.shape[2]))
                return coords
            coords = np.append(coords, currcoords[:, :, np.newaxis], axis=2)
        return coords

    teridx = []
    currter = 0
    topoend = False

    cryst1data = io.StringIO()
    topodata = io.StringIO()
    conectdata = io.StringIO()
    coorddata = io.StringIO()
    symmetrydata = io.StringIO()

    coords = None

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                cryst1data.write(line)
            if line.startswith('ATOM') or line.startswith('HETATM'):
                coorddata.write(line)
            if (line.startswith('ATOM') or line.startswith('HETATM')) and not topoend:
                topodata.write(line)
                teridx.append(str(currter))
            if line.startswith('TER'):
                currter += 1
            if (mode == 'pdb' and line.startswith('END')) or \
               (mode == 'pdbqt' and line.startswith('ENDMDL')):  # pdbqt should not stop reading at ENDROOT or ENDBRANCH
                topoend = True
            if line.startswith('CONECT'):
                conectdata.write(line)
            if line.startswith('MODEL'):
                coords = concatCoords(coords, coorddata)
                coorddata = io.StringIO()
            if line.startswith('REMARK 290   SMTRY'):  # TODO: Support BIOMT fields. It's a bit more complicated. Can't be done with pandas
                symmetrydata.write(line)

        cryst1data.seek(0)
        topodata.seek(0)
        conectdata.seek(0)
        symmetrydata.seek(0)

        coords = concatCoords(coords, coorddata)

        parsedbonds = read_fwf(conectdata, colspecs=bondcolspecs, names=bondnames, na_values=_NA_VALUES, keep_default_na=False)
        parsedcryst1 = read_fwf(cryst1data, colspecs=cryst1colspecs, names=cryst1names, na_values=_NA_VALUES, keep_default_na=False)
        enforceddtypes = {'resname': str, 'segid': str}
        parsedtopo = read_fwf(topodata, colspecs=topocolspecs, names=toponames, na_values=_NA_VALUES, keep_default_na=False, dtype=enforceddtypes, delimiter='\r\t') #, dtype=topodtypes)
        parsedsymmetry = read_fwf(symmetrydata, colspecs=symmetrycolspecs, names=symmetrynames, na_values=_NA_VALUES, keep_default_na=False)

    # Before stripping guess elements from atomname as the spaces are significant
    pdbGuessElementByName(parsedtopo)

    for field in topodtypes:
        if field in parsedtopo and topodtypes[field] == str and parsedtopo[field].dtype == object:
            parsedtopo[field] = parsedtopo[field].str.strip()

    # Fixing PDB format charges which can come after the number
    if parsedtopo.charge.dtype == 'object':
        parsedtopo.charge = parsedtopo.charge.str.strip()
        charges = np.zeros(len(parsedtopo.charge), dtype=np.float32)
        for i, c in enumerate(parsedtopo.charge):
            if not isinstance(c, str):
                continue
            if len(c) > 1:
                if c[1] == '-':
                    charges[i] = -1 * float(c[0])
                elif c[1] == '+':
                    charges[i] = float(c[0])
            elif len(c):
                charges[i] = float(c)
        parsedtopo.charge = charges

    # Fixing hexadecimal index and resids
    # Support for reading hexadecimal
    if parsedtopo.serial.dtype == 'object':
        logger.warning('Non-integer values were read from the PDB "serial" field. Dropping PDB values and assigning new ones.')
        if len(np.unique(parsedtopo.serial)) == len(parsedtopo.serial):
            parsedtopo.serial = sequenceID(parsedtopo.serial)
        else:
            parsedtopo.serial = np.arange(len(parsedtopo.serial))
    if parsedtopo.resid.dtype == 'object':
        logger.warning('Non-integer values were read from the PDB "resid" field. Dropping PDB values and assigning new ones.')
        parsedtopo.resid = sequenceID(parsedtopo.resid)

    if parsedtopo.insertion.dtype == np.float64 and not np.all(np.isnan(parsedtopo.insertion)):
        # Trying to minimize damage from resids overflowing into insertion field
        logger.warning('Integer values detected in "insertion" field. Your resids might be overflowing. Take care')
        parsedtopo.insertion = [str(int(x)) if not np.isnan(x) else '' for x in parsedtopo.insertion]

    if len(parsedtopo) > 99999:
        logger.warning('Reading PDB file with more than 99999 atoms. Bond information can be wrong.')

    crystalinfo = {}
    if len(parsedcryst1):
        crystalinfo = parsedcryst1.iloc[0].to_dict()
        if isinstance(crystalinfo['sGroup'], str) or not np.isnan(crystalinfo['sGroup']):
            crystalinfo['sGroup'] = crystalinfo['sGroup'].split()
    if len(parsedsymmetry):
        numcopies = int(len(parsedsymmetry)/3)
        crystalinfo['numcopies'] = numcopies
        crystalinfo['rotations'] = parsedsymmetry[['rot1', 'rot2', 'rot3']].as_matrix().reshape((numcopies, 3, 3))
        crystalinfo['translations'] = parsedsymmetry['trans'].as_matrix().reshape((numcopies, 3))

    topo = Topology(parsedtopo)

    # Bond formatting part
    # TODO: Speed this up. This is the slowest part for large PDB files. From 700ms to 7s
    serials = parsedtopo.serial.as_matrix()
    # if isinstance(serials[0], str) and np.any(serials == '*****'):
    #     logger.info('Non-integer serials were read. For safety we will discard all bond information and serials will be assigned automatically.')
    #     topo.serial = np.arange(1, len(serials)+1, dtype=np.int)
    if np.max(parsedbonds.max()) > np.max(serials):
        logger.info('Bond indexes in PDB file exceed atom indexes. For safety we will discard all bond information.')
    else:
        mapserials = np.empty(np.max(serials)+1)
        mapserials[:] = np.NAN
        mapserials[serials] = list(range(len(serials)))
        for i in range(len(parsedbonds)):
            row = parsedbonds.loc[i].tolist()
            for b in range(1, 5):
                if not np.isnan(row[b]):
                    topo.bonds.append([int(row[0]), int(row[b])])
        topo.bonds = np.array(topo.bonds, dtype=np.uint32)
        if topo.bonds.size != 0:
            mappedbonds = mapserials[topo.bonds[:]]
            wrongidx, _ = np.where(np.isnan(mappedbonds))  # Some PDBs have bonds to non-existing serials... go figure
            if len(wrongidx):
                logger.info('Discarding {} bonds to non-existing indexes in the PDB file.'.format(len(wrongidx)))
            mappedbonds = np.delete(mappedbonds, wrongidx, axis=0)
            topo.bonds = np.array(mappedbonds, dtype=np.uint32)

    if len(topo.segid) and np.all(np.array(topo.segid) == '') and currter != 0:  # If no segid was read, use the TER rows to define segments
        topo.segid = teridx

    if tempfile:
        os.unlink(filename)

    topo.crystalinfo = crystalinfo
    traj = Trajectory(coords=coords)
    return MolFactory.construct(topo, traj, filename, frame)


def PDBQTread(filename, frame=None, topoloc=None):
    return PDBread(filename, mode='pdbqt', frame=frame, topoloc=topoloc)


def PRMTOPread(filename, frame=None, topoloc=None):
    with open(filename, 'r') as f:
        topo = Topology()
        uqresnames = []
        residx = []
        bondsidx = []
        angleidx = []
        dihedidx = []
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
            elif line.startswith('%FLAG ANGLES_INC_HYDROGEN') or line.startswith('%FLAG ANGLES_WITHOUT_HYDROGEN'):
                section = 'angles'
            elif line.startswith('%FLAG DIHEDRALS_INC_HYDROGEN') or line.startswith('%FLAG DIHEDRALS_WITHOUT_HYDROGEN'):
                section = 'dihedrals'
            elif line.startswith('%FLAG BOX_DIMENSIONS'):
                section = 'box'
            elif line.startswith('%FLAG AMBER_ATOM_TYPE'):
                section = 'amberatomtype'
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
            elif section == 'angles':
                fieldlen = 8
                angleidx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                             if len(line[i:i + fieldlen].strip()) != 0]
            elif section == 'dihedrals':
                fieldlen = 8
                dihedidx += [int(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                             if len(line[i:i + fieldlen].strip()) != 0]
            elif section == 'amberatomtype':
                fieldlen = 4
                topo.atomtype += [line[i:i + fieldlen].strip() for i in range(0, len(line), fieldlen)
                                  if len(line[i:i + fieldlen].strip()) != 0]


    if len(topo.name) == 0:
        raise FormatError('No atoms read in PRMTOP file. Trying a different reader.')
    # Replicating unique resnames according to their start and end indeces
    residx.append(len(topo.name)+1)

    """
    NOTE: the atom numbers in the following arrays that describe bonds, angles, and dihedrals are coordinate array 
    indexes for runtime speed. The true atom number equals the absolute value of the number divided by three, plus one. 
    In the case of the dihedrals, if the fourth atom is negative, this implies that the dihedral is an improper. If the 
    third atom is negative, this implies that the end group interations are to be ignored. End group interactions are 
    ignored, for example, in dihedrals of various ring systems (to prevent double counting of 1-4 interactions) and 
    in multiterm dihedrals.
    """

    for i in range(len(residx) - 1):
        numresatoms = residx[i+1] - residx[i]
        topo.resname += [uqresnames[i]] * numresatoms
        topo.resid += [i+1] * numresatoms

    # Processing bond triplets
    for i in range(0, len(bondsidx), 3):
        topo.bonds.append([int(bondsidx[i] / 3), int(bondsidx[i+1] / 3)])

    # Processing angle quads
    for i in range(0, len(angleidx), 4):
        topo.angles.append([int(angleidx[i] / 3), int(angleidx[i + 1] / 3), int(angleidx[i + 2] / 3)])

    # Processing dihedral quints
    for i in range(0, len(dihedidx), 5):
        atoms = [int(dihedidx[i] / 3), int(dihedidx[i + 1] / 3), abs(int(dihedidx[i + 2] / 3)), int(dihedidx[i + 3] / 3)]
        if atoms[3] >= 0:
            topo.dihedrals.append(atoms)
        else:
            atoms[3] = abs(atoms[3])
            topo.impropers.append(atoms)
    return MolFactory.construct(topo, None, filename, frame)


def PSFread(filename, frame=None, topoloc=None):
    import re
    residinsertion = re.compile('(\d+)([a-zA-Z])')

    topo = Topology()

    with open(filename, 'r') as f:
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
    return MolFactory.construct(topo, None, filename, frame)


def XTCread(filename, frame=None, topoloc=None):
    """ Reads XTC file

    Parameters
    ----------
    filename : str
        Path of xtc file.
    frame : list
        A list of integer frames which we want to read from the file. If None will read all.

    Returns
    -------
    coords : nd.array
    box : nd.array
    boxangles : nd.array
    step : nd.array
    time : nd.array
    """
    givenframes = frame
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

    coords = None
    if givenframes is None:  # Read the whole XTC file at once
        retval = lib['libxtc'].xtc_read(
            ct.c_char_p(filename.encode("ascii")),
            natoms,
            nframes, deltat, deltastep)
        if not retval:
            raise IOError('XTC file {} possibly corrupt.'.format(filename))
        nframes = nframes[0]
        frames = range(nframes)
        coords = np.zeros((natoms[0], 3, nframes), dtype=np.float32)
    else:
        if not isinstance(givenframes, list) and not isinstance(givenframes, np.ndarray):
            givenframes = [givenframes]
        nframes = len(givenframes)
        frames = givenframes

    step = np.zeros(nframes, dtype=np.uint64)
    time = np.zeros(nframes, dtype=np.float32)
    box = np.zeros((3, nframes), dtype=np.float32)
    boxangles = np.zeros((3, nframes), dtype=np.float32)

    for i, f in enumerate(frames):
        if givenframes is not None:  # If frames were given, read specific frame
            retval = lib['libxtc'].xtc_read_frame(
                ct.c_char_p(filename.encode("ascii")),
                natoms,
                ct.c_int(f))
            if not retval:
                raise IOError('XTC file {} possibly corrupt.'.format(filename))
            if coords is None:
                coords = np.zeros((natoms[0], 3, nframes), dtype=np.float32)
            fidx = 0
        else:
            fidx = f

        step[i] = retval[fidx].step
        time[i] = retval[fidx].time
        box[:, i] = retval[fidx].box
        coords[:, :, i] = np.ctypeslib.as_array(retval[fidx].pos, shape=(natoms[0], 3))

        if givenframes is not None:
            lib['libc'].free(retval[0].pos)
            lib['libc'].free(retval)

    if givenframes is None:
        for f in range(len(frames)):
            lib['libc'].free(retval[f].pos)
        lib['libc'].free(retval)

    if np.size(coords, 2) == 0:
        raise NameError('Malformed XTC file. No frames read from: {}'.format(filename))
    if np.size(coords, 0) == 0:
        raise NameError('Malformed XTC file. No atoms read from: {}'.format(filename))

    coords *= 10.  # Convert from nm to Angstrom
    box *= 10.  # Convert from nm to Angstrom
    nframes = coords.shape[2]
    if len(step) != nframes or np.sum(step) == 0:
        step = np.arange(nframes)
    if len(time) != nframes or np.sum(time) == 0:
        time = np.zeros(nframes, dtype=np.float32)
    return MolFactory.construct(None, Trajectory(coords=coords, box=box, boxangles=boxangles, step=step, time=time), filename, frame)


def CRDread(filename, frame=None, topoloc=None):
    #default_name
    #  7196
    #  -7.0046035  10.4479194  20.8320000  -7.3970000   9.4310000  20.8320000
    #  -7.0486898   8.9066002  21.7218220  -7.0486899   8.9065995  19.9421780

    with open(filename, 'r') as f:
        lines = f.readlines()

        if lines[0].startswith('*'):
            raise FormatError('CRDread failed. Trying other readers.')

        coords = []
        fieldlen = 12
        k = 0
        for line in lines[2:]:  # skip first 2 lines
            coords += [float(line[i:i + fieldlen].strip()) for i in range(0, len(line), fieldlen)
                       if len(line[i:i + fieldlen].strip()) != 0]

    coords = np.vstack([coords[i:i + 3] for i in range(0, len(coords), 3)])[:, :, np.newaxis]
    return MolFactory.construct(None, Trajectory(coords=coords), filename, frame)


def CRDCARDread(filename, frame=None, topoloc=None):
    """ https://www.charmmtutorial.org/index.php/CHARMM:The_Basics
        title = * WATER
        title = *  DATE:     4/10/07      4:25:51      CREATED BY USER: USER
        title = *
        Number of atoms (NATOM)       = 6
        Atom number (ATOMNO)          = 1 (just an exmaple)
        Residue number (RESNO)        = 1
        Residue name (RESName)        = TIP3
        Atom type (TYPE)              = OH2
        Coordinate (X)                = -1.30910
        Coordinate (Y)                = -0.25601
        Coordinate (Z)                = -0.24045
        Segment ID (SEGID)            = W
        Residue ID (RESID)            = 1
        Atom weight (Weighting)       = 0.00000

        now what that looks like...

        * WATER
        *  DATE:     4/10/07      4:25:51      CREATED BY USER: USER
        *
            6
            1    1 TIP3 OH2   -1.30910  -0.25601  -0.24045 W    1      0.00000
            2    1 TIP3 H1    -1.85344   0.07163   0.52275 W    1      0.00000
            3    1 TIP3 H2    -1.70410   0.16529  -1.04499 W    1      0.00000
            4    2 TIP3 OH2    1.37293   0.05498   0.10603 W    2      0.00000
            5    2 TIP3 H1     1.65858  -0.85643   0.10318 W    2      0.00000
            6    2 TIP3 H2     0.40780  -0.02508  -0.02820 W    2      0.00000
    """
    coords = []
    topo = Topology()
    with open(filename, 'r') as f:
        lines = f.readlines()

        if not lines[0].startswith('*'):
            raise FormatError('CRDCARDread failed. Trying other readers.')

        i = 0
        while lines[i].startswith('*'):
            i += 1

        for line in lines[i+1:]:
            pieces = line.split()
            topo.resname.append(pieces[2])
            topo.name.append(pieces[3])
            coords.append([float(x) for x in pieces[4:7]])
            topo.segid.append(pieces[7])
            topo.resid.append(int(pieces[8]))
    coords = np.vstack(coords)[:, :, np.newaxis]
    return MolFactory.construct(topo, Trajectory(coords=coords), filename, frame)


def BINCOORread(filename, frame=None, topoloc=None):
    import struct
    with open(filename, 'rb') as f:
        dat = f.read(4)
        fmt = 'i'
        natoms = struct.unpack(fmt, dat)[0]
        dat = f.read(natoms * 3 * 8)
        fmt = 'd' * (natoms * 3)
        coords = struct.unpack(fmt, dat)
        coords = np.array(coords, dtype=np.float32).reshape((natoms, 3, 1))
    return MolFactory.construct(None, Trajectory(coords=coords), filename, frame)


def MDTRAJread(filename, frame=None, topoloc=None):
    import mdtraj as md
    traj = md.load(filename, top=topoloc)
    coords = np.swapaxes(np.swapaxes(traj.xyz, 0, 1), 1, 2) * 10
    step = traj.time / traj.timestep
    time = traj.time * 1000  # need to go from picoseconds to femtoseconds

    if traj.unitcell_lengths is None:
        box = None
    else:
        box = traj.unitcell_lengths.T.copy() * 10

    if traj.unitcell_angles is None:
        boxangles = None
    else:
        boxangles = traj.unitcell_angles.T.copy()
    traj = Trajectory(coords=coords.copy(), box=box, boxangles=boxangles, step=step, time=time)  # Copying coords needed to fix MDtraj stride
    return MolFactory.construct(None, traj, filename, frame)


def MDTRAJTOPOread(filename, frame=None, topoloc=None):
    translate = {'serial': 'serial', 'name': 'name', 'element': 'element', 'resSeq': 'resid', 'resName': 'resname',
                 'chainID': 'chain', 'segmentID': 'segid'}
    import mdtraj as md
    mdstruct = md.load(filename)
    topology = mdstruct.topology
    table, bonds = topology.to_dataframe()

    topo = Topology()
    for k in table.keys():
        topo.__dict__[translate[k]] = table[k].tolist()

    coords = np.array(mdstruct.xyz.swapaxes(0, 1).swapaxes(1, 2) * 10, dtype=np.float32)
    topo.bonds = bonds
    return MolFactory.construct(topo, Trajectory(coords=coords), filename, frame)


def GROTOPread(filename, frame=None, topoloc=None):
    # Reader for GROMACS .top file format:
    # http://manual.gromacs.org/online/top.html
    topo = Topology()
    section = None
    atmidx = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(';') or line.startswith('#') or len(line.strip()) == 0:
                continue
            if not line.startswith('[') and section == 'atoms':
                pieces = line.split()
                atmidx.append(int(pieces[0]))
                topo.resid.append(pieces[2])
                topo.resname.append(pieces[3])
                topo.name.append(pieces[4])
                topo.charge.append(pieces[6])
                topo.element.append(pieces[1])
            if not line.startswith('[') and section == 'bonds':
                pieces = line.split()
                topo.bonds.append([int(pieces[0]), int(pieces[1])])

            if '[ atoms ]' in line:
                section = 'atoms'
            elif '[ bonds ]' in line:
                section = 'bonds'
            elif line.startswith('['):
                section = None

    if section is None and len(topo.name) == 0:
        raise FormatError('No atoms read in GROTOP file. Trying a different reader.')

    atmidx = np.array(atmidx)
    atommapping = np.ones(np.max(atmidx) + 1) * -1
    atommapping[atmidx] = np.arange(len(atmidx))
    for i in range(len(topo.bonds)):
        topo.bonds[i][0] = atommapping[topo.bonds[i][0]]

    return MolFactory.construct(topo, None, filename, frame)

def PDBXMMCIFread(filename, frame=None, topoloc=None):
    from htmd.molecule.pdbx.reader.PdbxReader import PdbxReader
    myDataList = []
    ifh = open(filename, "r")
    pRd = PdbxReader(ifh)
    pRd.read(myDataList)
    ifh.close()

    # Taken from http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html#ATOMP
    atom_site_mapping = {'group_PDB': ('record', str),
                         'id': ('serial', int),
                         'auth_atom_id': ('name', str),
                         'label_alt_id': ('altloc', str),
                         'auth_comp_id': ('resname', str),
                         'auth_asym_id': ('chain', str),
                         'auth_seq_id': ('resid', int),
                         'pdbx_PDB_ins_code': ('insertion', str),
                         'label_entity_id': ('segid', str),
                         'type_symbol': ('element', str),
                         'occupancy': ('occupancy', float),
                         'B_iso_or_equiv': ('beta', float),
                         'pdbx_formal_charge': ('charge', float)}

    cryst1_mapping = {'length_a': ('a', float),
                      'length_b': ('b', float),
                      'length_c': ('c', float),
                      'angle_alpha': ('alpha', float),
                      'angle_beta': ('beta', float),
                      'angle_gamma': ('gamma', float),
                      'space_group_name_H-M': ('sGroup', str),
                      'Z_PDB': ('z', int)}

    topo = Topology()

    if len(myDataList) > 1:
        logger.warning('Multiple Data objects in mmCIF. Please report this issue to the HTMD issue tracker')

    dataObj = myDataList[0]

    def fixDefault(val, dtype):
        if val == '?':
            if dtype == float or dtype == int:
                val = 0
            if dtype == str:
                val = ''
        return val

    # Parsing CRYST1 data
    cryst = dataObj.getObj('cell')
    if cryst is not None and cryst.getRowCount() == 1:
        row = cryst.getRow(0)

        crystalinfo = {}
        for source_field, target in cryst1_mapping.items():
            target_field, dtype = target
            val = dtype(fixDefault(row[cryst.getAttributeIndex(source_field)], dtype))
            crystalinfo[target_field] = val

        if isinstance(crystalinfo['sGroup'], str) or not np.isnan(crystalinfo['sGroup']):
            crystalinfo['sGroup'] = crystalinfo['sGroup'].split()
        topo.crystalinfo = crystalinfo

    # Parsing ATOM and HETATM data
    allcoords = []
    coords = []
    currmodel = -1
    firstmodel = None
    atom_site = dataObj.getObj('atom_site')
    for i in range(atom_site.getRowCount()):
        row = atom_site.getRow(i)
        modelid = row[atom_site.getAttributeIndex('pdbx_PDB_model_num')]
        # On a new model, restart coords and append the old ones
        if currmodel != -1 and currmodel != modelid:
            currmodel = modelid
            allcoords.append(np.array(coords, dtype=np.float32))
            coords = []

        if currmodel == -1:
            currmodel = modelid
            firstmodel = modelid

        if currmodel == firstmodel:
            for source_field, target in atom_site_mapping.items():
                target_field, dtype = target
                val = row[atom_site.getAttributeIndex(source_field)]
                val = dtype(fixDefault(val, dtype))
                if source_field == 'label_alt_id' and val == '.':  # Atoms without altloc seem to be stored with a dot
                    val = ''
                topo.__dict__[target_field].append(val)

        coords.append([row[atom_site.getAttributeIndex('Cartn_x')],
                       row[atom_site.getAttributeIndex('Cartn_y')],
                       row[atom_site.getAttributeIndex('Cartn_z')]])

    if len(coords) != 0:
        allcoords.append(np.array(coords, dtype=np.float32))

    allcoords = np.stack(allcoords, axis=2)

    return MolFactory.construct(topo, Trajectory(coords=allcoords), filename, frame)




# Register here all readers with their extensions
_TOPOLOGY_READERS = {'prmtop': PRMTOPread,
                     'prm': PRMTOPread,
                     'psf': PSFread,
                     'mae': MAEread,
                     'mol2': MOL2read,
                     'gjf': GJFread,
                     'xyz': XYZread,
                     'pdb': PDBread,
                     'ent': PDBread,
                     'pdbqt': PDBQTread,
                     'top': [GROTOPread, PRMTOPread],
                     'crd': CRDCARDread,
                     'cif': PDBXMMCIFread}

from mdtraj.core.trajectory import _TOPOLOGY_EXTS as _MDTRAJ_TOPOLOGY_EXTS
_MDTRAJ_TOPOLOGY_EXTS = [x[1:] for x in _MDTRAJ_TOPOLOGY_EXTS]  # Removing the initial dot
for ext in _MDTRAJ_TOPOLOGY_EXTS:
    if ext not in _TOPOLOGY_READERS:
        _TOPOLOGY_READERS[ext] = MDTRAJTOPOread

_TRAJECTORY_READERS = {'xtc': XTCread}

_COORDINATE_READERS = {'crd': CRDread,
                       'coor': BINCOORread}

_MDTRAJ_TRAJECTORY_EXTS = ('dcd', 'binpos', 'trr', 'nc', 'h5', 'lh5', 'netcdf')
for ext in _MDTRAJ_TRAJECTORY_EXTS:
    if ext not in _TRAJECTORY_READERS:
        _TRAJECTORY_READERS[ext] = MDTRAJread

from htmd.util import ensurelist
_ALL_READERS = {}
for k in _TOPOLOGY_READERS:
    if k not in _ALL_READERS:
        _ALL_READERS[k] = []
    _ALL_READERS[k] += ensurelist(_TOPOLOGY_READERS[k])

for k in _TRAJECTORY_READERS:
    if k not in _ALL_READERS:
        _ALL_READERS[k] = []
    _ALL_READERS[k] += ensurelist(_TRAJECTORY_READERS[k])

for k in _COORDINATE_READERS:
    if k not in _ALL_READERS:
        _ALL_READERS[k] = []
    _ALL_READERS[k] += ensurelist(_COORDINATE_READERS[k])


from htmd.molecule.molecule import Molecule
from glob import glob
from natsort import natsorted
from unittest import TestCase
class TestReaders(TestCase):
    def testfolder(self, subname=None):
        from htmd.home import home
        if subname is None:
            return home(dataDir='molecule-readers')
        else:
            return os.path.join(home(dataDir='molecule-readers'), subname)

    def test_pdb(self):
        from htmd.home import home
        for f in glob(os.path.join(self.testfolder(), '*.pdb')):
            mol = Molecule(f)
        for f in glob(os.path.join(home(dataDir='pdb/'), '*.pdb')):
            mol = Molecule(f)

    def test_psf(self):
        mol = Molecule(os.path.join(self.testfolder('4RWS'), 'structure.psf'))

    def test_xtc(self):
        mol = Molecule(os.path.join(self.testfolder('4RWS'), 'traj.xtc'))

    def test_combine_topo_traj(self):
        testfolder = self.testfolder('4RWS')
        mol = Molecule(os.path.join(testfolder, 'structure.psf'))
        mol.read(os.path.join(testfolder, 'traj.xtc'))

    def test_prmtop(self):
        testfolder = self.testfolder('3AM6')
        mol = Molecule(os.path.join(testfolder, 'structure.prmtop'))

    def test_crd(self):
        testfolder = self.testfolder('3AM6')
        mol = Molecule(os.path.join(testfolder, 'structure.crd'))

    def test_mol2(self):
        testfolder = self.testfolder('3L5E')
        mol = Molecule(os.path.join(testfolder, 'protein.mol2'))
        mol = Molecule(os.path.join(testfolder, 'ligand.mol2'))

    def test_mae(self):
        for f in glob(os.path.join(self.testfolder(), '*.mae')):
            mol = Molecule(f)

    def test_append_trajectories(self):
        testfolder = self.testfolder('CMYBKIX')
        mol = Molecule(os.path.join(testfolder, 'filtered.pdb'))
        mol.read(natsorted(glob(os.path.join(testfolder, '*.xtc'))))

    def test_missing_crystal_info(self):
        mol = Molecule(os.path.join(self.testfolder(), 'weird-cryst.pdb'))

    def test_dcd(self):
        from htmd.home import home
        mol = Molecule(os.path.join(home(dataDir='1kdx'), '1kdx_0.pdb'))
        mol.read(os.path.join(home(dataDir='1kdx'), '1kdx.dcd'))

    def test_dcd_frames(self):
        from htmd.home import home
        mol = Molecule(os.path.join(home(dataDir='1kdx'), '1kdx_0.pdb'))
        mol.read(os.path.join(home(dataDir='1kdx'), '1kdx.dcd'))
        tmpcoo = mol.coords.copy()
        mol.read([os.path.join(home(dataDir='1kdx'), '1kdx.dcd')], frames=[8])
        assert np.array_equal(tmpcoo[:, :, 8], np.squeeze(mol.coords)), 'Specific frame reading not working'

    def test_gromacs_top(self):
        mol = Molecule(os.path.join(self.testfolder(), 'gromacs.top'))

    def test_mmcif_single_frame(self):
        mol = Molecule(os.path.join(self.testfolder(), '1ffk.cif'))
        assert mol.numAtoms == 64281
        assert mol.numFrames == 1

    def test_mmcif_multi_frame(self):
        mol = Molecule(os.path.join(self.testfolder(), '1j8k.cif'))
        assert mol.numAtoms == 1402
        assert mol.numFrames == 20

    def test_multiple_file_fileloc(self):
        from htmd.home import home
        mol = Molecule(os.path.join(home(dataDir='adaptive'), 'input', 'e1s1_1', 'structure.pdb'))
        mol.read(glob(os.path.join(home(dataDir='adaptive'), 'data', '*', '*.xtc')))
        # Try to vstack fileloc. This will fail with wrong fileloc shape
        fileloc = np.vstack(mol.fileloc)
        assert fileloc.shape == (12, 2)
        print('Correct fileloc shape with multiple file reading.')

    def test_topo_overwriting(self):
        from htmd.home import home
        # Testing overwriting of topology fields
        mol = Molecule(os.path.join(home(dataDir='test-ffevaluate'), '1dihedral', 'mol.psf'))
        atomtypes = mol.atomtype.copy()
        charges = mol.charge.copy()
        coords = np.array([[[0.],
                            [0.],
                            [-0.17]],

                           [[0.007],
                            [1.21],
                            [0.523]],

                           [[0.],
                            [0.],
                            [-1.643]],

                           [[-0.741],
                            [-0.864],
                            [-2.296]]], dtype=np.float32)

        mol.read(os.path.join(home(dataDir='test-ffevaluate'), '1dihedral', 'mol.pdb'))
        assert np.array_equal(mol.atomtype, atomtypes)
        assert np.array_equal(mol.charge, charges)
        assert np.array_equal(mol.coords, coords)
        print('Merging of topology fields works')

    def test_integer_resnames(self):
        mol = Molecule(os.path.join(self.testfolder(), 'errors.pdb'))
        assert np.unique(mol.resname) == '007'

    def test_star_indexes(self):
        mol = Molecule(os.path.join(self.testfolder(), 'errors.pdb'))
        assert np.all(mol.serial == np.arange(mol.numAtoms))

    def test_pdb_element_guessing(self):
        mol = Molecule(os.path.join(self.testfolder(), 'errors.pdb'))
        refelem = np.array(['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N', 'H', 'H', 'C', 'Cl', 'Ca'], dtype=object)
        assert np.array_equal(mol.element, refelem)

    def test_pdb_charges(self):
        mol = Molecule(os.path.join(self.testfolder(), 'errors.pdb'))
        refcharge = np.array([-1.,  1., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.], dtype=np.float32)
        assert np.array_equal(mol.charge, refcharge)


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2)





