# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import requests
import numpy as np
from htmd.molecule.bincoor import *
from htmd.molecule.pdbparser import *
from htmd.molecule.prmtop import *
from htmd.molecule.psf import *
from htmd.molecule.vmdparser import *
from htmd.molecule.xtc import *
from htmd.molecule.wrap import *
from htmd.rotationmatrix import rotationmatrix
from htmd.vmdviewer import getCurrentViewer
from math import pi
from copy import deepcopy
from os import path
import logging
import re
logger = logging.getLogger(__name__)


class Molecule:
    """ Class to manipulate molecular structures.

    Molecule contains all the fields of a PDB and it is independent of any force field. It can contain multiple
    conformations and trajectories, however all operations are done on the current frame.

    Parameters
    ----------
    filename : str
            Optionally load a PDB file from the specified file. If there's no file and the value is four characters long
            assume it is a PDB accession code and try to download from the RCSB web server.
    name : str
        Give a name to the Molecule that will be used for visualization

    Examples
    --------
    >>> mol = Molecule( './test/data/dhfr/dhfr.pdb' )  # doctest: +SKIP
    >>> mol = Molecule( '3PTB', name='Trypsin' )
    >>> print(mol)                                     # doctest: +ELLIPSIS
    Molecule with 1701 atoms and 1 frames
    PDB field - altloc shape: (1701,)
    PDB field - beta shape: (1701,)
    ...

    Properties
    ----------
    PDB Fields

    record : np.ndarray
        PDB record field.
    serial : np.ndarray
        PDB serial field.
    name : np.ndarray
        PDB name field.
    altloc : np.ndarray
        PDB alternative location field.
    resname : np.ndarray
        PDB residue name field.
    chain : np.ndarray
        PDB chain field.
    resid : np.ndarray
        PDB residue ID field.
    insertion : np.ndarray
        PDB insertion code field.
    occupancy : np.ndarray
        PDB occupancy field.
    beta : np.ndarray
        PDB beta value field.
    segid : np.ndarray
        PDB segment ID field.
    element : np.ndarray
        PDB element field.
    charge : np.ndarray
        PDB charge field.
    bonds : np.ndarray
        PDB bonds information.
    ssbonds : np.ndarray
        PDB secondary structure bonds field.
    coords : np.ndarray
        PDB coordinates. 3D matrix. [number of atoms x 3 x number of frames] where the second dimension are the [x,y,z]

    Other Fields

    box : np.ndarray
        Box dimensions of the simulation.
    charge : np.ndarray
        Charges read from prmtop or psf files.
    masses : np.ndarray
        Masses read from prmtop or psf files.
    frame : int
        The current frame. Atomselections and get commands will be calculated on this frame.
    fileloc : list
        The location of the files used to read this Molecule
    time : list
        The time for each frame of the simulation
    step : list
        The step for each frame of the simulation
    reps : :class:`Representations` object
        A list of representations that is used when visualizing the molecule
    viewname : str
        The name used for the molecule in the viewer
    """
    _pdb_fields = {
        'record': object,
        'serial': numpy.int,
        'name': object,
        'altloc': object,
        'resname': object,
        'chain': object,
        'resid': numpy.int,
        'insertion': object,
        'coords': numpy.float32,
        'occupancy': numpy.float32,
        'beta': numpy.float32,
        'segid': object,
        'element': object,
        'charge': numpy.float32
    }

    def __init__(self, filename=None, name=None):
        self.bonds = []
        self.ssbonds = []
        self.box = None
        self.charge = []
        self.masses = None
        self.frame = 0
        self.fileloc = []
        self._append_fields = self._pdb_fields.copy()
        self._append_fields['masses'] = numpy.float32
        self.time = []
        self.step = []
        self.reps = Representations(self)
        self._tempreps = Representations(self)
        self.viewname = name

        for k in self._pdb_fields:
            if k == 'coords':
                self.__dict__[k] = np.zeros((0, 0, 0), dtype=self._pdb_fields[k])
            else:
                self.__dict__[k] = np.zeros(0, dtype=self._pdb_fields[k])
        if filename:
            self.read(filename)
            if isinstance(filename, str):
                self.topoloc = os.path.abspath(filename)
                if name is None and isinstance(filename, str):
                    self.viewname = filename
                    if path.isfile(filename):
                        self.viewname = path.basename(filename)
        #if( self.bonds is None ) or (not len(self.bonds)):
        #    self.guessBonds()

    def insert(self, mol, index):
        """Insert the contents of one molecule into another at a specific index.

        Parameters
        ----------
        mol   : :class:`Molecule`
                Molecule to be inserted
        index : integer
                The atom index at which the passed molecule will be inserted

        Example
        -------
        >>> mol.insert(Molecule('3PTB'), 158)
        """
        backup = self.copy()
        mol = mol.copy()  # Copy because I'll modify its bonds
        if len(mol.bonds) > 0:
            mol.bonds += index

        for k in mol._pdb_fields:
            if mol.__dict__[k] is not None and np.size(mol.__dict__[k]) != 0:
                numatoms = np.size(mol.__dict__[k])
                break

        try:
            if len(self.bonds) > 0:
                self.bonds[self.bonds >= index] += mol.numAtoms
                if len(mol.bonds) > 0:
                    self.bonds = np.append(self.bonds, mol.bonds, axis=0)
            elif len(mol.bonds) > 0:
                self.bonds = mol.bonds

            for k in self._append_fields:
                if k == 'serial':
                    continue
                if mol.__dict__[k] is None or np.size(mol.__dict__[k]) == 0:
                    self.__dict__[k] = np.insert(self.__dict__[k], index, np.zeros(numatoms, dtype=self.__dict__[k].dtype), axis=0)
                elif k == 'coords':
                    self.coords = np.insert(self.coords, index, np.atleast_3d(mol.coords), axis=0)
                else:
                    self.__dict__[k] = np.insert(self.__dict__[k], index, mol.__dict__[k], axis=0)
            self.serial = np.arange(1, self.numAtoms+1)
        except:
            self = backup
            raise NameError('Failed to insert molecule.')
        # TODO: Don't allow user to insert atoms inside a residue, only between (?)

    def remove(self, selection, _logger=True):
        """ Remove atoms from the Molecule

        Parameters
        ----------
        selection : str
            Atomselection string selecting the atoms we want to remove

        Returns
        -------
        removed : np.array
            The list of atoms removed

        Example
        -------
        >>> mol.remove('name CA')               # doctest: +ELLIPSIS
        array([   1,    9,   16,   20,   24,   36,   43,   49,   53,   58,...

        """
        sel = np.where(self.atomselect(selection))[0]
        self._removeBonds(sel)
        for k in self._append_fields:
            self.__dict__[k] = np.delete(self.__dict__[k], sel, axis=0)
        if _logger:
            logger.info('Removed {} atoms. {} atoms remaining in the molecule.'.format(len(sel), self.numAtoms))
        return sel

    def get(self, field, sel=None):
        """Retrieve a specific PDB field based on the selection

        Parameters
        ----------
        field : str
            The PDB field we want to get
        sel : str
            Atom selection string selecting which atoms we want to get the field from. Default all.

        Returns
        ------
        vals : np.ndarray
            Array of values of `field` for all atoms in the selection.

        Examples
        --------
        >>> mol.get('resname')
        array(['ILE', 'ILE', 'ILE', ..., 'ASN', 'ASN', 'ASN'], dtype=object)
        >>> mol.get('resname', sel='resid 158')
        array(['LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU',
               'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU'], dtype=object)

        """
        if field != 'index' and field not in self._pdb_fields:
            raise NameError("Invalid field '" + field + "'")
        s = self.atomselect(sel)
        if field == 'coords':
            return np.squeeze(self.coords[s, :, self.frame])
        elif field == 'index':
            return np.where(s)[0]
        else:
            return self.__dict__[field][s]

    def set(self, field, value, sel=None):
        """Set a specific PDB field based on the selection

        Parameters
        ----------
        field     : str
                    Which field of the Molecule to set
        value     : string or integer
                    All atoms that match the atom selection will have the PDB field `field` set to this scalar value
                     (Or 3-vector if setting the coordinates)
        sel       : str
                    Atom selection string

        Examples
        --------
        >>> mol.set('segid', 'P', sel='protein')
        """
        if field not in self._pdb_fields:
            raise NameError("Invalid field '" + field + "'")
        s = self.atomselect(sel)
        if field == 'coords':
            self.__dict__[field][s, :, self.frame] = value
        else:
            self.__dict__[field][s] = value

    def align(self, sel, refmol=None, refsel=None, frames=None):
        """ Align the molecule to a reference structure

        Parameters
        ----------
        sel : str
            Atom selection string
        refmol : :class:`Molecule`, optional
            Optionally pass a reference Molecule on which to align. If none is given it will align on the first frame of the same Molecule
        refsel : str, optional
            An atom selection for the reference Molecule if one is given. Default: same as `sel`
        frames : list or range
            A list of frames which to align. By default it will align all frames of the Molecule

        Examples
        --------
        >>> mol.align('protein')
        >>> mol.align('name CA', refmol=Molecule('3PTB'))
        """
        if refmol is None:
            refmol = self
        if refsel is None:
            refsel = sel
        if frames is None:
            frames = range(self.numFrames)
        # if not isinstance(refmol, Molecule):
        # raise NameError('Reference molecule has to be a Molecule object')
        sel = self.atomselect(sel)
        refsel = refmol.atomselect(refsel)
        if np.sum(sel) != np.sum(refsel):
            raise NameError('Cannot align molecules. The two selections produced different number of atoms')
        for f in frames:
            P = self.coords[sel, :, f]
            Q = refmol.coords[refsel, :, refmol.frame]
            all1 = self.coords[:, :, f]
            (rot, tmp) = _pp_measure_fit(P, Q)
            # Translating mol to 0,0,0
            centroidP = np.mean(P, 0)
            centroidQ = np.mean(Q, 0)
            all1 = all1 - centroidP
            # Rotating mol
            all1 = np.dot(all1, np.transpose(rot))
            # Translating to centroid of refmol
            all1 = all1 + centroidQ
            self.coords[:, :, f] = all1

    def append(self, mol, collisions=False, coldist=1.3):
        """ Append a molecule at the end of the current molecule

        Parameters
        ----------
        mol : :class:`Molecule`
            Target Molecule which to append to the end of the current Molecule
        collisions : bool
            If set to True it will remove residues of `mol` which collide with atoms of this Molecule object.
        coldist : float
            Collision distance in Angstrom between atoms of the two molecules. Anything closer will be considered a collision.

        Example
        -------
        >>> mol.append(Molecule('3PTB'))
        """
        if collisions:
            # Set different occupancy to separate atoms of mol1 and mol2
            occ1 = self.get('occupancy')
            occ2 = mol.get('occupancy')
            self.set('occupancy', 1)
            mol.set('occupancy', 2)

        backup = self.copy()
        mol = mol.copy()  # Copy because I'll modify its bonds
        if len(mol.bonds) > 0:
            mol.bonds += self.numAtoms

        try:
            if np.size(self.coords) != 0 and (np.size(self.coords, 2) != 1 or np.size(mol.coords, 2) != 1):
                raise NameError('Cannot concatenate molecules which contain multiple frames.')

            if len(self.bonds) > 0 and len(mol.bonds) > 0:
                self.bonds = np.append(self.bonds, mol.bonds, axis=0)
            elif len(mol.bonds) > 0:
                self.bonds = mol.bonds

            for k in self._append_fields:
                dtype = self._append_fields[k]
                if self.__dict__[k] is None or np.size(self.__dict__[k]) == 0:
                    self.__dict__[k] = np.array(mol.__dict__[k], dtype=dtype)
                elif k == 'coords':
                    self.coords = np.append(self.coords, mol.coords, axis=0)
                else:
                    self.__dict__[k] = np.append(self.__dict__[k], np.array(mol.__dict__[k], dtype=dtype))
            self.serial = np.arange(1, self.numAtoms+1)
        except:
            self = backup
            raise NameError('Failed to append molecule.')

        if collisions:
            _resolveCollisions(self, occ1, occ2, coldist)

    def atomselect(self, sel, indexes=False, strict=False):
        """ Select a set of atoms based on a selection text

        Parameters
        ----------
        sel : str
            Text selection, e.g. 'name CA'. See VMD atomselect for documentation.
        indexes : bool
            If True returns the indexes instead of a bitmap
        strict: bool
            If True it will raise an error if no atoms were selected.

        Return
        ------
        asel : np.ndarray
            Either a bitmap of selected atoms or their indexes

        Examples
        --------
        >>> a = mol.atomselect('resname MOL')
        """
        if sel is None or (isinstance(sel, str) and sel == 'all'):
            s = np.ones(self.numAtoms, dtype=bool)
        elif isinstance(sel, str):
            selc = self.coords[:, :, self.frame].copy()
            s = vmdselection(sel, selc, self.element, self.name, self.resname, self.resid,
                               chain=self.chain,
                               segname=self.segid, insert=self.insertion, altloc=self.altloc, beta=self.beta,
                               occupancy=self.occupancy)
            if np.sum(s) == 0 and strict:
                raise NameError('No atoms were selected with atom selection "{}".'.format(sel))
        else:
            s = sel

        if indexes:
            return np.where(s)[0]
        else:
            return s

    def copy(self):
        ''' Create a copy of the molecule object

        Returns
        -------
        newmol : :class:`Molecule`
            A copy of the object
        '''
        return deepcopy(self)

    def filter(self, sel, _logger=True):
        '''Removes all atoms not included in the atomselection

        Parameters
        ----------
        sel: str
            Atom selection text

        Examples
        --------
        >>> mol.filter('protein')
        '''
        s = self.atomselect(sel)
        if np.all(s):  # If all are selected do nothing
            return

        if not isinstance(s, np.ndarray) or s.dtype != bool:
            raise NameError('Filter can only work with string inputs or boolean arrays')
        self.remove(np.invert(s), _logger=_logger)

    def _removeBonds(self, idx):
        ''' Renumbers bonds after removing atoms and removes non-existent bonds

        Needs to be called before removing atoms!
        '''
        if len(self.bonds) == 0:
            return
        map = np.ones(self.numAtoms, dtype=int)
        map[idx] = -1
        map[map == 1] = np.arange(self.numAtoms - len(idx))
        bonds = np.array(self.bonds, dtype=np.int32)  # Have to store in temp because bonds is uint and can't accept -1 values
        bonds[:, 0] = map[self.bonds[:, 0]]
        bonds[:, 1] = map[self.bonds[:, 1]]
        remA = bonds[:, 0] == -1
        remB = bonds[:, 1] == -1
        stays = np.invert(remA | remB)
        # Delete bonds between non-existant atoms
        self.bonds = bonds[stays, :]

    def guessBonds(self):
        """ Tries to guess the bonds in the Molecule

        Can fail badly when non-bonded atoms are very close together. Use with extreme caution.
        """
        self.bonds = guessbonds(self.coords, self.element, self.name, self.resname, self.resid, self.chain, self.segid, self.insertion, self.altloc)

    def moveBy(self, vector, sel=None):
        '''Move a selection of molecule atoms by a given vector

        Parameters
        ----------
        vector: list
            3D coordinates to add to the Molecule coordinates
        sel: str
            Atomselection of atoms which we want to move

        Examples
        --------
        >>> mol.moveBy([3, 45 , -8])
        '''
        vector = np.array(vector)
        if np.size(vector) != 3:
            raise NameError('Move vector must be a 1x3 dimensional vector.')
        vector.shape = [1, 3]  # Forces it to be row vector

        s = self.atomselect(sel)
        for f in range(self.numFrames):
            self.coords[s, :, f] += vector

    def rotate(self, axis, angle, sel=None):
        """ Rotate molecule atoms around a given axis

        Parameters
        ----------
        axis : 3dim vector
            Axis of rotation
        angle : float
            Angle of rotation in radians
        sel :
            selection text

        Examples
        --------
        >>> mol.rotate([0, 1, 0], 1.57)
        """
        M = rotationmatrix(axis, angle)
        s = self.atomselect(sel, indexes=True)
        for a in s:
            self.coords[a, :, self.frame] = np.dot(M, self.coords[a, :, self.frame])

    def rotateBy(self, M, center=[0, 0, 0], sel='all'):
        """ Rotate a selection of atoms by a given rotation around a center

        Parameters
        ----------
        M : np.ndarray
            The rotation matrix
        center : list
            The rotation center
        """
        coords = self.get('coords', sel=sel)
        newcoords = coords - center
        newcoords = np.dot(newcoords, np.transpose(M)) + center
        self.set('coords', newcoords, sel=sel)

    def center(self, loc=[0, 0, 0], sel='all'):
        """ Moves the geometric center of the Molecule to a given location

        Parameters
        ----------
        loc : list, optional
            The location to which to move the geometric center
        sel : str
            An Atomselection string of the atoms whose geometric center we want to center on the `loc` position

        Examples
        --------
        >>> mol.center()
        >>> mol.center([10, 10, 10], 'name CA')
        """
        coords = self.get('coords', sel=sel)
        com = np.mean(coords, 0)
        self.moveBy(-com)

    def read(self, filename, type=None, skip=None, frames=None, append=False, bond=False):
        """ Read any supported file (pdb, psf, prmtop, prm, xtc, mol2, gjf, mae)

        Detects from the extension the file type and loads it into Molecule

        Parameters
        ----------
        filename : str
            Name of the file we want to read
        type : str, optional
            File type of the file. If None, it's automatically determined by the extension
        skip : int, optional
            If the file is a trajectory, skip every `skip` frames
        frames : list, optional
            If the file is a trajectory, read only the given frames
        append : bool, optional
            If the file is a trajectory, append the coordinates to the previous coordinates
        bond : bool, optional
            Guess the connectivity of the atoms
        """
        if isinstance(filename, list) or isinstance(filename, np.ndarray):
            firstfile = filename[0]
        else:
            firstfile = filename

        from htmd.simlist import Sim
        if isinstance(filename, Sim):
            self._readPDB(filename.molfile)
            self._readTraj(filename.trajectory)
            return

        if type is not None:
            type = type.lower()

        if (type is None and firstfile.endswith(".psf")) or type == "psf":
            con = PSFread(filename)
            self.charge = numpy.asarray(con.charges, dtype=np.float32)
            self.masses = numpy.asarray(con.masses, dtype=np.float32)
            self.bonds = numpy.asarray(con.bonds, dtype=np.int32)
        elif (type is None and (firstfile.endswith(".prm") or firstfile.endswith(".prmtop"))) or type == "prmtop" or type == "prm":
            con = PRMTOPread(filename)
            self.charge = numpy.asarray(con.charges, dtype=np.float32)
            # self.masses = numpy.asarray(con.masses, dtype=np.float32)  # No masses in PRMTOP
            self.bonds = numpy.asarray(con.bonds, dtype=np.int32)
        elif (type is None and firstfile.endswith(".pdb")) or type == "pdb":
            self._readPDB(filename)
        elif (type is None and firstfile.endswith(".pdbqt")) or type == "pdbqt":
            self._readPDB(filename, mode='pdbqt')
        elif (type is None and firstfile.endswith(".xtc")) or type == "xtc":
            self._readTraj(filename, skip=skip, frames=frames, append=append)
        elif (type is None and firstfile.endswith(".coor")) or type == "coor":
            self._readBinCoordinates(filename)
        elif len(firstfile) == 4:  # Could be a PDB id. Try to load it from the PDB website
            self._readPDB(filename)
        elif (type is None and firstfile.endswith(".xyz")) or type == "xyz":
            self._readXYZ(filename)
        elif (type is None and firstfile.endswith(".gjf")) or type == "gjf":
            self._readGJF(filename)
        elif (type is None and firstfile.endswith(".mae")) or type == "mae":
            self._readMae(filename)
        elif (type is None and firstfile.endswith(".mol2")) or type == "mol2":
            self._readMOL2(filename)
        else:
            try:
                self._readTraj(filename, skip=skip, frames=frames, append=append, mdtraj=True)
            except:
                raise ValueError("Unknown file type")

        if bond:
            self.guessBonds()

    def _readXYZ(self, filename):
        f = open(filename, "r")
        natoms = int(f.readline())
        for k in self._pdb_fields:
            self.__dict__[k] = numpy.zeros(natoms, dtype=self._pdb_fields[k])
        self.__dict__["coords"] = numpy.zeros((natoms, 3, 1), dtype=numpy.float32)

        f.readline()
        for i in range(natoms):
            s=f.readline().split()
            self.record[i] = "HETATM"
            self.serial[i] = i+1
            self.element[i] = s[0]
            self.name[i] = s[0]
            self.coords[i, 0, 0] = float(s[1])
            self.coords[i, 1, 0] = float(s[2])
            self.coords[i, 2, 0] = float(s[3])
            self.resname[i] = "MOL"

    def _readGJF(self, filename):
        f = open(filename, "r")
        l = f.readlines()
        start = -1
        end = -1
        c = 0
        for i in range(len(l)):
            if len(l[i].strip()) == 0 and c == 0: c = 1
            elif len(l[i].strip()) == 0 and c == 1:
                c = 2
                start = i+2
            elif len(l[i].strip()) == 0 and c == 2:
                c = 3
                end = i

        natoms = end-start
        if start == -1 or end == -1 or natoms == 0: raise ValueError( "Invalid GJF file" )

        nn = 0
        nf = self.numFrames
        if len(self.coords):
            # resize
            nn = self.numAtoms
            nf = self.numFrames
            if nn != natoms: raise ValueError("Mismatch in teh number of atoms")
            self.coords = numpy.append(self.coords, numpy.zeros((natoms, 3, 1), dtype=numpy.float32), axis=2)
            self.box = numpy.append(self.box, numpy.zeros((3, 1), dtype=numpy.float32), axis=1)
        else:
            for k in self._pdb_fields:
                self.__dict__[k] = numpy.zeros((natoms), dtype=self._pdb_fields[k])
            self.coords = numpy.zeros((natoms, 3, 1), dtype=numpy.float32)
            self.box = numpy.zeros((3, 1), dtype=numpy.float32)

        for idx in range(natoms):
            s = l[idx+start].split()
            self.record[idx] = "HETATM"
            self.serial[idx] = i+1
            self.element[idx] = s[0]
            self.name[idx] = s[0]
            self.coords[idx, 0, nf] = float(s[1])
            self.coords[idx, 1, nf] = float(s[2])
            self.coords[idx, 2, nf] = float(s[3])
            self.resname[idx] = "MOL"

    def _readPDB(self, filename, mode='pdb'):
        mol = []
        if os.path.isfile(filename):
            mol = PDBParser(filename, mode)
        elif len(filename) == 4:
            # Try loading it from the PDB website
            r = requests.get(
                "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" + filename)
            if r.status_code == 200:
                tempfile = string_to_tempfile(r.content.decode('ascii'), "pdb")
                mol = PDBParser(tempfile, mode)
                os.unlink(tempfile)
            else:
                raise NameError('Invalid PDB code')
        else:
            raise NameError('File {} not found'.format(filename))

        natoms = len(mol.record)
        for k in self._pdb_fields:
            self.__dict__[k] = numpy.asarray(mol.__dict__[k], dtype=self._pdb_fields[k])
            # Pad any short list
            if k is not "coords":
                if len(self.__dict__[k]) != natoms:
                    self.__dict__[k] = numpy.zeros(natoms, dtype=self.__dict__[k].dtype)

        self.coords = np.atleast_3d(np.array(self.coords, dtype=np.float32))
        self.bonds = np.array(mol.bonds)
        self.ssbonds = np.array(mol.ssbonds)
        self.box = np.array(mol.box)

        if self.masses is None or len(self.masses) == 0:
            self.masses = numpy.zeros(natoms, dtype=numpy.float32)
        if self.charge is None or len(self.charge) == 0:
            self.charge = mol.charge.copy()
        if self.charge is None or len(self.charge) == 0:
            self.charge = numpy.zeros(natoms, dtype=numpy.float32)

        self.fileloc.append([filename, 0])

    def _readMOL2(self, filename):
        f = open(filename, "r")
        l = f.readlines()
        f.close()
        start = None
        end = None
        for i in range(len(l)):
            if l[i].startswith("@<TRIPOS>ATOM"): start = i+1
            if l[i].startswith("@<TRIPOS>BOND"): end = i-1

        if not start or not end:
            raise ValueError("File cannot be read")

        natoms = end - start + 1

        for k in self._pdb_fields:
            self.__dict__[k] = numpy.zeros((natoms), dtype=self._pdb_fields[k])
        self.__dict__["coords"] = numpy.zeros((natoms, 3, 1), dtype=numpy.float32)

        for i in range(natoms):
            s=l[i+start].strip().split()
            self.record[i] = "HETATM"
            self.serial[i] = int(s[0])
            self.element[i] = re.sub("[0123456789]*", "", s[1])
            self.name[i] = s[1]
            self.coords[i, 0, 0] = float(s[2])
            self.coords[i, 1, 0] = float(s[3])
            self.coords[i, 2, 0] = float(s[4])
            self.resname[i] = "MOL"

    def _readMae(self, filename):
        datadict = _maestroparser(filename)
        natoms = len(datadict['record'])
        for k in self._pdb_fields:
            self.__dict__[k] = numpy.asarray(datadict[k], dtype=self._pdb_fields[k])
            # Pad any short list
            if k is not "coords":
                if len(self.__dict__[k]) != natoms:
                    self.__dict__[k] = numpy.zeros(natoms, dtype=self.__dict__[k].dtype)

        self.segid = self.segid.astype('str').astype('object')  # Making sure the segids are strings because np.zeros gives int objects
        self.coords = np.atleast_3d(self.coords)
        for h in datadict['het']:
            self.set('record', 'HETATM', sel='resname {}'.format(h))

        # self.serial = np.arange(1, natoms+1)
        self.masses = np.array(datadict['masses'])
        self.bonds = np.array(datadict['bonds']) - 1  # convert to 0 indexing

        fnamestr = os.path.splitext(os.path.basename(filename))[0]
        self.viewname = fnamestr
        self.fileloc = [[fnamestr, 0]]

        self.box = np.max(self.coords, axis=0) - np.min(self.coords, axis=0)

    def _readBinCoordinates(self, filename):
        self.coords = BINCOORread(filename)

    def _readMDtraj(self, filename):
        class Struct:
            def __init__(self):
                return
        import mdtraj as md
        traj = md.load(filename, top=self.topoloc)
        s = Struct()
        s.coords = np.swapaxes(np.swapaxes(traj.xyz, 0, 1), 1, 2) * 10
        if traj.timestep == 1:
            s.time = np.zeros(traj.time.shape, dtype=traj.time.dtype)
            s.step = np.zeros(traj.time.shape, dtype=traj.time.dtype)
        else:
            s.time = traj.time
            s.step = s.time / 25  # DO NOT TRUST THIS. I just guess that there are 25 simulation steps in each picosecond
        s.box = traj.unitcell_lengths.T * 10
        return s

    def _readTraj(self, filename, skip=None, frames=None, append=False, mdtraj=False):
        if not append:
            self.coords = []
            self.box = []
        else:
            logger.warning('Appending trajectories not well tested yet')

        # If a single filename is specified, turn it into an array so we can iterate
        if isinstance(filename, str):
            filename = [filename]
        if not isinstance(filename, np.ndarray):
            filename = np.array(filename)
        # print(len(filename), len(frames), type(frames))

        #from IPython.core.debugger import Tracer
        #Tracer()()
        if frames is not None:
            if not isinstance(frames, list) and not isinstance(frames, np.ndarray):
                frames = [frames]
            if len(filename) != len(frames):
                raise NameError('Number of trajectories (' + str(len(filename)) + ') does not match number of frames (' + str(len(frames)) + ') given as arguments')

        self.fileloc = []
        for i, f in enumerate(filename):
            if frames is None:
                if mdtraj:
                    traj = self._readMDtraj(f)
                else:
                    traj = XTCread(f)
                for j in range(np.size(traj.coords, 2)):
                    self.fileloc.append([f, j])
            else:
                if mdtraj:
                    traj = self._readMDtraj(f)
                    traj.coords = traj.coords[:, :, frames[i]]
                else:
                    traj = XTCread(f, frames[i])
                self.fileloc.append([f, int(frames[i])])

            if self.numAtoms != 0 and np.size(traj.coords, 0) != self.numAtoms:
                raise ValueError('Trajectory # of atoms ' + str(np.size(self.coords, 0)) + ' mismatch with # of already loaded atoms ' + str(self.numAtoms))

            if len(self.coords) > 0 and self.coords.shape[0] > 0 and (self.coords.shape[0] != traj.coords.shape[0]):
                raise ValueError("Trajectory # of atoms mismatch with already loaded coordinates")
            # print(np.shape(traj.box), np.shape(self.box))
            # TODO : check step correct increment
            if len(self.coords) == 0:
                self.coords = traj.coords
                self.box = traj.box
            else:
                self.coords = np.append(self.coords, traj.coords, 2)
                self.box = np.append(self.box, traj.box, 1)

        if skip is not None:
            self.coords = self.coords[:, :, ::skip]  # Might actually not free memory! Check numpy views
            self.box = self.box[:, ::skip]
            self.fileloc = self.fileloc[::skip]

        self.coords = np.atleast_3d(self.coords)
        self.step = traj.step
        self.time = traj.time
        if len(traj.time) < 2:
            #logger.info('Trajectory has broken framestep. Cannot read correctly, setting to 0.1ns.')
            self.fstep = 0.1
        else:
            self.fstep = (traj.time[1] - traj.time[0]) / 1E6  # convert femtoseconds to nanoseconds

    def view(self, sel=None, style=None, color=None, guessbonds=True, viewer=None, hold=False, name=None, viewerhandle=None):
        """ Visualizes the molecule in a molecular viewer

        Parameters
        ----------
        sel : str
            Atomselection string for a representation.
        style : str
            Representation style.
        color : str
            Coloring mode or color ID.
        guessbonds : bool
            Allow VMD to guess bonds for the molecule
        viewer : str ('vmd','notebook')
            Choose viewer backend. Default is taken from htmd.config
        hold : bool
            If set to True, it will not visualize the molecule but instead collect representations until set back to False.
        name : str, optional
            A name to give to the molecule in VMD
        viewerhandle : :class:`VMD <htmd.vmdviewer.VMD>` object, optional
            A specific viewer in which to visualize the molecule. If None it will use the current default viewer.
        """
        from htmd.util import tempname

        if sel is not None or style is not None or color is not None:
            self._tempreps.add(sel=sel, style=style, color=color)

        if hold:
            return

        # Write out PDB and XTC files
        pdb = tempname(suffix=".pdb")
        self.write(pdb)
        xtc = None
        if self.numFrames > 1:
            xtc = tempname(suffix=".xtc")
            self.write(xtc)

        # Call the specified backend
        if viewer is None:
            from htmd.config import _config
            viewer = _config['viewer']
        if viewer.lower() == 'notebook':
            return self._viewMDTraj(pdb, xtc)
        elif viewer.lower() == 'vmd':
            self._viewVMD(pdb, xtc, viewerhandle, name, guessbonds)
        elif viewer.lower() == 'ngl':
            return self._viewNGL(pdb, xtc, guessbonds)

        # Remove temporary files
        if xtc:
            os.remove(xtc)
        os.remove(pdb)

    def _viewVMD(self, pdb, xtc, vhandle, name, guessbonds):
        if name is None:
            name = self.viewname
        if vhandle is None:
            vhandle = getCurrentViewer()

        if guessbonds:
            vhandle.send("mol new " + pdb)
        else:
            vhandle.send("mol new " + pdb + " autobonds off")

        if name is not None:
            vhandle.send('mol rename top "' + name + '"')
        else:
            vhandle.send('mol rename top "Mol [molinfo top]: pdb"')

        if xtc:
            vhandle.send('animate delete all')
            vhandle.send('mol addfile ' + xtc + ' type xtc waitfor all')
            if name is None:
                vhandle.send('mol rename top "Mol [molinfo top]: pdb+xtc"')

        self._tempreps.append(self.reps)
        self._tempreps._repsVMD(vhandle)
        self._tempreps.remove()

    def _viewMDTraj(self, pdb, xtc):
        from mdtraj.html import TrajectoryView, TrajectorySliderView, enable_notebook
        import mdtraj
        enable_notebook()

        if xtc:
            t = mdtraj.load(xtc, top=pdb)
            widget = TrajectorySliderView(t)
        else:
            t = mdtraj.load(pdb)
            widget = TrajectoryView(t)
        return widget

    def _viewNGL(self, pdb, xtc, guessb):
        from nglview import Trajectory
        import nglview

        class TrajectoryStreamer(Trajectory):
            def __init__(self, coords):
                self.coords = coords
            def get_coordinates_list(self, index):
                return self.coords[:, :, index].flatten().tolist()
            def get_frame_count(self):
                return np.size(self.coords, 2)

        struc = nglview.FileStructure(pdb)
        struc.params['dontAutoBond'] = not guessb
        if xtc:
            traj = TrajectoryStreamer(self.coords)
            w = nglview.NGLWidget(struc, traj)
        else:
            w = nglview.NGLWidget(struc)

        self._tempreps.append(self.reps)
        self._tempreps._repsNGL(w)
        self._tempreps.remove()
        return w

    def mutateResidue(self, sel, newres):
        """ Mutates a residue by deleting it's sidechain and renaming it

        Parameters
        ----------
        sel : str
            Atomselection for the residue we want to mutate. The selection needs to include all atoms of the residue.
        newres : str
            The name of the new residue

        Examples
        --------
        >>> mol.mutateResidue('resid 158', 'ARG')
        """
        s = self.atomselect(sel, strict=True)
        # Changed the selection from "and sidechain" to "not backbone" to remove atoms like phosphates which are bonded
        # but not part of the sidechain. Changed again the selection to "name C CA N O" because "backbone" works for
        # both protein and nucleic acid backbones and it confuses phosphates of modified residues for nucleic backbones.
        removed = self.remove(sel + ' and not name C CA N O', _logger=False)
        s = np.delete(s, removed)
        self.set('resname', newres, sel=s)

    def wrap(self, wrapsel=None):
        """ Wraps coordinates of the molecule into the simulation box

        Parameters
        ----------
        wrapsel : str
            Selection of atoms on which to center the wrapping box

        Examples
        --------
        >>> mol.wrap()
        >>> mol.wrap('protein')
        """
        # TODO: selection is not used. WHY?
        if len(self.bonds) == 0:
            bonds = guessbonds(self.coords, self.element, self.name, self.resname, self.resid, self.chain, self.segid, self.insertion, self.altloc)
        else:
            bonds = np.append(self.bonds, guessbonds(self.coords, self.element, self.name, self.resname, self.resid, self.chain, self.segid, self.insertion, self.altloc), axis=0)
        '''# Duplicating bonds in reverse to see if it helps
        uqatms = np.unique(bonds)
        for a in uqatms:
            idx = bonds[bonds[:, 0] == a, 1]
            idx = np.unique(np.append(idx, bonds[bonds[:, 1] == a, 0]))
            for i in idx:
                if not np.any(np.all(bonds == [a, i], axis=1)):
                    bonds = np.append(bonds, [[a, i]], axis=0)
                if not np.any(np.all(bonds == [i, a], axis=1)):
                    bonds = np.append(bonds, [[i, a]], axis=0)'''
        self.coords = wrap(self.coords, bonds, self.box)

    def write(self, filename, sel=None, type=None):
        """ Writes any of the supported formats (pdb, coor, psf, xtc)

        Parameters
        ----------
        filename : str
            The filename of the file we want to write to disk
        sel : str, optional
            The atomselections of the atoms we want to write. If None it will write all atoms
        type : str, optional
            The filetype we want to write. By default, detected from the file extension
        """
        if type:
            type = type.lower()

        if type == "coor" or filename.endswith(".coor"):
            self._writeBinCoordinates(filename, sel)
        elif type == "pdb" or filename.endswith(".pdb"):
            self._writePDB(filename, sel)
        elif type == "xyz" or filename.endswith(".xyz"):
            self._writeXYZ(filename, sel)
        elif type == "psf" or filename.endswith(".psf"):
            self._writeConnectivity(filename, sel)
        elif type == "xtc" or filename.endswith(".xtc"):
            self._writeTraj(filename, sel)
        else:
            try:
                import mdtraj as md
                from htmd.util import tempname
                tmppdb = tempname(suffix='.pdb')
                tmpxtc = tempname(suffix='.xtc')
                self.write(tmppdb)
                self.write(tmpxtc)
                traj = md.load(tmpxtc, top=tmppdb)
                #traj.xyz = np.swapaxes(np.swapaxes(self.coords, 1, 2), 0, 1) / 10
                #traj.time = self.time
                #traj.unitcell_lengths = self.box.T / 10
                traj.save(filename)
            except:
                raise ValueError("Unknown file type")

    def _writeXYZ( self, filename, sel="all" ):
        src = self
        if sel is not None: 
          src = sel.copy()
          src.filter(sel, _logger=False)
        fh = open( filename, "w" )
        natoms = len(src.record)
        print( "%d\n" % (natoms), file=fh )
        for i in range(natoms):
          e = src.element[i].strip()
          if( not len(e) ):
             e = re.sub( "[1234567890]*", "", src.name[i] )
          print("%s   %f   %f    %f" % ( e, src.coords[i,0,0], src.coords[i,1,0], src.coords[i,2,0] ), file=fh )
        fh.close()

    def _writeBinCoordinates(self, filename, sel):
        if self.frame < 0 or self.frame > self.numFrames:
            raise NameError("frame out of range")
        mol = self.copy()
        mol.coords = mol.coords[:, :, self.frame]
        mol.coords = np.atleast_3d(mol.coords.reshape((mol.coords.shape[0], 3, 1)))
        if sel is not None: mol.filter(sel, _logger=False)
        # Bincoor is in angstrom
        BINCOORwrite(mol.coords, filename)

    def _writeConnectivity(self, filename, sel):
        src = self
        if sel is not None:
            src = self.copy()
            src.filter(sel, _logger=False)
        PSFwrite(src, filename)

    def _writePDB(self, filename, sel='all'):
        src = self
        if sel is not None and sel != 'all':
            src = self.copy()
            src.filter(sel, _logger=False)

        pdb = PDBParser()
        for k in self._pdb_fields:
            pdb.__dict__[k] = src.__dict__[k].copy()

        pdb.coords = np.atleast_3d(pdb.coords[:, :, self.frame]) # Writing out only current frame

        pdb.bonds = src.bonds
        pdb.ssbonds = src.ssbonds  # TODO: Is there such a thing in pdb format?
        pdb.box = self.box

        pdb.serial = np.arange(1, np.size(pdb.coords, 0)+1)
        pdb.writePDB(filename)

    def _writeTraj(self, filename, sel):
        # Write xtc
        src = self
        if sel is not None:
            src = self.copy()
            src.filter(sel, _logger=False)
        if np.size(src.box, 1) != self.numFrames:
            src.box = np.tile(src.box, (1, self.numFrames))
        XTCwrite(src.coords, src.box, filename, self.time, self.step)

    def empty(self, numAtoms):
        """ Creates an empty molecule of N atoms.

        Parameters
        ----------
        numAtoms : int
            Number of atoms to create in the molecule.

        Example
        -------
        >>> mol = Molecule()
        >>> mol.empty(100)
        """
        self.record = np.array(['ATOM'] * numAtoms, dtype=self._pdb_fields['record'])
        self.chain = np.array(['X'] * numAtoms, dtype=self._pdb_fields['chain'])
        self.segid = np.array(['X'] * numAtoms, dtype=self._pdb_fields['segid'])
        self.occupancy = np.array([0] * numAtoms, dtype=self._pdb_fields['occupancy'])
        self.beta = np.array([0] * numAtoms, dtype=self._pdb_fields['beta'])
        self.insertion = np.array([''] * numAtoms, dtype=self._pdb_fields['insertion'])
        self.element = np.array([''] * numAtoms, dtype=self._pdb_fields['element'])
        self.altloc = np.array([''] * numAtoms, dtype=self._pdb_fields['altloc'])
        self.name = np.array(['UNK'] * numAtoms, dtype=self._pdb_fields['name'])
        self.resname = np.array(['UNK'] * numAtoms, dtype=self._pdb_fields['resname'])
        self.resid = np.array([999] * numAtoms, dtype=self._pdb_fields['resid'])
        self.coords = np.zeros((numAtoms, 3, 1), dtype=self._pdb_fields['coords'])
        self.serial = np.arange(1, numAtoms+1)


    def sequence(self, oneletter=True):
        """ Return the AA sequence of the Molecule.

        Parameters
        ----------
        oneletter : bool
            Whether to return one-letter or three-letter AA codes. There should be only one atom per residue.

        Returns
        -------
        sequence : str
            The primary sequence as a string

        Examples
        --------
        >>> m=Molecule("3PTB"); m.filter("protein")
        >>> m.sequence()
        {'0': 'IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN'}
        """
        residues = {'ARG': 'R', 'AR0': 'R',
                    'HIS': 'H', 'HID': 'H', 'HIE': 'H',
                    'LYS': 'K', 'LSN': 'K', 'LYN': 'K',
                    'ASP': 'D', 'ASH': 'D',
                    'GLU': 'E', 'GLH': 'E',
                    'SER': 'S',
                    'THR': 'T',
                    'ASN': 'N',
                    'GLN': 'Q',
                    'CYS': 'C', 'CYX': 'C',
                    'SEC': 'U',
                    'GLY': 'G',
                    'PRO': 'P',
                    'ALA': 'A',
                    'VAL': 'V',
                    'ILE': 'I',
                    'LEU': 'L',
                    'MET': 'M',
                    'PHE': 'F',
                    'TYR': 'Y',
                    'TRP': 'W'}
        from htmd.builder.builder import sequenceID
        prot = self.atomselect('protein')
        segs = np.unique(self.segid[prot])
        increm = sequenceID((self.resid, self.insertion, self.chain))
        sequence = {}
        olsequence = {}
        for seg in segs:
            sequence[seg] = []
            olsequence[seg] = ''
            segatoms = self.atomselect('protein and segid {}'.format(seg))
            resnames = self.resname[segatoms]
            incremseg = increm[segatoms]
            for i in np.unique(incremseg):
                resname = np.unique(resnames[incremseg == i])
                if len(resname) != 1:
                    raise AssertionError('Something went wrong here.')
                sequence[seg].append(resname)
                olsequence[seg] += residues[resname[0]]
        if oneletter:
            return olsequence
        else:
            return sequence

    @property
    def numFrames(self):
        """ Number of coordinate frames in the molecule
        """
        return np.size(np.atleast_3d(self.coords), 2)

    @property
    def numAtoms(self):
        """ Number of atoms in the molecule
        """
        return np.size(self.coords, 0)

    @property
    def x(self):
        """Get the x coordinates at the current frame"""
        return self.coords[:,0,self.frame]

    @property
    def y(self):
        """Get the y coordinates at the current frame"""
        return self.coords[:,1,self.frame]

    @property
    def z(self):
        """Get the z coordinates at the current frame"""
        return self.coords[:,2,self.frame]

    def __str__(self):
        def formatstr(name, field):
            if isinstance(field, np.ndarray) or isinstance(field, list):
                rep = '{} shape: {}'.format(name, np.shape(field))
            elif field == 'reps':
                rep = '{}: {}'.format(name, len(self.reps.replist))
            else:
                rep = '{}: {}'.format(name, field)
            return rep

        rep = 'Molecule with ' + str(self.numAtoms) + ' atoms and ' + str(self.numFrames) + ' frames'
        for p in sorted(self._pdb_fields):
            rep += '\n'
            rep += 'PDB field - ' + formatstr(p, self.__dict__[p])
        for j in sorted(self.__dict__.keys() - self._pdb_fields.keys()):
            if j[0] == '_':
                continue
            rep += '\n'
            rep += formatstr(j, self.__dict__[j])

        return rep


def _resolveCollisions(mol, occ1, occ2, gap):
    from htmd.molecule.util import sequenceID

    s1 = mol.atomselect('occupancy 1')
    s2 = mol.atomselect('occupancy 2')
    # Give unique "residue" beta number to all resids
    beta = mol.get('beta')
    mol.set('beta', sequenceID(mol.resid))
    # Calculate overlapping atoms
    overlaps = mol.atomselect('(occupancy 2) and same beta as exwithin {} of (occupancy 1)'.format(gap))
    # Restore original beta and occupancy
    mol.set('beta', beta)
    mol.set('occupancy', occ1, s1)
    mol.set('occupancy', occ2, s2)

    logger.info('Removed {} residues from appended Molecule due to collisions.'.format(len(np.unique(mol.resid[overlaps]))))
    # Remove the overlaps
    mol.remove(overlaps, _logger=False)


def _pp_measure_fit(P, Q):
    """
    PP_MEASURE_FIT - molecule alignment function.
    For documentation see http://en.wikipedia.org/wiki/Kabsch_algorithm
    the Kabsch algorithm is a method for calculating the optimal
    rotation matrix that minimizes the RMSD (root mean squared deviation)
    between two paired sets of points
    """
    centroidP = np.mean(P, 0)
    centroidQ = np.mean(Q, 0)

    # Centering them on 0,0,0
    # Can also be done with translation matrix if it's faster
    P = P - centroidP
    Q = Q - centroidQ

    covariance = np.dot(np.transpose(P), Q)

    (V, S, W) = np.linalg.svd(covariance)  # Matlab svd returns the W transposed compared to numpy.svd
    W = np.transpose(W)

    E0 = np.sum(np.sum(P * P)) + np.sum(np.sum(Q * Q))
    RMSD = E0 - (2 * np.sum(S.ravel()))
    RMSD = np.sqrt(np.abs(RMSD / np.size(P, 0)))

    d = np.sign(np.linalg.det(W) * np.linalg.det(V))
    z = np.eye(3)
    z[2, 2] = d
    U = np.dot(np.dot(W, z), np.transpose(V))
    return U, RMSD


class Representations:
    """ Class that stores representations for Molecule.

    Examples
    --------
    >>> from htmd.molecule.molecule import Molecule
    >>> mol = Molecule('3PTB')
    >>> mol.reps.add('protein', 'NewCartoon')
    >>> print(mol.reps)                     # doctest: +NORMALIZE_WHITESPACE
    rep 0: sel='protein', style='NewCartoon', color='Name'
    >>> mol.view()
    >>> mol.reps.remove()
    """
    def __init__(self, mol):
        self.replist = []
        self._mol = mol
        return

    def append(self, reps):
        if not isinstance(reps, Representations):
            raise NameError('You can only append Representations objects.')
        self.replist += reps.replist

    def add(self, sel=None, style=None, color=None):
        """ Adds a new representation for Molecule.

        Parameters
        ----------
        sel : str
            Atom selection for the given representation.
        style : str
            Representation visual style.
        color : str
            Color style or ID.
        """
        self.replist.append(_Representation(sel, style, color))

    def remove(self, index=None):
        """ Removed one or all representations.

        Parameters
        ----------
        index : int
            The index of the representation to delete. If none is given it deletes all.
        """
        if index is None:
            self.replist = []
        else:
            del self.replist[index]

    def list(self):
        """ Lists all representations. Equivalent to using print.
        """
        print(self)

    def __str__(self):
        s = ''
        for i, r in enumerate(self.replist):
            s += ('rep {}: sel=\'{}\', style=\'{}\', color=\'{}\'\n'.format(i, r.sel, r.style, r.color))
        return s

    def _translateNGL(self, rep):
        styletrans = {'newcartoon': 'cartoon', 'licorice': 'hyperball', 'lines': 'line', 'vdw': 'spacefill', 'cpk': 'ball+stick'}
        colortrans = {'name': 'element', 'index': 'atomindex', 'chain': 'chainindex', 'secondary structure': 'sstruc', 'colorid': 'color'}
        hexcolors = {0: '#0000ff', 1: '#ff0000', 2: '#333333', 3: '#ff6600', 4: '#ffff00', 5: '#4c4d00', 6: '#b2b2cc', 7: '#33cc33', 8: '#ffffff', 9: '#ff3399', 10: '#33ccff'}
        try:
            selidx = '@' + ','.join(map(str, self._mol.atomselect(rep.sel, indexes=True)))
        except:
            return None
        if rep.style.lower() in styletrans:
            style = styletrans[rep.style.lower()]
        else:
            style = rep.style
        if isinstance(rep.color, int):
            color = hexcolors[rep.color]
        elif rep.color.lower() in colortrans:
            color = colortrans[rep.color.lower()]
        else:
            color = rep.color
        return _Representation(sel=selidx, style=style, color=color)

    def _repsVMD(self, viewer):
        colortrans = {'secondary structure': 'Structure'}
        if len(self.replist) > 0:
            viewer.send('mol delrep 0 top')
            for rep in self.replist:
                if isinstance(rep.color, str) and rep.color.lower() in colortrans:
                    color = colortrans[rep.color.lower()]
                else:
                    color = rep.color
                viewer.send('mol selection {}'.format(rep.sel))
                viewer.send('mol representation {}'.format(rep.style))
                if isinstance(rep.color, int):
                    viewer.send('mol color ColorID {}'.format(color))
                else:
                    viewer.send('mol color {}'.format(color))
                viewer.send('mol addrep top')

    def _repsNGL(self, viewer):
        if len(self.replist) > 0:
            reps = []
            for r in self.replist:
                r2 = self._translateNGL(r)
                if r2 is not None:
                    reps.append({"type": r2.style, "params": {"sele": r2.sel, "color": r2.color}})
            if reps != []:
                viewer.representations = reps


def _maestroparser(fname):
    """ Reads maestro files.

    Parameters
    ----------
    fname : str
        .mae file

    Returns
    -------
    mol : Molecule

    """
    section = None
    section_desc = False
    section_data = False

    data = {}
    data['serial'] = []
    data['record'] = []
    data['name'] = []
    data['resname'] = []
    data['resid'] = []
    data['chain'] = []
    data['segid'] = []
    data['occupancy'] = []
    data['beta'] = []
    data['insertion'] = []
    data['element'] = []
    data['altloc'] = []
    data['coords'] = []
    data['bonds'] = []
    data['charge'] = []
    data['masses'] = []
    data['het'] = []

    import csv
    with open(fname, newline='') as fp:
        reader = csv.reader(fp, delimiter=' ', quotechar='"', skipinitialspace=True)
        for row in reader:
            if len(row) == 0:
                continue

            if row[0][0:6] == 'm_atom':
                section = 'atoms'
                section_desc = True
                section_cols = []
            elif row[0][0:6] == 'm_bond':
                section = 'bonds'
                section_desc = True
                section_cols = []
            elif row[0][0:18] == 'm_PDB_het_residues':
                section = 'hetresidues'
                section_desc = True
                section_cols = []
            elif section_desc and row[0] == ':::':
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
                    data['record'].append('ATOM')
                    row = np.array(row)
                    row[row == '<>'] = 0
                    if 'i_pdb_PDB_serial' in section_dict:
                        data['serial'].append(row[section_dict['i_pdb_PDB_serial']])
                    if 's_m_pdb_atom_name' in section_dict:
                        data['name'].append(row[section_dict['s_m_pdb_atom_name']].strip())
                    if 's_m_pdb_residue_name' in section_dict:
                        data['resname'].append(row[section_dict['s_m_pdb_residue_name']].strip())
                    if 'i_m_residue_number' in section_dict:
                        data['resid'].append(int(row[section_dict['i_m_residue_number']]))
                    if 's_m_chain_name' in section_dict:
                        data['chain'].append(row[section_dict['s_m_chain_name']])
                    if 's_pdb_segment_id' in section_dict:
                        data['segid'].append(row[section_dict['s_pdb_segment_id']])
                    if 'r_m_pdb_occupancy' in section_dict:
                        data['occupancy'].append(float(row[section_dict['r_m_pdb_occupancy']]))
                    if 'r_m_pdb_tfactor' in section_dict:
                        data['beta'].append(float(row[section_dict['r_m_pdb_tfactor']]))
                    if 's_m_insertion_code' in section_dict:
                        data['insertion'].append(row[section_dict['s_m_insertion_code']].strip())
                    if '' in section_dict:
                        data['element'].append('')  # TODO: Read element
                    if '' in section_dict:
                        data['altloc'].append('')  # TODO: Read altloc. Quite complex actually. Won't bother.
                    if 'r_m_x_coord' in section_dict:
                        data['coords'].append([float(row[section_dict['r_m_x_coord']]), float(row[section_dict['r_m_y_coord']]), float(row[section_dict['r_m_z_coord']])])
                    data['masses'].append(0)

                # Reading the data of the bonds section
                if section == 'bonds' and section_data:
                    data['bonds'].append([int(row[section_dict['i_m_from']]), int(row[section_dict['i_m_to']])])

                # Reading the data of the hetero residue section
                if section == 'hetresidues' and section_data:
                    data['het'].append(row[section_dict['s_pdb_het_name']].strip())
    return data


class _Representation:
    """ Class that stores a representation for Molecule

    Parameters
    ----------
    sel : str
        Atom selection for the given representation.
    style : str
        Representation visual style.
    color : str
        Color style.
    colorid: int
        Color ID, if `color` was set to 'colorID'.

    Examples
    --------
    >>> r = _Representation(sel='protein', style='NewCartoon', color='Index')
    >>> r = _Representation(sel='resname MOL', style='Licorice')
    >>> r = _Representation(sel='ions', style='VDW', color=1)
    """
    def __init__(self, sel=None, style=None, color=None):
        if sel is not None:
            self.sel = sel
        else:
            self.sel = 'all'
        if style is not None:
            self.style = style
        else:
            self.style = 'Lines'
        if color is not None:
            self.color = color
        else:
            self.color = 'Name'


if __name__ == "__main__":

    mol = Molecule('3PTB')
    mol_backup = mol.copy()

    import doctest
    doctest.testmod(extraglobs={'mol': mol})

    # Oddly, if these are moved before doctests, 1. extraglobs don't work; and 2. test failures are not printed.

    mol = mol_backup
    a = mol.get('resid', sel='resname TRP')
    a = mol.get('coords')
    print(a.ndim)
    mol.write('/tmp/test.pdb')
    mol.write('/tmp/test.coor')
    mol.write('/tmp/test.xtc')
    mol.moveBy([1, 1, 1])
    mol.rotate([1, 0, 0], pi / 2)
    mol.align('name CA')

    # test rotate
