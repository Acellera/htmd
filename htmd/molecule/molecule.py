# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import requests
import numpy as np
from htmd.molecule.pdbparser import PDBParser
from htmd.molecule.vmdparser import guessbonds, vmdselection
from htmd.molecule.readers import XTCread, CRDread, BINCOORread, PRMTOPread, PSFread, MAEread, MOL2read, GJFread, XYZread
from htmd.molecule.writers import XTCwrite, PSFwrite, BINCOORwrite
from htmd.molecule.support import string_to_tempfile
from htmd.molecule.wrap import *
from htmd.rotationmatrix import rotationMatrix
from htmd.vmdviewer import getCurrentViewer
from math import pi
from copy import deepcopy
from os import path
import logging
import re

logger = logging.getLogger(__name__)


class TopologyInconsistencyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

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
    angles : np.ndarray
        Angle terms, valid only if PSF read and molecule unmodified
    dihedrals : np.ndarray
        Dihedral terms, valid only if PSF read and molecule unmodified
    impropers : np.ndarray
        Improper terms, valid only if PSF read and molecule unmodified
    atomtype : np.ndarray
        Atom types, valid only if PSF read and molecule unmodified

    """
    _pdb_fields = ['record', 'serial', 'name', 'altloc', 'resname', 'chain', 'resid', 'insertion', 'coords',
                   'occupancy', 'beta', 'segid', 'element', 'charge']
    _append_fields = _pdb_fields + ['masses']

    _dtypes = {
        'record': object,
        'serial': np.int,
        'name': object,
        'altloc': object,
        'resname': object,
        'chain': object,
        'resid': np.int,
        'insertion': object,
        'coords': np.float32,
        'occupancy': np.float32,
        'beta': np.float32,
        'segid': object,
        'element': object,
        'charge': np.float32,
        'bonds': np.uint32,
        'angles': np.uint32,
        'dihedrals': np.uint32,
        'impropers': np.uint32,
        'atomtype': object,
        'masses': np.float32,
        'box': np.float32
    }

    _dims = {
        'record': 0,
        'serial': 0,
        'name': 0,
        'altloc': 0,
        'resname': 0,
        'chain': 0,
        'resid': 0,
        'insertion': 0,
        'coords': (0, 3, 0),
        'occupancy': 0,
        'beta': 0,
        'segid': 0,
        'element': 0,
        'charge': 0,
        'bonds': (0, 2),
        'angles': (0, 3),
        'dihedrals': (0, 4),
        'impropers': (0, 4),
        'atomtype': 0,
        'masses': 0,
        'box': (3, 1)
    }

    def __init__(self, filename=None, name=None):
        for field in self._dtypes:
            self.__dict__[field] = np.empty(self._dims[field], dtype=self._dtypes[field])
        self.ssbonds = []
        self.frame = 0
        self.fileloc = []
        self.time = []
        self.step = []

        self.reps = Representations(self)
        self._tempreps = Representations(self)
        self.viewname = name

        if filename:
            self.read(filename)
            if isinstance(filename, str):
                self.topoloc = os.path.abspath(filename)
                if name is None and isinstance(filename, str):
                    self.viewname = filename
                    if path.isfile(filename):
                        self.viewname = path.basename(filename)

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
        >>> mol=tryp.copy()
        >>> mol.numAtoms
        1701
        >>> mol.insert(tryp, 0)
        >>> mol.numAtoms
        3402
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
                    self.__dict__[k] = np.insert(self.__dict__[k], index,
                                                 np.zeros(numatoms, dtype=self.__dict__[k].dtype), axis=0)
                elif k == 'coords':
                    self.coords = np.insert(self.coords, index, np.atleast_3d(mol.coords), axis=0)
                else:
                    self.__dict__[k] = np.insert(self.__dict__[k], index, mol.__dict__[k], axis=0)
            self.serial = np.arange(1, self.numAtoms + 1)
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
        >>> mol=tryp.copy()
        >>> mol.remove('name CA')               # doctest: +ELLIPSIS
        array([   1,    9,   16,   20,   24,   36,   43,   49,   53,   58,...
        """
        sel = self.atomselect(selection, indexes=True)
        self._removeBonds(sel)
        for k in self._append_fields:
            self.__dict__[k] = np.delete(self.__dict__[k], sel, axis=0)
        if _logger:
            logger.info('Removed {} atoms. {} atoms remaining in the molecule.'.format(len(sel), self.numAtoms))

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
        >>> mol=tryp.copy()
        >>> mol.get('resname')
        array(['ILE', 'ILE', 'ILE', ..., 'HOH', 'HOH', 'HOH'], dtype=object)
        >>> mol.get('resname', sel='resid 158')
        array(['LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU'], dtype=object)


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
        >>> mol=tryp.copy()
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
        >>> mol=tryp.copy()
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
        >>> mol=tryp.copy()
        >>> mol.filter("not resname BEN")
        >>> lig=tryp.copy()
        >>> lig.filter("resname BEN")
        >>> mol.append(lig)
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
                if self.__dict__[k] is None or np.size(self.__dict__[k]) == 0:
                    self.__dict__[k] = np.array(mol.__dict__[k], dtype=self._dtypes[k])
                elif k == 'coords':
                    self.coords = np.append(self.coords, mol.coords, axis=0)
                else:
                    self.__dict__[k] = np.append(self.__dict__[k], np.array(mol.__dict__[k], dtype=self._dtypes[k]))
            self.serial = np.arange(1, self.numAtoms + 1)
        except Exception as err:
            self = backup
            raise NameError('Failed to append molecule. "{}"'.format(err))

        if collisions:
            _resolveCollisions(self, occ1, occ2, coldist)

    def _getBonds(self, fileBonds=True, guessBonds=True):
        """ Returns an array of all bonds.

        Parameters
        ----------
        fileBonds : bool
            If True will use bonds read from files.
        guessBonds : bool
            If True will use guessed bonds.

        Returns
        -------
        bonds : np.ndarray
            An array of bonds
        """
        bonds = np.empty((0, 2), dtype=np.uint32)
        if fileBonds:
            if len(self.bonds) == 0:  # This is a patch for the other readers not returning correct empty dimensions
                self.bonds = np.empty((0, 2), dtype=numpy.uint32)
            bonds = numpy.vstack((bonds, self.bonds))
        if guessBonds:
            bonds = numpy.vstack((bonds, self._guessBonds()))
        return bonds

    def atomselect(self, sel, indexes=False, strict=False, fileBonds=True, guessBonds=True):
        """ Select a set of atoms based on a selection text

        Parameters
        ----------
        sel : str
            Text selection, e.g. 'name CA'. See VMD atomselect for documentation.
        indexes : bool
            If True returns the indexes instead of a bitmap
        strict: bool
            If True it will raise an error if no atoms were selected.
        fileBonds : bool
            If True will use bonds read from files.
        guessBonds : bool
            If True will use guessed bonds.

        Return
        ------
        asel : np.ndarray
            Either a bitmap of selected atoms or their indexes

        Examples
        --------
        >>> mol=tryp.copy()
        >>> mol.atomselect('resname MOL')
        array([False, False, False, ..., False, False, False], dtype=bool)
        """
        if sel is None or (isinstance(sel, str) and sel == 'all'):
            s = np.ones(self.numAtoms, dtype=bool)
        elif isinstance(sel, str):
            selc = self.coords[:, :, self.frame].copy()
            s = vmdselection(sel, selc, self.element, self.name, self.resname, self.resid,
                             chain=self.chain,
                             segname=self.segid, insert=self.insertion, altloc=self.altloc, beta=self.beta,
                             occupancy=self.occupancy, bonds=self._getBonds(fileBonds, guessBonds))
            if np.sum(s) == 0 and strict:
                raise NameError('No atoms were selected with atom selection "{}".'.format(sel))
        else:
            s = sel

        if indexes and s.dtype == bool:
            return np.array(np.where(s)[0], dtype=np.int32)
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
        >>> mol=tryp.copy()
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
        bonds = np.array(self.bonds,
                         dtype=np.int32)  # Have to store in temp because bonds is uint and can't accept -1 values
        bonds[:, 0] = map[self.bonds[:, 0]]
        bonds[:, 1] = map[self.bonds[:, 1]]
        remA = bonds[:, 0] == -1
        remB = bonds[:, 1] == -1
        stays = np.invert(remA | remB)
        # Delete bonds between non-existant atoms
        self.bonds = bonds[stays, :]

    def _guessBonds(self):
        """ Tries to guess the bonds in the Molecule

        Can fail badly when non-bonded atoms are very close together. Use with extreme caution.
        """
        framecoords = self.coords[:, :, self.frame].copy()
        return guessbonds(framecoords, self.element, self.name, self.resname, self.resid, self.chain, self.segid,
                                self.insertion, self.altloc)

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
        >>> mol=tryp.copy()
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
        >>> mol=tryp.copy()
        >>> mol.rotate([0, 1, 0], 1.57)
        """
        M = rotationMatrix(axis, angle)
        s = self.atomselect(sel, indexes=True)
        for a in s:
            self.coords[a, :, self.frame] = np.dot(M, self.coords[a, :, self.frame])

    def rotateBy(self, M, center=(0, 0, 0), sel='all'):
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

    def center(self, loc=(0, 0, 0), sel='all'):
        """ Moves the geometric center of the Molecule to a given location

        Parameters
        ----------
        loc : list, optional
            The location to which to move the geometric center
        sel : str
            An Atomselection string of the atoms whose geometric center we want to center on the `loc` position

        Examples
        --------
        >>> mol=tryp.copy()
        >>> mol.center()
        >>> mol.center([10, 10, 10], 'name CA')
        """
        coords = self.get('coords', sel=sel)
        com = np.mean(coords, 0)
        self.moveBy(-com)
        self.moveBy(loc)

    def read(self, filename, type=None, skip=None, frames=None, append=False):
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
        """
        from htmd.simlist import Sim
        if isinstance(filename, list) or isinstance(filename, np.ndarray):
            for f in filename:
                if len(f) != 4 and not os.path.exists(f):
                    raise FileNotFoundError('File {} was not found.'.format(f))
            firstfile = filename[0]
        else:
            if not isinstance(filename, Sim) and len(filename) != 4 and not os.path.exists(filename):
                raise FileNotFoundError('File {} was not found.'.format(filename))
            firstfile = filename

        if isinstance(filename, Sim):
            self._readPDB(filename.molfile)
            self._readTraj(filename.trajectory)
            return

        if type is not None:
            type = type.lower()

        if (type is None and firstfile.endswith(".psf")) or type == "psf":
            topo = PSFread(filename)
            self._readTopology(topo, filename)
        elif (type is None and (
            firstfile.endswith(".prm") or firstfile.endswith(".prmtop"))) or type == "prmtop" or type == "prm":
            topo = PRMTOPread(filename)
            self._readTopology(topo, filename)
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
            topo, coords = XYZread(filename)
            self._readTopology(topo, filename)
            self.coords = np.atleast_3d(np.array(coords, dtype=self._dtypes['coords']))
        elif (type is None and firstfile.endswith(".gjf")) or type == "gjf":
            topo, coords = GJFread(filename)
            self._readTopology(topo, filename)
            self.coords = np.atleast_3d(np.array(coords, dtype=self._dtypes['coords']))
        elif (type is None and firstfile.endswith(".mae")) or type == "mae":
            topo, coords = MAEread(filename)
            self._readTopology(topo, filename)
            self.coords = np.atleast_3d(np.array(coords, dtype=self._dtypes['coords']))
        elif (type is None and firstfile.endswith(".mol2")) or type == "mol2":
            topo, coords = MOL2read(filename)
            self._readTopology(topo, filename)
            self.coords = np.atleast_3d(np.array(coords, dtype=self._dtypes['coords']))
        elif (type is None and firstfile.endswith(".crd")) or type == "crd":
            self.coords = np.atleast_3d(np.array(CRDread(filename), dtype=np.float32))
        else:
            try:
                self._readTraj(filename, skip=skip, frames=frames, append=append, mdtraj=True)
            except:
                raise ValueError("Unknown file type")

    def _readPDB(self, filename, mode='pdb'):
        if os.path.isfile(filename):
            mol = PDBParser(filename, mode)
        elif len(filename) == 4:
            # Try loading it from the pdb data directory
            localpdb = os.path.join(htmd.home(dataDir="pdb"), filename.lower() + ".pdb")
            if os.path.isfile(localpdb):
                logger.info("Using local copy for {:s}: {:s}".format(filename, localpdb))
                mol = PDBParser(localpdb, mode)
            else:
                # or the PDB website
                logger.info("Attempting PDB query for {:s}".format(filename))
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
            self.__dict__[k] = numpy.asarray(mol.__dict__[k], dtype=self._dtypes[k])
            # Pad any short list
            if k is not "coords":
                if len(self.__dict__[k]) != natoms:
                    self.__dict__[k] = numpy.zeros(natoms, dtype=self.__dict__[k].dtype)

        self.coords = np.atleast_3d(np.array(self.coords, dtype=np.float32))
        self.bonds = np.array(np.vstack((self.bonds, mol.bonds)), dtype=np.uint32)
        self.ssbonds = np.array(mol.ssbonds, dtype=np.uint32)
        self.box = np.array(mol.box)

        if self.masses is None or len(self.masses) == 0:
            self.masses = numpy.zeros(natoms, dtype=numpy.float32)
        if self.charge is None or len(self.charge) == 0:
            self.charge = mol.charge.copy()
        if self.charge is None or len(self.charge) == 0:
            self.charge = numpy.zeros(natoms, dtype=numpy.float32)

        self.fileloc.append([filename, 0])

    def _readTopology(self, topo, filename, overwrite='all'):
        if isinstance(overwrite, str):
            overwrite = (overwrite, )

        # Checking number of atoms that were read in the topology file
        natoms = []
        for field in topo.atominfo:
            if len(topo.__dict__[field]) != 0:
                natoms.append(len(topo.__dict__[field]))
        natoms = np.unique(natoms)
        if len(natoms) != 1:
            raise TopologyInconsistencyError('Different number of atoms read from file {} for different fields: {}.'
                                             .format(filename, natoms))
        natoms = natoms[0]

        if self.numAtoms == 0:
            self.empty(natoms)

        for field in topo.__dict__:
            newfielddata = np.array(topo.__dict__[field], dtype=self._dtypes[field])
            if len(newfielddata) == 0:
                continue
            if overwrite[0] == 'all' or field in overwrite or len(self.__dict__[field]) == 0:
                self.__dict__[field] = newfielddata
            else:
                if np.shape(self.__dict__[field]) != np.shape(newfielddata):
                    raise TopologyInconsistencyError(
                        'Different number of atoms read from topology file for field {}'.format(field))
                if not np.array_equal(self.__dict__[field], newfielddata):
                    raise TopologyInconsistencyError(
                        'Different atom information read from topology file for field {}'.format(field))

        fnamestr = os.path.splitext(os.path.basename(filename))[0]
        self.viewname = fnamestr
        self.fileloc = [[fnamestr, 0]]

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

        # from IPython.core.debugger import Tracer
        # Tracer()()
        if frames is not None:
            if not isinstance(frames, list) and not isinstance(frames, np.ndarray):
                frames = [frames]
            if len(filename) != len(frames):
                raise NameError(
                    'Number of trajectories (' + str(len(filename)) + ') does not match number of frames (' + str(
                        len(frames)) + ') given as arguments')

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
                raise ValueError('Trajectory # of atoms ' + str(
                    np.size(self.coords, 0)) + ' mismatch with # of already loaded atoms ' + str(self.numAtoms))

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
            self.coords = np.array(self.coords[:, :, ::skip])  # np.array is required to make copy and thus free memory!
            self.box = np.array(self.box[:, ::skip])
            self.fileloc = self.fileloc[::skip]

        self.coords = np.atleast_3d(self.coords)
        self.step = traj.step
        self.time = traj.time
        if len(traj.time) < 2:
            # logger.info('Trajectory has broken framestep. Cannot read correctly, setting to 0.1ns.')
            self.fstep = 0.1
        else:
            self.fstep = (traj.time[1] - traj.time[0]) / 1E6  # convert femtoseconds to nanoseconds

    def view(self, sel=None, style=None, color=None, guessBonds=True, viewer=None, hold=False, name=None,
             viewerhandle=None):
        """ Visualizes the molecule in a molecular viewer

        Parameters
        ----------
        sel : str
            Atomselection string for a representation.
        style : str
            Representation style.
        color : str
            Coloring mode or color ID.
        guessBonds : bool
            Allow VMD to guess bonds for the molecule
        viewer : str ('vmd', 'webgl')
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

        oldbonds = None
        if guessBonds:
            oldbonds = self.bonds
            self.bonds = self._getBonds()

        # Write out PDB and XTC files
        psf = tempname(suffix=".psf")
        self.write(psf)

        if guessBonds:
            self.bonds = oldbonds

        xtc = tempname(suffix=".xtc")
        self.write(xtc)

        # Call the specified backend
        if viewer is None:
            from htmd.config import _config
            viewer = _config['viewer']
        if viewer.lower() == 'notebook':
            return self._viewMDTraj(psf, xtc)
        elif viewer.lower() == 'vmd':
            self._viewVMD(psf, xtc, viewerhandle, name, guessBonds)
        elif viewer.lower() == 'ngl' or viewer.lower() == 'webgl':
            return self._viewNGL(psf, self.coords, guessBonds)
        else:
            raise ValueError('Unknown viewer.')

        # Remove temporary files
        os.remove(xtc)
        os.remove(psf)

    def _viewVMD(self, psf, xtc, vhandle, name, guessbonds):
        if name is None:
            name = self.viewname
        if vhandle is None:
            vhandle = getCurrentViewer()

        if guessbonds:
            vhandle.send("mol new " + psf)
        else:
            vhandle.send("mol new " + psf + " autobonds off")
        vhandle.send('animate delete all')
        vhandle.send('mol addfile ' + xtc + ' type xtc waitfor all')

        if name is not None:
            vhandle.send('mol rename top "' + name + '"')
        else:
            vhandle.send('mol rename top "Mol [molinfo top]: psf+xtc"')

        self._tempreps.append(self.reps)
        self._tempreps._repsVMD(vhandle)
        self._tempreps.remove()

    def _viewMDTraj(self, psf, xtc):
        from mdtraj.html import TrajectoryView, TrajectorySliderView, enable_notebook
        import mdtraj
        enable_notebook()

        t = mdtraj.load(xtc, top=psf)
        if self.numFrames > 1:
            widget = TrajectorySliderView(t)
        else:
            widget = TrajectoryView(t)
        return widget

    def _viewNGL(self, psf, coords, guessb):
        from nglview import Trajectory
        import nglview

        class TrajectoryStreamer(Trajectory):
            def __init__(self, coords):
                self.coords = coords

            def get_coordinates_list(self, index):
                return self.coords[:, :, index].flatten().tolist()

            def get_frame_count(self):
                return np.size(self.coords, 2)

        struc = nglview.FileStructure(psf)
        struc.params['dontAutoBond'] = not guessb

        traj = TrajectoryStreamer(coords)
        w = nglview.NGLWidget(struc, traj)
        #else:
        #    w = nglview.NGLWidget(struc)

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
        >>> mol=tryp.copy()
        >>> mol.mutateResidue('resid 158', 'ARG')
        """
        s = self.atomselect(sel, strict=True)
        # Changed the selection from "and sidechain" to "not backbone" to remove atoms like phosphates which are bonded
        # but not part of the sidechain. Changed again the selection to "name C CA N O" because "backbone" works for
        # both protein and nucleic acid backbones and it confuses phosphates of modified residues for nucleic backbones.
        remidx = self.atomselect(sel + ' and not name C CA N O', indexes=True)
        self.remove(remidx, _logger=False)
        s = np.delete(s, remidx)
        self.set('resname', newres, sel=s)

    def wrap(self, wrapsel=None, fileBonds=True, guessBonds=True):
        """ Wraps coordinates of the molecule into the simulation box

        Parameters
        ----------
        wrapsel : str
            Selection of atoms on which to center the wrapping box

        Examples
        --------
        >>> mol=tryp.copy()
        >>> mol.wrap()
        >>> mol.wrap('protein')
        """
        # TODO: selection is not used. WHY?
        if wrapsel is not None:
            centersel = self.atomselect(wrapsel, indexes=True)
        else:
            centersel = None
        self.coords = wrap(self.coords, self._getBonds(fileBonds, guessBonds), self.box, centersel=centersel)

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
                # traj.xyz = np.swapaxes(np.swapaxes(self.coords, 1, 2), 0, 1) / 10
                # traj.time = self.time
                # traj.unitcell_lengths = self.box.T / 10
                traj.save(filename)
            except:
                raise ValueError("Unknown file type")

    def _writeXYZ(self, filename, sel="all"):
        src = self
        if sel is not None:
            src = sel.copy()
            src.filter(sel, _logger=False)
        fh = open(filename, "w")
        natoms = len(src.record)
        print("%d\n" % (natoms), file=fh)
        for i in range(natoms):
            e = src.element[i].strip()
            if (not len(e)):
                e = re.sub("[1234567890]*", "", src.name[i])
            print("%s   %f   %f    %f" % (e, src.coords[i, 0, 0], src.coords[i, 1, 0], src.coords[i, 2, 0]), file=fh)
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

        pdb.coords = np.atleast_3d(pdb.coords[:, :, self.frame])  # Writing out only current frame

        pdb.bonds = src.bonds
        pdb.ssbonds = src.ssbonds  # TODO: Is there such a thing in pdb format?
        if pdb.box is not None:
            pdb.box = np.atleast_2d(np.atleast_2d(self.box)[:, self.frame])

        pdb.serial = np.arange(1, np.size(pdb.coords, 0) + 1)
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
        >>> newmol = Molecule()
        >>> newmol.empty(100)
        """
        self.record = np.array(['ATOM'] * numAtoms, dtype=self._dtypes['record'])
        self.chain = np.array([''] * numAtoms, dtype=self._dtypes['chain'])
        self.segid = np.array([''] * numAtoms, dtype=self._dtypes['segid'])
        self.occupancy = np.array([0] * numAtoms, dtype=self._dtypes['occupancy'])
        self.beta = np.array([0] * numAtoms, dtype=self._dtypes['beta'])
        self.insertion = np.array([''] * numAtoms, dtype=self._dtypes['insertion'])
        self.element = np.array([''] * numAtoms, dtype=self._dtypes['element'])
        self.altloc = np.array([''] * numAtoms, dtype=self._dtypes['altloc'])
        self.name = np.array([''] * numAtoms, dtype=self._dtypes['name'])
        self.resname = np.array([''] * numAtoms, dtype=self._dtypes['resname'])
        self.resid = np.array([0] * numAtoms, dtype=self._dtypes['resid'])
        self.coords = np.zeros((numAtoms, 3, 1), dtype=self._dtypes['coords'])
        self.charge = np.array([0] * numAtoms, dtype=self._dtypes['charge'])
        self.serial = np.arange(1, numAtoms + 1)

        self.masses = np.array([0] * numAtoms, dtype=self._dtypes['masses'])
        self.box = np.zeros(self._dims['box'], dtype=np.float32)


    def sequence(self, oneletter=True):
        """ Return the AA sequence of the Molecule.

        Parameters
        ----------
        oneletter : bool
            Whether to return one-letter or three-letter AA codes. There should be only one atom per residue.

        Returns
        -------
        sequence : str
            The primary sequence as a dictionary segid - string (if oneletter is True) or segid - list of
            strings (otherwise).

        Examples
        --------
        >>> mol=tryp.copy()
        >>> mol.sequence()
        {'0': 'IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN'}
        >>> sh2 = Molecule("1LKK")
        >>> pYseq = sh2.sequence(oneletter=False)
        >>> pYseq['1']
        ['PTR', 'GLU', 'GLU', 'ILE']
        >>> pYseq = sh2.sequence(oneletter=True)
        >>> pYseq['1']
        '?EEI'

        """
        residueTable = {'ARG': 'R', 'AR0': 'R',
                        'HIS': 'H', 'HID': 'H', 'HIE': 'H', 'HSE': 'H', 'HSD': 'H',
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

        segSequences = {}

        # Iterate over segments
        for seg in segs:
            segSequences[seg] = []
            segatoms = self.atomselect('protein and segid {}'.format(seg))
            resnames = self.resname[segatoms]
            incremseg = increm[segatoms]
            for i in np.unique(incremseg):  # Iterate over residues
                resname = np.unique(resnames[incremseg == i])
                if len(resname) != 1:
                    raise AssertionError('Unexpected non-uniqueness of chain, resid, insertion in the sequence.')
                resname = resname[0]
                if oneletter:
                    rescode = residueTable.get(resname, "?")
                    if rescode == "?":
                        logger.warning("Cannot provide one-letter code for non-standard residue %s" % resname)
                else:
                    rescode = resname
                segSequences[seg].append(rescode)

        # Join single letters into strings
        if oneletter:
            segSequences = {k: "".join(segSequences[k]) for k in segSequences}

        return segSequences

    def dropFrames(self, keep='all', drop=None):
        """ Removes trajectory frames from the Molecule

        Parameters
        ----------
        keep : int or list of ints
            Index of frame, or list of frame indexes which we want to keep (and drop all others).
        drop : int or list of ints
            Index of frame, or list of frame indexes which we want to drop (and keep all others).
        """
        if keep != 'all' and drop is not None:
            raise RuntimeError('Cannot both drop and keep trajectories. Please use only one of the two arguments.')
        if keep != 'all':
            self.coords = np.array(np.atleast_3d(self.coords[:, :, keep]))  # Copy array. Slices are dangerous with C
            self.box = np.array(np.atleast_2d(self.box[:, keep]))
        if drop is not None:
            self.coords = np.delete(self.coords, drop, axis=2)
            self.box = np.delete(self.box, drop, axis=1)

    @property
    def numFrames(self):
        """ Number of coordinate frames in the molecule
        """
        return np.size(np.atleast_3d(self.coords), 2)

    @property
    def numAtoms(self):
        """ Number of atoms in the molecule
        """
        return len(self.record)

    @property
    def x(self):
        """Get the x coordinates at the current frame"""
        return self.coords[:, 0, self.frame]

    @property
    def y(self):
        """Get the y coordinates at the current frame"""
        return self.coords[:, 1, self.frame]

    @property
    def z(self):
        """Get the z coordinates at the current frame"""
        return self.coords[:, 2, self.frame]

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
        for j in sorted(self.__dict__.keys() - self._pdb_fields):
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

    logger.info(
        'Removed {} residues from appended Molecule due to collisions.'.format(len(np.unique(mol.resid[overlaps]))))
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
    >>> mol = tryp.copy()
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
        styletrans = {'newcartoon': 'cartoon', 'licorice': 'hyperball', 'lines': 'line', 'vdw': 'spacefill',
                      'cpk': 'ball+stick'}
        colortrans = {'name': 'element', 'index': 'atomindex', 'chain': 'chainindex', 'secondary structure': 'sstruc',
                      'colorid': 'color'}
        hexcolors = {0: '#0000ff', 1: '#ff0000', 2: '#333333', 3: '#ff6600', 4: '#ffff00', 5: '#4c4d00', 6: '#b2b2cc',
                     7: '#33cc33', 8: '#ffffff', 9: '#ff3399', 10: '#33ccff'}
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
    # Unfotunately, tests affect each other because only a shallow copy is done before each test, so
    # I do a 'copy' before each.
    import doctest
    from htmd.home import home

    m = Molecule(path.join(home(), 'data', 'building-protein-membrane', '3PTB_clean.pdb'))
    doctest.testmod(extraglobs={'tryp': m.copy()})

    # Oddly, if these are moved before doctests, 1. extraglobs don't work; and 2. test failures are not printed. May
    # have to do with the vmd console?
    a = m.get('resid', sel='resname TRP')
    a = m.get('coords')
    print(a.ndim)
    m.write('/tmp/test.pdb')
    m.write('/tmp/test.coor')
    m.write('/tmp/test.xtc')
    m.moveBy([1, 1, 1])
    m.rotate([1, 0, 0], pi / 2)
    m.align('name CA')
    m = Molecule('2OV5')
    m.filter('protein or water')

    # test rotate
