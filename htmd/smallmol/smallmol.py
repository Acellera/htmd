import os
from copy import deepcopy
import multiprocessing
import math
import numpy as np
from rdkit import Chem
from htmd.smallmol.util import get_rotationMatrix, rotate, openbabelConvert, _depictMol, depictMultipleMols
from htmd.smallmol.chemlab.periodictable import _hybridizations_IdxToType, _bondtypes_IdxToType,  _hybridizations_StringToType, _bondtypes_StringToType
import logging
logger = logging.getLogger(__name__)

print('SmallMol module is in beta version')


class SmallMol:
    """
    Class to manipulate small molecule structures

    SmallMol extract all the molecule and atoms data from an rdkit.Chem.Molecule into the fields of this object.
    The fields can be access directly as attributes.
    Atom fields: 'idx', 'atomname', 'charge','formalcharge', 'element',  'chiral', 'hybridization',
                 'neighbors', 'bondtypes', 'coords', '_chiraltags'
    Molecule fields: 'ligname', 'totalcharge', '_mol'

    mol, ignore_errors=False, force_reading=False, fixHs=True, removeHs=False
    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol  or filename or smile or htmd.smallmol.smallmol.SmallMol
        (i) Rdkit molecule or (ii) Location of molecule file (".pdb"/".mol2") or (iii) a smile string or iv) another
        SmallMol object
    ignore_errors: bool
        If True, errors will not be raised.
    force_reading: bool
        If True, and the mol provided is not accepted, the molecule will be initially converted into sdf
    fixHs: bool
        If True, the missing hydrogens are assigned, the others are correctly assinged into the graph of the molecule
    removeHs: bool
        If True, remove the hydrogens

    Examples
    --------
    >>> sm = SmallMol( 'CCO' )  # doctest: +SKIP
    >>> sm = SmallMol( 'htmd/data/test-smallmol/ligand.pdb', fixHs=False, removeHs=True )
    >>> sm = SmallMol( 'htmd/data/test-smallmol/benzamidine.mol2')
    >>> print(mol)                                     # doctest: +ELLIPSIS
    SmallMol with 18 atoms and
    Atom field - atomname shape: (18,)
    Atom field - bondtypes shape: (18,)
    ...

    .. rubric:: Methods
    .. autoautosummary:: htmd.smallmol.smallmol.SmallMol
       :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.smallmol.smallmol.SmallMol
       :attributes:

    Attributes
    ----------

    listProps : list
        List of the molecule property
    listPropsAtoms: list
        List of the atom property
    numConformers: int
        Number of confomrmers
    numAtoms: int
        Number of atoms

    """

    _atom_fields = ['idx', 'atomname', 'charge','formalcharge', 'element',  'chiral', 'hybridization',
                    'neighbors', 'bondtypes', 'coords', '_chiraltags']

    _mol_fields = ['ligname', 'totalcharge', '_mol']

    ### Field types
    ## Development Note:  'expliciths' removed.
    _atom_dtypes = {'idx': np.int,
                    'atomname': object,
                    'charge': np.float32,
                    'formalcharge': np.int,
                    'element': object,
                    'chiral': object,
                    'hybridization': int,
                    'neighbors': object,
                    'bondtypes': int,
                    'coords': np.float32,
                    '_chiraltags': object
                    }

    _mol_dtypes = {'_mol': object,
                    'ligname': object,
                   'totalcharge': np.int
                   }

    ### Field shapes
    ## Development Note:  'expliciths' removed.
    _atom_dims = {'idx': (0,),
                  'atomname': (0,),
                  'charge': (0,),
                  'formalcharge': (0,),
                  'element': (0,),
                  'chiral': (0,),
                  'neighbors': (0,),
                  'bondtypes': (0,),
                  'hybridization': (0,),
                  'coords': (0, 3, 1),
                  '_chiraltags':(0,)
                    }

    _mol_dims = {'_mol': (0,),
                 'ligname': (0,),
                 'totalcharge': (0,)
                   }

    def __init__(self, mol, ignore_errors=False, force_reading=False, fixHs=True, removeHs=False):


        self.mol_fields = self._mol_fields.copy()
        self.atom_fields = self._atom_fields.copy()

        # init the __dict__ with the fields
        for field in self._atom_dtypes:
            self.__dict__[field] = np.empty( self._atom_dims[field], dtype=self._atom_dtypes[field] )
        for field in self._mol_dtypes:
            self.__dict__[field] = None

        # load the input and store the rdkit Mol obj
        # self._mol is used for debugging. Later on will be removed
        self._mol, smallMol = self._initializeMolObj(mol, force_reading)
        if not ignore_errors and self._mol == None:
            if not force_reading:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj. "
                                 "Try by setting the force_reading option as True.".format(mol))
            else:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj.".format(mol))

        if removeHs:
            self._mol = Chem.RemoveHs(self._mol)

        if fixHs:
            self._mol = Chem.AddHs(self._mol, addCoords=True)

        # fill the molecule and atom properties
        self._readMol(smallMol)

    @property
    def listProps(self):
        """
        The list of properties of the molecule

        Returns
        -------
        fields: list
            The list of properties
        """
        _fields = list(self.mol_fields)
        return _fields

    @property
    def listPropsAtoms(self):
        """
        The list of properties of the atoms

        Returns
        -------
        fields: list
            The list of properties
        """
        _fields = list(self.atom_fields)
        return _fields

    @property
    def numAtoms(self):
        """
        Returns the number of atoms in the molecule

        Returns
        -------
        natoms: int
            The number of atoms in the molecule

        """
        # getting the number of atoms from the shape of element array
        natoms = self.element.shape[0]
        return natoms


    @property
    def numConformers(self):
        """
        Returns the number of conformers in the rdkit molecule object

        Returns
        -------
        nconfs: int
            The number of conformers in the molecule
        """
        # getting the number of conformers from the shape of coords array
        _nConformers = self.coords.shape[2]
        return _nConformers

    def _initializeMolObj(self, mol, force_reading):
        """
        Read the input and tries to convert it into a rdkit.Chem.rdchem.Mol obj

        Parameters
        ----------
        mol: str or rdkit.Chem.rdchem.Mol or htmd.smallmol.smallmol.SmallMol
            i) rdkit.Chem.rdchem.Mol ii) The path to the pdb/mol2 to load iii) The smile string iv) SmallMol object
        force_reading: bool
           If the mol provided is not accepted, the molecule will be initially converted into sdf

        Returns
        -------
        _mol: rdkit.Chem.Molecule object
            The rdkit molecule
        smallMol: htmd.smallmol.smallmol.SmallMol
            The smallMol object if SmallMol was passed
        """
        _mol = None
        smallMol = None
        if isinstance(mol, SmallMol):
            _mol = mol._mol
            smallMol = mol

        if isinstance(mol, Chem.Mol):
            _mol = mol

        elif isinstance(mol, str):
            if os.path.isfile(mol):
                name_suffix = os.path.splitext(mol)[-1]
                # load mol2 file
                if name_suffix == ".mol2":
                    _mol = Chem.MolFromMol2File(mol, removeHs=False)
                # load pdb file
                elif name_suffix == ".pdb":
                    _mol = Chem.MolFromPDBFile(mol, removeHs=False)
                # if the file failed to be loaded and 'force_reading' = True, file convert to sdf and than loaded
                if _mol == None and force_reading:
                    logger.warning('Reading {} with force_reading procedure'.format(mol))
                    sdf = openbabelConvert(mol, name_suffix, 'sdf')
                    _mol = Chem.SDMolSupplier(sdf, removeHs=False)[0]

                    os.remove(sdf)

            # assuming is a smile
            # TODO validate it. Implement smarts recognition
            else:
                _mol = Chem.MolFromSmiles(mol)

        return _mol, smallMol

    def _readMol(self, smallMol):
        """
        Set up the parameters of SmallMol object.

        Parameters
        ----------
        smallMol: htmd.smallmol.smallmol.SmallMol or None
             If is not None the parameters are taken from this object
        """
        # the fields
        # Development note explicitHs removed
        idxs = []
        atomnames = []
        charges = []
        formalcharges = []
        elements = []
        chirals =  []
        hybridizations = []
        neighbors = []
        bondtypes = []
        chiraltags = []

        _mol = self._mol

        if smallMol is not None:
            for f, v in smallMol.__dict__.items():
                self.__dict__[f] = v
            return

        # Development note explicitHs removed
        for a in _mol.GetAtoms():
            i = a.GetIdx()
            e = a.GetSymbol()
            chiraltags.append(a.GetChiralTag())
            idxs.append(i)
            atomnames.append('{}{}'.format(e,i))
            formalcharges.append(a.GetFormalCharge())
            elements.append(e)
            neighbors.append([na.GetIdx() for na in a.GetNeighbors()])
            hybridizations.append(int(a.GetHybridization()))
            bondtypes.append([int(b.GetBondType()) for b in a.GetBonds()])
            if a.HasProp('_TriposPartialCharge'):
                charges.append(a.GetPropsAsDict()['_TriposPartialCharge'])
            else:
                charges.append(0.000)
            if a.HasProp('_CIPCode'):
                chirals.append(a.GetPropsAsDict()['_CIPCode'])
            else:
                chirals.append('')

        if _mol.HasProp('_Name'):
            self.__dict__['ligname'] = _mol.GetProp('_Name')
        else:
            self.__dict__['ligname'] = 'UNK'

        for k, v in _mol.GetPropsAsDict().items():
            self.setProp(k,v, overwrite=True)

        # Development note explicitHs removed
        self.__dict__['totalcharge'] = sum(formalcharges)
        self.__dict__['idx'] = np.array(idxs)
        self.__dict__['atomname'] = np.array(atomnames)
        self.__dict__['charge'] = np.array(charges )
        self.__dict__['formalcharge'] = np.array(formalcharges)
        self.__dict__['element'] = np.array(elements)
        self.__dict__['chiral'] = np.array(chirals)
        self.__dict__['hybridization'] = np.array(hybridizations)
        self.__dict__['neighbors'] = np.array(neighbors)
        self.__dict__['bondtypes'] = np.array(bondtypes)
        self.__dict__['_chiraltags'] = np.array(chiraltags)

        if _mol.GetNumConformers() != 0:
            coords = _mol.GetConformer(0).GetPositions()
            self.__dict__['coords'] =  coords[:,:,np.newaxis]
        else:
            self.__dict__['coords'] = np.zeros((self.numAtoms,3,1))

    def copy(self):
        """
        Create a copy of the molecule object

        Returns
        -------
        newsmallmol : :class:`SmallMol`
            A copy of the object
        """
        newsmallmol = deepcopy(self)
        return newsmallmol

    def setProp(self, key, value, overwrite=False):
        """
        Sets the property value based on the key

        Parameters
        ----------
        key: str
            The name of the property to set
        value: int, str, float, object
            The value of the property
        overwrite: bool
            If True, the property will be overwritten if already exist

        Example
        -------
        >>> sm.setProp('Ki', 25) # doctest: +SKIP
        >>> sm.listProps
        ['ligname', 'totalcharge', '_mol', Ki]
        >>> sm.setProp('ligname', 'myMol', overwrite=True)
        >>> sm.getProp('ligname')
        'myMol'

        """

        if not isinstance(key, str):
             raise ValueError('Wrong type {} for key.  Should be {} '.format(type(key), type('string')))

        _props = self.listProps

        if not overwrite and key in _props:
            raise ValueError('The key passed {} already exists. Set "overwrite" as True to overwrite an existing key'.format(key))

        self.mol_fields.append(key)
        self.__dict__[key] = value

    def setPropsAtoms(self, key, values, aIdxs=None):
        """
        Sets the property for one or more atoms

        Parameters
        ----------
        key: str
            The name of the property to set
        values: list
            The list of the values for the property
        aIdxs: list
            The list of the atom index to which set up the values

        Example
        -------
        >>> sm.setPropsAtoms('myTag', ['a', 'b', 'a', 'c'], aIdxs=[0,2,5,3]) # doctest: +SKIP
        >>> sm.myTag
        array(['a', None, 'b', 'c', None, 'a', None, None, None, None, None, None,
       None, None, None, None, None, None], dtype=object)
        """
        class ValueLengthErrors(Exception):
            pass

        if not isinstance(key, str):
             raise ValueError('Wrong type {} for key.  Should be {} '.format(type(key), type('string')))

        if aIdxs is None and len(values) != self.numAtoms:
            raise ValueLengthErrors('The number of values passed {} do not match '
                                     'the number of atoms {}'.format(len(values), self.numAtoms))

        if aIdxs is not None and len(aIdxs) != len(values):
            raise ValueLengthErrors('The number of values passed {} do not match '
                                     'the number of idxs {}'.format(len(values), len(aIdxs)))

        _props = self.listPropsAtoms

        # creation of new property. Requires the creation of a new array
        if key not in _props:
            tmp_array = np.array([None]*self.numAtoms)
            self.__dict__[key] = tmp_array

        self.atom_fields.append(key)
        self.__dict__[key][aIdxs] = values

    def getProp(self, key):
        """
        Returns the  molecule or the atoms property of the key passed

        Parameters
        ----------
        key: str
            The name of the molecule property to retrieve

        Returns
        -------
        value: str, int, float, object
            The value of the molecule property
        """
        if key not in self.listProps:
            raise KeyError('The property passed {} does not exist'.format(key))

        _value = self.__dict__[key]

        return _value

    def get(self,  returnField, sel, convertType=True, invert=False):
        """
        Returns the property for the atom specified with the selection. The selection is another atom property

        Parameters
        ----------
        returnField: str
            The field of the atom to return
        sel: str
            The selection string. atom field name followed by spaced values for that field
        convertType: bool
            If True, and where possible the returnField is converted in rdkit object
        invert: bool
            If True, the selection is inverted

        Returns
        -------
        values: np.array
            The array of values for the property

        Example
        -------
        >>> sm.get('element', 'idx 0 1 7')
        array(['C', 'C', 'H'], dtype='<U1')
        >>> sm.get('hybridization', 'element N')
        array([rdkit.Chem.rdchem.HybridizationType.SP2, rdkit.Chem.rdchem.HybridizationType.SP2], dtype=object)
        >>> sm.get('hybridization', 'element N', convertType=False)
        array([3, 3])
        >>> sm.get('element', 'hybridization sp2')
        array(['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N'],  dtype='<U1')
        >>> sm.get('element', 'hybridization S')
        array(['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],  dtype='<U1')
        >>> sm.get('element', 'hybridization 1')
        array(['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],  dtype='<U1')
        >>> sm.get('atomobject', 'element N')
        array([<rdkit.Chem.rdchem.Atom object at 0x7faf616dd120>, <rdkit.Chem.rdchem.Atom object at 0x7faf616dd170>], dtype=object)
        """

        # get the field key and the value to grep
        key = sel.split()[0]
        selector = sel.split()[1:]

        if key not in self.listPropsAtoms:
            raise KeyError('The property passed {} does not exist'.format(key))
        if len(selector) == 0:
            raise ValueError('No selection was provided')

        # get the returnField and process exceptional field
        _arrayFrom = self.__dict__[key]
        if returnField == 'atomobject':
            _arrayTo = self.getAtoms()
        elif returnField == 'hybridization' and convertType:
            _arrayTo = np.array([ _hybridizations_IdxToType[v] for v in self.__dict__[returnField] ], dtype=object)
        elif returnField == 'bondtypes' and convertType:
            _arrayTo = np.array([ [  _bondtypes_IdxToType[v] for v in b ] for b in self.__dict__[returnField]], dtype=object)
        else:
            _arrayTo = self.__dict__[returnField]

        # special selector for hybridization: can be idx, or rdkit.Chem.rdchem.HybridizationType
        if key == 'hybridization':
            try:
                selector = [_hybridizations_StringToType[s.upper()] for s in selector]
            except:
                pass

        _dtype = self._atom_dtypes[key]
        if _dtype is not object:
            selector = [ _dtype(s) for s in selector ]
        idxs = np.concatenate([ np.where( _arrayFrom == s )[0] for s in selector ])
        if invert:
            idxs = np.array([ i for i in self.idx if i not in idxs ])
        idxs = np.sort(idxs)

        return _arrayTo[idxs]

    def isChiral(self, returnDetails=False):
        """
        Returns True if the molecule has at least one chiral atom. If returnDetails is set as True,
        a list of tuples with the atom idx and chiral type is returned.

        Parameters
        ----------
        returnDetails: bool
            If True, returns the chiral atoms and their chiral types

        Returns
        -------
        ischiral: bool
            True if the atom has at least a chiral atom
        details: list
            A list of tuple with the chiral atoms and their types

        Example
        -------
        >>> sm.isChiral()
        True
        >>> sm.isChiral(returnDetails=True)
        (True, [('C2', 'R')])
        """

        _chirals = self.chiral

        idxs = np.where(_chirals != '')[0]
        if len(idxs) == 0:
            return False
        if returnDetails:
            idxs =  idxs.astype(str)
            idxs_str = ' '.join(idxs)
            atomnames = self.get( 'atomname', 'idx {}'.format(idxs_str))
            chirals = self.get( 'chiral', 'idx {}'.format(idxs_str))

            _details = [(a, c) for a, c in zip(atomnames, chirals)]

            return True, _details

        return True

    def foundBondBetween(self, sel1, sel2, bondtype=None):
        """
        Returns True if at least a bond is found between the two selections. Is possible to check for specific bond type.
        A tuple is return like (True, [ [(0,1), rdkit.Chem.rdchem.BondType.SINGLE]] ])

        Parameters
        ----------
        returnDetails: bool
            If True, returns the chiral atoms and their chiral types

        Returns
        -------
        sel1: str
            The selection for the first set of atoms
        sel2: str
            The selection for the second set of atoms
        bondtype: str or int
            The bondtype as index or string

        Returns
        -------
        isbond: bool
            True if a bond was found
        details: list
            A list of lists with the index of atoms in the bond and its type
        """

        if  isinstance(bondtype, str):
            _btype = _bondtypes_StringToType[bondtype]
        else:
            _btype = bondtype

        atomIdx_sel1 = self.get('idx', sel1)
        neighbors_sel1 = self.get('neighbors', sel1)
        bondtypes_sel1 = self.get('bondtypes', sel1)

        atomIdx_sel2 = self.get('idx', sel2)

        founds = []
        for aIdx, neighbors_set, btype_set in zip(atomIdx_sel1, neighbors_sel1, bondtypes_sel1):
            for neighbor, btype in zip(neighbors_set, btype_set):
                if neighbor in atomIdx_sel2:
                    btype_str = _bondtypes_IdxToType[btype]
                    if _btype is not None and _btype == btype:
                        founds.append([(aIdx, neighbor), btype_str])
                    elif bondtype == None:
                        founds.append( [(aIdx, neighbor), btype_str]  )
        if len(founds) != 0:
            return True, founds

        return False

    def getAtoms(self):
        """
        Retuns an array with the rdkit.Chem.rdchem.Atom present in the molecule
        """

        _mol = self.toRdkitMol()
        return np.array([a for a in _mol.GetAtoms()])

    def getCoords(self, confId=0):
        """
        Returns molecule coordinates of the conformer id passed.

        Parameters
        ----------
        id: int
            The id of the conformer

        Returns
        -------
        coords: numpy.array
            A numpy array of the coords for the conformer

        """
        coords = self.coords[:,:,confId]

        return coords

    def _getAtomTypes(self):
        """
        Returns ndarray of shape (n_atoms x n_properties) molecule atom types,
        according to the following definitions and order:
            0. Hydrophibic
            1. Aromatic
            2. Acceptor
            3. Donor
            4. - Ionizable
            5. + Ionizable
            6. Metal (empty)
            7. Occupancy (No hydrogens)
        """
        from htmd.smallmol.chemlab.periodictable import atom_mapping
        from htmd.smallmol.util import factory
        n_atoms = self.numAtoms
        _mol = self.toRdkitMol()

        feats = factory.GetFeaturesForMol(_mol)
        properties = np.zeros((n_atoms, 8), dtype=bool)

        for feat in feats:
            fam = feat.GetFamily()
            if fam not in atom_mapping:  # Non relevant property
                continue
            properties[feat.GetAtomIds(), atom_mapping[fam]] = 1

        # Occupancy, ignoring hydrogens.
        properties[:, 7] = self.element != 'H'
        return properties

    def _getChannelRadii(self):
        """
        Multiplies atom types by each atom vdW radius.
        """
        from htmd.molecule.vdw import radiidict
        radii = np.vectorize(radiidict.__getitem__)(self.element) * self._getAtomTypes().T
        return radii.T.copy()

    def getCenter(self, confId=0, coords=None):
        """
        Returns geometrical center of molecule conformation

        Parameters
        ----------
        confId: int
            The conformer
        coords: np.array
            The coords for which you want the center
        """
        if coords is None:
            coords = self.getCoords(confId)
        return coords.mean(axis=0).astype(np.float32)

    def generateConformers(self, num_confs=400,  optimizemode='mmff', align=True, append=True):
        """
        Generates ligand conformers

        Parameters
        ----------
        num_confs: int
           Number of conformers to generate.
        optimizemode: str
            The optimizemode to use. Can be  'uff', 'mmff'
        align: bool
            If True, the conformer are aligned to the first one
        append: bool
            If False, the current conformers are deleted

        """
        from rdkit.Chem.AllChem import UFFOptimizeMolecule, MMFFOptimizeMolecule, EmbedMultipleConfs
        from rdkit.Chem.rdMolAlign import AlignMolConformers

        if not append:
            self.removeConformers()

        # get the rdkit mol and copied it.
        _mol = self.toRdkitMol()
        mol = deepcopy(_mol)
        # hydrogens are added for safety
        mol = Chem.AddHs(mol)

        # if coords exist these are stored as a conformer object insider the rdkit molecule one
        if self.numConformers != 0:
            conf = self.getConformers([0])[0]
            mol.AddConformer(conf, 0)

        # generating conformations
        ids = EmbedMultipleConfs(mol, clearConfs=False,numConfs=num_confs, pruneRmsThresh=1., maxAttempts=10000)
        if optimizemode not in ['uff', 'mmff']:
            raise ValueError('Unknown optimizemode. Should be  "uff", "mmff"')
        # optimizing conformations depends on the optimizemode passed
        for id in ids:
            if optimizemode == 'mmff':
                MMFFOptimizeMolecule(mol, confId=id)
            elif optimizemode == 'uff':
                UFFOptimizeMolecule(mol, confId=id)

        if align:
            AlignMolConformers(mol)

        for i in ids:
            conf = mol.GetConformer(i)
            coords = conf.GetPositions()
            coords = coords[:,:,np.newaxis]
            if self.coords.shape[0] == 0:

                self.coords = coords
            else:
                self.coords = np.concatenate((self.coords, coords), axis=2).astype(np.float32)

    def getVoxels(self, center=None, size=24, resolution=1., rotation=None,
                   displacement=None, dtype=np.float32, confId=0):
        """
        Computes molecule voxelization.

        Parameters
        ----------
        center: array-like
            Geometrical coordinates where descriptors will be computed.
        size: int
            Size of resulting descriptor array.
        resolution: float
            Grid resolution of resulting array.

        rotation : array-like of shape (3,)
            Prior to voxelization rotates the molecule around its center give the
            rotation angles in radians.
        displacement: array-like of shape (3,)
            Prior to voxelization displaces the molecule by provided (X, Y, Z) distance before
            returning the voxelized representation.
        dtype : numpy datatype
            returns array of the specified type.
        Returns
        -------
        voxels: array-like
            Computed descriptors.
        """
        from htmd.molecule.voxeldescriptors import _getOccupancyC, _getGridCenters
        from htmd.smallmol.util import array_cache

        coords = self.getCoords(confId)
        lig_center = self.getCenter(coords=coords)

        if center is None:
            center = lig_center

        if rotation is not None:
            rotation = list(rotation)
            matx = get_rotationMatrix([1, 0, 0], rotation[0])
            maty = get_rotationMatrix([0, 1, 0], rotation[1])
            matz = get_rotationMatrix([0, 0, 1], rotation[2])

            coords = rotate(coords, matx, center=lig_center)
            coords = rotate(coords, maty, center=lig_center)
            coords = rotate(coords, matz, center=lig_center)

        if displacement is not None:
            coords += np.asarray(displacement)

        multisigmas = self._getChannelRadii()
        #if (size, resolution) not in SmallMol.array_cache:
        if (size, resolution) not in array_cache:
            N = [size, size, size]
            bbm = (np.zeros(3) - float(size * resolution / 2))
            centers = _getGridCenters(bbm, N, resolution)

            # Cache the array
            #SmallMol.array_cache[(size, resolution)] = centers.reshape(size**3, 3)
            array_cache[(size, resolution)] = centers.reshape(size ** 3, 3)
            centers2D = centers + center
        else:
            #centers2D = SmallMol.array_cache[(size, resolution)] + center
            centers2D = array_cache[(size, resolution)] + center

        voxels = _getOccupancyC(coords.astype(np.float32), centers2D,
                                multisigmas).reshape(size, size, size, 8).astype(dtype)
        return voxels

    def _getConformer(self, id):
        """
        Returns the rdkit.Chem.rdchem.Conformer based on the id.

        Parameters
        ----------
        id: int
            The id of the conformer to return

        Returns
        -------
        conf: rdkit.Chem.rdchem.Conformer
            The rdkit conformer

        """

        from rdkit.Chem import Conformer
        from rdkit.Geometry.rdGeometry import Point3D
        conf = Conformer()

        atoms_coords = self.getCoords(id)
        for n, coord in enumerate(atoms_coords):
            p = Point3D(*coord.tolist())
            conf.SetAtomPosition(n, p)

        return conf

    def getConformers(self, ids=None):
        """
        Returns the conformer of the molecule depends on the id passed. If None, all the conformers are
        returned

        Parameters
        ----------
        ids: list
            The list of ids for the molecule conformers to return. If None all the conformers are returned

        Returns
        -------
        _conformer: list of rdkit.Chem.rdchem.Conformer
            The conformers
        """
        _nConformers = self.coords.shape[-1]

        if ids == None:
            return [ self._getConformer(i) for i in range(_nConformers) ]

        if  max(ids) >= _nConformers:
            raise IndexError("The ids list contains conformers ids {} that do not exist. Available conformers: {}".format(ids, _nConformers))

        return [self._getConformer(id) for id in ids]

    def writeConformers(self, savefolder='conformations', savename="molConf", filetype="sdf", savefolder_exist_ok=False,
                         merge=False, ids=None):
        """
        Writes the conformers in one or multiple files in 'pdb' or 'sdf' file formats.

        Paramters
        ---------
        savefolder: str
            The name of the folder where to write the files
        savename: str
            The basename of the conformer file
        filetype: str ('sdf', 'pdb')
            The filetype of the output
        savefolder_exist_ok: bool
            Set as True to overwrite the output folder
        merge: bool
            Set as True to save in a unique file
        ids: list
            A list of the conformer ids to save. If None, all are written

        Example
        -------
        >>> sm.generateConformers(num_confs=3)
        >>> sm.writeConformers(savefolder='Confs')
        >>> glob('Confs/*')
        ['Confs/molConf_0.sdf', 'Confs/molConf_1.sdf', 'Confs/molConf_2.sdf', 'Confs/molConf_3.sdf']
        >>> sm.writeConformers(savefolder='Confs', savefolder_exist_ok=True, merge=True)
        [ 'Confs/molConf_merge.sdf', ...]

        """

        os.makedirs(savefolder, exist_ok=savefolder_exist_ok)

        if ids is None:
            ids = range(self.numConformers)
        elif not isinstance(ids, list):
            raise ValueError("The ids argument should be a list of conformer ids")

        _mol = self.toRdkitMol(includeConformer=True)

        # Init the Writer depends on the filetype passed
        if filetype == 'pdb':
            chemwrite = Chem.PDBWriter
        elif filetype == "sdf":
            chemwrite = Chem.SDWriter
        else:
            raise ValueError("Unknown file format. Cannot save to format '{}'".format(filetype))

        # If merge is set as True a unique file is generated
        if merge:
            fname = os.path.join(savefolder, '{}_merge.{}'.format(savename, filetype))
            writer = chemwrite(fname)

        for id in ids:
            # If merge is set as False a file is created for each conformer
            if not merge:

                fname = os.path.join(savefolder, '{}_{}.{}'.format(savename, id, filetype))
                writer = chemwrite(fname)
            writer.write(_mol, confId=id)

    def removeConformers(self, ids=None):
        """
        Deletes the conformers based on ids

        Parameters
        ----------
        ids: list
            The list of conformer id to delete. If None, all are removed except the first one
        """

        _nConformers = self.numConformers

        if ids is None:
            ids = np.arange(_nConformers)
        elif not isinstance(ids, list):
            raise TypeError("The ids argument should be list of conformer ids")

        self.coords = np.delete(self.coords, ids, axis=2)

    def invertChirality(self, atom):
        """
        Inverts the chirality of a specific atom

        Parameters
        ----------
        atom: int
            The atom index of the atom to invert the chirality

        """

        # The coordinates are not modified. Thus if you need the 3D you need to generate a conformer.
        #TODO
        # option for activate coordinates fixing
        _chiral = self.chiral[atom]
        _chiraltag = self._chiraltags[atom]

        if _chiral == 'R':
            self.chiral[atom] = 'S'
        elif _chiral == 'S':
            self.chiral[atom] = 'R'

        if _chiraltag == 1:
            _chiraltag = _chiraltag + 1
        else:
            _chiraltag = _chiraltag - 1

        self._chiraltags[atom] = _chiraltag

    def toRdkitMol(self, includeConformer=False, _debug=False):
        """
        Returns the rdkit.Chem.rdchem.Mol object.

        Parameter
        ---------
        includeConformer: bool
            If True, also the conformers coordinates are returned

        Returns
        -------
        rdkitmol: rdkit.Chem.rdchem.Mol
            The rdkit Molecule

        """
        #Development note: _debug for debugging purpose. If True, the rdkit molecule are not sanitized

        from rdkit.Chem import RWMol
        from rdkit.Chem import Atom
        from htmd.smallmol.chemlab.periodictable import _hybridizations_IdxToType, _bondtypes_IdxToType, _chiral_type_Dict

        rw = RWMol()

        # Development note: removed explicitHs
        for n in range(self.numAtoms):
            a = Atom(self.element[n])
            a.SetFormalCharge(int(self.formalcharge[n]))
            a.SetNoImplicit(1)
            if self.chiral[n] != '':
                a.SetProp('_CIPCode', self.chiral[n])
            chiral_type = _chiral_type_Dict[self._chiraltags[n]]
            a.SetChiralTag(chiral_type)
            a.SetHybridization(_hybridizations_IdxToType[self.hybridization[n]])
            rw.AddAtom(a)

        # add bonds
        for aIdx, (idxs, bonds) in enumerate(zip(self.neighbors, self.bondtypes)):
            for n, bt in zip(idxs, bonds):
                # TODO clean this behaviuour
                try:
                    if isinstance(bt, int):
                        bt = _bondtypes_IdxToType[bt]
                    elif isinstance(bt, np.int64):
                        bt = _bondtypes_IdxToType[int(bt)]
                        n = int(n)
                    rw.AddBond(aIdx, n, bt)
                except:
                    pass
        mol = rw.GetMol()

        if includeConformer:
            for id, conf in enumerate(self.getConformers()):
                conf.SetId(id)
                mol.AddConformer(conf)

        if _debug:
            return mol

        Chem.SanitizeMol(mol)
        return mol

    def toMolecule(self, formalcharges=False, ids=None):
        """
        Return the htmd.molecule.molecule.Molecule

        Parameters
        ----------
        formalcharges: bool
            If True,the forrmal charges are used instead of partial ones
        ids: list
            The list of conformer ids to store in the htmd Molecule object- If None, all are returned

        Returns
        -------
        mol: htmd.molecule.molecule.Molecule
         The htmd Molecule object

        """
        from htmd.molecule.molecule import Molecule
        class NoConformerError(Exception):
            pass

        _nConformers = self.numConformers
        if _nConformers == 0:
            raise NoConformerError("No Conformers are found in the molecule. Generate at least one conformer.")

        if ids == None:
            ids = list(range(_nConformers))

        elif not isinstance(ids, list):
            raise ValueError('The argument ids should be a list of confomer ids')

        molHtmd = None
        for n in ids:
            coords = self.getCoords(confId=n)
            elements = self.element
            mol = Molecule()
            mol.empty(self.numAtoms)
            mol.resname[:] = self.getProp('ligname')[:3]
            mol.resid[:] = 1
            mol.name[:] = elements
            mol.element[:] = elements
            if formalcharges:
                mol.charge[:] = self.charge
            else:
                mol.charge[:] = self.formalcharge
            mol.coords[:, :, 0] = coords
            mol.viewname = self.getProp('ligname')
            mol.bonds, mol.bondtype = self.getBonds()
            if molHtmd == None:
                molHtmd = mol
            else:
                molHtmd.appendFrames(mol)
        return molHtmd

    def getBonds(self):
        """
        Returns the rdkit bonds object and their type as two np.array

        Returns
        -------
        bonds, bondstype: numpy.array, numpy.array
         An array of rdkit bonds objects and an array with their types
        """

        from rdkit.Chem.rdchem import BondType
        bonds = []
        bondtypes = []
        _mol = self.toRdkitMol()
        for bo in _mol.GetBonds():#self._mol.GetBonds():
            bonds.append([bo.GetBeginAtomIdx(), bo.GetEndAtomIdx()])
            if bo.GetBondType() == BondType.SINGLE:
                bondtypes.append('1')
            elif bo.GetBondType() == BondType.DOUBLE:
                bondtypes.append('2')
            elif bo.GetBondType() == BondType.TRIPLE:
                bondtypes.append('3')
            elif bo.GetBondType() == BondType.AROMATIC:
                bondtypes.append('ar')
        return np.vstack(bonds), np.array(bondtypes)

    def depict(self, sketch=False, filename=None, ipython=False, optimize=False, optimizemode='std', removeHs=True,
                     atomlabels=None, highlightAtoms=None):
        """
        Depicts the molecules. It is possible to save it into an svg file and also generates a jupiter-notebook rendering

        Parameters
        ----------
        sketch: bool
            Set to True for 2D depiction
        filename: str
            Set the filename for the svg file
        ipython: bool
            Set to True to return the jupiter-notebook rendering
        optimize: bool
            Set to True to optimize the conformation. Works only with 3D.
        optimizemode: ['std', 'mmff']
            Set the optimization mode for 3D conformation
        removeHs: bool
            Set to True to hide hydrogens in the depiction
        atomlabels: str
            Accept any combinations of the following pararemters as unique string '%a%i%c%*' a:atom name, i:atom index,
            c:atom formal charge (+/-), *:chiral (* if atom is chiral)
        highlightAtoms: list
            List of atom to highlight. It can be also a list of atom list, in this case different colors will be used

        Returns
        -------
            ipython_svg: SVG object if ipython is set to True

        Example
        -------
        >>> sm.depict(ipython=True, optimize=True, optimizemode='std')
        >>> sm.depict(ipython=True, sketch=True)
        >>> sm.depict(ipython=True, sketch=True)
        >>> sm.depict(ipython=True, sketch=True, atomlabels="%a%i%c")
        >>> ids = np.intersect1d( sm.get('idx', 'hybridization SP2'),  sm.get('idx', 'element C'))
        >>> sm.depict(ipython=True, sketch=True,highlightAtoms=ids.tolist(), removeHs=False)

        """
        from rdkit import Chem
        from rdkit.Chem.AllChem import Compute2DCoords, EmbedMolecule, MMFFOptimizeMolecule, ETKDG

        if sketch and optimize:
            raise ValueError('Impossible to use optmization in  2D sketch representation')

        if optimizemode not in ['std', 'mmff']:
            raise ValueError('Optimization mode {} not understood. Can be "std" or  "ff"'.format(optimizemode))

        #_mol = self._mol
        _mol = self.toRdkitMol(includeConformer=True)

        elements = self.element
        indexes = self.idx
        formalcharges = self.formalcharge
        chirals = self.chiral

        if sketch:
            Compute2DCoords(_mol)

        if removeHs:
            _mol = Chem.RemoveHs(_mol)
            elements = self.get('element', 'element H', invert=True)
            indexes = self.get('idx', 'element H', invert=True)
            formalcharges = self.get('formalcharge', 'element H', invert=True)
            chirals = self.get('chiral', 'element H', invert=True)

        _labelsFunc = ['a', 'i', 'c', '*']

        if atomlabels is not None:
            labels = atomlabels.split('%')[1:]
            formalcharges = ['' if c == 0 else "+" if c == 1 else "-" for c in formalcharges ]
            chirals = ['' if c == '' else '*'  for c in chirals ]
            values = [elements, indexes, formalcharges, chirals]

            idxs = [ _labelsFunc.index(l) for l in labels]
            labels_required = [values[i] for i in idxs]
            atomlabels = ["".join([str(i) for i in a]) for a in list(zip(*labels_required))]

        if optimize:
            if optimizemode == 'std':
                EmbedMolecule(_mol, ETKDG())
            elif optimizemode == 'mmff':
                MMFFOptimizeMolecule(_mol)


        return _depictMol(_mol, filename=filename, ipython=ipython,  atomlabels=atomlabels, highlightAtoms=highlightAtoms)

    def __repr__(self):
        return '<{}.{} object at {}>\n'.format(self.__class__.__module__, self.__class__.__name__, hex(id(self))) \
               + self.__str__()

    def __str__(self):
        def formatstr(name, field):
            if isinstance(field, np.ndarray) or isinstance(field, list):
                rep = '{} shape: {}'.format(name, np.shape(field))
            elif field == 'reps':
                rep = '{}: {}'.format(name, len(self.reps.replist))
            else:
                rep = '{}: {}'.format(name, field)
            return rep

        rep = 'SmallMol with ' + str(self.numAtoms)
        for p in sorted(self._atom_fields):
            if p.startswith('_'):
                continue
            rep += '\n'
            rep += 'Atom field - ' + formatstr(p, self.__dict__[p])
        for j in sorted(self.__dict__.keys() - list(SmallMol._atom_fields)):
            if j[0] == '_':
                continue
            rep += '\n'
            rep += formatstr(j, self.__dict__[j])

        return rep

# multiprocess does not play well with class methods.
def unwrap_self(arg, **kwarg):
    return SmallMolLib.vox_fun(arg[0], arg[1], **kwarg)

class SmallMolLib:
    """
    Class to manage ligands databases (sdf). Ligands are stored as htmd.smallmol.smallmol.SmallMol objects and
    fields type in the sdf are stored in a list

    Parameters
    ----------
    sdf_file: str
        The sdf file path
    removeHs: bool
        If True, the hydrogens of the molecules will be removed
    fixHs: bool
        If True, the hydrogens are added and optimized

    Example
    -------
    >>> lib = SmallMolLib('htmd/data/test-smallmol/fda_drugs_light.sdf')
    100/100 [00:00<00:00, 546.10it/s]
    >>> lib.numMols
    100

    .. rubric:: Methods
    .. autoautosummary:: htmd.smallmol.smallmol.SmallMolLib
       :methods:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.smallmol.smallmol.SmallMolLib
       :attributes:

    Attributes
    ----------
    numMols: int
        Number of SmallMol molecules

    """

    def __init__(self, sdf_file=None, removeHs=False, fixHs=True):  # , n_jobs=1
        self._sdffile = sdf_file if  self._isSdfFile(sdf_file) else None

        self._mols = np.array([])
        self.fields = []

        if sdf_file != None:
            self._initializeMolObjs(sdf_file, removeHs, fixHs)

    def _isSdfFile(self, sdf_file):
        """
        Returns True if the file exists
        sdf_file: str
            The sdf file
        """

        if sdf_file == None: return None

        if not os.path.isfile(sdf_file):
            raise FileNotFoundError('The sdf file {} does not exist'.format(sdf_file))

        sdf_ext = os.path.splitext(sdf_file)[-1]
        if sdf_ext != '.sdf':
            raise TypeError('The file extension {} is not valid. Should be .sdf'.format(sdf_ext))

        return True

    def _initializeMolObjs(self, sdf_file, removeHs, fixHs):
        """
        Processes and loads the molecules inside in the sdf file

        Parameters
        ----------
        sdf_file: str
            The sdf file
        removeHs: bool
            If True, the hydrogens are removed
        fixHs: bool
            If True,, the hydrogens are added and optmized
        """

        from tqdm import tqdm
        supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
        mols_failed = []

        for i, mol in enumerate(tqdm(supplier)):
            if mol is not None:
                m = SmallMol(mol, removeHs=removeHs, fixHs=fixHs)
                self.appendSmallMol(m, strictDirection=2)
            else:
                mols_failed.append(i)

        if len(mols_failed) != 0:
            logger.warning('The following entries were skipped beacause could not be loaded: {}.'.format(mols_failed))

    @property
    def numMols(self):
        """
        Returns the number of molecules
        """
        return len(self._mols)

    def getMols(self, ids=None):
        """
        Returns the SmallMol objects that corresponds ot the indexes of the list passed

        Parameters
        ----------
        ids: list
            The index list of the molecules to return

        Returns
        -------
        smallmollist: list
            The list of SmallMol objects

        Example
        -------
        >>> lib2 = lib.getMols([1,2,3])
        >>> len(lib2)
        3
        """

        if ids == None:
            return self._mols
        if not isinstance(ids, list):
            raise TypeError("The argument ids {} should be list".format(type(ids)))
        _mols = np.array(self._mols)

        return _mols[ids]

    def writeSdf(self, sdf_name, fields=None):
        """
        Writes an sdf file with molecules stored. Is it possible also to manage which field will be written

        Parameters
        ----------
        sdf_name: str
            The ouput sdf filename
        fields: list
            A list of the fields to write. If None all are saved
        """

        from rdkit.Chem import SDWriter

        writer = SDWriter(sdf_name)
        if fields is not None:
            if not isinstance(fields, list):
                raise TypeError("The fields argument {} should be a list".format(type(fields)))
            writer.SetProps(fields)

        for m in self._mols:
            writer.write(m._mol)

    def appendSmallLib(self, smallLib, strictField=False, strictDirection=1):
        """
        Merge two htmd.smallmol.smallmol.SmallMolLib objects

        Parameters
        ----------
        smallLib: htmd.smallmol.smallmol.SmallMolLib
            The new SmallMolLib to merge
        strictField: bool
            If True, the new SmallMolLib can be merged only if they have exactly the same fields
        strictDirection: int
            The valid options are 1 or 2 only. With 1 only the fields of the current SmallMolLib are added to the new one.
            With 2 also the fields of the new SmallMolLib are added into the current one.

        """

        ### origianl sdf_filename should i store it???
        from tqdm import tqdm

        for sm in tqdm(smallLib._mols):
            self.appendSmallMol(sm, strictField, strictDirection)

    def appendSmallMol(self, smallmol, strictField=False, strictDirection=1):
        """
        Adds a new htmd.smallmol.smallmol.SmallMol object in the current SmallMolLib object

        Paramters
        ---------
        smallmol: htmd.smallmol.smallmol.SmallMol
            The SmallMol object to add
        strictField: bool
            If True, the new SmallMolLib can be merged only if they have exactly the same fields
        strictDirection: int
            The valid options are 1 or 2 only. With 1 only the fields of the current SmallMolLib are added to the new one.
            With 2 also the fields of the new SmallMolLib are added into the current one.
        """

        #### check fields and in case as zero ?  the same for the ones present?
        class NoSameField(Exception):
            pass

        if strictDirection  not in [1,2]:
            raise  ValueError("The strictDirections should be 1 (add fields into new mol) or 2 (add fields also in the"
                              "database mols)")
        tmp_fields = smallmol.listProps

        if strictField:
            areSameField = set(self.fields) == set(tmp_fields)
            if not areSameField:
                raise NoSameField("The fields of the new molecule does not match the current database. Set strictField "
                                  "as False to skip this error")

        # TODO
        #improve the speed of the following part of the code.
        ### slow down??
        if strictDirection >= 1:
            old_fields = set(self.fields) - set(tmp_fields)
            for f in old_fields:
                smallmol.setProp(f, np.nan, True)
        if strictDirection == 2:
            new_field = set(tmp_fields) - set(self.fields)
            self.fields = self.fields + list(new_field)
            for f in new_field:
                for m in self._mols:
                    m.setProp(f, np.nan, True)

        self._mols = np.append(self._mols, smallmol)

    def removeMols(self, ids):
        """
        Removes the htmd.smallmol.smallmol.SmallMol object based on the indexes in the list

        Parameters
        ----------
        ids: list
            The list of molecules index to remove from the SmallMolLib
        """

        if not isinstance(ids, list):
            raise TypeError('The argument ids {} is not valid. Should be list'.format(type(ids)))
        _oldNumMols = self.numMols
        self._mols = np.delete(self._mols, ids)

        logger.warning("[num mols before deleting: {}]. The molecules {} were removed, now the number of "
                    "molecules are {} ".format(_oldNumMols, ids, self.numMols))

    def toDataFrame(self, fields=None, molAsImage=True, sketch=True):
        """
        Returns a pandas.DataFrame of the SmallMolLib object.

        Parameters
        ----------
        fields: list
            The list of fields to convert into a pandas DataFrame column
        molAsImage: bool
            If True, the rdkit.Chem.rdchem.Mol is converted into an image
        sketch: bool
            If True, the molecule are rendered to be 2D

        Returns
        -------
        dataframe: pandas.DataFrame
            The pandas DataFrame
        """
        from rdkit.Chem import PandasTools
        from rdkit.Chem.AllChem import Compute2DCoords
        import pandas as pd

        if fields is not None:
            if not isinstance(fields, list):
                raise TypeError('The argument fields passed {} should be a list '.format(type(fields)))
        else:
            firstFields = ['ligname', '_mol', 'totalcharge'] if molAsImage else ['ligname', 'totalcharge']
            fields = firstFields + list(set(self.fields) - set(firstFields))

        records = []
        indexes = []
        for i, m in enumerate(self._mols):
            row = dict((f, m.getProp(f))  for f in fields)
            if sketch:
                mm = deepcopy(m._mol)
                Compute2DCoords(mm)
                row['_mol'] = mm
            records.append(row)
            indexes.append(i)

        df = pd.DataFrame(records, columns=fields, index=indexes)
        if molAsImage:
            Chem.PandasTools.ChangeMoleculeRendering(df)
        return df

    def copy(self):
        """
        Returns a copy of the SmallMolLib object
        """
        return deepcopy(self)

    def voxFun(mol):
        return None

    def __len__(self):
        return len(self._mols)

    def __getitem__(self, item):
        return self._mols[item]

    def __iter__(self):
        _mols = self._mols
        for smallmol in _mols:
            yield smallmol

    def __str__(self):
        _mols = self._mols()

        return ('Stack of Small molecules.'
                '\n\tContains {} Molecules.'
                '\n\tSource file: "{}".').format(len(_mols), self._sdffile)

    def voxel_generator(self, batch_size=32, center=None, boxsize=24, resolution=1., n_jobs=1):
        """
        Batch voxel generator.abs

        Parameters
        ----------
        batch_size: int
            The size to yield each batch.
        center:
            Either a list of centers or a np.array of shape (3, ) contanining a single one.
            By default it chooses its molecule's geometrical center.
        boxsize: int
            Resulting size of voxelized array.
        resolution: float
            Resolution in Amstrong of the resulting array.
        n_jobs: int
            Number of threads to use during voxelization.
        """
        from htmd.smallmol.util import _getGridCenters

        # Cache the box centers
        if (boxsize, resolution) not in SmallMol.array_cache:
            bbm = (np.zeros(3) - float(boxsize * resolution / 2))
            SmallMol.array_cache[(boxsize, resolution)] = \
                _getGridCenters(bbm, [boxsize]*3, 1.).reshape(boxsize ** 3, 3)

        num_batches = math.ceil(self.__len__() / batch_size)

        # Setup voxelization
        def get_vox(mol, xcenter=None):
            if mol is None:
                return None
            return SmallMol.get_voxels(mol, center=xcenter, size=boxsize, resolution=resolution)
        SmallMolStack.vox_fun = get_vox

        # Generate batches of data:
        if n_jobs == 1:  
            for batch in range(num_batches):
                idx_mols = enumerate(self._mols[batch * batch_size: (batch + 1) * batch_size])
                yield [get_vox(mol, center[i+(batch_size*batch)] if isinstance(center, list) else center)
                       for i, mol in idx_mols]

        elif n_jobs > 1 and n_jobs <= multiprocessing.cpu_count(): 
            pool = multiprocessing.Pool(n_jobs)
            for batch in range(num_batches):
                idx_mols = enumerate(self._mols[batch * batch_size: (batch + 1) * batch_size])

                yield pool.map(unwrap_self, [[mol, center[i+(batch_size*batch)]
                                             if isinstance(center, list)
                                             else center]
                                             for i, mol in idx_mols])
            pool.close()
        else:
            raise ValueError("n_jobs needs to be a positive integer!")

    def depict(self, ids=None, sketch=False, filename=None, ipython=False, optimize=False, optimizemode='std',
               removeHs=True,  legends=None, highlightAtoms=None, mols_perrow=3):

        """
        Depicts the molecules into a grid. It is possible to save it into an svg file and also generates a
        jupiter-notebook rendering

        Parameters
        ----------
        ids: list
            The index of the molecules to depict
        sketch: bool
            Set to True for 2D depiction
        filename: str
            Set the filename for the svg file
        ipython: bool
            Set to True to return the jupiter-notebook rendering
        optimize: bool
            Set to True to optimize the conformation. Works only with 3D.
        optimizemode: ['std', 'mmff']
            Set the optimization mode for 3D conformation
        removeHs: bool
            Set to True to hide hydrogens in the depiction
        legends: str
            The name to used for each molecule. Can be 'names':the name of themselves; or 'items': a incremental id
        highlightAtoms: list
            A List of atom to highligh for each molecule. It can be also a list of atom list, in this case different
            colors will be used
        mols_perrow: int
            The number of molecules to depict per row of the grid

        Returns
        -------
            ipython_svg: SVG object if ipython is set to True

        """
        from rdkit.Chem.AllChem import Compute2DCoords, EmbedMolecule, MMFFOptimizeMolecule, ETKDG
        from rdkit.Chem import RemoveHs

        if sketch and optimize:
            raise ValueError('Impossible to use optmization in  2D sketch representation')

        if legends is not None and legends not in ['names', 'items']:
            raise ValueError('The "legends" should be "names" or "items"')



        _smallmols = self.getMols(ids)

        if ids is None:
            _mols = [ _m.toRdkitMol() for _m in self._mols ]
        else:
            _mols = [ _m.toRdkitMol() for _m in self.getMols(ids)]

        if highlightAtoms is not None:
            if len(highlightAtoms) != len(_mols):
                raise ValueError('The highlightAtoms {} should have the same length of the mols {}'.format(len(highlightAtoms), len(_mols)))

        if sketch:
            for _m in _mols: Compute2DCoords(_m)

        if removeHs:
            _mols = [ RemoveHs(_m) for _m in _mols ]

        # activate 3D coords optimization
        if optimize:
            if optimizemode == 'std':
                for _m in _mols: EmbedMolecule(_m)
            elif optimizemode == 'mmff':
                for _m in _mols: MMFFOptimizeMolecule(_m, ETKDG())

        legends_list = []
        if legends == 'names':
            legends_list = [ _m.getProp('ligname') for _m in _smallmols ]
        elif legends == 'items':
            legends_list = [ str(n+1) for n in range(len(_smallmols))]

        return depictMultipleMols(_mols, ipython=ipython, legends=legends_list, highlightAtoms=highlightAtoms,
                                         filename=filename, mols_perrow=mols_perrow)

SmallMolStack = SmallMolLib

if __name__ == '__main__':
    pass
