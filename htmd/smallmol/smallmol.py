import os
import multiprocessing
import math
import numpy as np

from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.rdchem import BondType


from htmd.molecule.voxeldescriptors import _getOccupancyC, _getGridCenters
from htmd.smallmol.util import get_rotationMatrix, rotate, InputToOutput, _depictMol, depictMultipleMols, convertToString, _highlight_colors
from copy import deepcopy
import logging

logger = logging.getLogger(__name__)


rdBase.DisableLog('rdApp.error')
atom_mapping = {"Hydrophobe": 0,
                "LumpedHydrophobe": 0,
                "Aromatic": 1,
                "Acceptor": 2,
                "Donor": 3,
                "PosIonizable": 4,
                "NegIonizable": 5}

_chiral_type = {ChiralType.CHI_TETRAHEDRAL_CW:'clockwise',
                ChiralType.CHI_TETRAHEDRAL_CCW:'anitclockwise'}

_chiral_type_Dict = ChiralType.values

_heteroatoms = ['N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
_halogens = ['F', 'Cl', 'Br', 'I']

_hybridizations = {'S': HybridizationType.S,
                   'SP': HybridizationType.SP,
                   'SP2': HybridizationType.SP2,
                   'SP3': HybridizationType.SP3}

_bondtypes = {'SINGLE': BondType.SINGLE,
              'DOUBLE': BondType.DOUBLE,
               'TRIPLE': BondType.TRIPLE}

# multiprocess does not play well with class methods.
def unwrap_self(arg, **kwarg):
    return SmallMolStack.vox_fun(arg[0], arg[1], **kwarg)


class SmallMol:
    """
    SmallMol class using RDkit for featurization.
    """

    array_cache = {(16, 1.): _getGridCenters(np.array([-8] * 3),
                                               [16, 16, 16], 1.).reshape(16**3, 3),
                   (24, 1.): _getGridCenters(np.array([-12] * 3),
                                               [24, 24, 24], 1.).reshape(24**3, 3)}

    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    ### Fields
    _atom_fields = ['idx', 'atomname', 'charge','formalcharge', 'element',  'chiral', 'hybridization',
                    'neighbors', 'bondtypes', 'coords', 'expliciths', '_chiraltags']

    _mol_fields = ['ligname', 'totalcharge', '_mol']

    # field types
    _atom_dtypes = {'idx': np.int,
                    'atomname': object,
                    'charge': np.float32,
                    'formalcharge': np.int,
                    'element': object,
                    'chiral': object,
                    'hybridization': object,
                    'neighbors': object,
                    'bondtypes': object,
                    'coords': np.float32,
                    'expliciths': np.int,
                    '_chiraltags': object
                    }

    _mol_dtypes = {'_mol': object,
                    'ligname': object,
                   'totalcharge': np.int
                   }

    # field shapes
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
                  'expliciths':(0,),
                  '_chiraltags':(0,)
                    }

    _mol_dims = {'_mol': (0,),
                 'ligname': (0,),
                 'totalcharge': (0,)
                   }

    def __init__(self, mol, ignore_errors=False, force_reading=False, fixHs=True, removeHs=False):
        """
        Initializes small molecule object

        Parameters
        ----------
        mol: rdkit Molecule object or string
            (i) Rdkit molecule or (ii) Location of molecule file (".pdb"/".mol2") or (iii) a smile string or iv) another
            SmallMol object
        ignore_errors: bool
            If True errors will not be raised.
        force_reading: bool
            If the mol provided is not accepted, the molecule will be initially converted into sdf
        fixHs: bool
            The missing hydrogens are assigned, the others are correctly assinged into the graph of the molecule
        removeHs: bool
            Set as True to remove the hydrogens
        """

        self.mol_fields = self._mol_fields.copy()
        self.atom_fields = self._atom_fields.copy()

        # init the __dict__ with the fields
        for field in self._atom_dtypes:
            self.__dict__[field] = np.empty( self._atom_dims[field], dtype=self._atom_dtypes[field] )
        for field in self._mol_dtypes:
            self.__dict__[field] = None

        # load the input and store the rdkit Mol obj
        self._mol, smallMol = self._initializeMolObj(mol, force_reading)
        if not ignore_errors and self._mol == None:
            if not force_reading:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj. Try by setting the force_reading option as True.".format(mol))
            else:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj.".format(mol))

        if removeHs:
            self._mol = Chem.RemoveHs(self._mol)

        if fixHs:
            self._mol = Chem.AddHs(self._mol, addCoords=True)

        # fill the molecule and atom properties
        self._readMol(smallMol)

    def _initializeMolObj(self, mol, force_reading):
        """
        Read the input and it try to convert it into a rdkit Molecule obj

        Parameters
        ----------
        mol: str or rdkit Molecule object
            i) rdkit Molecule Object ii) The path to the pdb/mol2 to load iii) The smile string iv) SmallMOl object
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
                if name_suffix == ".mol2":
                    _mol = Chem.MolFromMol2File(mol, removeHs=False)
                elif name_suffix == ".pdb":
                    _mol = Chem.MolFromPDBFile(mol, removeHs=False)
                if _mol == None and force_reading:
                    logger.warning('Reading {} with force_reading procedure'.format(mol))
                    sdf = InputToOutput(mol, name_suffix, 'sdf')
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

        idxs = []
        atomnames = []
        charges = []
        formalcharges = []
        elements = []
        chirals =  []
        hybridizations = []
        neighbors = []
        bondtypes = []
        expliciths = []
        chiraltags = []

        _mol = self._mol

        if smallMol is not None:
            for f, v in smallMol.__dict__.items():
                self.__dict__[f] = v
            return

        for a in _mol.GetAtoms():
            i = a.GetIdx()
            e = a.GetSymbol()
            chiraltags.append(a.GetChiralTag())
            idxs.append(i)
            atomnames.append('{}{}'.format(e,i))
            formalcharges.append(a.GetFormalCharge())
            elements.append(e)
            hybridizations.append(a.GetHybridization())
            neighbors.append([na.GetIdx() for na in a.GetNeighbors()])
            bondtypes.append([b.GetBondType() for b in a.GetBonds()])
            expliciths.append(a.GetNumExplicitHs())
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

        self.__dict__['totalcharge'] = sum(formalcharges)
        self.__dict__['idx'] = np.array(idxs)
        self.__dict__['atomname'] = np.array(atomnames)
        self.__dict__['charge'] = np.array(charges )
        self.__dict__['formalcharge'] = np.array(formalcharges)
        self.__dict__['element'] = np.array(elements)
        self.__dict__['chiral'] = np.array(chirals)
        self.__dict__['hybridization'] = np.array(hybridizations, dtype=object)
        self.__dict__['neighbors'] = np.array(neighbors)
        self.__dict__['bondtypes'] = np.array(bondtypes)
        self.__dict__['expliciths'] = np.array(expliciths)
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
        return deepcopy(self)

    @property
    def listProps(self):
        """
        The list of properties of the moleule

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
            If set as True the property will be overwritten if already exist

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
        """
        class ValuesLengthErrors(Exception):
            pass

        if not isinstance(key, str):
             raise ValueError('Wrong type {} for key.  Should be {} '.format(type(key), type('string')))

        if aIdxs is None and len(values) != self.numAtoms:
            raise ValuesLengthErrors('The number of values passed {} do not match '
                                     'the number of atoms {}'.format(len(values), self.numAtoms))

        if aIdxs is not None and len(aIdxs) != len(values):
            raise ValuesLengthErrors('The number of values passed {} do not match '
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
            The name of the property to retrieve

        Returns
        -------
        value: str, int, float, object
            The value of the property
        """
        if key not in self.listProps:
            raise KeyError('The property passed {} does not exist'.format(key))

        _value = self.__dict__[key]

        return _value

    def get(self, sel, returnField, invert=False):
        """
        Returns the property for the atom specified with the selection. The selection is another atom property

        Parameters
        ----------
        sel: str
            The selection string like: element H C; idx 1 2 3; formalcharge 1
        returnField: str
            The field of the atom to return

        Returns
        -------
        values: np.array
            The array of values for the property
        """

        key = sel.split()[0]
        selector = sel.split()[1:]

        if key not in self.listPropsAtoms:
            raise KeyError('The property passed {} does not exist'.format(key))
        if len(selector) == 0:
            raise ValueError('No selection was provided')

        if key == 'hybridization':
            selector = [_hybridizations[s.upper()] for s in selector]

        _arrayFrom = self.__dict__[key]
        if returnField == 'atomobject':
            _arrayTo = self.getAtoms()
        else:
            _arrayTo = self.__dict__[returnField]

        _dtype = self._atom_dtypes[key]
        if _dtype is not object:
            selector = [ _dtype(s) for s in selector ]
        idxs = np.concatenate([ np.where( _arrayFrom == s )[0] for s in selector ])
        if invert:
            idxs = np.array([ i for i in self.idx if i not in idxs ])
        idxs = np.sort(idxs)

        return _arrayTo[idxs]

    @property
    def heteroAtoms(self):
        sel_str = convertToString(_heteroatoms)

        atom_idxs = self.get('element {}'.format(sel_str), 'idx')

        return atom_idxs

    @property
    def halogenAtoms(self):
        sel_str = convertToString(_halogens)

        atom_idxs = self.get('element {}'.format(sel_str), 'idx')

        return atom_idxs

    def isChiral(self, returnDetails=False):
        """
        Returns if the molecule has at least one chiral atom. If returnDetails is set as True, a list of tuples with the
        atom idx and chiral type is returned.

        Parameters
        ----------
        returnDetails: bool (default=False)
            Set as True to return the chiral atoms and their chiral types

        Returns
        -------
        ischiral: bool
            True if the atom has at least a chiral atom
        details: list
            A list of tuple with the chiral atoms and their types
        """

        _chirals = self.chiral

        idxs = np.where(_chirals != '')[0]
        if len(idxs) == 0:
            return False
        if returnDetails:
            idxs =  idxs.astype(str)
            idxs_str = ' '.join(idxs)
            atomnames = self.get('idx {}'.format(idxs_str), 'atomname')
            chirals = self.get('idx {}'.format(idxs_str), 'chiral')

            _details = [(a, c) for a, c in zip(atomnames, chirals)]

            return True, _details

        return True

    def foundBondBetween(self, sel1, sel2, bondtype=None):

        _bondtypes_reversed = { v:k for k, v in _bondtypes.items()}

        atomIdx_sel1 = self.get(sel1, 'idx')
        atomIdx_sel2 = self.get(sel2, 'idx')

        neighbors_sel1 = self.get(sel1, 'neighbors')
        bondtypes_sel1 = self.get(sel1, 'bondtypes')

        founds = []
        for aIdx, neighbors_set, btype_set in zip(atomIdx_sel1, neighbors_sel1, bondtypes_sel1):
            for neighbor, btype in zip(neighbors_set, btype_set):
                if neighbor in atomIdx_sel2:
                    btype_str = _bondtypes_reversed[btype]
                    if bondtype is not None and _bondtypes[bondtype.upper()] == btype:
                        founds.append([(aIdx, neighbor), btype_str])
                    elif bondtype == None:
                        founds.append( [(aIdx, neighbor), btype_str]  )
        if len(founds) != 0:
            return True, founds

        return False

    def getAtoms(self):
        """
        Retuns all the rdkit.Chem.rdchem.Atom present in the molecule
        """

        #_mol = self._mol
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
        # n_atoms = self.numAtoms
        #
        # conformer = self.getConformers([confId])[0]
        # coords = [[corobj.x, corobj.y, corobj.z] for corobj in [conformer.GetAtomPosition(i) for i in range(n_atoms)]]
        # return np.array(coords, dtype=np.float32)
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
        n_atoms = self.numAtoms
        _mol = self.toRdkitMol()

        feats = SmallMol.factory.GetFeaturesForMol(_mol)
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
        Returns geometrical center of molecule.
        """
        if coords is None:
            coords = self.getCoords(confId)
        return coords.mean(axis=0).astype(np.float32)

    def generateConformers(self, num_confs=400,  optimizemode='mmff', align=True, append=True):
        """
        Generates ligand conformer and saves the results to a folder.


        Parameters
        ----------
        num_confs: int
           Number of conformers to generate.
        optimizemode: str, (default='mmff')
            The optimizemode to use. Can be  'uff', 'mmff'
        append: bool
            If set as False the current conformers are deleted, excepted for the first one.

        """
        from rdkit.Chem.AllChem import UFFOptimizeMolecule, MMFFOptimizeMolecule, EmbedMultipleConfs
        from rdkit.Chem.rdMolAlign import AlignMolConformers

        if not append:
            self.removeConformers()

        #_mol = self._mol
        _mol = self.toRdkitMol()
        mol = deepcopy(_mol)
        mol = Chem.AddHs(mol)
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
            #conf = mol.GetConformer(id)
            #_mol.AddConformer(conf, id+1)

        if align:
            AlignMolConformers(mol)
            #AlignMolConformers(_mol)

#        print(list(ids))
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
        if (size, resolution) not in SmallMol.array_cache:
            N = [size, size, size]
            bbm = (np.zeros(3) - float(size * resolution / 2))
            centers = _getGridCenters(bbm, N, resolution)

            # Cache the array
            SmallMol.array_cache[(size, resolution)] = centers.reshape(size**3, 3)
            centers2D = centers + center
        else:
            centers2D = SmallMol.array_cache[(size, resolution)] + center

        voxels = _getOccupancyC(coords.astype(np.float32), centers2D,
                                multisigmas).reshape(size, size, size, 8).astype(dtype)
        return voxels

    @property
    def numAtoms(self):
        """
        Returns the number of atoms in the molecule

        Returns
        -------
        natoms: int
            The number of atoms in the molecule

        """
        #_mol = self._mol
        return self.element.shape[0]
        #return _mol.GetNumAtoms()

    @property
    def numConformers(self):
        """
        Returns the number of conformers in the rdkit molecule object

        Returns
        -------
        nconfs: int
            The number of conformers in the molecule
        """
        #_mol = self._mol
        #_nConformers = _mol.GetNumConformers()
        _nConformers = self.coords.shape[-1]

        return _nConformers

    def _getConformer(self, id):
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
        Returns the conformer of the molecule depends on the id passed. If id is equal to -1, all the conformers are
        returned

        Parameters
        ----------
        ids: list (default=None)
            The list of ids for the molecule conformers to return. If None all the conformers are returned

        Returns
        -------
        _conformer: list of rdkit.Chem.rdchem.Conformer
            The conformer
        """
        #_mol = self._mol
        #_conformers = list(_mol.GetConformers())
        #_nConformers = len(_conformers)
        _nConformers = self.coords.shape[-1]

        if ids == None:
            return [ self._getConformer(i) for i in range(_nConformers) ]

        if  max(ids) >= _nConformers:
            raise IndexError("The ids list contains conformers ids {} that do not exist. Available conformers: {}".format(ids, _nConformers))

        return [self._getConformer(id) for id in ids]

    def writeConformers(self, savefolder='conformations', savename="molConf", filetype="sdf", savefolder_exist_ok=False,
                         merge=False, ids=None):
        """
        Writes conformers in one or multiple files in 'pdb' or 'sdf' file formats.

        Paramters
        ---------
        savefolder: str (default='conformations')
            The name of the folder where to write the files
        savename: str (default='molConf')
            The basename of the output file
        filetype: str ('sdf', 'pdb') (default='sdf')
        savefolder_exist_ok: bool (default=False)
            Set as True to overwrite the output folder
        merge: bool (default=False)
            Set as True to save in a unique file
        ids: list (default=None)
            A list of the conformer ids to save. If None, all are written

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

    def _addAtoms(self, elements):
        """
        elements: list of {'element':'H', 'attachTo':atomIdx, 'bondType':idx/rdkittype,
                           'hybrizidation':idx/rdkithybrid, 'formalcharge':int, 'chiral':''}
        """
        _mol_ref = self.copy()
        constraints = self.idx
        for e in elements:
            if self.expliciths[e['attachTo']] != 0:
                self.expliciths[e['attachTo']] -= 1
            n_atom = self.numAtoms

            self.__dict__['idx'] = np.append(self.idx, n_atom)
            self.__dict__['element'] = np.append(self.element, e['element'])
            self.__dict__['atomname'] = np.append(self.atomname, '{}{}'.format(e['element'], n_atom) )
            self.__dict__['charge'] = np.append(self.charge, 0.000)
            self.__dict__['formalcharge'] = np.append(self.formalcharge, e['formalcharge'])
            self.__dict__['chiral'] = np.append(self.chiral, e['chiral'])
            self.__dict__['hybridization'] = np.append(self.hybridization, e['hybridization'])
            coords = np.array([0,0,0])
            coords = coords.reshape(1,3,1)
            self.__dict__['coords'] = np.concatenate((self.coords, coords), axis=0)
            self.__dict__['neighbors'] = np.array(self.neighbors.tolist() + [[e['attachTo']]])
            self.neighbors[e['attachTo']].append(n_atom)
            self.__dict__['bondtypes'] = np.array(self.bondtypes.tolist() + [[e['bondType']]])
            self.bondtypes[e['attachTo']].append(e['bondType'])
            self.__dict__['expliciths'] = np.append(self.expliciths, 0)
            self.__dict__['_chiraltags'] = np.append(self._chiraltags, 0)

        self._alignMol(_mol_ref, constraints=constraints)

    # working aonly for hydrogens
    def _removeAtoms(self, ids):
        if len(ids) == 0:
            return

        atom_toKeep = np.delete(self.idx, ids)
        atom_mapRenumbering = {idx:n for n, idx in enumerate(atom_toKeep)}

        self.__dict__['idx'] = np.array( list(atom_mapRenumbering.values()) )
        self.__dict__['element'] = np.delete(self.element, ids)
        self.__dict__['atomname'] = np.array([ '{}{}'.format(e, i) for e, i in zip(self.element, self.idx) ])
        self.__dict__['charge'] = np.delete(self.charge, ids)
        self.__dict__['formalcharge'] = np.delete(self.formalcharge, ids)
        self.__dict__['chiral'] = np.delete(self.chiral, ids)
        self.__dict__['hybridization'] = np.delete(self.hybridization, ids)
        self.__dict__['coords'] = np.delete(self.coords, ids, axis=0 )


        incr = np.bincount( np.concatenate(self.neighbors[ids]) )
        ids_forexpliciths = np.arange(incr.shape[0])
        self.expliciths[ids_forexpliciths] += incr
        self.__dict__['expliciths'] = np.delete(self.expliciths, ids)
        self.__dict__['_chiraltags'] = np.delete(self._chiraltags, ids)

        tmp_neighbors = np.delete(self.neighbors, ids)
        tmp_bondtypes = np.delete(self.bondtypes, ids)
        new_neighbors = []
        new_bondtypes = []

        for n_set, bt_set in zip(tmp_neighbors, tmp_bondtypes):
            new_n_set = []
            new_bt_set = []
            for n, bt in zip(n_set, bt_set):
                if n not in ids:
                    if n in atom_mapRenumbering:
                        n = atom_mapRenumbering[n]
                    new_n_set.append(n)
                    new_bt_set.append(bt)
            new_neighbors.append(new_n_set)
            new_bondtypes.append(new_bt_set)
        self.__dict__['neighbors'] = np.array(new_neighbors)
        self.__dict__['bondtypes'] = np.array(new_bondtypes)

    def _alignMol(self, smallmolref, constraints):

        from rdkit.Chem.AllChem import  EmbedMolecule
        from rdkit.Chem.AllChem import UFFGetMoleculeForceField
        from rdkit.Chem import  rdMolAlign

        _refmol = smallmolref.toRdkitMol(includeConformer=True)
        _refconf = smallmolref.getConformers()[0]

        _mol = self.toRdkitMol(includeConformer=True)

        ff = UFFGetMoleculeForceField(_mol)
        for idx in constraints:
            ff.AddFixedPoint(int(idx))
        ff.Minimize()
        #_constraintsMap = {int(n):_refconf.GetAtomPosition(int(n)) for n in constraints }
        #print(EmbedMolecule(_mol, coordMap=_constraintsMap, numZeroFail=10))

        pyO3A = rdMolAlign.GetO3A(_mol, _refmol)
        pyO3A.Align()
        new_coords = _mol.GetConformer().GetPositions()
        new_coords = new_coords[:, :, np.newaxis]

        self.coords = new_coords

    def removeConformers(self, ids=None):
        """
        Deletes the conformers passed

        Parameters
        ----------
        ids: list (default=None)
            The list of conformer id to delete. If None, all are removed except the first one
        """
        # _mol = self._mol
        #
        # _conformers = self.getConformers(ids)
        #
        # if ids is None:
        #     _conformers = _conformers[1:]
        # elif not isinstance(ids, list):
        #     raise TypeError("The ids argument should be list of confermer ids")
        #
        # _conformerIDs = [c.GetId() for c in _conformers]
        #
        # for id in _conformerIDs:
        #     _mol.RemoveConformer(id)

        _nConformers = self.numConformers

        if ids is None:
            ids = np.arange(_nConformers)
        elif not isinstance(ids, list):
            raise TypeError("The ids argument should be list of confermer ids")

        self.coords = np.delete(self.coords, ids, axis=2)

    def invertChirality(self, atom):
        #TODO now is idx
        # atom --> idx.  atomname
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
        from rdkit.Chem import RWMol
        from rdkit.Chem import Atom

        rw = RWMol()

        # add atoms
        # print(self.element, len(self.element))
        # print(self.formalcharge, len(self.formalcharge))
        # print(self.chiral, len(self.chiral))
        # print(self.expliciths, len(self.expliciths))
        #for n, (element, formalcharge, chiral, h) in enumerate(zip(self.element, self.formalcharge, self.chiral, self.expliciths)):
        for n in range(self.numAtoms):

        #    print(n, h)
            a = Atom(self.element[n])
            a.SetFormalCharge(int(self.formalcharge[n]))
            a.SetNumExplicitHs(int(self.expliciths[n]))
            a.SetNoImplicit(1)
            if self.chiral[n] != '':
                a.SetProp('_CIPCode', self.chiral[n])
            chiral_type = _chiral_type_Dict[self._chiraltags[n]]
            a.SetChiralTag(chiral_type)
            #print(element, h)
           # print('>> ', n)
            #print(a, element)
            rw.AddAtom(a)

        # add bonds
        for aIdx, (idxs, bonds) in enumerate(zip(self.neighbors, self.bondtypes)):
            for n, bt in zip(idxs, bonds):
                try:
                    rw.AddBond(aIdx, n, bt)
                except:
                    pass

        mol = rw.GetMol()
        #Chem.SanitizeMol(mol)

        if includeConformer:
            for id, conf in enumerate(self.getConformers()):
                conf.SetId(id)
                mol.AddConformer(conf)

        if _debug:
            return mol
        #try:
        Chem.SanitizeMol(mol)
        #except:
        #    Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL^Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE)
        return mol

    def toMolecule(self, formalcharges=False, ids=None):
        """
        Return the htmd.molecule.molecule.Molecule from the rdkit.Molecule

        Parameters
        ----------
        formalcharges: bool
            Set as True if you want formal charges instead of partial ones
        ids: list
            The list of conformer ids to store in the htmd Molecule object

        Returns
        -------
        mol: htmd.molecule.molecule.Molecule
         The htmd Molecule object

        """
        from htmd.molecule.molecule import Molecule
        class NoConformerError(Exception):
            pass

        #_mol = self._mol
        _mol = self.toRdkitMol(includeConformer=True)
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

        from rdkit.Chem import rdchem
        bonds = []
        bondtypes = []
        _mol = self.toRdkitMol()
        for bo in _mol.GetBonds():#self._mol.GetBonds():
            bonds.append([bo.GetBeginAtomIdx(), bo.GetEndAtomIdx()])
            if bo.GetBondType() == rdchem.BondType.SINGLE:
                bondtypes.append('1')
            elif bo.GetBondType() == rdchem.BondType.DOUBLE:
                bondtypes.append('2')
            elif bo.GetBondType() == rdchem.BondType.TRIPLE:
                bondtypes.append('3')
            elif bo.GetBondType() == rdchem.BondType.AROMATIC:
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
        optimizemode: ['std', 'mmff'], default='std'
            Set the optimization mode for 3D conformation
        removeHs: bool, default=True
            Set to True to hide hydrogens in the depiction
        atomlabels: str
            Accept any combinations of the following pararemters as unique string '%a%i%c%*' a:atom name, i:atom index,
            c:atom formal charge (+/-), *:chiral (* if atom is chiral)
        highlightAtoms: list
            List of atom to highligh. It can be also a list of atom list, in this case different colors will be used

        Returns
        -------
            ipython_svg: SVG object if ipython is set to True

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
            elements = self.get('element H', 'element', invert=True)
            indexes = self.get('element H', 'idx', invert=True)
            formalcharges = self.get('element H', 'formalcharge', invert=True)
            chirals = self.get('element H', 'chiral', invert=True)

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
        #return _depictMol(_mol, sketch, filename, ipython, optimize, optimizemode, removeHs, atomlabels, highlightAtoms)

class SmallMolLib:
    """
    Collection of objects of class SmallMol.
    """

    def __init__(self, sdf_file=None, removeHs=False, fixHs=True):  # , n_jobs=1
        self._sdffile = sdf_file if  self._isSdfFile(sdf_file) else None

        self._mols = np.array([])
        self.fields = []

        if sdf_file != None:
            self._initializeMolObjs(sdf_file, removeHs, fixHs)

    def _isSdfFile(self, sdf_file):

        if sdf_file == None: return None

        if not os.path.isfile(sdf_file):
            raise FileNotFoundError('The sdf file {} does not exist'.format(sdf_file))

        sdf_ext = os.path.splitext(sdf_file)[-1]
        if sdf_ext != '.sdf':
            raise TypeError('The file extension {} is not valid. Should be .sdf'.format(sdf_ext))

        return True

    def _initializeMolObjs(self, sdf_file, removeHs, fixHs):
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
        return len(self._mols)

    def getMols(self, ids=None):

        if ids == None:
            return self._mols
        if not isinstance(ids, list):
            raise TypeError("The argument ids {} should be list".format(type(ids)))
        _mols = np.array(self._mols)

        return _mols[ids]

    def writeSdf(self, sdf_name, fields=None):
        from rdkit.Chem import SDWriter

        writer = SDWriter(sdf_name)
        if fields is not None:
            if not isinstance(fields, list):
                raise TypeError("The fields argument {} should be a list".format(type(fields)))
            writer.SetProps(fields)

        for m in self._mols:
            writer.write(m._mol)

    def appendSmallLib(self, smallLib, strictField=False, strictDirection=1):
        ### sdf_file should i store it???
        from tqdm import tqdm

        for sm in tqdm(smallLib._mols):
            self.appendSmallMol(sm, strictField, strictDirection)

    def appendSmallMol(self, smallmol, strictField=False, strictDirection=1):
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
        #self._mols.append(smallmol)

    def removeMols(self, ids):

        if not isinstance(ids, list):
            raise TypeError('The argument ids {} is not valid. Should be list'.format(type(ids)))
        _oldNumMols = self.numMols
        self._mols = np.delete(self._mols, ids)

        logger.warning("[num mols before deleting: {}]. The molecules {} were removed, now the number of "
                    "molecules are {} ".format(_oldNumMols, ids, self.numMols))

    def toDataFrame(self, fields=None, molAsImage=True, sketch=True):
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
        sketch: bool
            Set to True for 2D depiction
        filename: str
            Set the filename for the svg file
        ipython: bool
            Set to True to return the jupiter-notebook rendering
        optimize: bool
            Set to True to optimize the conformation. Works only with 3D.
        optimizemode: ['std', 'mmff'], default='std'
            Set the optimization mode for 3D conformation
        removeHs: bool, default=True
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

        # return depictMultipleMols(_mols, sketch, filename, ipython, optimize, optimizemode,
                            #    removeHs, legends_list, highlightAtoms, mols_perrow)



SmallMolStack = SmallMolLib

if __name__ == '__main__':
    pass
