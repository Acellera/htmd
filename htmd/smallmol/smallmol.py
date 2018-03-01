import os
import multiprocessing
import math
import numpy as np

from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.rdchem import ChiralType

from htmd.molecule.voxeldescriptors import _getOccupancyC, _getGridCenters
from htmd.smallmol.util import get_rotationMatrix, rotate, InputToOutput, _depictMol, depictMultipleMols
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

    def __init__(self, mol, ignore_errors=False, force_reading=False, fixHs=True, removeHs=False):
        """
        Initializes small molecule object

        Parameters
        ----------
        mol: rdkit Molecule object or string
            (i) Rdkit molecule or (ii) Location of molecule file (".pdb"/".mol2") or (iii) a smile string.
        ignore_errors: bool
            If True errors will not be raised.
        force_reading: bool
            If the mol provided is not accepted, the molecule will be initially converted into sdf
        fixHs: bool
            The missing hydrogens are assigned, the others are correctly assinged into the graph of the molecule
        removeHs: bool
            Set as True to remove the hydrogens
        """

        # load the input and store the rdkit Mol obj
        self._mol = self._initializeMolObj(mol, force_reading)
        if not ignore_errors and self._mol == None:
            if not force_reading:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj. Try by setting the force_reading option as True.".format(mol))
            else:
                raise ValueError("Unkown '{}' provided. Not a valid mol2,pdb,smile, rdkitMol obj.".format(mol))

        if removeHs:
            self._mol = Chem.RemoveHs(self._mol)

        if fixHs:
            self._mol = Chem.AddHs(self._mol, addCoords=True)

        #### STORE INFO INSIDE MOL and ATOMS
        # molecule name
        if not self._mol.HasProp('_Name'):
            self.set_name()

        # set initial SmallMol atom charges
        self._initCharges()

        # atom parameters when generated
        self._initParameters()

    def _initializeMolObj(self, mol, force_reading):
        """
        Read the input and it try to convert it into a rdkit Molecule obj

        Parameters
        ----------
        mol: str or rdkit Molecule object
            i) rdkit Molecule Object ii) The path to the pdb/mol2 to load iii) The smile string
        force_reading: bool
           If the mol provided is not accepted, the molecule will be initially converted into sdf

        Returns
        -------
        _mol: rdkit Molecule object

        """

        _mol = None
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

        return _mol

    def _initParameters(self):
        """
        Stores several properties inside the rdkit atom object that can be used for future purposed. See atom idx in
        case of fragmentation
        """
        _atoms = self.get_atoms()

        for a in _atoms:
            chiral_tag = a.GetChiralTag()
            atom_idx = a.GetIdx()
            self.setPropAtom(a, '_SmallMolChiralTag', chiral_tag)
            self.setPropAtom(a, '_SmallMolAtomIdx', atom_idx)

    def _initCharges(self):
        """
        Initialize the atom property "_SmallMolCharge". If the atoms already have it or one of the following
        ('_TriposPartialCharge', '_GasteigerCharge'), it will use these charges

        """

        # init charges property. If Tripos or Gaesteiger than these charges are set otherwise all are 0
        _PossibleChargeTypes = ['_SmallMolCharge', '_TriposPartialCharge', '_GasteigerCharge']

        #init for the first atom
        _atom = self.get_atom(0)
        _FoundChargeTypes = self.listPropsAtom(_atom)

        matched = set(_PossibleChargeTypes) & set(_FoundChargeTypes)

        returnWarning = False

        if len(matched) == 0:
            for a in self.get_atoms():
                self.setPropAtom(a, '_SmallMolCharge', 0.000)
        else:
            # for loop of possible to get key than use it
            k = [k for k in _PossibleChargeTypes if k in matched][0]
            if k != '_SmallMolCharge':
                for a in self.get_atoms():
                    try:
                        self.setPropAtom(a, '_SmallMolCharge', self.getPropAtom(a, k))
                    except:
                        self.setPropAtom(a, '_SmallMolCharge', 0.000)
                        returnWarning = True

        if returnWarning:
            logger.warning('Found atoms without a charge. The "_SmallMolCharge" for these atoms are set to 0.0.')


    def copy(self):
        """
        Create a copy of the molecule object

        Returns
        -------
        newsmallmol : :class:`SmallMol`
            A copy of the object
        """
        return deepcopy(self)

    def _chiralType(self, atom):
        """
        Returns the chiral type of the passed atom. Can be R, S or None.

        Parameters
        ----------
        atom: rdkit.Chem.rdchem.Atom
            The rdkit atom object

        Returns
        -------
        chiral: str
            The chiral type 'R' or 'S'. None is returned if the atom is not a chiral one
        """

        chiral = atom.GetChiralTag()

        if chiral in _chiral_type.keys():
            return self.getPropAtom(atom, '_CIPCode')
        else:
            return None

    def _set_prop(self, key, value, obj):
        """
        Set the value of a property. Common for atom and molecule

        Parameters
        ----------
        key: str
            The property key
        value: str,bool,int,float
            The property value
        obj: The rdkit.Chem.Mol or rdkit.Chem.rdchem.Atom
            The rdkit molecule or the rdkit atom object

        """

        if isinstance(value, bool):
            addprop = obj.SetBoolProp

        if isinstance(value, float):
            addprop = obj.SetDoubleProp

        elif isinstance(value, int):
            addprop = obj.SetIntProp

        elif isinstance(value, str):
            addprop = obj.SetProp

        addprop(key, value)

    def setProp(self, key, value, overwrite=False):
        """
        Set a molecule property. By default the overwrite options is disabled

        Parameters
        ----------
        key: str
            The key identifying the property
        value: str,bool,float,int
            The value of the property
        overwrite: bool
            Set as True to enable property overwriting

        """

        if not isinstance(key, str):
            raise ValueError('Wrong type {} for key.  Should be {} '.format(type(key), type('string')))

        _props = self.listProps()
        if not overwrite and key in _props:
            raise ValueError('The key passed {} already exists. Set "overwrite" as True to overwrite an existing key'.format(key))

        _mol = self.get_mol()

        self._set_prop(key, value, _mol)


    def setPropAtom(self, atom,  key, value, overwrite=False):
        """
        Set an atom property. By default the overwrite options is disabled

        Parameters
        ----------
        key: str
            The key identifying the property
        value: str,bool,float,int
            The value of the property
        overwrite: bool
            Set as True to enable property overwriting

        """

        if not isinstance(key, str):
            raise ValueError('Wrong type {} for key.  Should be {} '.format(type(key), type('string')))

        _props = list(atom.GetPropsAsDict().keys())
        if not overwrite and key in _props:
            raise ValueError('The key passed {} already exists. Set "overwrite" as True to overwrite an existing key'.format(key))

        self._set_prop(key, value, atom)

    def getPropAtom(self, atom, key):
        """
        Returns the atom property for the atom and key passed

        Parameters
        ----------
        atom: rdkit.Chem.rdhem.Atom
            The rdkit Atom object
        key: str
            The property key

        Returns
        -------
        prop:
            The property
        """

        if key not in self.listPropsAtom(atom):
            raise KeyError('The property "{}" was not found '.format(key))

        prop = atom.GetPropsAsDict()[key]

        return prop

    def getProp(self, key):
        """
        Returns the molecule property for the key passed

        Parameters
        ----------
        key: str
            The property key

        Returns
        -------
        prop:
            The property
        """

        if key not in self.listProps():
            raise KeyError('The property "{}" was not found '.format(key))

        _mol = self.get_mol()
        prop = _mol.GetPropsAsDict()[key]
        return  prop

    def listProps(self):
        """
        Returns a list of the molecule properties
        """
        _mol = self.get_mol()

        return list(_mol.GetPropsAsDict().keys())

    def listPropsAtom(self, atom):
        """
        Returns a list of the atom properties
        """

        return list(atom.GetPropsAsDict().keys())

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

        stereocenters = [ (a.GetIdx(), self._chiralType(a)) for a in self.get_atoms() if self._chiralType(a) is not None]

        if len(stereocenters) == 0:
            return False

        if returnDetails:
            return True, stereocenters
        return True

    def get_mol(self):
        """
        Returns the rdkit Molecule object
        """
        return self._mol

    def get_neighbours(self, atom, returnAsIdx=False):
        """
        Returns the atom neighbours of the atom passed. If returnAsIdx is set as True, the idx of the atoms are returned
        otherwise the rdkit Atom objects

        Parameters
        ----------
        atom: int or rdkit.Chem.rdchem.Atom
            The atom as index or rdkit Atom object
        returnAsIdx: bool
            Set as True if you want to retrieve the neighbours as atom indexes

        Returns
        -------
        neighbours: list
            List of atom neighbours. If returnAsIdx set as True: list of atom indexes, otherwise list of rdkit Atom objs
        """

        from rdkit.Chem.rdchem import Atom

        if isinstance(atom, int):
            _mol = self.get_mol()
            atom = _mol.GetAtomWithIdx(atom)

        if not isinstance(atom, Atom):
            raise ValueError('type {} not valid. Should be "int" or "rdkit.Chem.rdchem.Atom"'.format(type(atom)))

        neighbours = atom.GetNeighbors()

        if returnAsIdx:
            return [ a.GetIdx() for a in neighbours ]

        return  neighbours

    def get_element(self, atom):
        """
        Returns the element of the atom. The rkdit.Chem.rdchem.Atom or the index can be passed

        Parameters
        ----------
        atom: int or rdkit.Chem.rdchem.Atom
            The atom you want to retrieve the element

        Returns
        -------
        element: str
            The element of the atom
        """
        from rdkit.Chem.rdchem import Atom

        if isinstance(atom, int):
            _mol = self.get_mol()
            atom = _mol.GetAtomWithIdx(atom)
        if not isinstance(atom, Atom):
            raise ValueError('type {} not valid. Should be "int" or "rdkit.Chem.rdchem.Atom"'.format(type(atom)))
        element = atom.GetSymbol()
        return element

    def get_atom(self, sel, returnAll=False):
        """
        Return the rdkit.Chem.rdchem.Atom of the sel. The sel can be the idx of the atom or the element. If returnAll
        is set to True, all the incindences all returned as list

        Parameters
        ----------
        sel: int or str
            The selector can be the element or the atom index
        returnAll: bool
            Set to True if you want all the incidences. In this case a list will be returned

        Returns
        -------
        atoms: rdkit.Chem.rdchem.Atom or list
            A single or a list of rdkit.Chem.rdchem.Atom based on the returnAll value
        """

        _mol = self.get_mol()
        if isinstance(sel, int):
            atoms = [_mol.GetAtomWithIdx(sel)]

        elif isinstance(sel, str):
            atoms = [ a for a in self.get_atoms() if self.get_element(a) == sel]
        else:
            raise ValueError('type {} not valid. Should be "int" or "str"'.format(type(sel)))

        if len(atoms) == 0:
            return None

        if not returnAll:
            return atoms[0]
        else:
            return atoms

    def get_atoms(self):
        """
        Retuns all the rdkit.Chem.rdchem.Atom present in the molecule
        """

        _mol = self.get_mol()
        return _mol.GetAtoms()

    def get_coords(self, id=0):
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
        n_atoms = self.get_natoms()

        conformer = self.get_conformers()[id]
        coords = [[corobj.x, corobj.y, corobj.z] for corobj in [conformer.GetAtomPosition(i) for i in range(n_atoms)]]
        return np.array(coords, dtype=np.float32)

    def get_elements(self):
        """
        Returns molecule elements.

        Returns
        -------
        elements: numpy.array
            A numpy array with elements of the molecule atoms
        """
        return np.array([self.get_element(atom) for atom in self.get_atoms()])

    def _get_atom_types(self):
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
        n_atoms = self.get_natoms()
        _mol = self.get_mol()

        feats = SmallMol.factory.GetFeaturesForMol(_mol)
        properties = np.zeros((n_atoms, 8), dtype=bool)

        for feat in feats:
            fam = feat.GetFamily()
            if fam not in atom_mapping:  # Non relevant property
                continue
            properties[feat.GetAtomIds(), atom_mapping[fam]] = 1

        # Occupancy, ignoring hydrogens.
        properties[:, 7] = self.get_elements() != 'H'
        return properties

    def _get_channel_radii(self):
        """
        Multiplies atom types by each atom vdW radius.
        """
        from htmd.molecule.vdw import radiidict
        radii = np.vectorize(radiidict.__getitem__)(self.get_elements()) * self._get_atom_types().T
        return radii.T.copy()

    def get_center(self, coords=None):
        """
        Returns geometrical center of molecule.
        """
        if coords is None:
            coords = self.get_coords()
        return coords.mean(axis=0).astype(np.float32)

    def generate_conformers(self, num_confs=400,  optimizemode='mmff', append=True):
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

        if not append:
            self.remove_conformers()

        _mol = self.get_mol()
        mol = deepcopy(_mol)
        mol = Chem.AddHs(mol)

        # generating conformations
        ids = EmbedMultipleConfs(mol, numConfs=num_confs, pruneRmsThresh=1., maxAttempts=10000)

        if optimizemode not in ['uff', 'mmff']:
            raise ValueError('Unknown optimizemode. Should be  "uff", "mmff"')
        # optimizing conformations depends on the optimizemode passed
        for id in ids:
            if optimizemode == 'mmff':
                MMFFOptimizeMolecule(mol, confId=id)
            elif optimizemode == 'uff':
                UFFOptimizeMolecule(mol, confId=id)
            conf = mol.GetConformer(id)
            _mol.AddConformer(conf, id+1)

    def get_voxels(self, center=None, size=24, resolution=1., rotation=None,
                   displacement=None, dtype=np.float32):
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
        coords = self.get_coords()
        lig_center = self.get_center(coords=coords)

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

        multisigmas = self._get_channel_radii()
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

    def get_name(self):
        """
        Returns the molecule name

        Returns
        -------
        name: str
            The molecule name

        """

        return self._mol.GetProp('_Name')

    def set_name(self, name='UNK'):
        """
        Set the property '_Name' into the rdkit Mol object

        Parameters
        ----------
        name: str
            The molecule name
        """
        self._mol.SetProp('_Name', name)

    def get_natoms(self):
        """
        Returns the number of atoms in the molecule

        Returns
        -------
        natoms: int
            The number of atoms in the molecule

        """
        _mol = self.get_mol()
        return _mol.GetNumAtoms()

    def get_nconformers(self):
        """
        Returns the number of conformers in the rdkit molecule object

        Returns
        -------
        nconfs: int
            The number of conformers in the molecule
        """
        _mol = self.get_mol()
        _nConformers = _mol.GetNumConformers()

        return _nConformers

    def get_conformers(self, ids=None):
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
        _mol = self.get_mol()
        _conformers = list(_mol.GetConformers())
        _nConformers = len(_conformers)

        if ids == None:
            return _conformers

        if  max(ids) >= _nConformers:
            raise IndexError("The ids list contains conformers ids {} that do not exist. Available conformers: {}".format(ids, _nConformers))

        return [_conformers[id] for id in ids]

    def write_conformers(self, savefolder='conformations', savename="molConf", filetype="sdf", savefolder_exist_ok=False,
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
            ids = range(self.get_nconformers())
        elif not isinstance(ids, list):
            raise ValueError("The ids argument should be a list of conformer ids")

        _mol = self.get_mol()

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
                _id = self.get_conformers()[id].GetId()
                fname = os.path.join(savefolder, '{}_{}.{}'.format(savename, _id, filetype))
                writer = chemwrite(fname)
            writer.write(_mol, confId=_id)

    def remove_conformers(self, ids=None):
        """
        Deletes the conformers passed

        Parameters
        ----------
        ids: list (default=None)
            The list of conformer id to delete. If None, all are removed except the first one
        """
        _mol = self.get_mol()

        _conformers = self.get_conformers(ids)

        if ids is None:
            _conformers = _conformers[1:]
        elif not isinstance(ids, list):
            raise TypeError("The ids argument should be list of confermer ids")

        _conformerIDs = [c.GetId() for c in _conformers]

        for id in _conformerIDs:
            _mol.RemoveConformer(id)

    def to_molecule(self, formalcharges=False, ids=None):
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

        _mol = self.get_mol()
        _nConformers = self.get_nconformers()
        if _nConformers == 0:
            raise NoConformerError("No Conformers are found in the molecule. Generate at least one confomer.")

        if ids == None:
            ids = list(range(_nConformers))

        elif not isinstance(ids, list):
            raise ValueError('The argument ids should be a list of confomer ids')

        molHtmd = None
        for n in ids:
            coords = self.get_coords(id=n)
            elements = self.get_elements()
            mol = Molecule()
            mol.empty(self.get_natoms())
            mol.resname[:] = self.get_name()[:3]
            mol.resid[:] = 1
            mol.name[:] = elements
            mol.element[:] = elements
            mol.charge[:] = self.get_charges(formal=formalcharges)
            mol.coords[:, :, 0] = coords
            mol.viewname = self.get_name()
            mol.bonds, mol.bondtype = self.get_bonds()
            if molHtmd == None:
                molHtmd = mol
            else:
                molHtmd.appendFrames(mol)
        return molHtmd

    def get_bonds(self):
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
        for bo in self._mol.GetBonds():
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

    def get_charges(self, formal=False):
        """
        Return the charges of the atoms in the molecules.

        Paramters
        ---------
        formal: bool
            Set as True for the formal charges, otherwise the partial ones are retrieved

        Returns
        -------
        charges: numpy.array
            An array with the atom charges of the molecule

        """

        charges = []
        atoms = self.get_atoms()
        for a in atoms:
            if formal:
                charges.append(a.GetFormalCharge())
            else:
                props = a.GetPropsAsDict()
                charge = props['_TriposPartialCharge'] if '_TriposPartialCharge' in props.keys() else 0.000
                charges.append(charge)
        return np.array(charges)


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
        _mol = self.get_mol()

        _labelsFunc = ['a', 'i', 'c', '*']

        if atomlabels is not None:
            labels = atomlabels.split('%')[1:]
            names = self.get_elements().tolist()
            indexes = [ self.getPropAtom(a, '_SmallMolAtomIdx') for a in self.get_atoms() ]
            charges = ['+' if c > 0 else "-" if c < 0 else "" for c in  self.get_charges(formal=True).tolist() ]
            chirals = ["" if self._chiralType(a) is None else "*" for a in self.get_atoms() ]
            values = [names, indexes, charges, chirals]
            idxs = [ _labelsFunc.index(l) for l in labels]
            labels_required = [values[i] for i in idxs]
            atomlabels = ["".join([str(i) for i in a]) for a in list(zip(*labels_required))]

        return  _depictMol(_mol, sketch, filename, ipython, optimize, optimizemode, removeHs, atomlabels, highlightAtoms)



class SmallMolStack:
    """
    Collection of objects of class SmallMol.
    """

    def __init__(self, sdf_file=None, removeHs=False, fixHs=True):  # , n_jobs=1
        from tqdm import tqdm

        self._sdffile = sdf_file if  self._isSdfFile(sdf_file) else None

        self._mols = None

        if sdf_file != None:
            self._mols = self._initializeMolObjs(sdf_file, removeHs, fixHs)

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
        mols = []
        for i, mol in enumerate(tqdm(supplier)):
            if mol is not None:
                mols.append(SmallMol(mol, removeHs=removeHs, fixHs=fixHs))
            else:
                mols.append(None)

        mols = mols

        invalid_mols = self._get_invalid_indexes(mols)

        if len(invalid_mols) != 0:
            logger.warning('The following entries could not be loaded: {}. Use clean_invalids to remove them from the '
                           'pool'.format(invalid_mols))

        return np.array(mols)

    def clean_invalids(self):
        _mols = self.get_mols()


        self._mols =  [ m for m in _mols if m is not None]

    def get_mols(self):

        return self._mols

    def get_nmols(self):
        return len(self._mols)

    def vox_fun(mol):
        return None

    def __len__(self):
        return len(self._mols)

    def __getitem__(self, item):
        return self._mols[item]

    def _get_invalid_indexes(self, mols):
        """
        Returns indexes of invalid molecules
        """
        return [i for i, mol in enumerate(mols) if mol is None]

    def __iter__(self):

        _mols = self.get_mols()
        for smallmol in _mols:
            yield smallmol

    def __str__(self):
        _mols = self.get_mols()

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

    def depict(self, sketch=False, filename=None, ipython=False, optimize=False, optimizemode='std',
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

        if legends is not None and legends not in ['names', 'items']:

            raise ValueError('The "legends" should be "names" or "items"')

        legends_list = []
        if legends == 'names':
            legends_list = [ _m.get_name() for _m in self._mols ]
        elif legends == 'items':
            legends_list = [ str(n+1) for n in range(len(self._mols))]

        _mols = [ _m.get_mol() for _m in self.get_mols() ]

        if highlightAtoms is not None:
            if len(highlightAtoms) != len(_mols):
                raise ValueError('The highlightAtoms {} should have the same length of the mols {}'.format(len(highlightAtoms), len(_mols)))


        return depictMultipleMols(_mols, sketch, filename, ipython, optimize, optimizemode,
                                removeHs, legends_list, highlightAtoms, mols_perrow)



if __name__ == '__main__':
    pass
