import os
import multiprocessing
import math
import numpy as np

from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
from rdkit.Chem import ChemicalFeatures

from htmd.molecule.voxeldescriptors import _getOccupancyC, _getGridCenters
from htmd.smallmol.util import get_rotationMatrix, rotate
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

    def __init__(self, mol, ignore_errors=False, addHs=True):
        """
        Initializes small molecule object

        Parameters
        ----------
        mol: rdkit Molecule object or string
            (i) Rdkit molecule or (ii) Location of molecule file (".pdb"/".mol2") or (iii) a smile string.
        ignore_errors: bool
            If True errors will not be raised.
        """

        #  Determine how to load molecule
        # Process as Rdkit molecule
        if isinstance(mol, Chem.Mol):
            self._mol = mol
        elif mol is None and not ignore_errors:
            self._mol = mol

        # Process as string
        elif isinstance(mol, str):
            name_sufix = os.path.splitext(mol)[-1]
            if name_sufix == ".mol2":
                self._mol = Chem.MolFromMol2File(mol)
            elif name_sufix == ".pdb":
                self._mol = Chem.MolFromPDBFile(mol)

            # We assume any string is a valid smile
            # TODO: validate the strings
            else:
                self._mol = Chem.MolFromSmiles(mol)

        # Don't feed garbage!
        else:
            raise ValueError("Unkown file type: '{}'.".format(type(mol)))

        # Add hydrogens
        if addHs:
            self._mol = Chem.RemoveHs(self._mol)
            self._mol = Chem.AddHs(self._mol, addCoords=True)

    def get_coords(self):
        """
        Returns molecule coordinates.
        """
        n_atoms = self._mol.GetNumAtoms()
        conformer = self._mol.GetConformer()
        coords = [[corobj.x, corobj.y, corobj.z] for corobj in [conformer.GetAtomPosition(i) for i in range(n_atoms)]]
        return np.array(coords, dtype=np.float32)

    def get_elements(self):
        """
        Returns molecule elements.
        """
        return np.array([atom.GetSymbol() for atom in self._mol.GetAtoms()])

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
        n_atoms = self._mol.GetNumAtoms()

        feats = SmallMol.factory.GetFeaturesForMol(self._mol)
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

    def generate_conformers(self, savefolder, savename="molecule_conformers", filetype="pdb",
                            savefolder_exist_ok=False, num_confs=400):
        """
        Generates ligand conformer and saves the results to a folder.


        Parameters
        ----------
        savefolder: str
            Path to directory where the results will be saved
        savename: str
           Name of the generated files. example filename: <savename>_1.pdb
        filetype: str
           must be 'pdb' or 'mol2'
        savefolder_exist_ok: bool
           if false returns an error if savefolder already exsits
        Nconformers: int
           Number of conforer to generate.

        """
        from rdkit.Chem import AllChem
        os.makedirs(savefolder, exist_ok=savefolder_exist_ok)

        mol = deepcopy(self._mol)
        mol = Chem.AddHs(mol)
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, pruneRmsThresh=1., maxAttempts=10000)
        for id in ids:
            AllChem.UFFOptimizeMolecule(mol, confId=id)
        for index, id in enumerate(ids):
            if filetype == "pdb":
                chemwrite = Chem.PDBWriter
            elif filetype == "sdf":
                chemwrite = Chem.SDWriter
            else:
                raise ValueError("Unknown file format. Cannot save to format '{}'".format(filetype))
            writer = chemwrite(os.path.join(savefolder, '{}_{}.{}'.format(savename, index + 1, filetype)))
            writer.write(mol, confId=id)


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
        return self._mol.GetProp('_Name')

    def get_natoms(self):
        return self._mol.GetNumAtoms()

    def to_molecule(self):
        from htmd.molecule.molecule import Molecule
        coords = self.get_coords()
        elements = self.get_elements()
        mol = Molecule()
        mol.empty(self.get_natoms())
        mol.resname[:] = self.get_name()[:3]
        mol.resid[:] = 1
        mol.name[:] = elements
        mol.element[:] = elements
        mol.charge[:] = self.get_charges()
        mol.coords[:, :, 0] = coords
        mol.viewname = self.get_name()
        mol.bonds, mol.bondtype = self.get_bonds()
        return mol

    def get_bonds(self):
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

    def get_charges(self):
        charges = []
        for a in self._mol.GetAtoms():
            charges.append(a.GetFormalCharge())
        return np.array(charges)


class SmallMolStack:
    """
    Collection of objects of class SmallMol.
    """

    def __init__(self, sdf_file, removeHs=True, addHs=True):  # , n_jobs=1
        from tqdm import tqdm
        if addHs and not removeHs:
            raise AttributeError('To add hydrogens with the addHs option please also enable the removeHs option.')

        supplier = Chem.SDMolSupplier(sdf_file, removeHs=removeHs)
        self.filepath = sdf_file
        mm = []
        for x in supplier:
            if x is not None:
                mm.append(SmallMol(x, addHs=addHs))
            else:
                mm.append(None)
        self._mols = np.array(mm)
        self.n_invalid = len(self.get_invalid_indexes())
        if self.n_invalid > 0:
            print('We detected {} errors when reading entries in the sdf file. Please'\
                  ' run SmallMol.get_invalid_indexes() and remove them accordingly from'\
                  ' SmallMol._mols as they can not be featurized.'.format(self.n_invalid))

    def vox_fun(mol):
        return None

    def __len__(self):
        return len(self._mols)

    def __getitem__(self, item):
        return self._mols[item]

    def get_invalid_indexes(self):
        """
        Returns indexes of invalid molecules
        """
        return [i for i, mol in enumerate(self._mols) if mol is None]

    def remove_invalid_indexes(self):
        self._mols = [m for m in self._mols if m is not None]

    def __iter__(self):
        for smallmol in self._mols:
            yield smallmol

    def __str__(self):
        return ('Stack of Small molecules.'
                '\n\tContains {} Molecules.'
                '\n\tSource file: "{}".').format(len(self._mols), self.filepath)

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


if __name__ == '__main__':
    pass
