from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.chemlab.periodictable import PeriodicTable, _heavy_atoms, _hetero_atoms, _hybridizations_IdxToType
import numpy as np
from htmd.smallmol.util import _getPerpendicular, _normalizeVector, _getRotationMatrix
import rdkit

class Builder:
    """
    A class to modify a SmallMol object: adding and removing hydrogens; creation of ligands from scracth

    Parameters
    ----------
    smallmol: htmd.smallmol.smallmol.SmallMol
        The SmallMol object
    checkInitialConformer: bool
        If True, a conformation is generated if None are found. (designed espcially for smiles)

    Example
    -------
    >>> sm = SmallMol('CCC')
    >>> sm.coords[0], sm.numAtoms
    ( array([[ 0.],[ 0.],[ 0.]]), 11)
    >>> B = Builder(sm)
    >>> sm_new = B.getSmallMol()
    >>> sm_new.coords[0]
    array([[ 1.1701015 ],[-0.21571056],[-0.21615667]])
    >>> B.removeHydrogens()
    >>> B.getSmallMol().numAtoms
    3
    """

    def __init__(self, smallmol=None, checkInitialConformer=True):

        self.smallmol = self._initializeMol(smallmol) if smallmol is not None else smallmol
        self.periodictable = PeriodicTable()

        if checkInitialConformer and smallmol is not None:
            self._prepareConformer()


    def _initializeMol(self, smallmol):
        rmol = smallmol.toRdkitMol(includeConformer=True)
        rdkit.Chem.Kekulize(rmol)
        sm = SmallMol(rmol)
        return sm

    def loadMol(self, smallmol, checkInitialConformer=True):
        """
        Loads the SmallMol object. The additional argument checkInitialConformer ensure that one valid conformer exists

        Parameters
        ---------
        smallmol: htmd.smallmol.smallmol.SmallMol
            The SmallMol object
        checkInitialConformer: bool
            If True, a conformer is generated if only not valid are found in the smallmol object
        """

        if not isinstance(smallmol, SmallMol):
            raise ValueError('The argument tyep {} is not accepted. Only SmallMol object. '.format(smallmol))

        self.smallmol = smallmol.copy()

        if checkInitialConformer:
            self._prepareConformer()

    def _prepareConformer(self):
        """
        Generate a conformer and stores the coords in the SmallMol object
        """

        from rdkit.Chem import AllChem

        sm = self.smallmol
        mol = sm.toRdkitMol(includeConformer=True)

        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        conf = mol.GetConformer(0)
        coords = conf.GetPositions()
        coords = coords[:, :, np.newaxis]
        sm.coords = coords

    def addHydrogens(self, onlyExplicit=False):
        """
        Adds the hydrogens to fill up the valence of the atoms. With onlyExplicit is possible to add only the polar ones

        Parameters
        ----------
        onlyExplicit: bool
            If True, only the polar hydrogens will be added
        """

        pT = self.periodictable
        sm = self.smallmol

        # polar hydrogens are the ones bonded to hetero atoms (no carbon)
        if onlyExplicit:
            heavy_sel = ' '.join(_hetero_atoms)
        else:
            heavy_sel = ' '.join(_heavy_atoms)

        atomidxs = sm.get('idx', 'element {}'.format(heavy_sel))
        formalcharges = sm.get('formalcharge', 'element {}'.format(heavy_sel))

        missing_atoms = np.array([pT.getMissingValence(sm, aId, fch) for aId, fch in zip(atomidxs, formalcharges)])
        ids = np.where(missing_atoms > 0)
        heavyToFill = atomidxs[ids]
        numHydrogens = missing_atoms[ids]
        # print("missing hydrogens: ", numHydrogens)
        h_atoms = []

        for n_h, a in zip(numHydrogens, heavyToFill):
            for i in range(n_h):
                new_atom = pT.getAtom('H', a)
                h_atoms.append(new_atom)

        self._addAtoms(h_atoms)

    def removeHydrogens(self, removePolars=True, removeNonPolars=True):
        """
        Removes the hydrogens. It is possible to choose the type of hydrogens to remove: polar and not polar

        Parameters
        ----------
        removePolars: bool
            If True, the polar hydrogens will be removed
        removeNonPolars: bool
            If True, the non polar hydrogens will be removed
        """

        # TODO: Fix the the kekulization problem when removing hydrogen from NH in aromatic system.
        sm = self.smallmol
        hydrogensIdx = sm.get('idx', 'element H')

        neighborsIdx = np.concatenate(sm.neighbors[hydrogensIdx])
        neighbors = sm.element[neighborsIdx]
        hydrogensPolarIdx = hydrogensIdx[np.where(neighbors != 'C')]
        hydrogensNonPolarIdx = hydrogensIdx[np.where(neighbors == 'C')]

        hydrogensIdx = np.array([], dtype=int)
        if removePolars:
            hydrogensIdx = np.hstack((hydrogensIdx, hydrogensPolarIdx))
        if removeNonPolars:
            hydrogensIdx = np.hstack((hydrogensIdx, hydrogensNonPolarIdx))

        self._removeAtoms(hydrogensIdx)

    def _addAtoms(self, atoms):
        """
        Adds the atoms one by one from the list of atoms data passed

        Parameters
        ----------
        atoms: list
            List of the atom data to add
        """

        sm = self.smallmol

        # c = 1
        for a in atoms:
            n_atom = sm.numAtoms
            # print('placing {} to {}'.format(n_atom, a['attachTo']))
            sm.idx = np.append(sm.idx, n_atom)
            sm.element = np.append(sm.element, a['element'])
            sm.atomname = np.append(sm.atomname, '{}{}'.format(a['element'], n_atom))
            sm.charge = np.append(sm.charge, a['charge'])
            sm.formalcharge = np.append(sm.formalcharge, a['formalcharge'])
            sm.chiral = np.append(sm.chiral, a['chiral'])
            sm.hybridization = np.append(sm.hybridization, _hybridizations_IdxToType[a['hybridization']])
            sm.coords = np.concatenate((sm.coords, a['coords']), axis=0)
            sm.neighbors = np.array(sm.neighbors.tolist() + [[a['attachTo']]])
            sm.bondtypes = np.array(sm.bondtypes.tolist() + [[a['bondtype']]])
            if isinstance(sm.neighbors[a['attachTo']], np.ndarray):
                new_neighbors = sm.neighbors.tolist()
                new_neighbors[a['attachTo']].append(n_atom)
                sm.neighbors = np.array(new_neighbors)
                new_bonds = sm.bondtypes.tolist()
                new_bonds[a['attachTo']].append(a['bondtype'])
                sm.bondtypes = np.array(new_bonds)

            else:
                sm.neighbors[a['attachTo']].append(n_atom)
                sm.bondtypes[a['attachTo']].append(a['bondtype'])

            # TODO: based on chiral
            sm._chiraltags = np.append(sm._chiraltags, 0)
            # Development note: removed assignement explicitHs

            # coords
            sm.coords[n_atom] = self._computeAtomCoords(n_atom, a['attachTo'])

            # c += 1

    def _computeAtomCoords(self, atomIdx, attachToIdx):
        """
        Retruns the coordinates of the new atom based on the VSEPR theory.

        Parameters
        ----------
        atomIdx: int
            The index of the new atom
        attachToIdx: int
            The index of the atom you want to bond to

        Returns
        -------
        newcoords: np.array
            The np.array of the new atom coords
        """
        # TODO: different hybridization from S SP SP2 SP3: Like S, P

        sm = self.smallmol
        pT = self.periodictable

        atom_element = sm.element[atomIdx]
        heavy_element = sm.element[attachToIdx]
        heavy_coords = sm.coords[attachToIdx].reshape(3)
        heavy_hybridization = sm.hybridization[attachToIdx]
        heavy_neighbors = sm.neighbors[attachToIdx]
        heavy_neighbors_num = len(heavy_neighbors)
        bondLength = pT.getBondRadius(heavy_element) + pT.getBondRadius(atom_element)

        dirVect = np.array([0, 0, 0])
        if heavy_neighbors_num == 1:
            dirVect[2] = 1

        elif heavy_neighbors_num == 2:

            nbr = self._filterNeighbors(attachToIdx, atomIdx)[0]
            nbr_coords = sm.coords[nbr].reshape(3)
            nbrVect = nbr_coords - heavy_coords

            nbrVect = _normalizeVector(nbrVect)
            nbrVect = -nbrVect

            perpVect = _getPerpendicular(nbrVect)
            perpVect = _normalizeVector(perpVect)

            if heavy_hybridization == 4:
                # the rotation of 180 - 109.471 comes from the tetrahedral geometry of the VSEPR theory.
                # I did not directly write 70.529, just to keep track of the 109.471 value that is clear for a chemist
                rot_matrix = _getRotationMatrix(perpVect, (180-109.471), deg=True)
                dirVect = np.dot(rot_matrix, nbrVect)

            elif heavy_hybridization == 3:

                nbr_neighbors = sm.neighbors[nbr]
                nbr_neighbors_num = len(nbr_neighbors)

                if nbr_neighbors_num > 1:
                    nbr2 = [n for n in sm.neighbors[nbr] if n != attachToIdx][0]
                    idx = sm.neighbors[nbr].index(nbr2)
                    nbr2_nbr_btype = sm.bondtypes[nbr][idx]

                    if nbr2_nbr_btype in [12, 2]:
                        nbr2_coords = sm.coords[nbr2].reshape(3)
                        nbr2Vect = nbr2_coords - nbr_coords
                        nbr2Vect = _normalizeVector(nbr2Vect)
                        perpVect = np.cross(nbr2Vect, nbrVect)
                        perpVect = _normalizeVector(perpVect)

                rot_matrix = _getRotationMatrix(perpVect, 60, deg=True)
                dirVect = np.dot(rot_matrix, nbrVect)
            elif heavy_hybridization == 2:
                dirVect = nbrVect


        elif heavy_neighbors_num == 3:
            nbrs = self._filterNeighbors(attachToIdx, atomIdx)
            nbr1, nbr2 = nbrs
            nbr1_coords, nbr2_coords = sm.coords[nbr1].reshape(3), sm.coords[nbr2].reshape(3)
            nbr1Vect = heavy_coords - nbr1_coords
            nbr2Vect = heavy_coords - nbr2_coords

            nbr1Vect = _normalizeVector(nbr1Vect)
            nbr2Vect = _normalizeVector(nbr2Vect)

            dirVect = nbr1Vect + nbr2Vect
            dirVect = _normalizeVector(dirVect)

            if heavy_hybridization == 4:
                nbrPerp = np.cross(nbr1Vect, nbr2Vect)
                rotAxis = np.cross(nbrPerp, dirVect)
                rotAxis = _normalizeVector(rotAxis)

                rot_matrix = _getRotationMatrix(rotAxis, (109.471/2), deg=True)

                dirVect = np.dot(rot_matrix, dirVect)

        elif heavy_neighbors_num == 4:
            # print('last attachment sp3 procedure')
            nbrs = self._filterNeighbors(attachToIdx, atomIdx)
            nbr1, nbr2, nbr3 = nbrs

            nbr1_coords, nbr2_coords, nbr3_coords = sm.coords[nbr1].reshape(3), sm.coords[nbr2].reshape(3), \
                                                    sm.coords[nbr3].reshape(3)

            nbr1Vect = heavy_coords - nbr1_coords
            nbr2Vect = heavy_coords - nbr2_coords
            nbr3Vect = heavy_coords - nbr3_coords

            nbr1Vect = _normalizeVector(nbr1Vect)
            nbr2Vect = _normalizeVector(nbr2Vect)
            nbr3Vect = _normalizeVector(nbr3Vect)

            dirVect = nbr1Vect + nbr2Vect + nbr3Vect
            dirVect = _normalizeVector(dirVect)

        atom_coords = heavy_coords + dirVect * bondLength

        return atom_coords.reshape(3, 1)

    def _filterNeighbors(self, atomIdx, exclude):
        """
        Returns all the indexes of the neighbors atom except for the excluded one

        Parameters
        ----------
        atomIdx: int
         The index of the central atom
        exclude: int
            The index of the atom to exclude

        Returns
        --------
        neighbors: list
            A list of the neighbors

        """

        sm = self.smallmol

        neighbors = sm.neighbors[atomIdx]
        return [n for n in neighbors if n != exclude]

    def _removeAtoms(self, ids):
        """
        Remove the atoms based on the indexes provided

        Paraemters
        ----------
        ids: list
            A list of atom indexes to remove
        """

        sm = self.smallmol

        for k in sm._atom_fields:

            if k == 'coords':
                sm.__dict__[k] = np.delete(sm.__dict__[k], ids, axis=0)
            else:
                sm.__dict__[k] = np.delete(sm.__dict__[k], ids)

        new_neighbors = []
        new_bondstype = []
        for n_set, b_set in zip(sm.neighbors, sm.bondtypes):
            new_n_set = []
            new_b_set = []
            for n, b in zip(n_set, b_set):
                if n not in ids:
                    new_n_set.append(n)
                    new_b_set.append(b)
            new_neighbors.append(new_n_set)
            new_bondstype.append(new_b_set)
        sm.neighbors = np.array(new_neighbors)
        sm.bondtypes = np.array(new_bondstype)

        atom_mapper = {hIdx: n for n, hIdx in enumerate(sm.idx)}

        self._fixAtomNumber(atom_mapper)

    def _fixAtomNumber(self, atom_mapper):
        """
        Renumber the atoms after they are removed.

        Parameters
        ----------
        atom_mapper: dict
            A dictionary of the old index: new index
        """

        sm = self.smallmol

        for i, idx in enumerate(sm.idx):
            sm.idx[i] = atom_mapper[idx]

        for i, el in zip(sm.idx, sm.element):
            sm.atomname[i] = ''.join([el, str(i)])

        for n_set in sm.neighbors:
            for i, n in enumerate(n_set):
                n_set[i] = atom_mapper[n]

    def getSmallMol(self):
        """
        Returns the SmallMol object.

        Returns
        -------
        newsmallmol: htmd.smallmol.smallmol.SmallMol
            The SmallMol object
        """
        return self.smallmol
