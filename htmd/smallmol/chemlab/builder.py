from htmd.smallmol.smallmol import SmallMol
import numpy as np
from htmd.smallmol.chemlab.periodictable import PeriodicTable
from rdkit.Chem.rdchem import HybridizationType, BondType
import math


_heavy_atoms = ['C', 'N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']

_hetero_atoms = ['N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']


_hybridizations_IdxToType = HybridizationType.values
_bondtypes_IdxToType = BondType.values

class Builder:

    def __init__(self, smallmol=None, checkInitialConformer=True):

        self.smallmol = smallmol
        self.periodictable = PeriodicTable()

        if checkInitialConformer and smallmol is not None:
            self._prepareConformer()


    def loadMol(self, smallmol, checkInitialConformer=True):

        if not isinstance(smallmol, SmallMol):
            raise ValueError('The argument tyep {} is not accepted. Only SmallMol object. '.format(smallmol))

        self.smallmol = smallmol.copy()

        if checkInitialConformer:
            self._prepareConformer()

    def _prepareConformer(self):
        from rdkit.Chem import AllChem

        sm = self.smallmol
        mol = sm.toRdkitMol(includeConformer=True)

        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        conf = mol.GetConformer(0)
        coords = conf.GetPositions()
        coords = coords[:, :, np.newaxis]
        sm.coords = coords


    def addHydrogens(self, onlyExplicit=False):

        pT = self.periodictable
        sm = self.smallmol

        if onlyExplicit:
            heavy_sel = ' '.join(_hetero_atoms)
        else:
            heavy_sel = ' '.join(_heavy_atoms)

        atomidxs = sm.get('element {}'.format(heavy_sel), 'idx')
        # testing
        formalcharges = sm.get('element {}'.format(heavy_sel), 'formalcharge')

        # testing : add formalcharges part
        missing_atoms = np.array([pT.getMissingValence(aId, sm, fch) for aId, fch in zip(atomidxs, formalcharges)])
        ids = np.where(missing_atoms > 0)
        heavyToFill = atomidxs[ids]
        numHydrogens = missing_atoms[ids]
        print("missing hydrogens: ", numHydrogens)
        h_atoms = []

        for n_h, a in zip(numHydrogens, heavyToFill):
            for i in range(n_h):
                new_atom = pT.getAtom('H', a)
                h_atoms.append(new_atom)

        self._addAtoms(h_atoms)


    def removeHydrogens(self, removePolars=True, removeNonPolars=True):

        sm = self.smallmol
        hydrogensIdx =  sm.get('element H', 'idx')

        neighborsIdx = np.concatenate(sm.neighbors[hydrogensIdx])
        neighbors = sm.element[neighborsIdx]
        hydrogensPolarIdx = hydrogensIdx[np.where( neighbors != 'C' )]
        hydrogensNonPolarIdx = hydrogensIdx[np.where(neighbors == 'C')]


        hydrogensIdx = np.array([],dtype=np.int)
        if removePolars:
            hydrogensIdx = np.append( hydrogensIdx, hydrogensPolarIdx)
        if removeNonPolars:
            hydrogensIdx = np.append(hydrogensIdx, hydrogensNonPolarIdx)

        self._removeAtoms(hydrogensIdx)

    def _addAtoms(self, atoms):

        sm = self.smallmol

        #debug
        c = 1
        for a in atoms:

            n_atom = sm.numAtoms
           # print('placing {} to {}'.format(n_atom, a['attachTo']))
            sm.__dict__['idx'] = np.append(sm.idx, n_atom)
            sm.__dict__['element'] = np.append(sm.element, a['element'])
            sm.__dict__['atomname'] = np.append(sm.atomname, '{}{}'.format(a['element'], n_atom))
            sm.__dict__['charge'] = np.append(sm.charge, a['charge'])
            sm.__dict__['formalcharge'] = np.append(sm.formalcharge, a['formalcharge'])
            sm.__dict__['chiral'] = np.append(sm.chiral, a['chiral'])
            sm.__dict__['hybridization'] = np.append(sm.hybridization, _hybridizations_IdxToType[a['hybridization']])
            sm.__dict__['coords'] = np.concatenate((sm.coords, a['coords']), axis=0)
            sm.__dict__['neighbors'] = np.array(sm.neighbors.tolist() + [[a['attachTo']]])
            sm.__dict__['bondtypes'] = np.array(sm.bondtypes.tolist() + [[a['bondtype']]])
            if isinstance(sm.neighbors[a['attachTo']], np.ndarray):
                new_neighbors = sm.neighbors.tolist()
                new_neighbors[a['attachTo']].append(n_atom)
                sm.__dict__['neighbors'] = np.array(new_neighbors)
                new_bonds = sm.bondtypes.tolist()
                new_bonds[a['attachTo']].append(a['bondtype'])
                sm.__dict__['bondtypes'] = np.array(new_bonds)
                
            else:
                sm.neighbors[a['attachTo']].append(n_atom)
                sm.bondtypes[a['attachTo']].append(a['bondtype'])



            #TODO based on chiral
            sm.__dict__['_chiraltags'] = np.append(sm._chiraltags, 0)
            #### default
            sm.__dict__['expliciths'] = np.append(sm.expliciths, 0)

            # coords

            sm.__dict__['coords'][n_atom]  = self._getAtomCoords(n_atom, a['attachTo'])

            #### REMOVE IT ( For debugging)
            #if c == 7: break
          #  print()
            c += 1

    def _getAtomCoords(self, atomIdx, attachToIdx):
        sm =self.smallmol
        pT = self.periodictable

        atom_element = sm.element[atomIdx]
        heavy_element = sm.element[attachToIdx]
        heavy_coords = sm.coords[attachToIdx].reshape(3)
        heavy_hybridization = sm.hybridization[attachToIdx]
        heavy_neighbors = sm.neighbors[attachToIdx]
        heavy_neighbors_num = len(heavy_neighbors)
        bondLength = pT.getRadiusBond(heavy_element) + pT.getRadiusBond(atom_element)

      #  print(heavy_neighbors)
      #  print(heavy_hybridization)
#
        dirVect = np.array([0,0,0])
        if heavy_neighbors_num == 1:
            dirVect[2] = 1

        elif heavy_neighbors_num == 2:
            if isinstance(heavy_hybridization, int):
                heavy_hybridization = _hybridizations_IdxToType[heavy_hybridization]
            else:
                heavy_hybridization = heavy_hybridization

            nbr = self._filterNeighbors(attachToIdx, atomIdx)[0]
            nbr_coords = sm.coords[nbr].reshape(3)
            nbrVect = nbr_coords - heavy_coords
    #       print('neighbor vector: ', nbrVect)

            nbrVect = normalizeVector(nbrVect)
            nbrVect = -nbrVect
            #print('neighbor normalized and inversed: ', nbrVect)

            perpVect = getPerpendicular(nbrVect)
            #              print('perp: ', perpVect)
            perpVect = normalizeVector(perpVect)
            #             print('perp norm: ',perpVect)

            if heavy_hybridization == HybridizationType.SP3:
               # print('sp3 procedure')
                rot_matrix = getRotationMatrix(perpVect, (180-109.471), deg=True)
                dirVect = np.dot(rot_matrix, nbrVect)

            elif heavy_hybridization == HybridizationType.SP2:
              #  print('sp2 procedure')

                nbr_neighbors = sm.neighbors[nbr]
                nbr_neighbors_num = len(nbr_neighbors)

               # print("nbr neighbors: ", nbr_neighbors_num)

                if nbr_neighbors_num > 1:

                    nbr2 = [n for n in sm.neighbors[nbr] if n != attachToIdx][0]
                    idx = sm.neighbors[nbr].index(nbr2)
                    nbr2_nbr_btype = sm.bondtypes[nbr][idx]
                    nbr2_nbr_btype = _bondtypes_IdxToType[nbr2_nbr_btype]
                    if nbr2_nbr_btype in [BondType.AROMATIC, BondType.DOUBLE]:
                       # print('taking into account another atom')
                        nbr2_coords = sm.coords[nbr2].reshape(3)
                        nbr2Vect = nbr2_coords - nbr_coords
                        nbr2Vect = normalizeVector(nbr2Vect)
                        perpVect = np.cross(nbr2Vect, nbrVect)
                        perpVect = normalizeVector(perpVect)

                rot_matrix = getRotationMatrix(perpVect, 60, deg=True)
                dirVect = np.dot(rot_matrix, nbrVect)




        elif heavy_neighbors_num == 3:
            if isinstance(heavy_hybridization, int):
                heavy_hybridization = _hybridizations_IdxToType[heavy_hybridization]
            else:
                heavy_hybridization = heavy_hybridization
            nbrs = self._filterNeighbors(attachToIdx, atomIdx)
            nbr1, nbr2 = nbrs
            nbr1_coords, nbr2_coords = sm.coords[nbr1].reshape(3), sm.coords[nbr2].reshape(3)
            nbr1Vect = heavy_coords - nbr1_coords
            nbr2Vect = heavy_coords - nbr2_coords

            nbr1Vect = normalizeVector(nbr1Vect)
            nbr2Vect = normalizeVector(nbr2Vect)

            dirVect = nbr1Vect + nbr2Vect
            dirVect = normalizeVector(dirVect)

            if heavy_hybridization == HybridizationType.SP3:
                #print('second attachment sp3 procedure')
                nbrPerp = np.cross(nbr1Vect, nbr2Vect)

                rotAxis = np.cross(nbrPerp, dirVect)
                rotAxis = normalizeVector(rotAxis)

                rot_matrix = getRotationMatrix(rotAxis, (109.471/2), deg=True)

                dirVect = np.dot(rot_matrix, dirVect)


        elif heavy_neighbors_num == 4:
          #  print('last attachment sp3 procedure')
            nbrs = self._filterNeighbors(attachToIdx, atomIdx)
            nbr1, nbr2, nbr3 = nbrs

            nbr1_coords, nbr2_coords, nbr3_coords = sm.coords[nbr1].reshape(3), sm.coords[nbr2].reshape(3), sm.coords[nbr3].reshape(3)

            nbr1Vect = heavy_coords - nbr1_coords
            nbr2Vect = heavy_coords - nbr2_coords
            nbr3Vect = heavy_coords - nbr3_coords

            nbr1Vect = normalizeVector(nbr1Vect)
            nbr2Vect = normalizeVector(nbr2Vect)
            nbr3Vect = normalizeVector(nbr3Vect)

            dirVect = nbr1Vect + nbr2Vect + nbr3Vect
            dirVect = normalizeVector(dirVect)



        atom_coords = heavy_coords + dirVect * bondLength

        return  atom_coords.reshape(3,1)

    def _filterNeighbors(self, atomIdx, exclude):

        sm = self.smallmol

        neighbors = sm.neighbors[atomIdx]
        return [ n for n in neighbors if n != exclude]

    def _removeAtoms(self, ids):
        sm = self.smallmol

        for k in sm._atom_fields:

            if k == 'coords':
                sm.__dict__[k] = np.delete(sm.__dict__[k], ids, axis=0)
            else:
                sm.__dict__[k] = np.delete(sm.__dict__[k], ids)

        new_neighbors = []
        new_bondstype = []
        for n_set, b_set in zip(sm.__dict__['neighbors'], sm.__dict__['bondtypes']):
            new_n_set = []
            new_b_set = []
            for n, b in zip(n_set, b_set):
                if n not in ids:
                    new_n_set.append(n)
                    new_b_set.append(b)
            new_neighbors.append(new_n_set)
            new_bondstype.append(new_b_set)
        sm.__dict__['neighbors'] = np.array(new_neighbors)
        sm.__dict__['bondtypes'] = np.array(new_bondstype)

        atom_mapper = { hIdx:n for n, hIdx in enumerate(sm.__dict__['idx']) }

        self._fixAtomNumber(atom_mapper)

    def _fixAtomNumber(self, atom_mapper):

        sm = self.smallmol

        for i, idx in enumerate(sm.idx):
            sm.idx[i] = atom_mapper[idx]

        for i, el in zip(sm.idx, sm.element):
            sm.atomname[i] = ''.join([el, str(i)])

        for n_set in sm.neighbors:
            for i, n in enumerate(n_set):
                n_set[i] = atom_mapper[n]

    def getSmallMol(self):
        return self.smallmol

def normalizeVector(v):
    from math import sqrt
    l = sqrt(v[0] ** 2 + v[1] ** 2. + v[2] ** 2)

    if l == 0:
        print('warning normalized vector goes to 0')
        return np.array([0,0,0])
    return v / l

def getPerpendicular(v):
    V = [0,0,0]
    X, Y, Z = v[0], v[1], v[2]
    if X != 0:
        if Y != 0:
            V[0] = Y
            V[1] = -X
        elif Z != 0:
            V[2] = -X
            V[0] = Z
        else:
            V[1] = 1
    elif Y != 0:
        if Z != 0:
            V[2] = -Y
            V[1] = Z
        else:
            V[0] = 1
    elif Z !=0:
        V[0] = 1

    return np.array(V)

def getRotationMatrix(axis, theta, deg=False):
    if deg:
        theta = math.radians(theta)

    c, s = math.cos(theta), math.sin(theta)

    t = 1 - c
    X = axis[0]
    Y = axis[1]
    Z = axis[2]

    return np.array([ [t*X*X+c, t*X*Y-s*Z, t*X*Z+s*Y],
                      [t*X*Y+s*Z, t*Y*Y+c, t*Y*Z-s*X],
                      [t*X*Z-s*Y, t*Y*Z+s*X, t*Z*Z+c] ])


def getRotationMatrix2(axis, theta, deg=False):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    if deg:
        theta = math.radians(theta)

    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
