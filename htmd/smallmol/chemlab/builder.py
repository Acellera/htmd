from htmd.smallmol.smallmol import SmallMol
import numpy as np

class Builder:

    def __init__(self, smallmol=None):

        self.smallmol = smallmol


    def loadMol(self, smallmol):

        if not isinstance(smallmol, SmallMol):
            raise ValueError('The argument tyep {} is not accepted. Only SmallMol object. '.format(smallmol))

        self.smallmol = smallmol.copy()


    def removeHydrogens(self, removePolars=True, removeNonPolars=True):

        sm = self.smallmol
        hydrogensIdx =  sm.get('element H', 'idx')

        neighborsIdx = np.concatenate(sm.neighbors[hydrogensIdx])
        neighbors = sm.element[neighborsIdx]
        hydrogensPolarIdx = hydrogensIdx[np.where( neighbors != 'C' )]
        hydrogensNonPolarIdx = hydrogensIdx[np.where(neighbors == 'C')]


        hydrogensIdx = np.array([],dtype=np.int)
        if removeNonPolars:
            hydrogensIdx = np.append( hydrogensIdx, hydrogensPolarIdx)
        if removeNonPolars:
            hydrogensIdx = np.append(hydrogensIdx, hydrogensNonPolarIdx)

        self._removeAtoms(hydrogensIdx)


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


        sm.__dict__['neighbors'] = new_neighbors
        sm.__dict__['bondtypes'] = new_bondstype

        atom_mapper = { hIdx:n for n, hIdx in enumerate(sm.__dict__['idx']) }

        self._fixAtomNumber(atom_mapper)

    def _fixAtomNumber(self, atom_mapper):
        sm = self.smallmol

        for i, idx in enumerate(sm.idx):
            sm.idx[i] = atom_mapper[idx]

        for n_set in sm.neighbors:
            for i, n in enumerate(n_set):
                n_set[i] = atom_mapper[n]



    def getSmallMol(self):
        return self.smallmol



