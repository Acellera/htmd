from htmd.smallmol.smallmol import SmallMol
import numpy as np
from htmd.smallmol.chemlab.periodictable import PeriodicTable
from rdkit.Chem.rdchem import HybridizationType, BondType

_heavy_atoms = ['C', 'N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']

_hetero_atoms = ['N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']


_hybridizations_IdxToType = HybridizationType.values
_bondtypes_IdxToType = BondType.values

class Builder:

    def __init__(self, smallmol=None):

        self.smallmol = smallmol
        self.periodictable = PeriodicTable()


    def loadMol(self, smallmol):

        if not isinstance(smallmol, SmallMol):
            raise ValueError('The argument tyep {} is not accepted. Only SmallMol object. '.format(smallmol))

        self.smallmol = smallmol.copy()


    def addHydrogens(self, onlyExplicit=False):

        pT = self.periodictable
        sm = self.smallmol

        if onlyExplicit:
            heavy_sel = ' '.join(_hetero_atoms)
        else:
            heavy_sel = ' '.join(_heavy_atoms)

        atomidxs = sm.get('element {}'.format(heavy_sel), 'idx')

        missing_atoms = np.array([pT.getMissingValence(aId, sm) for aId in atomidxs])
        ids = np.where(missing_atoms > 0)
        heavyToFill = atomidxs[ids]
        numHydrogens = missing_atoms[ids]

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

        for a in atoms:
            n_atom = sm.numAtoms

            sm.__dict__['idx'] = np.append(sm.idx, n_atom)
            sm.__dict__['element'] = np.append(sm.element, a['element'])
            sm.__dict__['atomname'] = np.append(sm.atomname, '{}{}'.format(a['element'], n_atom))
            sm.__dict__['charge'] = np.append(sm.charge, a['charge'])
            sm.__dict__['formalcharge'] = np.append(sm.formalcharge, a['formalcharge'])
            sm.__dict__['chiral'] = np.append(sm.chiral, a['chiral'])
            sm.__dict__['hybridization'] = np.append(sm.hybridization, _hybridizations_IdxToType[a['hybridization']])
            sm.__dict__['coords'] = np.concatenate((sm.coords, a['coords']), axis=0)
            sm.__dict__['neighbors'] = np.array(sm.neighbors.tolist() + [[a['attachTo']]])
            sm.neighbors[a['attachTo']].append(n_atom)
            sm.__dict__['bondtypes'] = np.array(sm.bondtypes.tolist() + [[a['bondtype']]])
            sm.bondtypes[a['attachTo']].append(a['bondtype'])

            #TODO based on chiral
            sm.__dict__['_chiraltags'] = np.append(sm._chiraltags, 0)
            #### default
            sm.__dict__['expliciths'] = np.append(sm.expliciths, 0)


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
        #TODO
        # atomname --> element, idx
        sm = self.smallmol

        for i, idx in enumerate(sm.idx):
            sm.idx[i] = atom_mapper[idx]

        for n_set in sm.neighbors:
            for i, n in enumerate(n_set):
                n_set[i] = atom_mapper[n]



    def getSmallMol(self):
        return self.smallmol



