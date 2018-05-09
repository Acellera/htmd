from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.chemlab.periodictable import PeriodicTable
from htmd.smallmol.chemlab.builder import Builder
from htmd.smallmol.util import _ensurenestedlists, _flatnestedlists
import numpy as np
import rdkit
import copy
from itertools import combinations
import logging
logger = logging.getLogger(__name__)


class MoietyRecognition:


    def __init__(self, smallmol):

        if not isinstance(smallmol, SmallMol):
            raise ValueError('Not a SmallMol object. You should provide a valid SmallMol object')

        self.periodicTable = PeriodicTable()

        self.smallmol = smallmol
        self.moieties = []

    def run(self):

    # moieties can be of four general types:
    # 1: rings. This are moieties on their own. The only procedure is checking fused rings
    # 2: moieties with hetero atoms.
    # 3: with  halogen atoms
    # 4: without any hetero atoms


        pT = self.periodicTable
        smallmol = self.smallmol

        moietiesIdentified = self.moieties
        atoms_lasting = set(smallmol.idx)

    # 1 - rings
        mois, atomsMarked = self._getMoietiesRing()
        moietiesIdentified.extend(mois)

        atoms_placed = sorted(atomsMarked)
        atoms_lasting = set(sorted(atoms_lasting - set(atoms_placed)))

    # 2 - moieities with heteroatoms
        # 2.1 mark  heteroatoms
        # 2.2 mark atoms bonded to heteroatoms with = or # [NOT AROMATIC]
        # 2.3 carbon in C=C or C#C substructure [NOT AROMATIC]
        # 2.4 carbon in acetal. C sp3 attached with single bond to at least two N,O,S

    # 2.1 - mark heteroatoms # 2.2 - mark atoms = or # to hetero [NOT AROMATIC]

        heteroAtoms = [a for a in atoms_lasting if pT.isHeteroAtom(smallmol.element[a])]
        #
        mois, atomsMarked = self._getMoietiesHeteroConnected(heteroAtoms)
        moietiesIdentified.extend(mois)

        atoms_placed = sorted(atomsMarked)
        atoms_lasting = set(sorted(atoms_lasting - set(atoms_placed)))

    # # 2.3 - carbons in C=C or C#C substructure [NOT AROMATIC]

        #mois, atomsMarked =
    #     carbonsMarked = self._marksCarbons(atoms_lasting)
    #
    #     atoms_lasting = sorted(atoms_lasting - set(carbonsMarked))
    #
    # # 2.4 - carbon in acetal-like
    #     print(atoms_lasting)

    def _getMoietiesRing(self):
        smallmol = self.smallmol

        ring_atoms = self.getRingsAtoms(smallmol)

        mois = []
        for atoms in ring_atoms:
            moi = Moiety(smallmol, atoms)
            mois.append(moi)

        return mois, _flatnestedlists(ring_atoms)

    def _getMoietiesHeteroConnected(self, heteroatoms):

        smallmol = self.smallmol

        atomsConnected = []
        for nh, heteroA in enumerate(heteroatoms):
            atomsConnected.append([heteroA])

            for n, btype in enumerate(smallmol.bondtypes[heteroA]):
                #
                if btype >= 2:
                    atomsConnected[nh].append(smallmol.neighbors[heteroA][n])

        atomsConnected_cleaned = []
        while len(atomsConnected) != 0:
            i = atomsConnected.pop(0)
            atomsConnected_cleaned.append(i)
            merged = []
            for i2 in atomsConnected:
                if set(i) >= set(i2):
                    merged.append(i2)
            for i2 in merged:
                atomsConnected.remove(i2)

        mois = []
        for atoms in atomsConnected_cleaned:
            moi = Moiety(smallmol, atoms)
            mois.append(moi)

        return mois, _flatnestedlists(atomsConnected_cleaned)

    def _marksCarbons(self, atoms):

        smallmol = self.smallmol

        atoms_lasting = atoms

        carbons_lasting = [a for a in atoms_lasting if smallmol.element[a] == 'C']


        carbons_pairs = list(combinations(carbons_lasting, 2))

        carbonsSP2 = [[a1, a2] for a1, a2 in carbons_pairs if
                      smallmol.foundBondBetween('idx {}'.format(a1), 'idx {}'.format(a2), bondtype=2)]

        if len(carbonsSP2) != 0:
            carbonsSP2 = np.concatenate(carbonsSP2).tolist()

        carbonsSP3 = [[a1, a2] for a1, a2 in carbons_pairs if
                      smallmol.foundBondBetween('idx {}'.format(a1), 'idx {}'.format(a2), bondtype=3)]

        if len(carbonsSP3) != 0:
            carbonsSP3 = np.concatenate(carbonsSP3).tolist()

        carbonsMarked = carbonsSP2 + carbonsSP3

        return carbonsMarked


    def getRingsAtoms(self, sm):

        rmol = sm.toRdkitMol()
        _rings = rdkit.Chem.GetSymmSSSR(rmol)

        rings_atoms_tmp = [ list(r) for r in _rings ]

        atoms_list = copy.deepcopy(rings_atoms_tmp)

        rings = self._mergeFusedRings(atoms_list)

        rings_atoms = [ [rings_atoms_tmp[ri] for ri in nr] for nr in rings.values() ]

        return rings_atoms

    def _mergeFusedRings(self, atoms):

        rings = {}
        rings_assigned = []

        queue = list(atoms)

        n_ring = 0
        while len(queue) != 0:
            r = queue.pop(0)
            if n_ring not in rings_assigned:
                rings[n_ring] = [n_ring]
                rings_assigned.append(n_ring)
            connected = True

            while connected:
                connected = False
                for r2 in atoms:
                    r2_idx = atoms.index(r2)
                    if len(set(r) & set(r2) ) >= 2 and r2_idx not in rings_assigned:
                        connected = True
                        rings_assigned.append(r2_idx)
                        rings[n_ring].append(r2_idx)
                        r.extend(r2)
            n_ring += 1

        return rings

class Moiety:

    def __init__(self, parentsmallmol, atoms=None):


        if atoms is None:
            logging.warning("Moiety object instanciated without atoms")

        elif not isinstance(atoms, list):
            raise ValueError("The atoms need to be passed as a list")

        self.parentsmallmol = parentsmallmol.copy()
        self.atoms = _ensurenestedlists(atoms)
        self.atomsAttached = self._getBreakPoints()

        self.smallmol = self._createSmallMol()

    def _getBreakPoints(self):
        parentsmallmol = self.parentsmallmol
        atoms = self.atoms

        atoms_string = " ".join(_flatnestedlists(atoms, str))
        neighbors = parentsmallmol.get('neighbors', 'idx {}'.format(atoms_string))

        neighbors = [ n for a, ns in zip(atoms, neighbors) for n in ns
                      if parentsmallmol.element[n] != 'H']

        neighbors = [a for a in neighbors if a not in _flatnestedlists(atoms) ]

        return neighbors


    def _createSmallMol(self):

        sm = self.parentsmallmol
        atoms = self.atoms

        atoms_string = " ".join( _flatnestedlists(atoms, str) )
        b = Builder(sm)
        atomsToRemove = sm.get('idx', 'idx {}'.format(atoms_string), invert=True)
        b._removeAtoms(atomsToRemove)
        b.addHydrogens(onlyExplicit=True)

        return b.getSmallMol()
