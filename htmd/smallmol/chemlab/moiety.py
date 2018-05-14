from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.chemlab.periodictable import PeriodicTable
from htmd.smallmol.chemlab.builder import Builder
from htmd.smallmol.util import _ensurenestedlists, _flatnestedlists
import numpy as np
import rdkit
import os
from glob import glob
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
        # 2.0 halogens
        # 2.1 mark  heteroatoms
        # 2.2 mark atoms bonded to heteroatoms with = or # [NOT AROMATIC]
        # 2.3 carbon in C=C or C#C substructure [NOT AROMATIC]
        # 2.4 carbon in acetal. C sp3 attached with single bond to at least two N,O,S

    # 2.0 - halogens
        halogens = [ [a] for a in atoms_lasting if pT.isHalogenAtom(smallmol.element[a])]
        mois = [ Moiety(smallmol, a) for a in halogens ]

        moietiesIdentified.extend(mois)
        atoms_placed = sorted(_flatnestedlists(halogens))
        atoms_lasting = set(sorted(atoms_lasting - set(atoms_placed)))

    # 2.1 - mark heteroatoms # 2.2 - mark atoms = or # to hetero [NOT AROMATIC] # 2.4 acetal likes moieties
        heteroAtoms = [a for a in atoms_lasting if pT.isHeteroAtom(smallmol.element[a])]
        #
        mois, atomsMarked = self._getMoietiesHeteroConnected(heteroAtoms)
        moietiesIdentified.extend(mois)

        atoms_placed = sorted(atomsMarked)
        atoms_lasting = set(sorted(atoms_lasting - set(atoms_placed)))

    # # 2.3 - carbons in C=C or C#C substructure [NOT AROMATIC]

        carbonAtoms = [a for a in atoms_lasting if smallmol.element[a] == 'C']

        mois, atomsMarked = self._getMoietiesSatureCarbons(carbonAtoms)
        moietiesIdentified.extend(mois)

        atoms_placed = sorted(atomsMarked)
        atoms_lasting = set(sorted(atoms_lasting - set(atoms_placed)))

        print(atoms_lasting)


    def _getMoietiesSatureCarbons(self, carbonAtons):

        smallmol = self.smallmol

        carbons_pairs = list(combinations(carbonAtons, 2))
        carbonsSP2 = [[a1, a2] for a1, a2 in carbons_pairs if
                      smallmol.foundBondBetween('idx {}'.format(a1), 'idx {}'.format(a2), bondtype=2)]
        carbonsSP3 = [[a1, a2] for a1, a2 in carbons_pairs if
                      smallmol.foundBondBetween('idx {}'.format(a1), 'idx {}'.format(a2), bondtype=3)]
        carbonsMarked = carbonsSP2 + carbonsSP3

        mois = []
        for atoms in carbonsMarked:
            moi = Moiety(smallmol, atoms)
            mois.append(moi)

        return mois, _flatnestedlists(carbonsMarked)

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

        # 2.3 hetero connected to C= or C#
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
        otherMois = []
        for atoms in atomsConnected_cleaned:
            moi = Moiety(smallmol, atoms)
            if len(atoms) == 1:
                otherMois.append(moi)
            else:
                mois.append(moi)

        # 2.4 acetals

        otherMois, atomsMarked = self._getMoietiesAcetal(otherMois)

        mois.extend(otherMois)
        atomsConnected_cleaned.extend(atomsMarked)

        return mois, _flatnestedlists(atomsConnected_cleaned)

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

    def _getMoietiesAcetal(self, mois):

        lmois = copy.deepcopy(mois)

        mois_cleaned = []
        atomsMarked = []
        while len(lmois) != 0:
            moi = lmois.pop(0)
            moiLinks = moi.links
            for lmoi in lmois:
                lmoiLinks = lmoi.links
                common = set(moiLinks) & set(lmoiLinks)
                if len(common) != 0:
                    atomsMarked.append(common)
                    moi.mergeMoiety(lmoi)
                    lmois.remove(lmoi)
                    break
            mois_cleaned.append(moi)

        return  mois_cleaned, atomsMarked


    def depict(self, showLinks=True, useDummies=True, showLabels=True, molspercol=3):

        from tempfile import NamedTemporaryFile
        import matplotlib.pyplot as plt
        from PIL import Image

        tmpdir = NamedTemporaryFile().name
        os.mkdir(tmpdir)

        mois = self.moieties

        for n, moi in enumerate(mois):
            im = moi.depict(showLinks, useDummies, showLabels)
            fname = os.path.join(tmpdir, '%03d' % n)
            f = open(fname + '.svg', 'w')
            f.write(im.data)
            f.close()
            os.system('convert {}.svg {}.png'.format(fname, fname))
        images = [Image.open(im) for im in sorted(glob(tmpdir + '/*.png'))]

        columns = molspercol
        fig = plt.figure(figsize=(40, 40))

        plt.axis('off')
        for i, image in enumerate(images):
            ax = plt.subplot(len(images) / columns + 1, columns, i + 1)
            ax.axis('off')
            ax.set_facecolor('white')
            plt.imshow(image)

class Moiety:

    def __init__(self, parentsmallmol, fragmentatoms=None):

        if fragmentatoms is None:
            logging.warning("Moiety object instanciated without atoms")

        elif not isinstance(fragmentatoms, list):
            raise ValueError("The atoms need to be passed as a list")

        self.parentsmallmol = parentsmallmol.copy()
        self.fragments = _ensurenestedlists(fragmentatoms)
        self.atoms = np.unique(sorted(_flatnestedlists(self.fragments)))
        self.links = []

        self.moismallmol = self._prepareSmallMol()

    def _prepareSmallMol(self):

        # Qui si deve sistemare gli azoto aromatici.
        # 1 idrogeni legati ad essi
        # 2 attenzione che devi iterare sugli atomi del moiety e non su tutta la molecola, se no potresti prenderne
        # da altri moities.
        # 3 ricorda di rimappare sull'originale così da poter trovare gli idrogeni corretti
        # Infine, aggiungi gli idrogeni non con la funzione addHydrogens ma con addatom in modo specifico
        parentsmallmol = self.parentsmallmol
        fragments = self.fragments
        atoms = self.atoms

        atoms_string = " ".join([str(a) for a in atoms])
        atomsToRemove = parentsmallmol.get('idx', 'idx {}'.format(atoms_string), invert=True)
        b = Builder(parentsmallmol)
        b._removeAtoms(atomsToRemove)
        nitrogenAromatic = [ a for a in atoms if parentsmallmol.element[a] == 'N' and parentsmallmol.isaromatic[a] ]
        #nhydrogens = [  for n in nitrogenAromatic parentsmallmol.neighbors[n]  ]
        print(nhydrogens)
        if len(nitrogenAromatic) != 0:
            b.addHydrogens(onlyExplicit=True)

        sm = b.getSmallMol()

        self._setBreakPoints(parentsmallmol, sm)

        return sm

    def _setBreakPoints(self, parentsmallmol, sm):

        atoms = self.atoms

        neighbors = []
        for a in atoms:
            ns = [ na for na in parentsmallmol.neighbors[a] if parentsmallmol.element[na] != 'H' and na not in atoms]
            neighbors.append(ns)
        sm.__dict__['links'] = np.array(neighbors)

        self.links = _flatnestedlists(neighbors)

    def depict(self, showLinks=True, useDummies=True, showLabels=True):
        from htmd.smallmol.util import _depictMol

        pT = PeriodicTable()
        parentsmallmol = self.parentsmallmol
        sm = copy.deepcopy(self.moismallmol)
        atoms = self.atoms
        print(sm.element)
        print(sm.idx)

        elements = parentsmallmol.element[atoms]
        print(elements)
        indexes = parentsmallmol.idx[atoms]
        formalcharges = parentsmallmol.formalcharge[atoms]
        if showLinks:
            b = Builder(sm, checkInitialConformer=False)
            for n, links in enumerate(sm.links):
                if len(links)  != 0:
                    #el = ('Du', '' )if useDummies else (parentsmallmol.element[l[0]], parentsmallmol.element[l[0]])
                    for l in links:
                        el = ('Du', '') if useDummies else (parentsmallmol.element[l], parentsmallmol.element[l])
                        indexes = np.append(indexes, parentsmallmol.idx[l])
                        elements = np.append(elements, el[1])
                        formalcharges = np.append(formalcharges, 0)
                        A = pT.getAtom(el[0], n)
                        b._addAtoms([A])
            sm = b.getSmallMol()
        _mol = sm.toRdkitMol(includeConformer=True)

        rdkit.Chem.AllChem.Compute2DCoords(_mol)

        if showLabels:
            formalcharges = ['' if c == 0 else "+" if c == 1 else "-" for c in formalcharges]
            values = [elements, indexes, formalcharges]
            atomlabels = ["".join([str(i) for i in a]) for a in list(zip(*values))]

        rdkit.Chem.Kekulize(_mol)
        svg = _depictMol(_mol, filename=None, ipython=True, atomlabels=atomlabels,
                          highlightAtoms=None)

        return svg

    def mergeMoiety(self, moi):

        atoms = self.atoms
        links = self.links
        parentsmallmol = self.parentsmallmol

        common_links = list(set(moi.links) & set(links))
        newAtoms = _flatnestedlists([moi.atoms, self.atoms, common_links])

        self.atoms = sorted(newAtoms)
        self.fragments = _ensurenestedlists(self.atoms)
        
        self.moismallmol = self._prepareSmallMol()