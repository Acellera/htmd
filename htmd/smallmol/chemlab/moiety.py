from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.chemlab.periodictable import PeriodicTable, _halogen_atoms
from htmd.smallmol.chemlab.builder import Builder
from htmd.smallmol.util import _ensurenestedlists, _flatnestedlists
import numpy as np
import rdkit
import os
from glob import glob
import copy
from itertools import combinations
import logging
import yaml
from htmd.home import home
import os

logger = logging.getLogger(__name__)

yTypes = os.path.join(home(), 'smallmol', 'chemlab', 'moietiesFiles', 'moietiesTypes.yaml')

moiTypes = yaml.load(open(yTypes, 'r'))


class MoietyRecognition:


    def __init__(self, smallmol):

        if not isinstance(smallmol, SmallMol):
            raise ValueError('Not a SmallMol object. You should provide a valid SmallMol object')

        self.periodicTable = PeriodicTable()

        self.smallmol = smallmol
        self.moieties = []
        self.moietiesType = []

    def run(self, classify=True, mode='moieities'):

    # moieties can be of four general types:
    # 1: rings. This are moieties on their own. The only procedure is checking fused rings
    # 2: moieties with hetero atoms.
    # 3: with  halogen atoms
    # 4: without any hetero atoms

    # modes: blocks, moieties


        pT = self.periodicTable
        smallmol = self.smallmol

        moietiesIdentified = []
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

        #print(atoms_lasting)
    # merge all connected moieties (important for amide, esther and so on)

        self.moieties = self._mergeConnectedMoieties(moietiesIdentified, mode)

        if classify:
            self._classifyMoieities()

    def _classifyMoieities(self):

        for moi in self.moieties:
            mtype, morder = moi._getMoietyType()
            print(mtype, morder, moi.atoms)

    def _mergeConnectedMoieties(self, moieties, mode):

        # mode blocks or moieties

        confirmedMoieties = []
        tocheckMoieties = []

        for moi in moieties:
            if moi.isring and mode == 'blocks':
                confirmedMoieties.append(moi)
            else:
                tocheckMoieties.append(moi)

        moiN = 0
        while len(tocheckMoieties) != 0:
            moi = tocheckMoieties.pop(0)
            #print("moi1 ", moi.atoms)
            merged = False

            for moi2 in tocheckMoieties:
                if moi2.isring:
                    continue
                #print('moi2 ', moi2.atoms)
                matched1to2 = [alink for alink in moi.links if alink in moi2.atoms]
                matched2to1 = [alink for alink in moi2.links if alink in moi.atoms]
                #print(matched1to2, matched2to1)

                if len(matched1to2) != 0 and len(matched2to1) != 0:
                    #print(moi2.atoms, moi.mergedTo, moi2.mergedTo)
                    if moi.mergedTo is not None:
                        moi.mergedTo.mergeMoiety(moi2)
                        moi2.mergedTo = moi.mergedTo
                    elif moi2.mergedTo is not None:
                        moi2.mergedTo.mergeMoiety(moi)
                        moi.mergedTo = moi2.mergedTo

                    else:
                        moi.mergeMoiety(moi2)
                        moi2.mergedTo = moi
                    #tocheckMoieties.remove(moi2)
            if moi.mergedTo is None:
                confirmedMoieties.append(moi)
            moiN+= 1

        return confirmedMoieties

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
            moi = Moiety(smallmol, atoms, isring=True)
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
                #if set(i) >= set(i2):
                if len(set(i) & set(i2)) != 0:
                    atms = np.unique(np.concatenate((atomsConnected_cleaned[-1], i2)))
                    atomsConnected_cleaned[-1] = atms.tolist()
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


    def depict(self, ipython=False, filename=None, showLinks=True, showConnectorAs='dummies', showLabels=True, molspercol=3):

        from tempfile import NamedTemporaryFile
        from PIL import Image

        tmpdir = NamedTemporaryFile().name
        os.mkdir(tmpdir)

        mois = self.moieties

        for n, moi in enumerate(mois):
            im = moi.depict(True, showLinks=showLinks, showConnectorAs=showConnectorAs, showLabels=showLabels)
            fname = os.path.join(tmpdir, '%03d' % n)
            f = open(fname + '.svg', 'w')
            f.write(im.data)
            f.close()
            os.system('convert {}.svg {}.png'.format(fname, fname))

        images = [im for im in sorted(glob(tmpdir + '/*.png'))]
        im_wsize = int(1000 / molspercol)

        new_im = Image.new('RGB', (1000, 1000), 'white')

        col = 0
        row = 0
        for i, elem in enumerate(images):
            im = Image.open(elem)
            wpercent = im_wsize/float(im.size[0])
            im_hsize = int((im.size[1] * wpercent))
            im.thumbnail((im_wsize, im_hsize), Image.ANTIALIAS)
            new_im.paste(im, (im_wsize * col, im_hsize * row))
            col += 1
            if col == molspercol:
                row += 1
                col = 0

        if filename is not None:
            fname, extension = os.path.splitext(filename)
            if extension != '.png':
                filename = fname + '.png'
            new_im.save(filename)

        if ipython:
            return new_im

class Moiety:

    def __init__(self, parentsmallmol, fragmentatoms=None, isring=False):

        if fragmentatoms is None:
            logging.warning("Moiety object instanciated without atoms")

        elif not isinstance(fragmentatoms, list):
            raise ValueError("The atoms need to be passed as a list")

        self.parentsmallmol = parentsmallmol.copy()
        self.fragments = _ensurenestedlists(fragmentatoms)
        self.atoms = np.unique(sorted(_flatnestedlists(self.fragments)))
        self.links = self._setBreakPoints()
        self.isring = isring
        self.mergedTo = None
        #self.moismallmol = self._prepareSmallMol()
        #self._setBreakPoints(self.parentsmallmol, self.moismallmol)

    def getMoiSmallmol(self, includeLinks=False):

        parentsmallmol = self.parentsmallmol
        fragments = self.fragments
        atoms = self.atoms if not includeLinks else self.atoms.tolist() + self.links

        atoms_string = " ".join([str(a) for a in atoms])
        atomsToRemove = parentsmallmol.get('idx', 'idx {}'.format(atoms_string), invert=True)
        b = Builder(parentsmallmol, checkInitialConformer=False)
        b._removeAtoms(atomsToRemove)

        if self.isring:
            self._fixAromaticNitrogen(b)

        sm = b.getSmallMol()

        return sm

    def _fixAromaticNitrogen(self, builder):
        pT = PeriodicTable()

        parentsmallmol = self.parentsmallmol
        atoms = self.atoms

        nitrogenAromatic = [a for a in atoms if parentsmallmol.element[a] == 'N' and parentsmallmol.isaromatic[a]]
        nonAromaticBonded = [n for n in nitrogenAromatic if len(parentsmallmol.neighbors[n]) != sum(parentsmallmol.isaromatic[parentsmallmol.neighbors[n]])]

        if len(nonAromaticBonded) == 0:
            return

        nitrogensWithHydrogen = np.where(self.atoms == nonAromaticBonded[0])[0]

        for n in nitrogensWithHydrogen:
            atom = [na for na in parentsmallmol.neighbors[atoms[n]] if na not in atoms][0]
            element = parentsmallmol.element[atom]
            self.atoms = np.append(atoms, atom)
            a = pT.getAtom(element, n)
            builder._addAtoms([a])

    def _setBreakPoints(self):

        atoms = self.atoms
        parentsmallmol = self.parentsmallmol

        neighbors = []
        for a in atoms:
            ns = [ na for na in parentsmallmol.neighbors[a] if parentsmallmol.element[na] != 'H' and na not in atoms]
            neighbors.append(ns)

        return _flatnestedlists(neighbors)

    def depict(self, ipython=False, filename=None, showLinks=True, showConnectorAs='dummies',  showLabels=True):
        # showConnectorAs dummies, atoms, groups

        from htmd.smallmol.util import _depictMol

        showconnectChoices = ['dummies', 'atoms', 'groups']
        if showConnectorAs not in showconnectChoices:
            raise ValueError('The showConnectAs argument {} is not a valid one. Should be {}'.format(showConnectorAs, showconnectChoices))

        pT = PeriodicTable()
        parentsmallmol = self.parentsmallmol
        # sm = self.getMoiSmallmol()
        atoms = self.atoms.tolist()
        links = self.links
        elements = parentsmallmol.element[atoms] if not showLinks else parentsmallmol.element[atoms + links]
        indexes = parentsmallmol.idx[atoms] if not showLinks else parentsmallmol.element[atoms + links]
        formalcharges = parentsmallmol.formalcharge[atoms] if not showLinks else parentsmallmol.element[atoms + links]
        if showLinks:
            sm = self.getMoiSmallmol(includeLinks=True)

        _mol = sm.toRdkitMol(includeConformer=True)

        rdkit.Chem.AllChem.Compute2DCoords(_mol)

        if showLabels:
            formalcharges = ['' if c == 0 else "+" if c == 1 else "-" for c in formalcharges]
            values = [elements, indexes, formalcharges]
            atomlabels = ["".join([str(i) for i in a]) for a in list(zip(*values))]

        rdkit.Chem.Kekulize(_mol)
        svg = _depictMol(_mol, filename=filename, ipython=ipython, atomlabels=atomlabels,
                          highlightAtoms=None)

        return svg

    def _getAtomLinkFormat(self, format, linkAtom):

        parentsmallmol = self.parentsmallmol

        atomName_Label = [None, None]

        if format == 'dummies':
            atomName_Label = ('Du', '')

        elif format == 'atoms':

            atomName_Label = (parentsmallmol.element[linkAtom], parentsmallmol.element[linkAtom])

        elif format == 'groups':
            atomName_Label = ('Du', '') if not parentsmallmol.isaromatic[linkAtom] else (parentsmallmol.element[linkAtom], 'Ar')

        return atomName_Label

    def mergeMoiety(self, moi):

        atoms = self.atoms
        links = self.links

        common_links = list(set(moi.links) & set(links))
        newAtoms = _flatnestedlists([moi.atoms, atoms, common_links])


        self.atoms = np.unique(sorted(newAtoms))
        self.fragments = _ensurenestedlists(self.atoms)
        self.links = self._setBreakPoints()

    def _getMoietyOrder(self, moismallmol, moiType):

        order = None
        element = None
        if moiType in ['amine', 'amide']:
            #print('itself')
            element = 'N'

        elif moiType in ['alchol', 'thiol'] :
            #print('first atom link')
            element = 'C'
        else:
            return order

        idx = np.where( moismallmol.element == element)[0]
        nbrs = moismallmol.element[moismallmol.neighbors[idx][0]]
        order = sum([ 1 for nbr in nbrs if nbr != 'H' ])
        return  order

    def _getMoietyType(self):

        atoms = self.atoms.tolist()
        maintype = 'ring' if self.isring else 'notring'
        n_atoms = len(atoms)
        elements = self.parentsmallmol.element[atoms]

        sm = self.getMoiSmallmol(includeLinks=True)
        rmol = sm.toRdkitMol()
        #sma = rdkit.Chem.MolToSmarts(rmol)
        moiType = None
        moiOrder = None
        try:
            elements = [ el for el in elements.tolist() if el  not in _halogen_atoms ]
            elementsString= "".join(sorted(elements))
            #print("atoms: ", elementsString)
            #print(moiTypes[maintype])
            list_moieties = moiTypes[maintype]['{}atoms'.format(n_atoms)][elementsString]
            #print("moieties ", list_moieties)
            moiType = None

            for mt, s in list_moieties.items():
                match = rmol.HasSubstructMatch(rdkit.Chem.MolFromSmarts(s))
                if match:
                    moiType = mt
                    break
            moiOrder = self._getMoietyOrder(sm, moiType)
            #print(sma)
            #print(moiType)
        except:
            pass
            #print('Unknown')

        return moiType, moiOrder

#
#
# from htmd.smallmol.smallmol import SmallMolLib
#
# lib = SmallMolLib('/shared/alberto/Projects/REPOS/Testing_repos/SmallMol_implementations/MoietyRecognition/MyTest.sdf')
# for n in range(lib.numMols):
#     sm = lib._mols[n]
#     mf = MoietyRecognition(sm)
#     mf.run()