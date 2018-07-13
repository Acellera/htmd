from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.chemlab.moiety import MoietyRecognition
from htmd.smallmol.chemlab.pkapredictorFiles.pKaTable import PKaTable
from htmd.smallmol.chemlab.periodictable import PeriodicTable
from htmd.smallmol.util import getAtomDistances, _highlight_colors
from htmd.smallmol.chemlab.periodictable import _atoms_positions
import numpy as np
from math import log10

import logging
logger = logging.getLogger(__name__)


class pKaPredictor:
    distanceLimit = {'low': _atoms_positions[2],
                     'medium-low': _atoms_positions[3],
                     'medium-high': _atoms_positions[4],
                     'high': _atoms_positions[5],
                     'very-high': _atoms_positions[6]}

    def __init__(self, smallmol, accuracy='medium-high'):



        if not isinstance(smallmol, SmallMol):
            raise ValueError('Not a SmallMol object. You should provide a valid SmallMol object')

        self.smallmol = self._prepareSmallmol(smallmol)

        self.periodictable = PeriodicTable()
        self.pKaTable = PKaTable()

        self._accuracy = self.distanceLimit[accuracy]

    def accuracy(self, value):

        self._accuracy = self.distanceLimit[value]

    def _prepareSmallmol(self, smallmol):

        if hasattr(smallmol, 'moieties'):
            return smallmol

        logger.info('The molecule passed has not already moieties. The MoietyRecognition feature will be called ')
        mf = MoietyRecognition(smallmol)
        mf.run()
        return smallmol

    def run(self):
        sm = self.smallmol
        moieties = sm.moieties

        pka_moieities, other_moieties = self._inspectMoieties(moieties)

        # TODO accuracy level force moieties merging (eg CX3 type where X is alogen)


        for moi in pka_moieities:
            if len(moi.nestedMoieties) != 0:
                print("Ionizable Moiety")
                pKas_nestedMoieties = []
                for moiNested in moi.nestedMoieties:
                    moiTypeNested = moiNested.moiType
                    if moiTypeNested is None:
                        continue
                    method = 'pKa' + moiTypeNested.title()
                    if not hasattr(self, method):
                        continue
                    pKaMoiNested = getattr(self, method)(moiNested)
                    moiNested.pKa = pKaMoiNested
                    pKas_nestedMoieties.append(pKaMoiNested)
                if len(pKas_nestedMoieties) != 0:
                    self._applyStatisticalFactor(moi.nestedMoieties)
                    moi.pKa = round(min([ mn.pKa for mn in moi.nestedMoieties  ]), 2)
            else:
                #print(moi.moiType, moi.moiOrder)
                method = 'pKa' + moi.moiType.title()
                if not hasattr(self, method):
                    continue
                pKaMoi = getattr(self, method)(moi)
                moi.pKa = round(pKaMoi, 2)

        self._applyStatisticalFactor(moieties)

        return
        #return pka_moieities, other_moieties
        moiTypes = [moi.moiType for moi in moieties]

        # TODO accuracy level force moieties merging (eg CX3 type where X is alogen)

        moieties_analyzed = []

        # for alifatic strain amine
        for i, moiT in enumerate(moiTypes):
           # print(moiT)
            if moiT is None:
                continue
            method = 'pKa' + moiT.title()
            if not hasattr(self, method):
                continue
            moi = moieties[i]
            pKaMoi = getattr(self, method)(moi)
            moi.pKa = pKaMoi
            moieties_analyzed.append(moi)
        return

    def _applyStatisticalFactor(self, moieties):
        print('Trying')

        mois = moieties.copy()

        mois_counts = []

        n = 0

        while len(mois) != 0:
            curr_moi = mois.pop(0)
            curr_Type = curr_moi.moiType
            curr_Order = curr_moi.moiOrder
            curr_pKa = curr_moi.pKa

            mois_counts.append(1)
            for m in moieties:
                if m == curr_moi:
                    continue
                m_Type = m.moiType
                m_Order = m.moiOrder
                m_pKa =  m.pKa

                if m_Type == curr_Type and m_Order == curr_Order and curr_pKa == m_pKa:
                    mois_counts[n] += 1

            n += 1

        for moi, mcount in zip(moieties, mois_counts):
            if moi.pKa is None:
                continue
            if mois_counts != 1:
                moi.pKa += round(log10(mcount), 2)


    def _inspectMoieties(self, moieties):
        sm = self.smallmol

        pka_moieties = []
        other_moieties = []
        for moi in moieties:
            moi.nestedMoieties = []
            toUse = False
            isIonazable = self.pKaTable.isMoietyIonazible(moi.moiType)
            if isIonazable:
                toUse = True
            else:
                # amines ?
                Ns = [ a for a in moi.atoms if sm.element[a] == 'N']
                for n in Ns:
                    moi.addNestedMoiety([n])
                    toUse = True

            if toUse:
                pka_moieties.append(moi)
            else:
                other_moieties.append(moi)

        return pka_moieties, other_moieties


    def depictpKas(self):
        sm = self.smallmol
        moieties = sm.moieties

        highlightatoms = []
        label_pKa = []
        colors = []
        for moi in moieties:
            if moi.pKa is not None:
                highlightatoms.append(moi.atoms.tolist())
                label_pKa.append(moi.pKa)

        for n, moi in enumerate(label_pKa):
            color = _highlight_colors[n % len(_highlight_colors)]
            color = '#%02x%02x%02x' % tuple([ int(c*255) for c in color])
            colors.append(color)


        svg =  sm.depict(ipython=True, highlightAtoms=highlightatoms)

        svg_data = svg.data
        for n, (pKa, color) in enumerate(zip(label_pKa, colors)):
            y = 20 + 20*n
            text =  '<text style="font-size:14px;font-style:normal;font-weight:bold;' \
                    'fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;' \
                    'fill:{}" x="10" y="{}"><tspan>pKa1= {}</tspan></text>'.format(color, y, pKa)
            svg_data = svg_data.replace('</svg>', '{}</svg>'.format(text))
        svg.data = svg_data


        return svg


    def pKaAmine(self, moi, aIdx=None, moieties=None):

        sm = self.smallmol

        validPositions = self._getValidPositions()


        if aIdx is None:
            N_atomIdx = moi.atoms[0]
            query = moi.moiType.title() + "-" + str(moi.moiOrder)
        else:
            N_atomIdx = aIdx

            order = moi._getMoietyOrder(moi.getMoiSmallmol(includeLinks=True), 'amine')
            query = 'Amine' + "-" + str(order)
        # pKa_base - pKa_delta_Methyls - sum(pKa_delta_moieties)

        pKa_base = self.pKaTable.search(query)

        if pKa_base is None:
            raise ValueError('Not available pKa for {}-{} from which start the prediction'.format(moi.moiType, moi.moiOrder))

        print("Base pKa for {} --> {}".format(query, pKa_base))

        pKa_delta_NMethylation = self._deltaPKaNMethylation(N_atomIdx, moi)

        tmp_otherMois = [ m for m in sm.moieties if m != moi]
        otherMois = []
        for m in tmp_otherMois:
            if len(m.nestedMoieties) == 0:
                otherMois.append(m)
            else:
                for mn in m.nestedMoieties:
                    if mn == moi:
                        continue
                    otherMois.append(mn)
        atoms = moi.atoms

        pKa_delta_mois = []
        for omoi in otherMois:
            if omoi.isring and not omoi.isaromatic and not omoi.isNested:
                continue
            oatoms = omoi.atoms

            print(N_atomIdx, oatoms)
            minDist, minPath  = [ getAtomDistances(sm, N_atomIdx, oatoms, mode='min', returnPath=True) for a in atoms ][0]
            print('\n', omoi.moiType, minDist, minPath)

            omoi_sm = omoi.getMoiSmallmol(includeLinks=True, dropAtoms=minPath[:-1])
            #sma = omoi_sm.toSmarts()
            rmol = omoi_sm.toRdkitMol()

            res = self.pKaTable.search(rmol, 'pKaDelta_AmineSubstituents.csv')
            if res is None:
                continue

            position = self.converDistanceToPosition(minDist, minPath, weightbonds=True)
            if isinstance(position, int):
                continue
            if position not in validPositions:
                continue

            table_position = position if position == 'α' else 'β'
            pKas_delta_moi = [res['{}-H-deltaPKa'.format(table_position)],
                              res['{}-CH3-deltaPKa'.format(table_position)],
                              res['{}-R-deltaPKa'.format(table_position)]]


            isnans = np.isnan(pKas_delta_moi)
            #print(pKas_delta_moi)
            #print(">>> ", isnans)
            if sum(isnans) == 0:
                # get subtype
                #print(">>> ", omoi.atoms, omoi.links, omoi.linksTypes)
                olinksType = [ omoi.linksTypes[n]  for n, olink in enumerate(omoi.links) if olink not in minPath ]
                #print(">>> >>>> ", olinksType)
                subType = 0
                if 'R' in olinksType:
                    subType = 2
                elif 'CH3' in olinksType:
                    subType = 1

                pKa_delta_moi = pKas_delta_moi[subType]

            else:
                idx = np.where(isnans == False)[0][0]

                pKa_delta_moi = pKas_delta_moi[idx]

            if position not in ['α', 'β']:
                delta_position = _atoms_positions.index(position) - _atoms_positions.index(table_position)
                pKa_delta_moi = pKa_delta_moi * 0.4**delta_position
            print("pKas available: ", pKas_delta_moi)
            print("pKa delta found: ", pKa_delta_moi)
            pKa_delta_mois.append(pKa_delta_moi)

        pKa = pKa_base - pKa_delta_NMethylation - sum(pKa_delta_mois)

        print("The pKa is ", pKa, '. This is the result of:')
        print('- pKa base ',pKa_base)
        print('- pKa Methylation ', pKa_delta_NMethylation)
        print('- pKa Moieities ', pKa_delta_mois)

        return round(pKa,2)

    def _deltaPKaNMethylation(self, Nidx, moi):
        sm = self.smallmol
        atoms = moi.atoms

        N_atomIdx = Nidx
        N_nbrsIdx = [aIdx for aIdx in sm.neighbors[N_atomIdx] if aIdx not in atoms]
        N_nMethyls = sum([self._isMethyl(nbrIdx, atoms) for nbrIdx in N_nbrsIdx])

        pKa_delta_Methyls = 0.2 * N_nMethyls

        return  pKa_delta_Methyls

    def _isMethyl(self, aIdx, atoms):
        sm = self.smallmol
        a_element = sm.element[aIdx]

        if a_element == 'H':
            return False
        a_nbrsElement = [sm.element[nbrIdx]  for nbrIdx in  sm.neighbors[aIdx] if nbrIdx not in atoms ]
        n_hydrogens = sum([ 1  for el in a_nbrsElement if el == 'H' ])

        if n_hydrogens == 3:
            return  True

        return False

    def converDistanceToPosition(self, distance, path, weightbonds=False, ):


        sm = self.smallmol
        position = distance - 1

        if weightbonds:
            double_triple_bonds = []
            for i in range(len(path) - 1):
                a = path[i]
                n = path[i + 1]
                nbr_idx = sm.neighbors[a].index(n)
                btype = sm.bondtypes[a][nbr_idx]
                if btype in [2, 3]:
                    double_triple_bonds.append(1)
            position = position - sum(double_triple_bonds)

        try:
            return _atoms_positions[position]
        except IndexError:
            return _atoms_positions



    def _getValidPositions(self):

        idx_l = 1
        idx_r =  _atoms_positions.index(self._accuracy) + 1

        validPositions = _atoms_positions[idx_l:idx_r ]

        return validPositions