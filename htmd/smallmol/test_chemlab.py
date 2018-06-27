import os
import unittest
from tempfile import NamedTemporaryFile
from glob import glob
from htmd.home import home
from htmd.molecule.molecule import Molecule
from htmd.smallmol.smallmol import SmallMol, SmallMolLib
from htmd.smallmol.chemlab.periodictable import *
from htmd.smallmol.chemlab.builder import *
from htmd.smallmol.util import calculateAngle
import rdkit
from rdkit.Chem import MolFromSmiles
from collections import Counter
from itertools import combinations

BOND_RADIUS_C = 0.77
BOND_RADIUS_H = 0.33

BENZAMIDINE_N_CHARGED_IDX = 12
BENZAMIDINE_N_CHARGED_MISSING_VALENCE = 0

BENZAMIDINE_ELEMENTS = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'N', 'N', 'H', 'H', 'H', 'H']
BENZAMIDINE_ELEMENTS_NOHS = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N']
BENZAMIDINE_ELEMENTS_NO_NONPOLARHS = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N', 'H', 'H', 'H', 'H']
BENZAMIDINE_ELEMENTS_NO_POLARHS = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'N', 'N']

ETANOLAMMINE_SMILE = 'CCN'
ETANOLAMMINE_N_MISSING_VALENCE = 2
ETANOLAMMINE_ELEMENTS_ALL =  ['C', 'C', 'N', 'H', 'H', 'H', 'H', 'H', 'H', 'H']
ETANOLAMMINE_ELEMENTS_POLARHS =['C', 'C', 'N', 'H', 'H']

INFO_H = {'element': 'H', 'charge': 0, 'formalcharge': 0, 'chiral': '', 'hybridization': 1, 'bondtype': 1,
          'coords': np.array([[[0], [0], [0]]]), 'attachTo': ''}


METHANE_SMILE = 'C'
METHANE_ANGLE = 109.47
ETHENE_SMILE = 'C=C'
ETHENE_ANGLE = 120.0
ETHINE_SMILE = 'C#C'
ETHINE_ANGLE = 180.0

class TestSmallMol(unittest.TestCase):

    def setUp(self):
        self.dataDir= home('test-smallmol')

    def _getangle(self, smallmol, centerIdx, atom1Idx, atom2Idx):

        centercoords = smallmol.coords[centerIdx]
        a1coords = smallmol.coords[atom1Idx]
        a2coords = smallmol.coords[atom2Idx]

        return calculateAngle(centercoords, a1coords, a2coords, deg=True)

    def _getCombinations(self, smallmol, centers):
        neighbors = [smallmol.neighbors[c] for c in centers]
        neighbors_combination = [list(combinations(ns, 2)) for ns in neighbors]

        atoms_set_angle = [ [ca, c[0], c[1]]  for ca, comb in zip(centers, neighbors_combination) for c in comb ]

        return atoms_set_angle

    def test_01_getBondRadius(self):
        PT = PeriodicTable()
        bond_radius_C = PT.getBondRadius('C')
        bond_radius_H = PT.getBondRadius('H')

        self.assertEqual(bond_radius_C, BOND_RADIUS_C, msg="The bond radius for the C atom is not as expected")
        self.assertEqual(bond_radius_H, BOND_RADIUS_H, msg="The bond radius for the H atom is not as expected")

    def test_02_getMissingValence(self):

        mol2file =  os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)

        PT = PeriodicTable()

        atom_formalcharge = sm.formalcharge[BENZAMIDINE_N_CHARGED_IDX]

        omv = PT.getMissingValence(sm, BENZAMIDINE_N_CHARGED_IDX, formalcharge=atom_formalcharge)

        self.assertEqual(omv, BENZAMIDINE_N_CHARGED_MISSING_VALENCE, msg="The missing valence for the N charged is not "
                                                                         "as expected")
        sm = SmallMol(ETANOLAMMINE_SMILE, removeHs=False, fixHs=False)

        n_idx = np.where(sm.element == 'N')[0][0]
        n_formalcharge = sm.formalcharge[n_idx]

        omv = PT.getMissingValence(sm, n_idx, formalcharge=n_formalcharge)

        self.assertEqual(omv, ETANOLAMMINE_N_MISSING_VALENCE, msg="The missing valence for the N is not "
                                                                         "as expected")
    def test_03_getAtomInfo(self):

        def checkDictionary(dict1, dict2):
            match = True

            for k in info_H:
                if k == 'coords':
                    _match = np.array_equal(info_H[k], INFO_H[k])
                else:
                    _match = info_H[k] == INFO_H[k]
                if _match == False:
                    match = _match
            return match

        PT = PeriodicTable()

        info_H = PT.getAtom('H')

        self.assertTrue(checkDictionary(info_H, INFO_H), msg="The atom information for H are not the expected one")

        # passing coords
        info_H = PT.getAtom('H', coords=np.array([1,2,3]))
        self.assertFalse(checkDictionary(info_H, INFO_H), msg="The atom information for H are not the expected one")

    def test_04_Builder(self):

        sm = SmallMol(ETANOLAMMINE_SMILE, removeHs=True, fixHs=False)
        sm_c = sm.copy()

        b = Builder(sm_c)
        sm_b = b.getSmallMol()

        expectedFalse = np.array_equal(sm_b.coords, sm.coords)

        self.assertIsInstance(b, Builder, msg="The molecule was not creted in the Builder")
        self.assertFalse(expectedFalse, msg="A initial conformation was not created")

    def test_05_Builder_load(self):
        sm = SmallMol(ETANOLAMMINE_SMILE, removeHs=True, fixHs=False)
        sm_c = sm.copy()

        b = Builder()
        b.loadMol(sm_c)
        sm_b = b.getSmallMol()

        expectedFalse = np.array_equal(sm_b.coords, sm.coords)

        self.assertIsInstance(b, Builder, msg="The molecule was not creted in the Builder")
        self.assertFalse(expectedFalse, msg="A initial conformation was not created")

    def test_06_getSmallMol(self):
        sm = SmallMol(ETANOLAMMINE_SMILE, removeHs=True, fixHs=False)
        sm_c = sm.copy()

        b = Builder(sm_c)
        sm_b = b.getSmallMol()

        self.assertIsInstance(sm_b, SmallMol, msg="The molecule was not loaded correctly")

    def test_07_addHydrogens(self):
        sm = SmallMol(ETANOLAMMINE_SMILE, removeHs=True, fixHs=False)
        sm_c = sm.copy()

        b = Builder(sm_c)
        b.addHydrogens()

        sm_b = b.getSmallMol()
        sm_b_elements = sm_b.element.tolist()

        self.assertListEqual(sm_b_elements, ETANOLAMMINE_ELEMENTS_ALL)

        sm_c = sm.copy()

        b = Builder(sm_c)
        b.addHydrogens(onlyExplicit=True)
        sm_b = b.getSmallMol()
        sm_b_elements = sm_b.element.tolist()

        self.assertListEqual(sm_b_elements, ETANOLAMMINE_ELEMENTS_POLARHS)

    def test_08_removeHydrogens(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')

        sm = SmallMol(mol2file)

        sm_c = sm.copy()

        b = Builder(sm_c)
        b.removeHydrogens()

        sm_b = b.getSmallMol()

        sm_b_element = sm_b.element.tolist()
        self.assertListEqual(sm_b_element, BENZAMIDINE_ELEMENTS_NOHS)

        sm_c = sm.copy()

        b = Builder(sm_c)
        b.removeHydrogens(removeNonPolars=False)

        sm_b = b.getSmallMol()

        sm_b_element = sm_b.element.tolist()
        self.assertListEqual(sm_b_element, BENZAMIDINE_ELEMENTS_NO_POLARHS)

        sm_c = sm.copy()

        b = Builder(sm_c)
        b.removeHydrogens(removePolars=False)

        sm_b = b.getSmallMol()

        sm_b_element = sm_b.element.tolist()
        self.assertListEqual(sm_b_element, BENZAMIDINE_ELEMENTS_NO_NONPOLARHS)

    def  test_09_anglesMETHANE(self):

        sm = SmallMol(METHANE_SMILE, fixHs=False, removeHs=False)

        b = Builder(sm)
        b.addHydrogens()

        sm_b = b.getSmallMol()

        centercoords = sm_b.coords[0]

        setone = [sm_b.coords[1], sm_b.coords[2]]
        settwo = [sm_b.coords[2], sm_b.coords[3]]
        setthree = [sm_b.coords[3], sm_b.coords[4]]
        setfour = [sm_b.coords[4], sm_b.coords[1]]

        sets = [setone, settwo, setthree, setfour]

        angles = []
        for s in sets:
            a1coords = s[0]
            a2coords = s[1]
            angle = calculateAngle(centercoords, a1coords, a2coords, deg=True)
            angles.append(angle)


        angle_sample = angles[0]
        equal_angles = Counter(angles)[METHANE_ANGLE]

        self.assertEqual(angle_sample, METHANE_ANGLE, msg="The computed angle is not the expected one. Probably atoms "
                                           "coords not correctly predict")
        self.assertEqual(equal_angles, 4, msg="Not all the angle have the same value.")

    def test_10_angleETHENE(self):
        sm = SmallMol(ETHENE_SMILE, fixHs=False, removeHs=False)

        b = Builder(sm)
        b.addHydrogens()

        sm_b = b.getSmallMol()

        centercoords1 = sm_b.coords[0]
        centercoords2 = sm_b.coords[1]

        set1one = [sm_b.coords[1], sm_b.coords[2]]
        set1two = [sm_b.coords[2], sm_b.coords[3]]
        set1three = [sm_b.coords[3], sm_b.coords[1]]

        set2one = [sm_b.coords[0], sm_b.coords[4]]
        set2two = [sm_b.coords[4], sm_b.coords[5]]
        set2three = [sm_b.coords[5], sm_b.coords[0]]

        sets1 = [set1two, set1two, set1three]
        sets2 = [set2two, set2two, set2three]

        angles1 = []
        for s in sets1:
            a1coords = s[0]
            a2coords = s[1]
            angle = calculateAngle(centercoords1, a1coords, a2coords, deg=True)
            angles1.append(angle)

        angles2 = []
        for s in sets2:
            a1coords = s[0]
            a2coords = s[1]
            angle = calculateAngle(centercoords2, a1coords, a2coords, deg=True)
            angles2.append(angle)

        angle1_sample = angles1[0]
        equal1_angles = Counter(angles1)[ETHENE_ANGLE]

        angle2_sample = angles2[0]
        equal2_angles = Counter(angles2)[ETHENE_ANGLE]

        self.assertEqual(angle1_sample, ETHENE_ANGLE, msg="The computed angle is not the expected one. Probably atoms "
                                           "coords not correctly predict")
        self.assertEqual(equal1_angles, 3, msg="Not all the angle have the same value.")

        self.assertEqual(angle2_sample, ETHENE_ANGLE, msg="The computed angle is not the expected one. Probably atoms "
                                                          "coords not correctly predict")
        self.assertEqual(equal2_angles, 3, msg="Not all the angle have the same value.")

    def test_11_angleETHINE(self):
        sm = SmallMol(ETHINE_SMILE, removeHs=False, fixHs=False)

        b = Builder(sm)
        b.addHydrogens()

        sm_b = b.getSmallMol()

        centercoords1 = sm_b.coords[0]
        centercoords2 = sm_b.coords[1]

        set1one = [sm_b.coords[1], sm_b.coords[2]]
        set2one = [sm_b.coords[0], sm_b.coords[3]]

        sets1 = [set1one]
        sets2 = [set2one]

        angles1 = []
        for s in sets1:
            a1coords = s[0]
            a2coords = s[1]
            angle = calculateAngle(centercoords1, a1coords, a2coords, deg=True)
            angles1.append(angle)

        angles2 = []
        for s in sets2:
            a1coords = s[0]
            a2coords = s[1]

            angle = calculateAngle(centercoords2, a1coords, a2coords, deg=True)
            angles2.append(angle)

        angle1_sample = angles1[0]
        equal1_angles = Counter(angles1)[ETHINE_ANGLE]

        angle2_sample = angles2[0]
        equal2_angles = Counter(angles2)[ETHINE_ANGLE]

        self.assertEqual(angle1_sample, ETHINE_ANGLE, msg="The computed angle is not the expected one. Probably atoms "
                                                          "coords not correctly predict")
        self.assertEqual(equal1_angles, 1, msg="Not all the angle have the same value.")

        self.assertEqual(angle2_sample, ETHINE_ANGLE, msg="The computed angle is not the expected one. Probably atoms "
                                                          "coords not correctly predict")
        self.assertEqual(equal2_angles, 1, msg="Not all the angle have the same value.")



if __name__ == '__main__':
    tloader = unittest.loader.TestLoader()
    tloader.sortTestMethodsUsing = None
    unittest.main(verbosity=2,testLoader=tloader)