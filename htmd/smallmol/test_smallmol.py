import os
import unittest
from tempfile import NamedTemporaryFile
from glob import glob
from htmd.home import home
from htmdmol.molecule import Molecule
from htmd.smallmol.smallmol import SmallMol, SmallMolLib
import rdkit
from rdkit.Chem import MolFromSmiles

BENZAMIDINE_N_ATOMS = 18
BENZAMIDINE_N_HEAVYATOMS = 9
BENZAMIDINE_BONDTYPES = ['ar', 'ar', '1', 'ar', '1', 'ar', '1', 'ar', '1', 'ar', '1', '1', '2', '1', '1', '1', '1', '1']
BENZAMIDINE_BOND_ATOMS = [[0, 1], [0, 5], [0, 6], [1, 2], [1, 7], [2, 3], [2, 8], [3, 4], [3, 9], [4, 5], [4, 10],
                          [5, 11], [6, 12], [6, 13], [12, 16], [12, 17], [13, 14], [13, 15]]

LIGAND_N_ATOMS = 64
LIGAND_N_HEAVYATOMS = 35

SMILE_SMI = 'c1ccccc1O'
SMILE_N_ATOMS = 13

SDF_N_MOLS = 100

PHENOL_ELEMENT_IDX_1 = 'C'
PHENOL_ELEMENT_NEIGHBORS_OX = [5, 12]
PHENOL_BTYPES_OX = [1, 1]

CHIRAL_SMI = 'C[C@H](Cl)F'
CHIRAL_DETAILS = [('C1', 'S')]

FOUNDBOND_SMI = 'C=CN'

class TestSmallMol(unittest.TestCase):

    def setUp(self):
        self.dataDir= home('test-smallmol')

    def test_01_loadMol2file(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, BENZAMIDINE_N_ATOMS, 'Atoms not correctly loaded. '
                                                        'Expected: {}; Now: {}'.format(BENZAMIDINE_N_ATOMS, n_atoms))

    def test_02_loadPdbfile(self):
        pdbfile = os.path.join(self.dataDir, 'ligand.pdb')
        sm = SmallMol(pdbfile)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, LIGAND_N_ATOMS, 'Atoms not correctly loaded. '
                                                        'Expected: {}; Now: {}'.format(LIGAND_N_ATOMS, n_atoms))

    def test_03_loadSmile(self):
        smi = SMILE_SMI
        sm = SmallMol(smi)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, SMILE_N_ATOMS, 'Atoms not correctly loaded. '
                                                  'Expected: {}; Now: {}'.format(SMILE_N_ATOMS, n_atoms))


    def test_04_loadSdffile(self):
        sdffile =  os.path.join(self.dataDir, 'fda_drugs_light.sdf')
        lib = SmallMolLib(sdffile)
        n_mols = lib.numMols
        self.assertEqual(n_mols, SDF_N_MOLS, 'Molecules not correctly loaded. '
                                                  'Expected: {}; Now: {}'.format(SDF_N_MOLS, n_mols))

    def test_05_getAtoms(self):
        smi = SMILE_SMI
        sm = SmallMol(smi)
        element_idx_1 = sm.get('element', 'idx 1')[0]
        neighbors_element_O = sm.get('neighbors', 'element O')[0]
        btypes_element_O = sm.get('bondtypes', 'element O', convertType=False)[0]

        self.assertEqual(element_idx_1, PHENOL_ELEMENT_IDX_1, 'Element of the first atom does not correspond'
                                                              'Expect: {}; Now: {}'.format(element_idx_1, PHENOL_ELEMENT_IDX_1))
        self.assertListEqual(neighbors_element_O, PHENOL_ELEMENT_NEIGHBORS_OX,  'Neighbors atoms of the oxygen atom do not correspond'
                           'Expected: {}; Now: {}'.format(PHENOL_ELEMENT_NEIGHBORS_OX, neighbors_element_O))

        self.assertListEqual(btypes_element_O, PHENOL_BTYPES_OX, 'Bondtypes of the oxygen atom do not correspond:'
                                                                 'Expeected: {}; Now: {}'.format(btypes_element_O, PHENOL_BTYPES_OX))

    def test_06_isChiral(self):
        smi = CHIRAL_SMI
        sm = SmallMol(smi)
        ischiral, details = sm.isChiral(returnDetails=True)
        self.assertListEqual(details, CHIRAL_DETAILS, 'chiral atom does not match.'
                                                      'Expected: {}; Now: {}'.format(CHIRAL_DETAILS, details))

    def test_07_foundBond(self):
        smi = FOUNDBOND_SMI
        sm = SmallMol(smi)
        isbond_0_N = sm.foundBondBetween('idx 0', 'element N')
        isbond_0_1_single = sm.foundBondBetween('idx 0', 'idx 1', bondtype=1)
        isbond_0_1_double, _ = sm.foundBondBetween('idx 0', 'idx 1', bondtype=2)


        self.assertFalse(isbond_0_N, 'Bond between atom 0 and any nitrogens should not be present')
        self.assertFalse(isbond_0_1_single, 'Bond between atom 0 1 should not be single')
        self.assertTrue(isbond_0_1_double, 'Bond between atom 0 1 should  be double')


    def test_08_generateConformers(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        current_conformer = sm.numConformers
        sm.generateConformers(num_confs=10, append=False)
        n_conformers = sm.numConformers

        self.assertGreater(n_conformers, current_conformer, 'The generation of conforemr should provide at least the '
                                                            'same amount of conformer')

    def test_09_writeGenerateAndWriteConformers(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        sm.generateConformers(num_confs=10, append=False)
        tmpdir = NamedTemporaryFile().name
        sm.writeConformers(savefolder=tmpdir)
        direxists = os.path.isdir(tmpdir)
        n_files = len(glob(os.path.join(tmpdir, '*.sdf')))
        self.assertTrue(direxists, 'The directory where to store the conformations where not created')
        self.assertGreater(n_files, 1, 'None conformations were written. At least one should be present')


    def test_10_removeGenerateConformer(self):
        molsmile = SMILE_SMI
        sm = SmallMol(molsmile)
        sm.generateConformers(num_confs=10, append=False)
        n_confs = sm.numConformers
        sm.removeConformers([0])
        n_confs_del = sm.numConformers
        sm.removeConformers()
        n_confs_zero = sm.numConformers

        self.assertEqual(n_confs_del, n_confs - 1, "The number of conformations after the deletion was not reduced of "
                                                   "exactly one unit")
        self.assertEqual(n_confs_zero, 0, "The number of conformations after the deletion was not reduced to 0")

    def test_11_invertChirality(self):
        molsmile = CHIRAL_SMI
        sm = SmallMol(molsmile)
        aname = CHIRAL_DETAILS[0][0]
        chiral = CHIRAL_DETAILS[0][1]
        aidx = sm.get('idx', 'atomname {}'.format(aname))[0]

        sm.invertChirality(aidx)
        newchiral = sm.isChiral(returnDetails=True)[1][0][-1]

        self.assertNotEqual(chiral, newchiral, msg="The chirality was not formally changed")

        sm.generateConformers(num_confs=1, append=False)
        m = sm.toMolecule()
        fname = NamedTemporaryFile().name + '.mol2'
        m.write(fname)

        sm2 = SmallMol(fname)
        newchiral_confirm = sm2.isChiral(returnDetails=True)[1][0][-1]

        self.assertEqual(newchiral, newchiral_confirm, msg="The chirality was not structurally changed")

    def test_12_convertToRdkit(self):
        smimol = SMILE_SMI
        sm = SmallMol(smimol, removeHs=True, fixHs=False)
        mrd = MolFromSmiles(smimol)
        mrd_natom = mrd.GetNumAtoms()

        sm_rd = sm.toRdkitMol(includeConformer=True)
        sm_rd_natoms = sm_rd.GetNumAtoms()


        self.assertIsInstance(sm_rd, rdkit.Chem.rdchem.Mol, msg="The conversion of the SmallMol object into the rdkit"
                                                                "Mol one get wrong")
        self.assertEqual(sm_rd_natoms, mrd_natom, msg="NUmber of atoms different. The handle and convertion of the "
                                                      "SmallMol object into the rdkit Mol one probably get wrong")


    def test_13_convertToMolecule(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        mol = Molecule(mol2file)
        mol_elements = mol.element.tolist()

        sm_htmd = sm.toMolecule()
        sm_htmd_element = sm_htmd.element.tolist()

        self.assertListEqual(sm_htmd_element, mol_elements, msg="The elements found are different. The handle and convertion"
                                                            "of the SmallMol object into the htmd.Molecule one, probably "
                                                            "get wrong")

    def test_14_getBonds(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)

        bonds, bondstype = sm.getBonds()

        self.assertListEqual(bonds.tolist(), BENZAMIDINE_BOND_ATOMS, msg="The atoms in bonds are not the same of the reference")

        self.assertListEqual(bondstype.tolist(), BENZAMIDINE_BONDTYPES, msg="The bonds type are not the same of the reference")

    def test_15_depict(self):
        import IPython
        refimg = os.path.join(self.dataDir, 'benzamidine.svg')
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')

        sm = SmallMol(mol2file)

        img_name = NamedTemporaryFile().name + '.svg'
        sm.depict(sketch=True, filename=img_name)
        _img = sm.depict(sketch=True, ipython=True)

        refimg_size = os.path.getsize(refimg)
        sm_img_size = os.path.getsize(img_name)

        self.assertIsInstance(_img, IPython.core.display.SVG, msg="The object is not an IPython image as expected")
        self.assertEqual(sm_img_size, refimg_size, msg="The svg image does not have the same size of the reference")






if __name__ == '__main__':
    tloader = unittest.loader.TestLoader()
    tloader.sortTestMethodsUsing = None
    unittest.main(verbosity=2,testLoader=tloader)