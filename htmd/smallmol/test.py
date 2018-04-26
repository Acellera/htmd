import os
import unittest
from tempfile import NamedTemporaryFile
from glob import glob
from htmd.home import home
from htmd.smallmol.smallmol import SmallMol, SmallMolLib

BENZAMIDINE_N_ATOMS = 18
BENZAMIDINE_N_HEAVYATOMS = 9

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

    def test_loadMol2file(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, BENZAMIDINE_N_ATOMS, 'Atoms not correctly loaded. '
                                                        'Expected: {}; Now: {}'.format(BENZAMIDINE_N_ATOMS, n_atoms))

    def test_loadPdbfile(self):
        pdbfile = os.path.join(self.dataDir, 'ligand.pdb')
        sm = SmallMol(pdbfile)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, LIGAND_N_ATOMS, 'Atoms not correctly loaded. '
                                                        'Expected: {}; Now: {}'.format(LIGAND_N_ATOMS, n_atoms))

    def test_loadSmile(self):
        smi = SMILE_SMI
        sm = SmallMol(smi)
        n_atoms = sm.numAtoms
        self.assertEqual(n_atoms, SMILE_N_ATOMS, 'Atoms not correctly loaded. '
                                                  'Expected: {}; Now: {}'.format(SMILE_N_ATOMS, n_atoms))


    def test_loadSdffile(self):
        sdffile =  os.path.join(self.dataDir, 'fda_drugs_light.sdf')
        lib = SmallMolLib(sdffile)
        n_mols = lib.numMols
        self.assertEqual(n_mols, SDF_N_MOLS, 'Molecules not correctly loaded. '
                                                  'Expected: {}; Now: {}'.format(SDF_N_MOLS, n_mols))

    def test_getAtoms(self):
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

    def test_isChiral(self):
        smi = CHIRAL_SMI
        sm = SmallMol(smi)
        ischiral, details = sm.isChiral(returnDetails=True)
        self.assertListEqual(details, CHIRAL_DETAILS, 'chiral atom does not match.'
                                                      'Expected: {}; Now: {}'.format(CHIRAL_DETAILS, details))

    def test_foundBond(self):
        smi = FOUNDBOND_SMI
        sm = SmallMol(smi)
        isbond_0_N = sm.foundBondBetween('idx 0', 'element N')
        isbond_0_1_single = sm.foundBondBetween('idx 0', 'idx 1', bondtype=1)
        isbond_0_1_double, _ = sm.foundBondBetween('idx 0', 'idx 1', bondtype=2)


        self.assertFalse(isbond_0_N, 'Bond between atom 0 and any nitrogens should not be present')
        self.assertFalse(isbond_0_1_single, 'Bond between atom 0 1 should not be single')
        self.assertTrue(isbond_0_1_double, 'Bond between atom 0 1 should  be double')


    def test_generateConformers(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        current_conformer = sm.numConformers
        sm.generateConformers(num_confs=10, append=False)
        n_conformers = sm.numConformers

        self.assertGreater(n_conformers, current_conformer, 'The generation of conforemr should provide at least the '
                                                            'same amount of conformer')

    def test_writeGenerateAndWriteConformers(self):
        mol2file = os.path.join(self.dataDir, 'benzamidine.mol2')
        sm = SmallMol(mol2file)
        sm.generateConformers(num_confs=10, append=False)
        tmpdir = NamedTemporaryFile().name
        sm.writeConformers(savefolder=tmpdir)
        direxists = os.path.isdir(tmpdir)
        n_files = len(glob(os.path.join(tmpdir, '*.sdf')))
        self.assertTrue(direxists, 'The directory where to store the conformations where not created')
        self.assertGreater(n_files, 1, 'None conformations were written. At least one should be present')


    def test_removeGenerateConformer(self):
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

    def test_invertChirality(self):
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







if __name__ == '__main__':
    unittest.main(verbosity=2)