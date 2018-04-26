import os
import unittest
from tempfile import NamedTemporaryFile
from glob import glob
from htmd.home import home
from htmd.molecule.molecule import Molecule
from htmd.smallmol.smallmol import SmallMol, SmallMolLib
from htmd.smallmol.chemlab.periodictable import *
from htmd.smallmol.chemlab.builder import *
import rdkit
from rdkit.Chem import MolFromSmiles


class TestSmallMol(unittest.TestCase):

    def setUp(self):
        self.dataDir= home('test-smallmol')




if __name__ == '__main__':
    tloader = unittest.loader.TestLoader()
    tloader.sortTestMethodsUsing = None
    unittest.main(verbosity=2,testLoader=tloader)