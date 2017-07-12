# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import shutil
import subprocess
import os
from tempfile import TemporaryDirectory
from htmd.parameterization.ff import RTF, PRM, AmberRTF, AmberPRM
from enum import Enum


class FFTypeMethod(Enum):
    CHARMM = 1
    AMBER = 2
    CGenFF_2b6 = 1000
    GAFF = 1001
    GAFF2 = 1002


class FFType:
    """
    Class to assign atom types and force field parameters for a given molecule.

    The assigment can be done:
      1. For CHARMM CGenFF_2b6 with MATCH (method = FFTypeMethod.CGenFF_2b6);
      2. For AMBER GAFF with antechamber (method = FFTypeMethod.GAFF);
      3. For AMBER GAFF2 with antechamber (method = FFTypeMethod.GAFF2);

    Parameters
    ----------
    mol : FFMolecule
        Molecule to use for the assigment
    method : FFTypeMethod
        Assigment method
    tmpDir: str
        Directory for temporary files. If None, a directory is created and
        deleted automatically.
    """

    def __init__(self, mol, method=FFTypeMethod.CGenFF_2b6, tmpDir=None):

        # Find the executables
        if method == FFTypeMethod.GAFF or method == FFTypeMethod.GAFF2:
            antechamber_binary = shutil.which("antechamber")
            if not antechamber_binary:
                raise RuntimeError("antechamber executable not found")

            parmchk2_binary = shutil.which("parmchk2")
            if not parmchk2_binary:
                raise RuntimeError("parmchk2 executable not found")

        elif method == FFTypeMethod.CGenFF_2b6:
            match_binary = shutil.which("match-typer")
            if not match_binary:
                raise RuntimeError("match-typer executable not found")

        else:
            raise ValueError('method')

        # Create a temporary directory
        with TemporaryDirectory() as tmpdir:

            # HACK to keep the files
            tmpdir = tmpdir if tmpDir is None else tmpDir

            if method == FFTypeMethod.GAFF or method == FFTypeMethod.GAFF2:

                # Write the molecule to a file
                mol.write(os.path.join(tmpdir, 'mol.mol2'))

                # Run antechamber
                if method == FFTypeMethod.GAFF:
                    atomtype = "gaff"
                elif method == FFTypeMethod.GAFF2:
                    atomtype = "gaff2"
                else:
                    raise ValueError('method')
                returncode = subprocess.call([antechamber_binary,
                                              '-at', atomtype,
                                              '-nc', str(mol.netcharge),
                                              '-fi', 'mol2',
                                              '-i', 'mol.mol2',
                                              '-fo', 'prepi',
                                              '-o', 'mol.prepi'], cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"antechamber" failed')

                # Run parmchk2
                returncode = subprocess.call([parmchk2_binary,
                                              '-f', 'prepi',
                                              '-i', 'mol.prepi',
                                              '-o', 'mol.frcmod',
                                              '-a', 'Y'], cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"parmchk2" failed')

                # Read the results
                self._rtf = AmberRTF(mol, os.path.join(tmpdir, 'mol.prepi'),
                                          os.path.join(tmpdir, 'mol.frcmod'))
                self._prm = AmberPRM(os.path.join(tmpdir, 'mol.prepi'),
                                     os.path.join(tmpdir, 'mol.frcmod'))

            elif method == FFTypeMethod.CGenFF_2b6:

                # Write the molecule to a file
                mol.write(os.path.join(tmpdir, 'mol.pdb'))

                # Run match-type
                returncode = subprocess.call([match_binary,
                                              '-charge', str(mol.netcharge),
                                              '-forcefield', 'top_all36_cgenff_new',
                                              'mol.pdb'], cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"match-typer" failed')

                # Read the results
                self._rtf = RTF(os.path.join(tmpdir, 'mol.rtf'))
                self._prm = PRM(os.path.join(tmpdir, 'mol.prm'))

            else:
                raise ValueError('method')

if __name__ == '__main__':

    import sys
    import re
    from htmd.home import home
    from htmd.parameterization.ffmolecule import FFMolecule

    # BUG: MATCH does not work on Mac!
    if 'TRAVIS_OS_NAME' in os.environ:
        if os.environ['TRAVIS_OS_NAME'] == 'osx':
            sys.exit(0)

    molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
    refDir = home('test-fftype/benzamidine')
    mol = FFMolecule(molFile)

    with TemporaryDirectory() as tmpDir:

        ff = FFType(mol, method=FFTypeMethod.CGenFF_2b6)
        ff._rtf.write(os.path.join(tmpDir, 'cgenff.rtf'))
        ff._prm.write(os.path.join(tmpDir, 'cgenff.prm'))

        ff = FFType(mol, method=FFTypeMethod.GAFF)
        ff._prm.writeFrcmod(ff._rtf, os.path.join(tmpDir, 'gaff.frcmod'))

        ff = FFType(mol, method=FFTypeMethod.GAFF2)
        ff._prm.writeFrcmod(ff._rtf, os.path.join(tmpDir, 'gaff2.frcmod'))

        # TODO: remove this condition when htmd-data is updated.
        if not os.path.isdir(refDir):
             sys.exit(0)

        for testFile in os.listdir(refDir):
            print(testFile)
            with open(os.path.join(refDir, testFile)) as refFile:
                with open(os.path.join(tmpDir, testFile)) as tmpFile:
                    # The line order for antichamber is unstable, so the files are sorted.
                    # Also, it get rids of HTMD version string
                    refData = sorted([line for line in refFile.readlines() if not re.search('HTMD', line)])
                    tmpData = sorted([line for line in tmpFile.readlines() if not re.search('HTMD', line)])
                    assert refData == tmpData

    with TemporaryDirectory() as tmpDir:

        ff = FFType(mol, method=FFTypeMethod.CGenFF_2b6, tmpDir=tmpDir)
        assert sorted(os.listdir(tmpDir)) == ['mol.pdb', 'mol.prm', 'mol.rtf', 'top_mol.rtf']

    with TemporaryDirectory() as tmpDir:

        ff = FFType(mol, method=FFTypeMethod.GAFF2, tmpDir=tmpDir)
        assert sorted(os.listdir(tmpDir)) == ['ANTECHAMBER.FRCMOD', 'ANTECHAMBER_AC.AC', 'ANTECHAMBER_AC.AC0',
                                              'ANTECHAMBER_BOND_TYPE.AC', 'ANTECHAMBER_BOND_TYPE.AC0',
                                              'ANTECHAMBER_PREP.AC', 'ANTECHAMBER_PREP.AC0', 'ATOMTYPE.INF',
                                              'NEWPDB.PDB', 'PREP.INF', 'mol.frcmod', 'mol.mol2', 'mol.prepi']

