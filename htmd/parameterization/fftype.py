# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import shutil
import subprocess
import os
from tempfile import TemporaryDirectory
from htmd.parameterization.readers import readPREPI, readRTF
from enum import Enum
import numpy as np
import parmed
import logging

logger = logging.getLogger(__name__)


class FFTypeMethod(Enum):
    NONE = 0
    CHARMM = 1
    AMBER = 2
    CGenFF_2b6 = 1000
    GAFF = 1001
    GAFF2 = 1002


def fftype(mol, rtfFile=None, prmFile=None, method=FFTypeMethod.CGenFF_2b6, acCharges=None, tmpDir=None, netcharge=None):
    """
    Function to assign atom types and force field parameters for a given molecule.

    The assigment can be done:
      1. For CHARMM CGenFF_2b6 with MATCH (method = FFTypeMethod.CGenFF_2b6);
      2. For AMBER GAFF with antechamber (method = FFTypeMethod.GAFF);
      3. For AMBER GAFF2 with antechamber (method = FFTypeMethod.GAFF2);

    Parameters
    ----------
    mol : FFMolecule
        Molecule to use for the assigment
    rtfFile : str
        Path to a RTF file from which to read the topology
    prmFile : str
        Path to a PRM file from which to read the parameters
    method : FFTypeMethod
        Assigment method
    acCharges : str
        Optionally assign charges with antechamber. Check `antechamber -L` for available options. Caution: This will
        overwrite any charges defined in the mol2 file.
    tmpDir: str
        Directory for temporary files. If None, a directory is created and
        deleted automatically.
    netcharge : float
        The net charge of the molecule.

    Returns
    -------
    prm : :class:`ParameterSet<parmed.parameters.ParameterSet>` object
        Returns a parmed ParameterSet object with the parameters.
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The modified Molecule object with the matching atom types for the ParameterSet
    """
    if netcharge is None:
        netcharge = np.sum(mol.charge)
    netcharge = int(round(netcharge))

    if rtfFile and prmFile:
        logger.info('Reading FF parameters from {} and {}'.format(rtfFile, prmFile))
        prm = parmed.charmm.CharmmParameterSet(rtfFile, prmFile)
        names, elements, atomtypes, charges, masses, impropers = readRTF(rtfFile)
        # addParmedResidue(prm, names, elements, atomtypes, charges, impropers)
    else:
        logger.info('Assigned atom types with {}'.format(method.name))
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
                cmd = [antechamber_binary,
                       '-at', atomtype,
                       '-nc', str(netcharge),
                       '-fi', 'mol2',
                       '-i', 'mol.mol2',
                       '-fo', 'prepi',
                       '-o', 'mol.prepi']
                if acCharges is not None:
                    cmd += ['-c', acCharges]
                returncode = subprocess.call(cmd, cwd=tmpdir)
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
                prm = parmed.amber.AmberParameterSet(os.path.join(tmpdir, 'mol.frcmod'))
                names, elements, atomtypes, charges, masses, impropers = readPREPI(mol, os.path.join(tmpdir, 'mol.prepi'))
                # addParmedResidue(prm, names, elements, atomtypes, charges, impropers)
            elif method == FFTypeMethod.CGenFF_2b6:

                # Write the molecule to a file
                mol.write(os.path.join(tmpdir, 'mol.pdb'))

                # Run match-type
                returncode = subprocess.call([match_binary,
                                              '-charge', str(netcharge),
                                              '-forcefield', 'top_all36_cgenff_new',
                                              'mol.pdb'], cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"match-typer" failed')

                prm = parmed.charmm.CharmmParameterSet(os.path.join(tmpdir, 'mol.rtf'), os.path.join(tmpdir, 'mol.prm'))
                names, elements, atomtypes, charges, masses, impropers = readRTF(os.path.join(tmpdir, 'mol.rtf'))
                # addParmedResidue(prm, names, elements, atomtypes, charges, impropers)
            else:
                raise ValueError('Invalide method {}'.format(method))

    # Substituting values from the read-in topology
    mol.name = names
    mol.element = elements
    mol.atomtype = atomtypes
    mol.charge = charges
    mol.impropers = impropers
    if len(mol.masses) == 0:
        mol.masses = masses

    return prm, mol


if __name__ == '__main__':

    import sys
    import re
    from htmd.home import home
    from htmd.molecule.molecule import Molecule
    from htmd.parameterization.writers import writeRTF, writePRM, writeFRCMOD
    from htmd.parameterization.util import canonicalizeAtomNames, getEquivalentsAndDihedrals

    # BUG: MATCH does not work on Mac!
    if 'TRAVIS_OS_NAME' in os.environ:
        if os.environ['TRAVIS_OS_NAME'] == 'osx':
            sys.exit(0)

    molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
    refDir = home(dataDir='test-fftype/benzamidine')
    mol = Molecule(molFile)
    mol = canonicalizeAtomNames(mol)
    mol, equivalents, all_dihedrals = getEquivalentsAndDihedrals(mol)

    with TemporaryDirectory() as tmpDir:

        parameters, mol = fftype(mol, method=FFTypeMethod.CGenFF_2b6)
        writeRTF(mol, parameters, 0, os.path.join(tmpDir, 'cgenff.rtf'))
        writePRM(mol, parameters, os.path.join(tmpDir, 'cgenff.prm'))

        parameters, mol = fftype(mol, method=FFTypeMethod.GAFF)
        writeFRCMOD(mol, parameters, os.path.join(tmpDir, 'gaff.frcmod'))

        parameters, mol = fftype(mol, method=FFTypeMethod.GAFF2)
        writeFRCMOD(mol, parameters, os.path.join(tmpDir, 'gaff2.frcmod'))

        for testFile in os.listdir(refDir):
            print(testFile)
            with open(os.path.join(refDir, testFile)) as refFile:
                with open(os.path.join(tmpDir, testFile)) as tmpFile:
                    # The line order for antechamber is unstable, so the files are sorted.
                    # Also, it get rids of HTMD version string
                    refData = sorted([line for line in refFile.readlines() if not re.search('HTMD', line)])
                    tmpData = sorted([line for line in tmpFile.readlines() if not re.search('HTMD', line)])
                    assert refData == tmpData

    with TemporaryDirectory() as tmpDir:
        parameters, mol = fftype(mol, method=FFTypeMethod.CGenFF_2b6, tmpDir=tmpDir)
        assert sorted(os.listdir(tmpDir)) == ['mol.pdb', 'mol.prm', 'mol.rtf', 'top_mol.rtf']

    with TemporaryDirectory() as tmpDir:
        parameters, mol = fftype(mol, method=FFTypeMethod.GAFF2, tmpDir=tmpDir)
        assert sorted(os.listdir(tmpDir)) == ['ANTECHAMBER.FRCMOD', 'ANTECHAMBER_AC.AC', 'ANTECHAMBER_AC.AC0',
                                              'ANTECHAMBER_BOND_TYPE.AC', 'ANTECHAMBER_BOND_TYPE.AC0',
                                              'ANTECHAMBER_PREP.AC', 'ANTECHAMBER_PREP.AC0', 'ATOMTYPE.INF',
                                              'NEWPDB.PDB', 'PREP.INF', 'mol.frcmod', 'mol.mol2', 'mol.prepi']

