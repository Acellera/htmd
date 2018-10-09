# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import shutil
import subprocess
import os
from tempfile import TemporaryDirectory
from htmd.parameterization.readers import readPREPI, readFRCMOD, readRTF
import numpy as np
import parmed
import logging
import unittest

logger = logging.getLogger(__name__)


fftypemethods = ('CGenFF_2b6', 'GAFF', 'GAFF2')


def listFftypemethods():
    print('\n'.join(fftypemethods))
    return


def fftype(mol, rtfFile=None, prmFile=None, method='GAFF2', acCharges=None, tmpDir=None, netcharge=None):
    """
    Function to assign atom types and force field parameters for a given molecule.

    The assignment can be done:
      1. For CHARMM CGenFF_2b6 with MATCH (method = 'CGenFF_2b6');
      2. For AMBER GAFF with antechamber (method = 'GAFF');
      3. For AMBER GAFF2 with antechamber (method = 'GAFF2');

    Parameters
    ----------
    mol : Molecule
        Molecule to use for the assignment
    rtfFile : str
        Path to a RTF file from which to read the topology
    prmFile : str
        Path to a PRM file from which to read the parameters
    method : str
        Atomtyping assignment method.
        Use :func:`fftype.listFftypemethods <htmd.parameterization.fftype.listFftypemethods>` to get a list of available
        methods.
        Default: :func:`fftype.defaultFftypemethod <htmd.parameterization.fftype.defaultFftypemethod>`
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
    prm : :class:`ParameterSet <parmed.parameters.ParameterSet>` object
        Returns a parmed ParameterSet object with the parameters.
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The modified Molecule object with the matching atom types for the ParameterSet
    """

    if method not in fftypemethods:
        raise ValueError('Invalid method {}. Available methods {}'.format(method, ','.join(fftypemethods)))

    if netcharge is None:
        netcharge = int(round(np.sum(mol.charge)))
        logger.warning('Using atomic charges from molecule object to calculate net charge')

    prm = names = elements = atomtypes = charges = impropers = masses = None
    if rtfFile and prmFile:
        logger.info('Reading FF parameters from {} and {}'.format(rtfFile, prmFile))
        prm = parmed.charmm.CharmmParameterSet(rtfFile, prmFile)
        names, elements, atomtypes, charges, masses, impropers = readRTF(rtfFile)
        # addParmedResidue(prm, names, elements, atomtypes, charges, impropers)
    else:
        logger.info('Assigning atom types with {}'.format(method))
        # Find the executables
        antechamber_binary = None
        parmchk2_binary = None
        match_binary = None
        if method in ('GAFF', 'GAFF2'):
            antechamber_binary = shutil.which('antechamber')
            if not antechamber_binary:
                raise RuntimeError('antechamber executable not found')
            parmchk2_binary = shutil.which('parmchk2')
            if not parmchk2_binary:
                raise RuntimeError('parmchk2 executable not found')
        elif method == 'CGenFF_2b6':
            match_binary = shutil.which('match-typer')
            if not match_binary:
                raise RuntimeError('match-typer executable not found')
        else:
            raise ValueError('method {} is not implemented'.format(method))

        # Create a temporary directory
        with TemporaryDirectory() as tmpdir:
            # HACK to keep the files
            tmpdir = tmpdir if tmpDir is None else tmpDir

            if method in ('GAFF', 'GAFF2'):
                from htmd.molecule.molecule import Molecule
                # Write the molecule to a file
                mol.write(os.path.join(tmpdir, 'mol.mol2'))

                # Run antechamber

                atomtype = method.lower()

                cmd = [antechamber_binary,
                       '-at', atomtype,
                       '-nc', str(netcharge),
                       '-fi', 'mol2',
                       '-i', 'mol.mol2',
                       '-fo', 'prepi',
                       '-o', 'mol.prepi']
                if acCharges is not None and acCharges != 'None':
                    cmd += ['-c', acCharges]
                returncode = subprocess.call(cmd, cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"antechamber" failed')

                # Run parmchk2
                returncode = subprocess.call([parmchk2_binary,
                                              '-f', 'prepi',
                                              '-s', atomtype,
                                              '-i', 'mol.prepi',
                                              '-o', 'mol.frcmod',
                                              '-a', 'Y'], cwd=tmpdir)
                if returncode != 0:
                    raise RuntimeError('"parmchk2" failed')

                # Check if antechamber did changes in atom names (and suggest the user to fix the names)
                acmol = Molecule(os.path.join(tmpdir, 'NEWPDB.PDB'), type='pdb')
                acmol.name = np.array([n.upper() for n in acmol.name]).astype(np.object)
                if len(np.setdiff1d(mol.name, acmol.name)) != 0 or len(np.setdiff1d(acmol.name, mol.name)) != 0:
                    raise RuntimeError('Initial atom names {} were changed by antechamber to {}. This probably means '
                                       'that the start of the atom name does not match element symbol. Please check '
                                       'the molecule.'.format(','.join(np.setdiff1d(mol.name, acmol.name)),
                                                              ','.join(np.setdiff1d(acmol.name, mol.name))
                                                              )
                                       )

                # Read the results
                prm = parmed.amber.AmberParameterSet(os.path.join(tmpdir, 'mol.frcmod'))
                names, atomtypes, charges, impropers = readPREPI(mol, os.path.join(tmpdir, 'mol.prepi'))
                masses, elements = readFRCMOD(atomtypes, os.path.join(tmpdir, 'mol.frcmod'))
            elif method == 'CGenFF_2b6':

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
                raise ValueError('Invalid method {}'.format(method))

    # Substituting values from the read-in topology
    mol = mol.copy()
    mol.name = names
    mol.element = elements
    mol.atomtype = atomtypes
    if acCharges is not None:
        mol.charge = charges
    mol.impropers = impropers
    if np.sum(mol.masses) == 0:
        mol.masses = masses

    return prm, mol


class TestFftype(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestFftype, self).__init__(*args, **kwargs)
        from htmd.home import home
        from htmd.molecule.molecule import Molecule
        from htmd.parameterization.util import canonicalizeAtomNames, getEquivalentsAndDihedrals
        from yaml import load as yamlload
        from pickle import load as pickleload
        self.refDir = home(dataDir='test-fftype')
        self.Molecule = Molecule
        self.canonicalizeAtomNames = canonicalizeAtomNames
        self.getEquivalentsAndDihedrals = getEquivalentsAndDihedrals
        self.yamlload = yamlload
        self.pickleload = pickleload

    def assertListAlmostEqual(self, list1, list2, places=7):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def setUp(self):
        pass

    def _init_mol(self, molName, ffTypeMethod, chargetuple):
        molFile = os.path.join(self.refDir, '{}.mol2'.format(molName))
        if chargetuple == 'None':
            acCharges = None
            netcharge = None
        else:
            acCharges = chargetuple[0]
            netcharge = chargetuple[1]

        testMolecule = self.Molecule(molFile)
        testMolecule = self.canonicalizeAtomNames(testMolecule, ffTypeMethod)
        testMolecule, _, _ = self.getEquivalentsAndDihedrals(testMolecule)

        from tempfile import mkdtemp
        tmpDir = mkdtemp(suffix='fftype')
        testParameters, testMolecule = fftype(testMolecule,
                                              method=ffTypeMethod,
                                              acCharges=acCharges,
                                              netcharge=netcharge,
                                              tmpDir=tmpDir)
        testIntermediaryFiles = sorted(os.listdir(tmpDir))

        self.testMolecule = testMolecule
        self.testParameters = testParameters
        self.testIntermediaryFiles = testIntermediaryFiles

    def _generate_references(self, name, method):
        import numbers
        import numpy
        from pickle import dump as pickledump

        def mapping(value):
            if isinstance(value, str):
                return '\'{}\''.format(value)
            elif isinstance(value, numbers.Real) or isinstance(value, numbers.Integral):
                return str(value)
            elif isinstance(value, numpy.ndarray):
                return '[{}]'.format(', '.join(map(mapping, list(value))))
            else:
                raise Exception('No mapping for type {}'.format(type(value)))

        print('Copy these to mol_props.yaml')
        for prop in ['name', 'element', 'atomtype', 'charge', 'impropers']:
            print('{}: {}'.format(prop if prop.endswith('s') else '{}s'.format(prop),
                                  '[{}]'.format(', '.join(map(mapping, getattr(self.testMolecule, prop))))))
        print('\nVerify these. They are already written through pickle')
        with open(os.path.join(self.refDir, name, method, 'params.p'), 'wb') as outfile:
            for i in self.testParameters.__dict__:
                print(i, getattr(self.testParameters, i))
            pickledump(self.testParameters, outfile)

        print('\nCopy these to intermediary_files.yaml')
        print('[{}]'.format(', '.join('\'{}\''.format(i) for i in self.testIntermediaryFiles)))

    def _test_mol_props(self, names, elements, atomtypes, charges, impropers):
        self.assertEqual(list(self.testMolecule.name), names)
        self.assertEqual(list(self.testMolecule.element), elements)
        self.assertEqual(list(self.testMolecule.atomtype), atomtypes)
        self.assertListAlmostEqual(list(self.testMolecule.charge), charges)
        if len(impropers) != 0:
            for test, ref in zip(list(self.testMolecule.impropers), impropers):
                self.assertEqual(list(test), ref)
        else:
            self.assertEqual(list(self.testMolecule.impropers), impropers)

    def _test_params(self, params):
        for i in self.testParameters.__dict__:
            self.assertEqual(getattr(self.testParameters, i), getattr(params, i))

    def _test_intermediary_files(self, files):
        self.assertEqual(self.testIntermediaryFiles, files)

    def test_mol(self):

        toTest = (
            ('ethanolamine', 'GAFF2', 'None'),
            ('1o5b_ligand', 'GAFF2', ('gas', 1)),
            ('1z95_ligand', 'GAFF2', ('gas', 0)),
            ('4gr0_ligand', 'GAFF2', ('gas', -1))
        )

        for name, method, chargetuple in toTest:
            with self.subTest(name=name, method=method, chargetuple=chargetuple):
                refDir = os.path.join(self.refDir, name, method)
                self._init_mol(name, method, chargetuple)
                # self._generate_references(name, method)
                with open(os.path.join(refDir, 'mol_props.yaml'), 'r') as infile:
                    refProps = self.yamlload(infile)
                self._test_mol_props(**refProps)
                with open(os.path.join(refDir, 'params.p'), 'rb') as infile:
                    refParams = self.pickleload(infile)
                self._test_params(refParams)
                with open(os.path.join(refDir, 'intermediary_files.yaml'), 'r') as infile:
                    refFiles = self.yamlload(infile)
                self._test_intermediary_files(refFiles)

    def test_broken_atomname(self):
        molFile = os.path.join(self.refDir, '{}.mol2'.format('ethanolamine_wrongnames'))

        testMolecule = self.Molecule(molFile)
        testMolecule = self.canonicalizeAtomNames(testMolecule, 'GAFF2')
        testMolecule, _, _ = self.getEquivalentsAndDihedrals(testMolecule)

        from tempfile import mkdtemp
        tmpDir = mkdtemp(suffix='fftype')
        with self.assertRaises(RuntimeError) as cm:
            _, _ = fftype(testMolecule,
                          method='GAFF2',
                          acCharges=None,
                          netcharge=None,
                          tmpDir=tmpDir)

        self.assertEqual(cm.exception.args[0], 'Initial atom names BR1 were changed by antechamber to CR1. '
                                               'This probably means that the start of the atom name does not '
                                               'match element symbol. Please check the molecule.')


if __name__ == '__main__':

    import sys

    unittest.main(verbosity=2)

    import re
    from htmd.home import home
    from htmd.molecule.molecule import Molecule
    from htmd.parameterization.writers import writeRTF, writePRM
    from htmd.parameterization.util import canonicalizeAtomNames, getEquivalentsAndDihedrals

    # BUG: MATCH does not work on Mac!
    if 'TRAVIS_OS_NAME' in os.environ:
        if os.environ['TRAVIS_OS_NAME'] == 'osx':
            sys.exit(0)

    molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
    refDir = home(dataDir='test-fftype/benzamidine')
    mol = Molecule(molFile)

    with TemporaryDirectory() as tmpDir:

        mol = canonicalizeAtomNames(mol, 'CGenFF_2b6')
        mol, equivalents, all_dihedrals = getEquivalentsAndDihedrals(mol)
        parameters, mol = fftype(mol, method='CGenFF_2b6')
        writeRTF(mol, parameters, 0, os.path.join(tmpDir, 'cgenff.rtf'))
        writePRM(mol, parameters, os.path.join(tmpDir, 'cgenff.prm'))

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
        parameters, mol = fftype(mol, method='CGenFF_2b6', tmpDir=tmpDir)
        assert sorted(os.listdir(tmpDir)) == ['mol.pdb', 'mol.prm', 'mol.rtf', 'top_mol.rtf']
