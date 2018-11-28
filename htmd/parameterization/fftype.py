# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import subprocess
import os
import re
import unittest
from tempfile import TemporaryDirectory, TemporaryFile

import numpy as np

logger = logging.getLogger(__name__)


fftypemethods = ('CGenFF_2b6', 'GAFF', 'GAFF2')


def _canonicalizeAtomNames(mol):
    """
    This fixes up the atom naming and reside name to be consistent.
    NB this scheme matches what MATCH does.
    Don't change it or the naming will be inconsistent with the RTF.
    """

    mol = mol.copy()

    mol.segid[:] = 'L'
    logger.debug('Rename segment to %s' % mol.segid[0])
    mol.resname[:] = 'MOL'
    logger.debug('Rename residue to %s' % mol.resname[0])

    sufices = {}
    for i in range(mol.numAtoms):
        name = mol.element[i].upper()
        sufices[name] = sufices.get(name, 0) + 1
        name += str(sufices[name])

        logger.debug('Rename atom %d: %-4s --> %-4s' % (i, mol.name[i], name))
        mol.name[i] = name

    return mol


def fftype(mol, rtfFile=None, prmFile=None, method='GAFF2', acCharges=None, tmpDir=None, netcharge=None):
    """
    Assing atom types and force field parameters for a given molecule.
    Additionally, atom masses and improper dihedral are set.
    Optionally, atom charges can be set if `acCharges` is set (see below).

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
        Optionally assign charges with antechamber. Check `antechamber -L` for available options.
        Note: only works for GAFF and GAFF2.
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

    import parmed

    if method not in fftypemethods:
        raise ValueError('Invalid method {}. Available methods {}'.format(method, ','.join(fftypemethods)))

    if method == 'CGenFF_2b6' and acCharges:
        raise ValueError('acCharges')

    if netcharge is None:
        netcharge = int(round(np.sum(mol.charge)))
        logger.warning('Molecular charge is set to {} by adding up the atomic charges'.format(netcharge))

    if rtfFile and prmFile:

        from htmd.parameterization.readers import readRTF

        logger.info('Reading FF parameters from {} and {}'.format(rtfFile, prmFile))
        prm = parmed.charmm.CharmmParameterSet(rtfFile, prmFile)
        names, elements, atomtypes, charges, masses, impropers = readRTF(rtfFile)

    else:
        logger.info('Assigning atom types with {}'.format(method))

        renamed_mol = _canonicalizeAtomNames(mol)

        # Create a temporary directory
        with TemporaryDirectory() as tmpdir:

            # HACK to keep the files
            tmpdir = tmpdir if tmpDir is None else tmpDir
            logger.debug('Temporary directory: {}'.format(tmpdir))

            if method in ('GAFF', 'GAFF2'):

                from htmd.molecule.molecule import Molecule
                from htmd.parameterization.readers import readPREPI, readFRCMOD

                # Write the molecule to a file
                renamed_mol.write(os.path.join(tmpdir, 'mol.mol2'))

                atomtype = method.lower()

                # Set arguments
                cmd = ['antechamber',
                       '-at', atomtype,
                       '-nc', str(netcharge),
                       '-fi', 'mol2',
                       '-i', 'mol.mol2',
                       '-fo', 'prepi',
                       '-o', 'mol.prepi']
                if acCharges is not None:
                    cmd += ['-c', acCharges]

                # Run antechamber
                with TemporaryFile() as stream:
                    if subprocess.call(cmd, cwd=tmpdir, stdout=stream, stderr=stream) != 0:
                        raise RuntimeError('"antechamber" failed')
                    stream.seek(0)
                    for line in stream.readlines():
                        logger.debug(line)

                # Set arguments
                cmd = ['parmchk2',
                       '-f', 'prepi',
                       '-s', atomtype,
                       '-i', 'mol.prepi',
                       '-o', 'mol.frcmod',
                       '-a', 'Y']

                # Run parmchk2
                with TemporaryFile() as stream:
                    if subprocess.call(cmd, cwd=tmpdir, stdout=stream, stderr=stream) != 0:
                        raise RuntimeError('"parmchk2" failed')
                    stream.seek(0)
                    for line in stream.readlines():
                        logger.debug(line)

                # Check if antechamber did changes in atom names (and suggest the user to fix the names)
                acmol = Molecule(os.path.join(tmpdir, 'NEWPDB.PDB'), type='pdb')
                acmol.name = np.array([n.upper() for n in acmol.name]).astype(np.object)
                changed_mol_acmol = np.setdiff1d(renamed_mol.name, acmol.name)
                changed_acmol_mol = np.setdiff1d(acmol.name, renamed_mol.name)
                if len(changed_mol_acmol) != 0 or len(changed_acmol_mol) != 0:
                    raise RuntimeError('Initial atom names {} were changed by antechamber to {}. '
                                       'This probably means that the start of the atom name does not match '
                                       'element symbol. '
                                       'Please check the molecule.'
                                       ''.format(','.join(changed_mol_acmol), ','.join(changed_acmol_mol)))

                # Read the results
                prm = parmed.amber.AmberParameterSet(os.path.join(tmpdir, 'mol.frcmod'))
                names, atomtypes, charges, impropers = readPREPI(renamed_mol, os.path.join(tmpdir, 'mol.prepi'))
                masses, elements = readFRCMOD(atomtypes, os.path.join(tmpdir, 'mol.frcmod'))

            elif method == 'CGenFF_2b6':

                from htmd.parameterization.readers import readRTF

                # Write the molecule to a file
                renamed_mol.write(os.path.join(tmpdir, 'mol.pdb'))

                # Set arguments
                cmd = ['match-typer',
                       '-charge', str(netcharge),
                       '-forcefield', 'top_all36_cgenff_new',
                       'mol.pdb']

                # Run match-type
                with TemporaryFile() as stream:
                    if subprocess.call(cmd, cwd=tmpdir, stdout=stream, stderr=stream) != 0:
                        raise RuntimeError('"match-typer" failed')
                    stream.seek(0)
                    for line in stream.readlines():
                        logger.debug(line)

                prm = parmed.charmm.CharmmParameterSet(os.path.join(tmpdir, 'mol.rtf'), os.path.join(tmpdir, 'mol.prm'))
                names, elements, atomtypes, charges, masses, impropers = readRTF(os.path.join(tmpdir, 'mol.rtf'))

            else:
                raise ValueError('Invalid method {}'.format(method))

        assert np.all(renamed_mol.name == names)

    assert np.all(mol.element == elements)

    mol = mol.copy()
    mol.atomtype = atomtypes
    mol.masses = masses
    mol.impropers = impropers
    if acCharges is not None:
        mol.charge = charges

    return prm, mol


class TestFftypeGAFF(unittest.TestCase):

    def setUp(self):
        from htmd.home import home
        self.refDir = home(dataDir='test-fftype')

    def assertListAlmostEqual(self, list1, list2, places=7):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def _init_mol(self, molName, ffTypeMethod, chargetuple):

        from htmd.molecule.molecule import Molecule

        molFile = os.path.join(self.refDir, '{}.mol2'.format(molName))

        if chargetuple == 'None':
            acCharges = None
            netcharge = None
        else:
            acCharges = chargetuple[0]
            netcharge = chargetuple[1]

        mol = Molecule(molFile)

        with TemporaryDirectory() as tmpDir:
            self.testParameters, self.testMolecule = fftype(mol,
                                                            method=ffTypeMethod,
                                                            acCharges=acCharges,
                                                            netcharge=netcharge,
                                                            tmpDir=tmpDir)
            self.testIntermediaryFiles = sorted(os.listdir(tmpDir))

    def _generate_references(self, name, method):

        import numbers
        import pickle

        def mapping(value):
            if isinstance(value, str):
                return '\'{}\''.format(value)
            elif isinstance(value, numbers.Real) or isinstance(value, numbers.Integral):
                return str(value)
            elif isinstance(value, np.ndarray):
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
            pickle.dump(self.testParameters, outfile)

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

        import pickle
        import yaml

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

                with open(os.path.join(refDir, 'mol_props.yaml')) as infile:
                    refProps = yaml.load(infile)
                self._test_mol_props(**refProps)

                with open(os.path.join(refDir, 'params.p'), 'rb') as infile:
                    refParams = pickle.load(infile)
                self._test_params(refParams)

                with open(os.path.join(refDir, 'intermediary_files.yaml')) as infile:
                    refFiles = yaml.load(infile)
                self._test_intermediary_files(refFiles)

    def test_broken_atomname(self):

        from htmd.molecule.molecule import Molecule

        molFile = os.path.join(self.refDir, 'ethanolamine_wrongnames.mol2')

        mol = Molecule(molFile, guessNE=['bonds'], guess=[])

        with self.assertRaises(RuntimeError) as cm:
            fftype(mol, method='GAFF2')

        self.assertEqual(cm.exception.args[0], 'Initial atom names BR1 were changed by antechamber to CR1. '
                                               'This probably means that the start of the atom name does not '
                                               'match element symbol. Please check the molecule.')


class TestFftypeCGenFF(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None

        from htmd.home import home
        from htmd.molecule.molecule import Molecule

        molFile = os.path.join(home('building-protein-ligand'), 'benzamidine.mol2')
        self.mol = Molecule(molFile, guessNE=['bonds'], guess=['angles', 'dihedrals'])

    def test_rtf_prm(self):

        from htmd.home import home
        from htmd.parameterization.writers import writeRTF, writePRM

        refDir = home(dataDir='test-fftype/benzamidine')
        with TemporaryDirectory() as resDir:
            parameters, mol = fftype(self.mol, method='CGenFF_2b6')
            writeRTF(mol, parameters, 0, os.path.join(resDir, 'cgenff.rtf'))
            writePRM(mol, parameters, os.path.join(resDir, 'cgenff.prm'))

            for testFile in ('cgenff.rtf', 'cgenff.prm'):
                with self.subTest(testFile=testFile):
                    # Get rid of the first linw with HTMD version string
                    with open(os.path.join(refDir, testFile)) as refFile:
                        refData = refFile.readlines()[1:]
                    with open(os.path.join(resDir, testFile)) as resFile:
                        resData = resFile.readlines()[1:]
                    self.assertListEqual(refData, resData, msg=testFile)

    def test_tmp_files(self):

        with TemporaryDirectory() as tmpDir:
            _, _ = fftype(self.mol, method='CGenFF_2b6', tmpDir=tmpDir)
            self.assertListEqual(sorted(os.listdir(tmpDir)), ['mol.pdb', 'mol.prm', 'mol.rtf', 'top_mol.rtf'])


if __name__ == '__main__':

    unittest.main(verbosity=2)