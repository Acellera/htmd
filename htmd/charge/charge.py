# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import os
import subprocess
from tempfile import TemporaryDirectory, TemporaryFile

import numpy as np

from htmdmol.molecule import Molecule

logger = logging.getLogger(__name__)


def fitGasteigerCharges(mol):
    """
    Fit Gasteiger atomic charges

    Parameters
    ----------
    mol: Molecule
        Molecule to fit the charges

    Return
    ------
    results: Molecule
        Copy of the molecule with the charges set

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmdmol.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = 0

    >>> new_mol = fitGasteigerCharges(mol)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.411509...,  0.205754...,  0.205754...], dtype=float32)
    """

    from rdkit.Chem.rdmolfiles import MolFromMol2File
    from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    # Set atom types to elements, overwise rdkit refuse to read a MOL2 file
    htmd_mol = mol.copy()
    htmd_mol.atomtype = htmd_mol.element

    # Convert Molecule to rdkit Mol
    with TemporaryDirectory() as tmpDir:
        filename = os.path.join(tmpDir, 'mol.mol2')
        htmd_mol.write(filename)
        rdkit_mol = MolFromMol2File(filename, removeHs=False)
    assert mol.numAtoms == rdkit_mol.GetNumAtoms()

    # Compute and store Gasteiger charges
    ComputeGasteigerCharges(rdkit_mol, throwOnParamFailure=True)
    mol = mol.copy()
    mol.charge[:] = [atom.GetDoubleProp('_GasteigerCharge') for atom in rdkit_mol.GetAtoms()]

    return mol


def fitChargesWithAntechamber(mol, type='gas', molCharge=None):
    """
    Fit atomic charges with Antechamber

    Parameters
    ----------
    mol: Molecule
        Molecule to fit the charges
    type: str
        Charge type
    molCharge: int
        Molecular charge

    Return
    ------
    results: Molecule
        Copy of the molecule with the charges set

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmdmol.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = 0

    >>> new_mol = fitChargesWithAntechamber(mol)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.411518...,  0.205759...,  0.205759...], dtype=float32)

    >>> new_mol = fitChargesWithAntechamber(mol, type='gas')
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.411518...,  0.205759...,  0.205759...], dtype=float32)

    >>> new_mol = fitChargesWithAntechamber(mol, type='gas', molCharge=10)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.411518...,  0.205759...,  0.205759...], dtype=float32)

    >>> new_mol = fitChargesWithAntechamber(mol, type='bcc')
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.785...,  0.39...,  0.39...], dtype=float32)

    >>> new_mol = fitChargesWithAntechamber(mol, type='bcc', molCharge=0)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.785...,  0.39...,  0.39...], dtype=float32)
    """

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    if type not in ('gas', 'bcc'):
        raise ValueError('"type" has to be "gas" or "bcc"')

    if molCharge is None:
        molCharge = int(round(np.sum(mol.charge)))
        logger.info('Using partial atomic charges to calculate molecular charge')

    mol = mol.copy()
    mol.charge[:] = molCharge/mol.numAtoms

    with TemporaryDirectory() as tmpDir:
        old_name = os.path.join(tmpDir, 'old.mol2')
        new_name = os.path.join(tmpDir, 'new.mol2')

        mol.write(old_name)

        cmd = ['antechamber',
               '-fi', 'mol2', '-i', old_name,
               '-nc', str(molCharge),
               '-c', type,
               '-fo', 'mol2', '-o', new_name]

        with TemporaryFile() as stream:
            if subprocess.call(cmd, cwd=tmpDir, stdout=stream, stderr=stream) != 0:
                raise RuntimeError('"antechamber" failed')
            stream.seek(0)
            for line in stream.readlines():
                logger.debug(line)

        mol.charge[:] = Molecule(new_name).charge

    return mol


def fitESPCharges(mol, qm, outdir, apply_bounds=True, restraint_factor=0, fixed=()):
    """
    Fit ESP atomic charges

    Parameters
    ----------
    mol: Molecule
        Molecule to fit the charges
    qm: Psi4
        Psi4 instance for QM calculations
    outdir: str
        Output directory for the QM calculations
    apply_bounds: boolean
        Apply bounds to atomic charges
    restraint_factor: float
        Restraint factor for heavy elements
    fixed : list of ints
        List of fixed atom indices

    Return
    ------
    mol: Molecule
        Copy of the molecule with the charges set
    extra: dict
        extra['qm_dipole']
            QM dipole moments
        extra['esp_loss']
            ESP fitting loss

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmdmol.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = 0

    >>> from tempfile import TemporaryDirectory
    >>> from htmd.qm import Psi4
    ffevaluate module is in beta version

    >>> np.random.seed(20181113)
    >>> with TemporaryDirectory() as tmpDir:
    ...     new_mol, extra = fitESPCharges(mol, Psi4(), tmpDir)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.3908...,  0.1954...,  0.1954...], dtype=float32)
    """

    from htmd.qm.base import QMBase
    from htmd.charge.esp import MoleculeGrid, ESP

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    if not issubclass(type(qm), QMBase):
        raise ValueError('"qm" has to be instance of {}'.format(QMBase))

    apply_bounds = bool(apply_bounds)
    restraint_factor = float(restraint_factor)

    # Get ESP points
    point_file = os.path.join(outdir, "00000", "grid.dat")
    if os.path.exists(point_file):
        # Load a point file if one exists from a previous job
        logger.info('Reusing ESP grid from %s' % point_file)
        esp_points = np.loadtxt(point_file)
    else:
        # Generate ESP points
        logger.info('Generating ESP grid')
        esp_points = MoleculeGrid(mol).getPoints()

    # Run QM simulation
    qm.molecule = mol
    qm.esp_points = esp_points
    qm.optimize = False
    qm.restrained_dihedrals = None
    qm.directory = outdir
    qm_results = qm.run()
    if qm_results[0].errored:
        raise RuntimeError('QM calculation failed! Check logs at {}'.format(outdir))

    # Safeguard QM code from changing coordinates :)
    assert np.all(np.isclose(mol.coords, qm_results[0].coords, atol=1e-6))

    # Fit ESP charges
    esp = ESP()
    esp.molecule = mol
    esp.qm_results = qm_results
    esp.apply_bounds = apply_bounds
    esp.restraint_factor = restraint_factor
    esp.fixed = fixed
    esp_result = esp.run()

    # Update the charges
    mol = mol.copy()
    mol.charge[:] = esp_result['charges']
    extra = {'qm_dipole': qm_results[0].dipole,
             'esp_loss': esp_result['loss'],
             'esp_rmsd': esp_result['RMSD']}

    return mol, extra

def symmetrizeCharges(mol):
    """
    Average the charges of equivalent atoms

    Parameters
    ----------
    mol: Molecule
        Molecule to symmetrize the charges

    Return
    ------
    results: Molecule
        Copy of the molecule with the symmetrized charges

    Examples
    --------
    >>> from htmd.home import home
    >>> from htmdmol.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = [0.5, -0.5, 0.0]

    >>> new_mol = symmetrizeCharges(mol)
    >>> assert new_mol is not mol
    >>> new_mol.charge
    array([ 0.5 , -0.25, -0.25], dtype=float32)
    """

    from htmd.parameterization.detect import detectEquivalentAtoms

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    mol = mol.copy()
    molCharge = np.sum(mol.charge)

    equivalent_groups, _, _ = detectEquivalentAtoms(mol)
    for atoms in equivalent_groups:
        atoms = list(atoms)
        mol.charge[atoms] = np.mean(mol.charge[atoms])

    assert np.isclose(np.sum(mol.charge), molCharge)

    return mol


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)