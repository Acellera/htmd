# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import os
from tempfile import TemporaryDirectory

import numpy as np

from htmd.molecule.molecule import Molecule

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
    >>> from htmd.molecule.molecule import Molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = 0

    >>> new_mol = fitGasteigerCharges(mol)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.411509...,  0.205754...,  0.205754...], dtype=float32)
    """

    from rdkit.Chem.rdmolfiles import MolFromMol2File
    from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges

    from htmd.parameterization.util import guessElementFromName

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    # Guess and set elements, overwise rdkit refuse to read a MOL2 file
    htmd_mol = mol.copy()
    htmd_mol.atomtype[:] = [guessElementFromName(name) for name in htmd_mol.name]

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


def fitESPCharges(mol, qm, outdir, fixed=()):
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
    >>> from htmd.molecule.molecule import Molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O.mol2')
    >>> mol = Molecule(molFile)
    >>> mol.charge[:] = 0

    >>> from tempfile import TemporaryDirectory
    >>> from htmd.qm import Psi4
    ffevaluate module is in beta version
    >>> with TemporaryDirectory() as tmpDir:
    ...     new_mol, extra = fitESPCharges(mol, Psi4(), tmpDir)
    >>> assert new_mol is not mol
    >>> new_mol.charge # doctest: +ELLIPSIS
    array([-0.394059...,  0.197029...,  0.197029...], dtype=float32)
    """

    from htmd.qm import Psi4
    from htmd.charge.esp import ESP

    if not isinstance(mol, Molecule):
        raise TypeError('"mol" has to be instance of {}'.format(Molecule))
    if mol.numFrames != 1:
        raise ValueError('"mol" can have just one frame, but it has {}'.format(mol.numFrames))

    if not isinstance(qm, Psi4):
        raise ValueError('"qm" has to be instance of {}'.format(Psi4))

    # Get ESP points
    point_file = os.path.join(outdir, "00000", "grid.dat")
    if os.path.exists(point_file):
        # Load a point file if one exists from a previous job
        logger.info('Reusing ESP grid from %s' % point_file)
        esp_points = np.loadtxt(point_file)
    else:
        # Generate ESP points
        logger.info('Generating ESP grid')
        esp_points = ESP.generate_points(mol)[0]

    # Run QM simulation
    qm.molecule = mol
    qm.esp_points = esp_points
    qm.optimize = False
    qm.restrained_dihedrals = None
    qm.directory = outdir
    qm_results = qm.run()
    if qm_results[0].errored:
        raise RuntimeError('\nQM calculation failed! Check logs at %s\n' % espDir)

    # Safeguard QM code from changing coordinates :)
    assert np.all(np.isclose(mol.coords, qm_results[0].coords, atol=1e-6))

    # Fit ESP charges
    esp = ESP()
    esp.molecule = mol
    esp.qm_results = qm_results
    esp.fixed = fixed
    esp_result = esp.run()

    # Update the charges
    mol = mol.copy()
    mol.charge[:] = esp_result['charges']
    extra = {'qm_dipole': qm_results[0].dipole,
             'esp_loss': esp_result['loss']}

    return mol, extra


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)