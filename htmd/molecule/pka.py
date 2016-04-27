import logging
import tempfile
import propka
import propka.lib
import propka.molecular_container

from htmd.molecule.molecule import Molecule
from htmd.proteinpreparation.residuedata import ResidueData


logger = logging.getLogger(__name__)


def pka(mol, pH=7.0):
    """Compute pKa values via propKa 3.1

    Parameters
    ----------
    mol : Molecule
        The molecule to be analyzed
    pH : float
        The pH at which to do the analysis

    Returns
    -------
    rd : ResidueData
        An object (see), containing resid, chain, pKa arrays.
        Properties unrelated to pKa are unset.

    Examples
    --------
    >>> m=Molecule("3ptb")
    >>> rd = pka(m,pH=7.0)
    >>> rd.pKa[rd.resid == 189].round(2)
    array([ 4.95])
    """

    return _pka_backend(mol, pH)


def _pka_backend(m, pH):
    """Internal function - may allow switching between internal and external propka if necessary """
    return _pka_backend_internal(m, pH)


def _pka_backend_internal(mol, pH):
    tmpmol = mol.copy()
    tmpmol.filter("noh")
    pka_pdb = tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb")
    tmpmol.write(pka_pdb.name)

    pka_options, _ = propka.lib.loadOptions('--pH', pH)  # add '--no-print' when released
    pka_molecule = propka.molecular_container.Molecular_container(pka_pdb.name, pka_options)
    pka_pdb.close()

    # calculating pKa values for ionizable residues -
    pka_molecule.calculate_pka()

    rd = ResidueData()
    rd._setPKAs(pka_molecule)
    return rd


# A test method
if __name__ == "__main__":
    import doctest
    from htmd.molecule.molecule import Molecule
    doctest.testmod()

