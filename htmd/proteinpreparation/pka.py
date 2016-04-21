import logging
import tempfile
import numpy as np
import propka
import propka.lib

# If necessary: http://stackoverflow.com/questions/16981921/relative-imports-in-python-3
from htmd.proteinpreparation.residuedata import ResidueData


from htmd.molecule.molecule import Molecule
import propka.molecular_container


logger = logging.getLogger(__name__)


def pka(mol, pH=7.0, pka_options=None):
    """Compute pKa values via propKa3.1

    Parameters
    ----------

    Examples
    --------
    >> m=Molecule("3ptb")
    >> rd = pka(m,pH=7.0)

    """

    tmpmol = mol.copy()
    tmpmol.filter("noh")
    pka_pdb = tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb")
    tmpmol.write(pka_pdb.name)

    if not pka_options:
        pka_options, _ = propka.lib.loadOptions()
    pka_options.pH = pH

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

    # doctest.testmod()

    m = Molecule("3ptb")
    rd = pka(m, pH=7.0)
    pass
