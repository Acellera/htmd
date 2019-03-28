# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import tempfile

import propka
import propka.lib
import propka.molecular_container

from htmd.builder.preparationdata import PreparationData

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
    rd : DataFrame
        A pandas dataframe containing resid, chain, pKa arrays.
        Properties unrelated to pKa are unset.

    See Also
    --------
    The ResidueData object.


    Examples
    --------
    >>> m=Molecule("3ptb")
    >>> rd = pka(m)
    >>> rd.pKa[rd.resid==189]
    39    4.949079
    Name: pKa, dtype: float64
    """

    return _pka_backend_internal(mol, pH)



def _pka_backend_internal(mol, pH):
    tmpmol = mol.copy()
    tmpmol.filter("noh")
    pka_pdb = tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb")
    tmpmol.write(pka_pdb.name)

    pka_options, _ = propka.lib.loadOptions('--pH', pH, '--quiet')
    pka_molecule = propka.molecular_container.Molecular_container(pka_pdb.name, pka_options)
    pka_pdb.close()

    # calculating pKa values for ionizable residues -
    pka_molecule.calculate_pka()

    rd = PreparationData()
    rd._importPKAs(pka_molecule)
    return rd.data



def _pka_backend_external(molecule, propka=None, pH=7.0):
    '''PKA - Predict per-residue pKa, using propka

    Parameters
    ----------
    molecule : A Molecule object
    propka   : Path of propka program. Default: find automatically
    pH       : pH at which to perform calculation. Default: 7.0

    Return
    ------

    List of dicts containin the keys 'resname', 'resid', 'chainid', 'pka'

    Example
    -------

    pkas = pka( Molecule('4DFR'), pH=6.5 )
    '''

    from tempfile import mkstemp
    import os
    from shutil import which
    from subprocess import check_output

    try:
        propka = which("propka31", mode=os.X_OK)
    except:
        raise NameError("Cannot find 'propka31' in PATH");

    if not os.access(propka, os.X_OK):
        raise NameError("'propka' not found");

    try:
        pH = float(pH)
    except:
        raise NameError("pH value invalid")
    if (pH <= 0. or pH >= 14.):
        raise NameError("pH value invalid")

    fh, fnp = mkstemp()
    fn = fnp + ".pdb"
    molecule.write(fn)

    dd = os.getcwd()
    try:
        os.chdir(os.path.dirname(fn))
        op = check_output([propka, fn, "-o", str(pH)])
        op = op.decode("ascii").split("\n");
    except:
        raise NameError("Failed to execute PropKa")
    os.chdir(dd)

    ret = []
    # print(op)
    atpka = False
    for line in op:
        if not atpka:
            if "model-pKa" in line:
                atpka = True
                continue
            else:
                continue
            if "-" in line:
                atpka = False
                continue
        print(line)
        s = line.split()
        if len(s) >= 4:
            pka = {"resname": s[0], "resid": s[1], "chainid": s[2], "pka": float(s[3])}
            ret.append(pka)
    os.unlink(fn)
    os.unlink(fnp + ".pka")
    os.unlink(fnp + ".propka_input")
    return ret




# A test method
if __name__ == "__main__":
    import doctest
    from moleculekit.molecule import Molecule
    doctest.testmod()

