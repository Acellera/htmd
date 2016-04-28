import logging
import tempfile
import numpy as np
import random

# If necessary: http://stackoverflow.com/questions/16981921/relative-imports-in-python-3
from htmd.proteinpreparation.residuedata import ResidueData
from htmd.proteinpreparation.pdb2pqr.src.pdbParser import readPDB
from htmd.proteinpreparation.pdb2pqr.main import runPDB2PQR

from htmd.molecule.molecule import Molecule

# Tried to make runs reproducible, but does not work
random.seed(2016)

logger = logging.getLogger(__name__)


def _createMolecule(name, resname, chain, resid, insertion, coords, segid, elements):
    numAtoms = len(name)
    mol = Molecule()
    mol.empty(numAtoms)

    mol.name = np.array(name, dtype=mol._pdb_fields['name'])
    mol.resname = np.array(resname, dtype=mol._pdb_fields['resname'])
    mol.chain = np.array(chain, dtype=mol._pdb_fields['chain'])
    mol.resid = np.array(resid, dtype=mol._pdb_fields['resid'])
    mol.insertion = np.array(insertion, dtype=mol._pdb_fields['insertion'])
    mol.coords = np.array(np.atleast_3d(np.vstack(coords)), dtype=mol._pdb_fields['coords'])
    mol.segid = np.array(segid, dtype=mol._pdb_fields['segid'])
    mol.element = np.array(elements, dtype=mol._pdb_fields['element'])
    return mol


def prepareProtein(mol_in,
                   pH=7.0,
                   verbose=0,
                   returnDetails=False,
                   keep=None):
    """A system preparation wizard for HTMD. 

    Returns a Molecule object, where residues have been renamed to follow
    internal conventions on protonation (below). Coordinates are changed to
    optimize the H-bonding network. This is very preliminar and should be
    roughly equivalent to mdweb and Maestro's wizard.

    The following residue names are used in the returned molecule:

        ASH 	Neutral ASP
        CYX 	SS-bonded CYS
        CYM 	Negative CYS
        GLH 	Neutral GLU
        HIP 	Positive HIS
        HID 	Neutral HIS, proton HD1 present
        HIE 	Neutral HIS, proton HE2 present
        LYN 	Neutral LYS
        TYM 	Negative TYR
        AR0     Neutral ARG


    Features
    --------
     - assign protonation states via propKa
     - flip residues to optimize H-bonding network
     - debump collisions
     - fill-in missing atoms, e.g. hydrogen atoms


    Parameters
    ----------
    mol_in : htmd.Molecule
        the object to be optimized
    pH : float
        pH to decide titration
    verbose : int
        verbosity
    returnDetails : bool
        whether to return just the prepared Molecule (False, default) or a molecule *and* a ResidueInfo
        object including computed properties
    keep : bool
        TODO


    Returns
    -------
    mol_out : Molecule
        the molecule titrated and optimized. The molecule object contains an additional attribute,
    resdata_out : ResidueData
        a table of residues with the corresponding protonation states, pKas, and other information


    Examples
    --------
    >> tryp = Molecule('3PTB')
    >> tryp_op = prepareProtein(tryp, pH=1.0)
    >> tryp_op.write('3PTB-opt-ph1.pdb')


    Unsupported/To Do/To Check
    --------------------------
     - ligands
     - termini
     - force residues
     - multiple chains
     - nucleic acids
     - reporting in machine-readable form
     - coupled titrating residues
     - Disulfide bridge detection (implemented but unused)

    """

    oldLoggingLevel = logger.level
    if verbose:
        logger.setLevel(logging.DEBUG)
    logger.info("Starting.")

    # We could transform the molecule into an internal object, but for now I prefer to rely on the strange
    # internal parser to avoid hidden quirks.
    tmpin = tempfile.NamedTemporaryFile(suffix=".pdb", mode="w+")
    logger.info("Temporary file is " + tmpin.name)
    mol_in.write(tmpin.name)  # Not sure this is sound unix

    pdblist, errlist = readPDB(tmpin)
    if len(pdblist) == 0 and len(errlist) == 0:
        raise Exception('Internal error in preparing input to pdb2pqr')

    # We could set additional options here
    import propka.lib
    propka_opts, dummy = propka.lib.loadOptions()
    propka_opts.verbose = verbose

    # Note on naming. The behavior of PDB2PQR is controlled by two parameters, ff and ffout. My understanding is
    # that the ff parameter sets which residues are SUPPORTED by the underlying FF, PLUS the charge and radii.
    # The ffout parameter sets the naming scheme. Therefore, I want ff to be as general as possible, which turns out
    # to be "parse". Then I pick a convenient ffout.

    # Relying on defaults
    header, lines, missedLigands, pdb2pqr_protein = runPDB2PQR(pdblist,
                                                               ph=pH, verbose=verbose,
                                                               ff="parse", ffout="amber",
                                                               propkaOptions=propka_opts)
    tmpin.close()

    # Diagnostics
    for missedligand in missedLigands:
        logger.warning("The following residue has not been optimized: " + missedligand)

    # Here I parse the returned protein object and recreate a Molecule,
    # because I need to access the properties.
    logger.info("Building Molecule object.")
    mol_out = Molecule()

    name = []
    resid = []
    chain = []
    insertion = []
    coords = []
    resname = []
    segids = []
    elements = []

    resdata_out = ResidueData()

    for residue in pdb2pqr_protein.residues:
        if 'ffname' in residue.__dict__:
            curr_resname = residue.ffname
            if len(curr_resname) >= 4:
                curr_resname = curr_resname[-3:]
                logger.info("Residue %s has internal name %s, replacing with %s" %
                            (residue, residue.ffname, curr_resname))
        else:
            curr_resname = residue.name

        resdata_out._setProtonation(residue, curr_resname)

        if 'patches' in residue.__dict__:
            for patch in residue.patches:
                resdata_out._appendPatches(residue, patch)
                if patch != "PEPTIDE":
                    logger.info("Residue %s has patch %s set" % (residue, patch))

        for atom in residue.atoms:
            name.append(atom.name)
            resid.append(residue.resSeq)
            chain.append(residue.chainID)
            insertion.append(residue.iCode)
            coords.append([atom.x, atom.y, atom.z])
            resname.append(curr_resname)
            segids.append(atom.segID)
            elements.append(atom.element)

    mol_out = _createMolecule(name, resname, chain, resid, insertion, coords, segids, elements)

    resdata_out._importPKAs(pdb2pqr_protein.pka_molecule)
    resdata_out.pdb2pqr_protein = pdb2pqr_protein
    resdata_out.pka_protein = pdb2pqr_protein.pka_molecule
    resdata_out.pka_dict = pdb2pqr_protein.pkadic
    resdata_out.missedLigands = missedLigands

    logger.info("Returning.")
    logger.setLevel(oldLoggingLevel)

    if returnDetails:
        return mol_out, resdata_out
    else:
        return mol_out


# A test method
if __name__ == "__main__":
    tryp = Molecule('3PTB')

    tryp_op = prepareProtein(tryp, pH=1.0)
    tryp_op.write('proteinpreparation-test-main-ph-1.pdb')

    tryp_op = prepareProtein(tryp, pH=14.0)
    tryp_op.write('proteinpreparation-test-main-ph-14.pdb')

    tryp_op, prepData = prepareProtein(tryp, returnDetails=True)
    tryp_op.write('proteinpreparation-test-main-ph-7.pdb')
    print(prepData)

    mol = Molecule("1r1j")
    mo, prepData = prepareProtein(mol, returnDetails=True)
    his = prepData.resname == "HIS"  # Has to be checked better, due to Zn++
    list(zip(prepData.protonation[his], prepData.resid[his]))
    pass
