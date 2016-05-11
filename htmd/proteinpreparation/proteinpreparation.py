import logging
import tempfile
import numpy as np
import random

# If necessary: http://stackoverflow.com/questions/16981921/relative-imports-in-python-3
from htmd.proteinpreparation.residuedata import ResidueData
from htmd.proteinpreparation.pdb2pqr.src.pdbParser import readPDB
from htmd.proteinpreparation.pdb2pqr.main import runPDB2PQR

from htmd.molecule.molecule import Molecule

logger = logging.getLogger(__name__)




def _fillMolecule(name, resname, chain, resid, insertion, coords, segid, elements):
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
                   hydrophobicThickness=None,
                   keep=None):
    """A system preparation wizard for HTMD.

    Returns a Molecule object, where residues have been renamed to follow
    internal conventions on protonation (below). Coordinates are changed to
    optimize the H-bonding network. This should be roughly equivalent to mdweb and Maestro's
    preparation wizard.

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

    If hydrophobicThickness is set to a positive value 2*h, a warning is produced for titratable residues
    having -h<z<h and are buried in the protein by less than 75%. The list of such residues can be accessed setting
    returnDetails to True. Note that the heuristic for the detection of membrane-exposed residues is very crude;
    the "buried fraction" computation (from propka) is approximate; also, in the presence of cavities,
    residues may be solvent-exposed independently from their z location.


    Notes
    -----
    In case of problems, exclude water and other dummy atoms.


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
    hydrophobicThickness : float
        the thickness of the membrane in which the protein is embedded, or None if globular protein.
        Used to provide a warning about membrane-exposed residues.
    keep : bool
        TODO


    Returns
    -------
    mol_out : Molecule
        the molecule titrated and optimized. The molecule object contains an additional attribute,
    resData : ResidueData
        a table of residues with the corresponding protonation states, pKas, and other information


    Examples
    --------
    >>> tryp = Molecule('3PTB')

    >>> tryp_op = prepareProtein(tryp, pH=1.0)
    >>> tryp_op.write('proteinpreparation-test-main-ph-1.pdb')

    >>> tryp_op = prepareProtein(tryp, pH=14.0)
    >>> tryp_op.write('proteinpreparation-test-main-ph-14.pdb')

    >>> tryp_op, prepData = prepareProtein(tryp, returnDetails=True)
    >>> tryp_op.write('proteinpreparation-test-main-ph-7.pdb')
    >>> prepData
    ResidueData object about 287 residues. Please find the full info in the .data property.
      resname  resid insertion chain       pKa protonation    buried    patches
    0     ILE     16               A       NaN         ILE       NaN    [NTERM]
    1     VAL     17               A       NaN         VAL       NaN  [PEPTIDE]
    2     GLY     18               A       NaN         GLY       NaN  [PEPTIDE]
    3     GLY     19               A       NaN         GLY       NaN  [PEPTIDE]
    4     TYR     20               A  9.590845         TYR  0.146429  [PEPTIDE]
     . . .
    >>> prepData.data.to_excel("/tmp/tryp-report.xlsx")

    >>> mol = Molecule("1r1j")
    >>> mo, prepData = prepareProtein(mol, returnDetails=True)
    >>> prepData.missedLigands
    ['NAG', 'ZN', 'OIR']

    >>> his = prepData.data.resname == "HIS"
    >>> prepData.data[his][["resid","insertion","chain","resname","protonation"]]
         resid insertion chain resname protonation
    160    214               A     HIS         HID
    163    217               A     HIS         HID
    383    437               A     HIS         HID
    529    583               A     HIS         HID
    533    587               A     HIS         HIP
    583    637               A     HIS         HID
    627    681               A     HIS         HID
    657    711               A     HIS         HIP
    679    733               A     HIS         HID

    >>> mor = Molecule(os.path.join(home(dataDir="mor"), "4dkl.pdb"))
    >>> mor.filter("protein and noh")
    >>> mor_opt, mor_data = prepareProtein(mor, returnDetails=True, hydrophobicThickness=32.0)
    >>> exposedRes = mor_data.data.membraneExposed
    >>> mor_data.data[exposedRes].to_excel("/tmp/mor_exposed_residues.xlsx")

    See Also
    --------
    The ResidueData object.


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
    propka_opts, dummy = propka.lib.loadOptions('--quiet')
    propka_opts.verbosity = verbose
    propka_opts.verbose = verbose  # Will be removed in future propKas

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

    name = []
    resid = []
    chain = []
    insertion = []
    coords = []
    resname = []
    segids = []
    elements = []

    resData = ResidueData()

    for residue in pdb2pqr_protein.residues:
        if 'ffname' in residue.__dict__:
            curr_resname = residue.ffname
            if len(curr_resname) >= 4:
                curr_resname = curr_resname[-3:]
                logger.info("Residue %s has internal name %s, replacing with %s" %
                            (residue, residue.ffname, curr_resname))
        else:
            curr_resname = residue.name

        resData._setProtonationState(residue, curr_resname)

        if 'patches' in residue.__dict__:
            for patch in residue.patches:
                resData._appendPatches(residue, patch)
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

    mol_out = _fillMolecule(name, resname, chain, resid, insertion, coords, segids, elements)

    resData._importPKAs(pdb2pqr_protein.pka_molecule)
    resData.pdb2pqr_protein = pdb2pqr_protein
    resData.pka_dict = pdb2pqr_protein.pkadic
    resData.missedLigands = missedLigands

    if hydrophobicThickness:
        resData._setMembraneExposure(hydrophobicThickness)


    logger.info("Returning.")
    logger.setLevel(oldLoggingLevel)

    if returnDetails:
        return mol_out, resData
    else:
        return mol_out


# A test method
if __name__ == "__main__":
    from htmd import home
    import os


    import doctest
    doctest.testmod()
