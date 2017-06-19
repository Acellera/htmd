# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import tempfile

import numpy as np
import propka.lib
from pdb2pqr.main import runPDB2PQR
from pdb2pqr.src.pdb import readPDB

from htmd.builder.preparationdata import PreparationData
from htmd.molecule.molecule import Molecule

logger = logging.getLogger(__name__)


def _selToHoldList(mol, sel):
    if sel:
        tx = mol.copy()
        tx.filter(sel)
        tx.filter("name CA")
        ret = list(zip(tx.resid, tx.chain, tx.insertion))
    else:
        ret = None
    return ret


def _fillMolecule(name, resname, chain, resid, insertion, coords, segid, element,
                  occupancy, beta, charge, record):
    numAtoms = len(name)
    mol = Molecule()
    mol.empty(numAtoms)

    mol.name = np.array(name, dtype=mol._dtypes['name'])
    mol.resname = np.array(resname, dtype=mol._dtypes['resname'])
    mol.chain = np.array(chain, dtype=mol._dtypes['chain'])
    mol.resid = np.array(resid, dtype=mol._dtypes['resid'])
    mol.insertion = np.array(insertion, dtype=mol._dtypes['insertion'])
    mol.coords = np.array(np.atleast_3d(np.vstack(coords)), dtype=mol._dtypes['coords'])
    mol.segid = np.array(segid, dtype=mol._dtypes['segid'])
    mol.element = np.array(element, dtype=mol._dtypes['element'])
    mol.occupancy = np.array(occupancy, dtype=mol._dtypes['occupancy'])
    mol.beta = np.array(beta, dtype=mol._dtypes['beta'])
    # mol.charge = np.array(charge, dtype=mol._dtypes['charge'])
    # mol.record = np.array(record, dtype=mol._dtypes['record'])
    return mol


def _fixupWaterNames(mol):
    """Rename WAT / OW HW HW atoms as O H1 H2"""
    mol.set("name", "O", sel="resname WAT and name OW")
    mol.set("name", "H1", sel="resname WAT and name HW and serial % 2 == 0")
    mol.set("name", "H2", sel="resname WAT and name HW and serial % 2 == 1")


def _warnIfContainsDUM(mol):
    """Warn if any DUM atom is there"""
    if any(mol.atomselect("resname DUM")):
        logger.warning("OPM's DUM residues must be filtered out before preparation. Continuing, but crash likely.")


def _buildResAndMol(pdb2pqr_protein):
    # Here I parse the returned protein object and recreate a Molecule,
    # because I need to access the properties.
    logger.debug("Building Molecule object.")

    name = []
    resid = []
    chain = []
    insertion = []
    coords = []
    resname = []
    segid = []
    element = []
    occupancy = []
    beta = []
    record = []
    charge = []

    prepData = PreparationData()

    for i, residue in enumerate(pdb2pqr_protein.residues):
        # if 'ffname' in residue.__dict__:
        if getattr(residue, 'ffname', None):
            curr_resname = residue.ffname
            if len(curr_resname) >= 4:
                curr_resname = curr_resname[-3:]
                logger.debug("Residue %s has internal name %s, replacing with %s" %
                             (residue, residue.ffname, curr_resname))
        else:
            curr_resname = residue.name

        prepData._setProtonationState(residue, curr_resname)

        # Removed because not really useful
        # if getattr(residue, 'patches', None):
        #     for patch in residue.patches:
        #         prepData._appendPatches(residue, patch)
        #         if patch != "PEPTIDE":
        #             logger.debug("Residue %s has patch %s set" % (residue, patch))

        if getattr(residue, 'wasFlipped', 'UNDEF') != 'UNDEF':
            prepData._setFlipped(residue, residue.wasFlipped)

        prepData._set(residue, 'pdb2pqr_idx', i)

        for atom in residue.atoms:
            name.append(atom.name)
            resid.append(residue.resSeq)
            chain.append(residue.chainID)
            insertion.append(residue.iCode)
            coords.append([atom.x, atom.y, atom.z])
            resname.append(curr_resname)
            segid.append(atom.segID)
            # Fixup element fields for added H (routines.addHydrogens)
            elt = "H" if atom.added and atom.name.startswith("H") else atom.element
            element.append(elt)
            occupancy.append(0.0 if atom.added else atom.occupancy)
            beta.append(99.0 if atom.added else atom.tempFactor)
            charge.append(atom.charge)
            record.append(atom.type)
            if atom.added:
                logger.debug("Coordinates of atom {:s} in residue {:s} were guessed".format(residue.__str__(),atom.name))
                prepData._set(residue, 'guessedAtoms', atom.name, append=True)


    mol_out = _fillMolecule(name, resname, chain, resid, insertion, coords, segid, element,
                            occupancy, beta, charge, record)
    # mol_out.set("element", " ")
    # Re-calculating elements
    mol_out.element[:] = ''
    mol_out.element = mol_out._guessMissingElements()

    prepData._importPKAs(pdb2pqr_protein.pka_protein)

    return mol_out, prepData


def proteinPrepare(mol_in,
                   pH=7.0,
                   verbose=0,
                   returnDetails=False,
                   hydrophobicThickness=None,
                   holdSelection=None):
    """A system preparation wizard for HTMD.

    Returns a Molecule object, where residues have been renamed to follow
    internal conventions on protonation (below). Coordinates are changed to
    optimize the H-bonding network. This should be roughly equivalent to mdweb and Maestro's
    preparation wizard.

    The following residue names are used in the returned molecule:

    === ===============================
    ASH Neutral ASP
    CYX SS-bonded CYS
    CYM Negative CYS
    GLH Neutral GLU
    HIP Positive HIS
    HID Neutral HIS, proton HD1 present
    HIE Neutral HIS, proton HE2 present
    LYN Neutral LYS
    TYM Negative TYR
    AR0 Neutral ARG
    === ===============================

    ========= ======= =========
    Charge +1 Neutral Charge -1
    ========= ======= =========
    -         ASH     ASP
    -         CYS     CYM
    -         GLH     GLU
    HIP       HID/HIE -
    LYS       LYN     -
    -         TYR     TYM
    ARG       AR0     -
    ========= ======= =========

    A detailed table about the residues modified is returned (as a second return value) when
    returnDetails is True (see PreparationData object).

    If hydrophobicThickness is set to a positive value 2*h, a warning is produced for titratable residues
    having -h<z<h and are buried in the protein by less than 75%. Note that the heuristic for the
    detection of membrane-exposed residues is very crude; the "buried fraction" computation
    (from propKa) is approximate; also, in the presence of cavities,
    residues may be solvent-exposed independently from their z location.


    Notes
    -----
    In case of problems, exclude water and other dummy atoms.

    Features:
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
    holdSelection : str
        (Untested) Atom selection to be excluded from optimization.
        Only the carbon-alpha atom will be considered for the corresponding residue.


    Returns
    -------
    mol_out : Molecule
        the molecule titrated and optimized. The molecule object contains an additional attribute,
    resData : PreparationData
        a table of residues with the corresponding protonation states, pKas, and other information


    Examples
    --------
    >>> tryp = Molecule('3PTB')

    >>> tryp_op, prepData = proteinPrepare(tryp, returnDetails=True)
    >>> tryp_op.write('proteinpreparation-test-main-ph-7.pdb')
    >>> prepData.data.to_excel("/tmp/tryp-report.xlsx")
    >>> prepData                                                        # doctest: +NORMALIZE_WHITESPACE
    PreparationData object about 290 residues.
    Unparametrized residue names: CA, BEN
    Please find the full info in the .data property, e.g.:
      resname  resid insertion chain       pKa protonation flipped     buried
    0     ILE     16               A       NaN         ILE     NaN        NaN
    1     VAL     17               A       NaN         VAL     NaN        NaN
    2     GLY     18               A       NaN         GLY     NaN        NaN
    3     GLY     19               A       NaN         GLY     NaN        NaN
    4     TYR     20               A  9.590845         TYR     NaN  14.642857
     . . .
    >>> x_HIE91_ND1 = tryp_op.get("coords","resid 91 and  name ND1")
    >>> x_SER93_H =   tryp_op.get("coords","resid 93 and  name H")
    >>> len(x_SER93_H) == 3
    True
    >>> np.linalg.norm(x_HIE91_ND1-x_SER93_H) < 3
    True

    >>> tryp_op = proteinPrepare(tryp, pH=1.0)
    >>> tryp_op.write('proteinpreparation-test-main-ph-1.pdb')

    >>> tryp_op = proteinPrepare(tryp, pH=14.0)
    >>> tryp_op.write('proteinpreparation-test-main-ph-14.pdb')

    >>> mol = Molecule("1r1j")
    >>> mo, prepData = proteinPrepare(mol, returnDetails=True)
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

    >>> mor = Molecule("4dkl")
    >>> mor.filter("protein and noh")
    >>> mor_opt, mor_data = proteinPrepare(mor, returnDetails=True,
    ...                                    hydrophobicThickness=32.0)
    >>> exposedRes = mor_data.data.membraneExposed
    >>> mor_data.data[exposedRes].to_excel("/tmp/mor_exposed_residues.xlsx")

    >>> im=Molecule("4bkj")
    >>> imo,imd=proteinPrepare(im,returnDetails=True)
    >>> imd.data.to_excel("/tmp/imatinib_report.xlsx")

    See Also
    --------
    :class:`htmd.builder.PreparationData.PreparationData`

    Notes
    -----
    Unsupported/To Do/To Check:
     - ligands
     - termini
     - multiple chains
     - nucleic acids
     - coupled titrating residues
     - Disulfide bridge detection (implemented but unused)
    """

    oldLoggingLevel = logger.level
    if verbose:
        logger.setLevel(logging.DEBUG)
    logger.debug("Starting.")

    _warnIfContainsDUM(mol_in)

    # We could transform the molecule into an internal object, but for
    # now I prefer to rely on the strange internal parser to avoid
    # hidden quirks.
    tmpin = tempfile.NamedTemporaryFile(suffix=".pdb", mode="w+")
    logger.debug("Temporary file is " + tmpin.name)
    mol_in.write(tmpin.name)  # Not sure this is good unix

    pdblist, errlist = readPDB(tmpin)
    if len(pdblist) == 0 and len(errlist) == 0:
        raise Exception('Internal error in preparing input to pdb2pqr')

    # An ugly hack to silence non-prefixed logging messages
    for h in propka.lib.logger.handlers:
        if h.formatter._fmt == '%(message)s':
            propka.lib.logger.removeHandler(h)

    propka_opts, dummy = propka.lib.loadOptions('--quiet')
    propka_opts.verbosity = verbose
    propka_opts.verbose = verbose  # Will be removed in future propKas

    # Note on naming. The behavior of PDB2PQR is controlled by two
    # parameters, ff and ffout. My understanding is that the ff
    # parameter sets which residues are SUPPORTED by the underlying
    # FF, PLUS the charge and radii.  The ffout parameter sets the
    # naming scheme. Therefore, I want ff to be as general as
    # possible, which turns out to be "parse". Then I pick a
    # convenient ffout.

    # Hold list (None -> None)
    hlist = _selToHoldList(mol_in, holdSelection)
    if hlist:
        logger.warning("The holdSelection option is untested and deprecated. Please use reprepare()")

    # Relying on defaults
    pqr_res = runPDB2PQR(pdblist,
                         ph=pH, verbose=verbose,
                         ff="parse", ffout="amber",
                         ph_calc_method="propka31",
                         ph_calc_options=propka_opts,
                         holdList=hlist)
    try:
        header, pqr, missedLigands, pdb2pqr_protein, pdb2pqr_routines = \
            pqr_res['header'], pqr_res['lines'], pqr_res['missedligands'], pqr_res['protein'], pqr_res['routines']
    except:
        logger.error("Problem calling pdb2pqr. Make sure you have htmd-pdb2pqr >= 2.1.2a9")
        raise

    tmpin.close()

    # Diagnostics
    for missedligand in missedLigands:
        logger.warning("The following residue has not been optimized: " + missedligand)

    mol_out, resData = _buildResAndMol(pdb2pqr_protein)

    mol_out.box = mol_in.box
    _fixupWaterNames(mol_out)

    # Misc. info
    resData.header = header
    resData.pqr = pqr

    # Return residue information
    resData.pdb2pqr_protein = pdb2pqr_protein
    resData.pdb2pqr_routines = pdb2pqr_routines
    resData.missedLigands = missedLigands

    # Store un-reprepared info
    resData.data['default_protonation'] = resData.data['protonation']

    resData._listNonStandardResidues()
    resData._warnIfpKCloseTopH(pH)
    resData.warnIfTerminiSuspect()

    if hydrophobicThickness:
        resData._setMembraneExposureAndWarn(hydrophobicThickness)

    logger.debug("Returning.")
    logger.setLevel(oldLoggingLevel)

    if returnDetails:
        return mol_out, resData
    else:
        return mol_out



# Reproducibility test
# rm mol-test-*; for i in `seq 9`; do py ./proteinpreparation.py ./1r1j.pdb > mol-test-$i.log ; cp ./mol-test.pdb mol-test-$i.pdb; cp mol-test.csv mol-test-$i.csv ; done


# A test method
if __name__ == "__main__":
    import sys
    import htmd
    import os

    # No arguments - quick travis test
    if len(sys.argv) == 1:
        tryp_op, prepData = proteinPrepare(Molecule("3PTB"), returnDetails=True)
        d = prepData.data
        assert d.protonation[d.resid == 40].iloc[0] == 'HIE'
        assert d.protonation[d.resid == 57].iloc[0] == 'HIP'
        assert d.protonation[d.resid == 91].iloc[0] == 'HID'

    # Long test
    elif sys.argv[1] == "long-test":
        pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
        for pdb in pdbids:
            mol = Molecule(pdb)
            mol.filter("protein")
            mol_op, prepData = proteinPrepare(mol, returnDetails=True)
            mol_op.write("./{}-prepared.pdb".format(pdb))
            prepData.data.to_csv("./{}-prepared.csv".format(pdb), float_format="%.2f")

            compareDir = htmd.home(dataDir=os.path.join('test-proteinprepare', pdb))
            htmd.util.assertSameAsReferenceDir(compareDir)
        import doctest
        doctest.testmod()

    # Stand-alone executable: prepare the PDB given as argument
    else:
        mol = Molecule(sys.argv[1])
        mol.filter("protein")
        mol_op, prepData = proteinPrepare(mol, returnDetails=True)
        mol_op.write("./prepared.pdb")
        prepData.data.to_excel("./prepared.xlsx")
        prepData.data.to_csv("./prepared.csv")

        # x_HIS91_ND1 = tryp_op.get("coords","resid 91 and  name ND1")
        # x_SER93_H =   tryp_op.get("coords","resid 93 and  name H")
        # assert len(x_SER93_H) == 3
        # assert np.linalg.norm(x_HIS91_ND1-x_SER93_H) > 2
        # assert tryp_op.get("resname","resid 91 and  name CA") == "HIE"

        # prepData.data.loc[d.resid == 40, 'new_protonation'] = 'HIP'
        # mHIP40, pHIP40 = prepData.reprepare()
        # mHIP40.write("./mol-test-hip40.pdb")
        # pHIP40.data.to_excel("./mol-test-hip40.xlsx")


"""
    # Code to regenerate reference files. Run with PYTHONHASHSEED=1
    from htmd import *
    pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
    for p in pdbids:
        preparedInputDir = home(dataDir=os.path.join('test-proteinprepare', p))
        m=Molecule(p)
        m.filter("protein")
        mp,dp = proteinPrepare(m, returnDetails=True)
        inFile = os.path.join(preparedInputDir, "{}-prepared".format(p))
        mp.write(inFile+".pdb")
        dp.data.to_csv(inFile+".csv",float_format="%.2f")
"""
