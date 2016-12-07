# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging

import numpy as np
import pandas as pd

from htmd.molecule.molecule import Molecule

logger = logging.getLogger(__name__)


def prettyPrintResidue(r):
    rs = "{:4s} {:4d}{:1s} {:1s}".format(r.resname, r.resid, r.insertion, r.chain)
    return rs


# Define a type for holding information on residues decisions
class PreparationData:
    """Results of the system preparation and optimization steps.

    Contains the results of an optimization operation, notably, for each residue name, id, and chain, the
    corresponding pKa and protonation state.

    The most important properties are accessible via the .data property, a pandas DataFrame. As such, they
    can be subset, converted, and saved in several ways (see examples below).

    Examples
    --------
    >>> tryp = Molecule("3PTB")
    >>> tryp_op, ri = proteinPrepare(tryp, returnDetails=True)
    >>> ri                                  # doctest: +NORMALIZE_WHITESPACE
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
    >>> "%.2f" % ri.data.pKa[ri.data.resid==189]
    '4.95'
    >>> ri.data.patches[ri.data.resid==57]
    39    [PEPTIDE, HIP]
    Name: patches, dtype: object
    >>> ri.data.to_csv("/tmp/report.csv")

    Attributes
    ----------
    data : :class:`DataFrame <pandas.core.frame.DataFrame>` object
        A pandas DataFrame with these columns: resname "Residue name, as per the original PDB", resid "Residue ID",
        insertion "Insertion code (resid suffix)", chain "Chain", pKa "pKa value computed by propKa",
        "protonation" Forcefield-independent protonation code, flipped "Whether the residue was flipped during the
        optimization", buried "Fraction of residue which is buried", membraneExposed "Whether residue is exposed to
        membrane", patches "Additional information (may change)" (etc)
    missedLigands : str
        List of ligands residue names which were not optimized
    header : str
        Messages and warnings from PDB2PQR
    propkaContainer : propka.molecular_container.Molecular_container
        Detailed information returned by propKa 3.1.
    """

    # Important- all must be listed or "set_value" will silently ignore them
    _columns = ['resname', 'resid', 'insertion', 'chain',
                'pKa', 'protonation', 'flipped', 'patches',
                'buried', 'z', 'membraneExposed',
                'pka_group_id',
                'pka_residue_type', 'pka_type', 'pka_charge',
                'pka_atom_name', 'pka_atom_sybyl_type']

    # Columns printed by the __str__ method
    _printColumns = ['resname', 'resid', 'insertion', 'chain',
                     'pKa', 'protonation', 'flipped', 'buried']

    def __init__(self):
        self.propkaContainer = None
        self.thickness = None
        self.missedLigands = []

        self.data = pd.DataFrame(columns=self._columns)
        self.data.resid = self.data.resid.astype(int)
        self.data.pKa = self.data.pKa.astype(float)
        self.data.buried = self.data.buried.astype(float)
        self.data.z = self.data.z.astype(float)
        self.data.pka_group_id = self.data.pka_group_id.astype(float)  # should be int, but NaN not allowed
        # self.data.flipped = self.data.flipped.astype(float)             #  should be bool, but NaN not allowed

    def __str__(self):
        r = "PreparationData object about {:d} residues.\n".format(len(self.data))
        if len(self.missedLigands) > 0:
            r += "Unparametrized residue names: " + ", ".join(self.missedLigands) + "\n"
        r += "Please find the full info in the .data property, e.g.: \n".format(len(self.data))
        r += str(self.data[self._printColumns].head())
        r += "\n . . ."
        return r

    def __repr__(self):
        return self.__str__()

    def _findRes(self, a_resname, a_resid, a_icode, a_chain, forceAppend=False):
        icode_pad = "{:1.1s}".format(a_icode)  # Pad and truncate to 1 char
        chain_pad = "{:1.1s}".format(a_chain)
        # Identity check should ignore residue name (self.data.resname == a_resname). Not doing this because
        # the N+ and C- resnames are indeed duplicated
        mask = (self.data.chain == chain_pad) & (self.data.resid == a_resid) & \
               (self.data.insertion == icode_pad) & (self.data.resname == a_resname)
        if sum(mask) == 0 or forceAppend:
            self.data = self.data.append({'resname': a_resname,
                                          'resid': a_resid,
                                          'insertion': icode_pad,
                                          'chain': chain_pad,
                                          'patches': []}, ignore_index=True)
            pos = len(self.data) - 1
        elif sum(mask) == 1:
            pos = np.argwhere(mask)
            pos = int(pos)
        else:
            assert False, "More than one resid matched: either duplicated chain-residue-icode, or internal error (please report if the latter)."
        return pos

    # Generic setter in the pandas table. Maybe one should use actual indices instead.
    def _set(self, residue, key, val):
        pos = self._findRes(residue.name, residue.resSeq, residue.iCode, residue.chainID)
        self.data.set_value(pos, key, val)

    # residue is e.g. pdb2pqr.src.aa.ILE
    def _setProtonationState(self, residue, state):
        # logger.debug("_setProtonationState %s %s" % (residue, state))
        self._set(residue, 'protonation', state)

    def _setFlipped(self, residue, state):
        logger.debug("_setFlipped %s %s" % (residue, state))
        self._set(residue, 'flipped', state)

    def _appendPatches(self, residue, patch):
        # logger.debug("_appendPatches %s %s" % (residue, patch))
        pos = self._findRes(residue.name, residue.resSeq, residue.iCode, residue.chainID)
        self.data.patches[pos].append(patch)

    def _importPKAs(self, pkaCont):
        logger.debug("Called _importPKAs")
        self.propkaContainer = pkaCont
        for i, grp in enumerate(self.propkaContainer.conformations['AVR'].groups):
            forceAppend = False
            # This is the key
            # Other places for the resname: grp.type  -  grp.atom.resName  grp.residue_type
            resname = grp.atom.resName
            if grp.residue_type in ['N+', 'C-']:  # Separate info about termini
                resname = grp.residue_type
                forceAppend = True
            elif grp.atom.sybyl_assigned:  # A ligand - a hack to allow multiple groups overriding key
                forceAppend = True
            resid = grp.atom.resNumb
            chain = grp.atom.chainID
            icode = grp.atom.icode
            pos = self._findRes(resname, resid, icode, chain, forceAppend)
            self.data.set_value(pos, 'pKa', grp.pka_value)
            self.data.set_value(pos, 'buried', grp.buried * 100.0)
            self.data.set_value(pos, 'z', grp.atom.z)
            self.data.set_value(pos, 'pka_group_id', i)
            self.data.set_value(pos, 'pka_residue_type', grp.residue_type)
            self.data.set_value(pos, 'pka_type', grp.type)
            self.data.set_value(pos, 'pka_charge', grp.charge)
            self.data.set_value(pos, 'pka_atom_name', grp.atom.name)
            self.data.set_value(pos, 'pka_atom_sybyl_type', grp.atom.sybyl_type)

    def _setMembraneExposureAndWarn(self, thickness, maxBuried=75.0):
        self.thickness = thickness
        ht = thickness / 2.0
        inSlab = (self.data.z > -ht) & (self.data.z < ht)
        notBuried = self.data.buried < maxBuried
        inSlabNotBuried = inSlab & notBuried
        self.data.membraneExposed = inSlabNotBuried
        if np.any(inSlabNotBuried):
            # dl = self._prettyPrintResidues(inSlabNotBuried)
            logger.warning(
                ("Predictions for {:d} residues may be incorrect because they are " +
                 "exposed to the membrane ({:.1f}<z<{:.2f} and buried<{:.1f}%).").format(
                    sum(inSlabNotBuried), -ht, ht, maxBuried))

    def _listNonStandardResidues(self):
        changed = self.data.resname != self.data.protonation
        cl = []
        for i, cr in self.data[changed].iterrows():
            if cr.resname in ['N+', 'C-'] or cr.protonation in ['WAT'] or type(cr.protonation) == float:
                continue
            cl.append("{:s} ({:s})".format(prettyPrintResidue(cr), cr.protonation))
        if cl:
            logger.info("The following residues are in a non-standard state: " + ", ".join(cl))

    def _warnIfpKCloseTopH(self, ph, tol=1.0):
        # Looks like NaN < 5 is False today
        dubious = abs(self.data.pKa - ph) < tol
        nd = sum(dubious)
        if nd > 1:
            logger.warning(
                "Dubious protonation state: the pKa of {:d} residues is within {:.1f} units of pH {:.1f}."
                    .format(nd, tol, ph))
            for i, dr in self.data[dubious].iterrows():
                drs = prettyPrintResidue(dr)
                logger.warning("Dubious protonation state:    {:s} (pKa={:5.2f})".format(drs, dr.pKa))

    def reprepare(self):
        """Repeat the system preparation, after changin the .data table.

        You should only modify the value of the .data.new_protonation column on the basis of the values
        in the .data.resid, .data.insertion, .data.chain attributes. Any other change will be ignored.

        Returns
        -------
        mol_out : Molecule
            the molecule titrated and optimized. The molecule object contains an additional attribute,
        resData : ResidueData
            a table of residues with the corresponding protonation states, pKas, and other information

        Examples
        --------
        mol, prepData = proteinPrepare(Molecule("3PTB"), returnDetails=True)
        d = prepData.data
        d.loc[d.resid == 40, 'new_protonation'] = 'HIP'
        mHIP40, pHIP40 = prepData.reprepare()

        """

        from pdb2pqr.src.hydrogens import hydrogenRoutines
        from pdb2pqr.src.forcefield import Forcefield
        from pdb2pqr.src.definitions import Definition
        from htmd.builder.preparation import _buildResAndMol

        d = self.data
        routines = self.pdb2pqr_routines
        p = routines.protein

        neutraln = neutralc = False
        assign_only = clean = False
        debump = opt = True


        # Code lifted from resinter.py
        routines.removeHydrogens()
        for index, oldResidue in enumerate(p.getResidues()):
            chain = p.chainmap[oldResidue.chainID]
            chainIndex = chain.residues.index(oldResidue)

            d_idx = d.pdb2pqr_idx == index
            if sum(d_idx) != 1:
                logger.warning("Residue {:s} appears {:d} times in data table".format(str(oldResidue), sum(d_idx)))
                continue

            newResidueName = d.new_protonation[d_idx].iloc[0]
            if pd.isnull(newResidueName):
                # newResidueName = d.protonation[d_idx].iloc[0]
                continue
            logger.debug("Replacing {} with {}".format(oldResidue, newResidueName))

            # Create the replacement residue
            residueAtoms = oldResidue.atoms
            newResidue = routines.protein.createResidue(residueAtoms, newResidueName)
            # Make sure our names are cleaned up for output.
            newResidue.renameResidue(newResidueName)
            # Drop it in
            p.residues[index] = newResidue
            chain.residues[chainIndex] = newResidue
            # Run the meaty bits of PDB2PQR
        routines.setTermini(neutraln, neutralc)
        routines.updateBonds()

        if not clean and not assign_only:
            routines.updateSSbridges()
            if debump:
                routines.debumpProtein()
            routines.addHydrogens()
            hydRoutines = hydrogenRoutines(routines)
            if debump:
                routines.debumpProtein()
            if opt:
                hydRoutines.setOptimizeableHydrogens()
                hydRoutines.initializeFullOptimization()
                hydRoutines.optimizeHydrogens()
            else:
                hydRoutines.initializeWaterOptimization()
                hydRoutines.optimizeHydrogens()
            # Special for GLH/ASH, since both conformations were added
            hydRoutines.cleanup()


        ff = "parse"
        ffout = "amber"
        usernames = userff = None

        routines.setStates()
        mydef = Definition()
        myForcefield = Forcefield(ff, mydef, userff, usernames)
        hitlist, misslist = routines.applyForcefield(myForcefield)
        # reslist, charge = routines.getCharge() # <--- ?

        # Copied from runPDB2PQR = ?
        if not ffout is None:
            if ffout != ff:
                myNameScheme = Forcefield(ffout, mydef, userff)
            else:
                myNameScheme = myForcefield
                routines.applyNameScheme(myNameScheme)


        newMol, newResData = _buildResAndMol(p)
        return newMol, newResData


if __name__ == "__main__":
    from htmd.builder.preparation import proteinPrepare

    import doctest

    doctest.testmod()
