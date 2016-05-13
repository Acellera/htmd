import numpy as np
import pandas as pd
import logging

from htmd.molecule.molecule import Molecule

logger = logging.getLogger(__name__)


# Define a type for holding information on residues decisions
class ResidueData:
    """Results of the system preparation and optimization steps.

    Contains the results of an optimization operation, notably, for each residue name, id, and chain, the
    corresponding pKa and protonation state.

    The most important properties are accessible via the .data property, a pandas DataFrame. As such, they
    can be subset, converted, and saved in several ways (see examples below).

    Examples
    --------
    >>> tryp = Molecule("3PTB")
    >>> tryp_op, ri = prepareProtein(tryp, returnDetails=True)
    >>> ri
    ResidueData object about 287 residues. Please find the full info in the .data property.
      resname  resid insertion chain       pKa protonation    buried    patches
    0     ILE     16               A  7.413075         ILE  0.839286    [NTERM]
    1     VAL     17               A       NaN         VAL       NaN  [PEPTIDE]
    2     GLY     18               A       NaN         GLY       NaN  [PEPTIDE]
    3     GLY     19               A       NaN         GLY       NaN  [PEPTIDE]
    4     TYR     20               A  9.590845         TYR  0.146429  [PEPTIDE]
     . . .
    >>> "%.2f" % ri.data.pKa[ri.data.resid==189]
    '4.95'
    >>> ri.data.patches[ri.data.resid==57]
    39    [PEPTIDE, HIP]
    Name: patches, dtype: object
    >>> ri.data.to_csv("/tmp/report.csv")

    Properties
    ----------
    .data
        A pandas DataFrame with these columns
        resname : str
            Residue name, as per the original PDB
        resid : int
            Residue ID
        insertion : str
            Insertion code (resid suffix)
        chain : str
            Chain
        pKa : float
            pKa value computed by propKa
        protonation : str
            Forcefield-independent protonation code
        buried : float
            Fraction of residue which is buried
        membraneExposed: bool
            Whether residue is exposed to membrane
        patches : str[]
            Additional information (may change)
        (etc)

     .missedLigands : str
            List of ligands residue names which were not optimized

     .propkaContainer : propka.molecular_container.Molecular_container
            Detailed information returned by propKa 3.1. See e.g.
                propkaContainer.conformations['AVR'].groups[4].__dict__
                propkaContainer.conformations['AVR'].groups[4].atom.__dict__
    """

    propkaContainer = None
    thickness = None

    # Important- all must be listed or "set_value" will silently ignore them
    _columns = ['resname', 'resid', 'insertion',
                'chain', 'pKa', 'protonation', 'buried',
                'patches', 'z', 'membraneExposed',
                'pka_group_id',
                'pka_residue_type', 'pka_type', 'pka_charge',
                'pka_atom_name', 'pka_atom_sybyl_type']

    # Columns printed by the __str__ method
    _printColumns = ['resname', 'resid', 'insertion',
                'chain', 'pKa', 'protonation', 'buried',
                'patches']

    def __init__(self):
        self.data = pd.DataFrame(columns=self._columns)
        self.data.resid = self.data.resid.astype(int)
        self.data.pKa = self.data.pKa.astype(float)
        self.data.buried = self.data.buried.astype(float)
        self.data.z = self.data.z.astype(float)
        self.data.pka_group_id = self.data.pka_group_id.astype(int)

    def __str__(self):
        r="ResidueData object about {:d} residues. Please find the full info in the .data property.\n".format(len(self.data))
        r+=str(self.data[self._printColumns].head())
        r+="\n . . ."
        return(r)

    def __repr__(self):
        return self.__str__()

    def _findRes(self, a_resname, a_resid, a_icode, a_chain, forceAppend=False):
        icode_pad = "{:1.1s}".format(a_icode)  # Pad and truncate to 1 char
        chain_pad = "{:1.1s}".format(a_chain)
        # Identity check should ignore residue name (self.resname == a_resname)
        mask = (self.data.resname == a_resname) & (self.data.resid == a_resid) & \
               (self.data.insertion == icode_pad) & (self.data.chain == chain_pad)
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
            assert False, "More than one resid matched (internal error, please report)"
        return pos

    # residue is e.g. pdb2pqr.src.aa.ILE
    def _setProtonationState(self, residue, state):
        logger.debug("_setProtonationState %s %s" % (residue, state))
        pos = self._findRes(residue.name, residue.resSeq, residue.iCode, residue.chainID)
        self.data.set_value(pos, 'protonation', state)

    def _appendPatches(self, residue, patch):
        logger.debug("_appendPatches %s %s" % (residue, patch))
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
            if grp.residue_type in ['N+','C-']: # Separate info about termini
                resname = grp.residue_type
                forceAppend = True
            elif grp.atom.sybyl_assigned:       # A ligand - a hack to allow multiple groups overriding key
                forceAppend = True
            resid = grp.atom.resNumb
            chain = grp.atom.chainID
            icode = grp.atom.icode
            pos = self._findRes(resname, resid, icode, chain, forceAppend)
            self.data.set_value(pos, 'pKa', grp.pka_value)
            self.data.set_value(pos, 'buried', grp.buried)
            self.data.set_value(pos, 'z', grp.atom.z)
            self.data.set_value(pos, 'pka_group_id', i)
            self.data.set_value(pos, 'pka_residue_type', grp.residue_type)
            self.data.set_value(pos, 'pka_type', grp.type)
            self.data.set_value(pos, 'pka_charge', grp.charge)
            self.data.set_value(pos, 'pka_atom_name', grp.atom.name)
            self.data.set_value(pos, 'pka_atom_sybyl_type', grp.atom.sybyl_type)

    def _setMembraneExposure(self, thickness, maxBuried=.75):
        self.thickness = thickness
        ht = thickness / 2.0
        inSlab = (self.data.z > -ht) & (self.data.z < ht)
        notBuried = self.data.buried < maxBuried
        inSlabNotBuried = inSlab & notBuried
        self.data.membraneExposed = inSlabNotBuried
        if np.any(inSlabNotBuried):
            logger.warning(
                ("Predictions for {:d} residues may be incorrect because they are "+
                 "exposed to the membrane ({:.1f}<z<{:.2f} and buried<{:.1f}%).").format(
                    sum(inSlabNotBuried), -ht, ht, maxBuried * 100.0))


if __name__ == "__main__":
    from htmd.proteinpreparation.proteinpreparation import prepareProtein

    import doctest
    doctest.testmod()
