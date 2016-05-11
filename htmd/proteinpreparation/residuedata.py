import numpy as np
import pandas as pd
import logging
import textwrap

from htmd.molecule.molecule import Molecule

logger = logging.getLogger(__name__)


# Define a type for holding information on residues decisions
class ResidueData:
    """ Class to hold results of the system preparation and optimization steps.

    ResiduesInfo contains the results of an optimization operation, notably, for each residue name, id, and chain, the
    corresponding pKa and protonation state.

    Examples
    --------
    >>> tryp = Molecule("3PTB")
    >>> tryp_op, ri = prepareProtein(tryp, returnDetails=True)
    >>> ri
    ResidueData object with 290 residues. Please find the full info in the .data property.
      resname  resid insertion chain      pKa protonation    buried    patches
    0     ILE     16               A      NaN         ILE       NaN    [NTERM]
    1     VAL     17               A      NaN         VAL       NaN  [PEPTIDE]
    2     GLY     18               A      NaN         GLY       NaN  [PEPTIDE]
    3     GLY     19               A      NaN         GLY       NaN  [PEPTIDE]
    4     TYR     20               A  9.59085         TYR  0.146429  [PEPTIDE]
     . . .
    >>> ri.data.pKa[ri.data.resid==189].iat[0]
    4.9490792910341349
    >>> ri.data.patches[ri.data.resid==57]
    39    [PEPTIDE, HIP]
    Name: patches, dtype: object

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

     .missedLigands : str
            List of ligands residue names which were not optimized

     .propkaContainer : propka.molecular_container.Molecular_container
            Detailed information returned by propKa 3.1. See e.g.
                propkaContainer.conformations['AVR'].groups[4].__dict__
                propkaContainer.conformations['AVR'].groups[4].atom.__dict__
    """

    propkaContainer = None
    thickness = None
    _columns = ['resname', 'resid', 'insertion',
                'chain', 'pKa', 'protonation', 'buried',
                'patches', 'z', 'membraneExposed']
    _printColumns = ['resname', 'resid', 'insertion',
                'chain', 'pKa', 'protonation', 'buried',
                'patches']

    def __init__(self):
        self.data = pd.DataFrame(columns=self._columns)
        self.data.resid = self.data.resid.astype(int)
        self.data.pKa = self.data.pKa.astype(float)
        self.data.buried = self.data.buried.astype(float)
        self.data.z = self.data.z.astype(float)

    def __str__(self):
        r="ResidueData object about {:d} residues. Please find the full info in the .data property.\n".format(len(self.data))
        r+=str(self.data[self._printColumns].head())
        r+="\n . . ."
        return(r)

    def __repr__(self):
        return self.__str__()

    def _findRes(self, a_resname, a_resid, a_icode, a_chain):
        icode_pad = "{:1.1s}".format(a_icode)  # Pad and truncate to 1 char
        chain_pad = "{:1.1s}".format(a_chain)
        # Identity check should ignore residue name (self.resname == a_resname)
        mask = (self.data.resname == a_resname) & (self.data.resid == a_resid) & \
               (self.data.insertion == icode_pad) & (self.data.chain == chain_pad)
        assert (sum(mask) <= 1), "More than one resid matched"
        if sum(mask) == 0:
            self.data = self.data.append({'resname': a_resname,
                                          'resid': a_resid,
                                          'insertion': icode_pad,
                                          'chain': chain_pad,
                                          'patches': []}, ignore_index=True)
            pos = len(self.data) - 1
        else:
            pos = np.argwhere(mask)
            pos = int(pos)
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
        self.propkaContainer = pkaCont
        for grp in self.propkaContainer.conformations['AVR'].groups:
            pka = grp.pka_value
            buried = grp.buried
            resname = grp.residue_type
            # Other places for the resname: grp.type  -  grp.atom.resName
            resid = grp.atom.resNumb
            chain = grp.atom.chainID
            icode = grp.atom.icode
            z = grp.atom.z
            pos = self._findRes(resname, resid, icode, chain)
            self.data.set_value(pos, 'pKa', pka)
            self.data.set_value(pos, 'buried', buried)
            self.data.set_value(pos, 'z', z)

    def _setMembraneExposure(self, thickness, maxBuried=.75):
        self.thickness = thickness
        ht = thickness / 2.0
        inSlab = (self.data.z < -ht) | (self.data.z > ht)
        notBuried = self.data.buried < maxBuried
        inSlabNotBuried = inSlab & notBuried
        self.data.membraneExposed = inSlabNotBuried
        if np.any(inSlabNotBuried):
            logger.warning(
                "Predictions for {:d} residues may be incorrect because they are exposed to the membrane ({:.1f}<z<{:.2f} and buried<{:.1f}%).".format(
                    sum(inSlabNotBuried), -ht, ht, maxBuried * 100.0))


if __name__ == "__main__":
    from htmd.proteinpreparation.proteinpreparation import prepareProtein
    import doctest

    doctest.testmod()
