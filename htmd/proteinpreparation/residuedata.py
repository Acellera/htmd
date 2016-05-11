import numpy as np
import logging

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
    >>> ri                                                      # doctest: +ELLIPSIS
    ILE    16  A : pKa   nan, buried  nan, state=ILE, patches=NTERM
    VAL    17  A : pKa   nan, buried  nan, state=VAL, patches=PEPTIDE
    GLY    18  A : pKa   nan, buried  nan, state=GLY, patches=PEPTIDE
    GLY    19  A : pKa   nan, buried  nan, state=GLY, patches=PEPTIDE
    TYR    20  A : pKa  9.59, buried 0.15, state=TYR, patches=PEPTIDE
    ...
    >>> ri.pKa[ri.resid == 189]
    array([ 4.94907929])
    >>> ri.patches[ri.resid == 57]
    array([['PEPTIDE', 'HIP']], dtype=object)


    Properties
    ----------
    List-like objects, one per residue:
        resname : np.ndarray(str)
            Residue name, as per the original PDB
        resid : np.ndarray(int)
            Residue ID
        insertion : np.ndarray(str)
            Insertion code (resid suffix)
        chain : np.ndarray(str)
            Chain
        pKa : np.ndarray(np.float32)
            pKa value computed by propKa
        protonation : np.ndarray(str)
            Forcefield-independent protonation code
        buried : np.ndarray(np.float32)
            Fraction of residue's surface which is buried
        membraneExposed: np.ndarray(np.bool)
            Whether residue is exposed to membrane
        patches : np.ndarray of lists
            Additional information (may change)

    Other objects:
        missedLigands : str
            List of ligands residue names which were not optimized
        propkaContainer : propka.molecular_container.Molecular_container
            Detailed information returned by propKa 3.1. See e.g.
                propkaContainer.conformations['AVR'].groups[4].__dict__
                propkaContainer.conformations['AVR'].groups[4].atom.__dict__
    """

    _residuesinfo_fields = {
        'resname': object,
        'resid': np.int,
        'insertion': object,
        'chain': object,
        'pKa': np.float32,
        'protonation': object,
        'buried': np.float32,
        'membraneExposed': np.bool,
        'z': np.float32,
        'patches': object
    }

    def __init__(self):
        for k in self._residuesinfo_fields:
            self.__dict__[k] = np.empty(0, dtype=self._residuesinfo_fields[k])

    def __str__(self):
        n = len(self.resid)
        r = ""
        for i in range(n):
            r += "{:4s} {:4d}{:1s} {:1s} : pKa {:5.2f}, buried {:4.2f}, state={:3s}, patches={:s}\n".format(
                self.resname[i],
                self.resid[i],
                self.insertion[i],
                self.chain[i],
                self.pKa[i],
                self.buried[i],
                self.protonation[i],
                ",".join(self.patches[i]))
        return r

    def __repr__(self):
        return self.__str__()

    def listResidues(self, sel):
        """List a subset of the residues in the object.

        Parameters
        ----------
        sel : np.array(bool)
            A vector of boolean values

        Returns
        -------
        rls : str
            A list of residues, as strings
        """
        rl = []
        for v in range(len(sel)):
            if sel[v]:
                rl.append("{:4s} {:4d}{:1s} {:1s}".format(self.resname[v], self.resid[v],
                                                          self.insertion[v], self.chain[v]))
        return rl

    def _findRes(self, a_resname, a_resid, a_icode, a_chain):
        icode_pad = "{:1.1s}".format(a_icode)  # Pad and truncate to 1 char
        chain_pad = "{:1.1s}".format(a_chain)
        # Identity check should ignore residue name (self.resname == a_resname)
        mask = (self.resname == a_resname) & (self.resid == a_resid) & \
               (self.insertion == icode_pad) & (self.chain == chain_pad)
        assert (sum(mask) <= 1), "More than one resid matched"
        if sum(mask) == 0:
            # Growing non-mutables is not the best way, but let's be consistent with Molecule
            self.resname = np.append(self.resname, a_resname)
            self.resid = np.append(self.resid, a_resid)
            self.insertion = np.append(self.insertion, icode_pad)
            self.chain = np.append(self.chain, chain_pad)
            self.protonation = np.append(self.protonation, "UNK")
            self.pKa = np.append(self.pKa, np.NaN)
            self.buried = np.append(self.buried, np.NaN)
            self.membraneExposed = np.append(self.membraneExposed, False)
            self.z = np.append(self.z, np.NaN)
            self.patches = np.append(self.patches, "")
            self.patches[-1] = []
            pos = len(self.resid) - 1
        else:
            pos = np.argwhere(mask)
            pos = int(pos)
        return pos

    # residue is e.g. pdb2pqr.src.aa.ILE
    def _setProtonationState(self, residue, state):
        logger.debug("_setProtonationState %s %s" % (residue, state))
        pos = self._findRes(residue.name, residue.resSeq, residue.iCode, residue.chainID)
        self.protonation[pos] = state

    def _appendPatches(self, residue, patch):
        pos = self._findRes(residue.name, residue.resSeq, residue.iCode, residue.chainID)
        self.patches[pos].append(patch)

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
            self.pKa[pos] = pka
            self.buried[pos] = buried
            self.z[pos] = z

    def _checkMembraneExposure(self, thickness, maxBuried=.75):
        import textwrap
        ht = thickness / 2.0
        old_settings = np.seterr(invalid='ignore')
        inSlab = np.logical_or(self.z < -ht, self.z > ht)
        notBuried = self.buried < maxBuried
        inSlabNotBuried = np.logical_and(inSlab, notBuried)
        np.seterr(**old_settings)
        self.membraneExposed[inSlabNotBuried] = True
        if np.any(inSlabNotBuried):
            rl = self.listResidues(inSlabNotBuried)
            rls = ", ".join(rl)
            rlsw = textwrap.wrap(rls)
            logger.warning(
                "Predictions for these residues may be incorrect because they are exposed to the membrane ({:.2f}<z<{:.2f} and buried<{:.2f}%).".format(
                    -ht, ht, maxBuried*100))
            logger.warning("\n".join(rlsw))




if __name__ == "__main__":
    from htmd.proteinpreparation.proteinpreparation import prepareProtein
    import doctest

    doctest.testmod()
