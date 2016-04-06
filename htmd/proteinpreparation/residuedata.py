import numpy as np
import logging


logger = logging.getLogger(__name__)


# Define a type for holding information on residues decisions
class ResidueData:
    """ Class to hold results of the system preparation and optimization steps.

    ResiduesInfo contains the results of an optimization operation, notably, for each residue name, id, and chain, the
    corresponding pKa and protonation state.

    Examples
    --------
    >> tryp = Molecule("3PTB")
    >> tryp_op, info = systemPreparation(tryp, returnDetails=True)
    >> ri=info.residuesInfo
    >> ri.pKa[ri.resid==189]

    Properties
    ----------

    resid : np.ndarray
        Residue ID
    resname : np.ndarray
        Residue name, as per the original PDB
    chain : np.ndarray
        Chain
    pKa : np.ndarray
        pKa value computed by propKa
    protonation : np.ndarray
        Forcefield-independent protonation code
    patches : np.ndarray
        Additional information (may change)

    """
    _residuesinfo_fields = {
        'resid': np.int,
        'resname': object,
        'chain': object,
        'pKa': np.float32,
        'protonation': object,
        'patches': object
    }

    def __init__(self):
        for k in self._residuesinfo_fields:
            self.__dict__[k] = np.zeros(0, dtype=self._residuesinfo_fields[k])

    def __str__(self):
        n = len(self.resid)
        r = ""
        for i in range(n):
            r += "%4s %4d %1s : pKa=%f, state=%s, patches=%s\n" % (
                self.resname[i],
                self.resid[i],
                self.chain[i],
                self.pKa[i],
                self.protonation[i],
                self.patches[i])
        return r

    def __repr__(self):
        return self.__str__()

    def _findRes(self, a_resid, a_resname, a_chain):
        mask = (self.resid == a_resid) & (self.resname == a_resname) & (self.chain == a_chain)
        assert (sum(mask) <= 1)
        if sum(mask) == 0:
            self.resid = np.append(self.resid, a_resid)
            self.resname = np.append(self.resname, a_resname)
            self.chain = np.append(self.chain, a_chain)
            self.protonation = np.append(self.protonation, "UNK")
            self.pKa = np.append(self.pKa, np.NaN)
            self.patches = np.append(self.patches, "")
            pos = len(self.resid) - 1
        else:
            pos = np.argwhere(mask)
        return pos

    # residue is e.g. pdb2pqr.src.aa.ILE
    def setProtonation(self, residue, protonation):
        logger.debug("setProtonation %s %s" % (residue, protonation))
        pos = self._findRes(residue.resSeq, residue.name, residue.chainID)
        self.protonation[pos] = protonation

    # TODO this should actually append to a list
    def appendPatches(self, residue, patch):
        pos = self._findRes(residue.resSeq, residue.name, residue.chainID)
        self.patches[pos] += patch + "/"

    def setPKAs(self, pka_molecule):
        for grp in pka_molecule.conformations['AVR'].groups:
            resname = grp.residue_type
            resid = grp.atom.resNumb
            chain = grp.atom.chainID
            pka = grp.pka_value
            pos = self._findRes(resid, resname, chain)
            self.pKa[pos] = pka
