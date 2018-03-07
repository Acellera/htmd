from rdkit import Chem
import logging
from htmd.smallmol.util import  convertToString

logger = logging.getLogger(__name__)

def getMaximumCommonSubstructure(smallmol_list, returnAtomIdxs=False):
    from rdkit.Chem import rdFMCS

    rdkitMols_list = [ sm.toRdkitMol() for sm in smallmol_list]

    mcs = rdFMCS.FindMCS(rdkitMols_list)

    logger.info('MCS found a substructure of {} atoms and {} bonds'.format(mcs.numAtoms, mcs.numBonds))

    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    if not returnAtomIdxs:
        return mcs_mol

    atoms_mcs_list = []
    atoms_no_mcs_list = []
    for sm, m in zip(smallmol_list, rdkitMols_list):
        match = m.GetSubstructMatch(mcs_mol)
        sel_str = convertToString(match)

        atoms_mcs = sm.get('idx {}'.format(sel_str), 'idx')
        atoms_no_mcs = sm.get('idx {}'.format(sel_str), 'idx', invert=True)

        atoms_mcs_list.append(atoms_mcs.tolist())
        atoms_no_mcs_list.append(atoms_no_mcs.tolist())

    return mcs_mol, atoms_mcs_list, atoms_no_mcs_list


