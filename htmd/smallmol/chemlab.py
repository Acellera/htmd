from rdkit import Chem
import numpy as np
import logging
from htmd.smallmol.util import  convertToString

logger = logging.getLogger(__name__)

def getMaximumCommonSubstructure(smallmol_list, removeHs=True, returnAtomIdxs=False):
    from rdkit.Chem import rdFMCS

    smallmol_list = [sm.copy() for sm in smallmol_list]

    if  removeHs:
        for sm in smallmol_list:
            sm._removeAtoms(sm.get('element H', 'idx'))

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


def cluster(smallmol_list, method, removeHs=True ):
    from sklearn.cluster import DBSCAN  # clusterizzatore
    import sys
    this_module = sys.modules[__name__]
    from rdkit.Chem.Fingerprints import FingerprintMols  # calcola path fingerprints
    from rdkit.Chem.AtomPairs import Pairs  # Atom pairs
    from rdkit.Chem.AtomPairs import Torsions  # Topological Torsions
    from rdkit.Chem import AllChem  # calcola circular fingerprints

    _methods = [ 'maccs','pathFingerprints','atomFingerprints','torsionsFingerprints',
                'circularFingerprints', 'shape', 'mcs']
    # _methods_compute = {'maccs' : MACCSkeys.GenMACCSKeys,
    #                     'pathFingerprints' : FingerprintMols.FingerprintMol,
    #                     'atomFingerprints' : Pairs.GetAtomPairFingerprint,
    #                     'torsionsFingerprints':Torsions.GetTopologicalTorsionFingerprint,
    #                     'circularFingerprints':AllChem.GetMorganFingerprint,
    #                 }

    smallmol_list = [sm.copy() for sm in smallmol_list]

    if removeHs:
        for sm in smallmol_list:
            sm._removeAtoms(sm.get('element H', 'idx'))

    rdkitMols_list = [sm.toRdkitMol() for sm in smallmol_list]

    clustmethod = getattr(this_module, '_{}Clustering'.format(method))

    matrix_distance = clustmethod(rdkitMols_list)

    db = DBSCAN(eps=0.2, min_samples=0, metric='precomputed' ).fit(matrix_distance)

    return db.labels_


def _maccsClustering(rdkit_mols):
    from rdkit.Chem import MACCSkeys  # calcola MACCS keys

    fps = [MACCSkeys.GenMACCSKeys(m) for m in rdkit_mols ]

    tanimoto_matrix = TanimotoDistances(fps)

    return tanimoto_matrix

def TanimotoDistances(fps):
    from rdkit import DataStructs  # fingerprint similarity

    tanimoto_matrix = []

    for fp1 in fps:
        tanimoto_row = []
        for fp2 in fps:
            tani = 1 - DataStructs.FingerprintSimilarity(fp1, fp2)
            tanimoto_row.append(tani)
        tanimoto_matrix.append(tanimoto_row)
    tanimoto_matrix = np.array(tanimoto_matrix)

    return tanimoto_matrix





