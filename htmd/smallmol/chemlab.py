from rdkit import Chem
import numpy as np
import logging
from htmd.smallmol.util import  convertToString
from tqdm import tqdm
from htmd.parallelprogress import ParallelExecutor, delayed

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

def cluster(smallmol_list, method, distThresholds=0.2, returnDetails=True, removeHs=True ):
    from sklearn.cluster import DBSCAN  # clusterizzatore

    import sys
    this_module = sys.modules[__name__]

    _methods = ['maccs', 'pathFingerprints', 'atomsFingerprints', 'torsionsFingerprints',
                'circularFingerprints', 'shape', 'mcs']

    if method not in _methods:
        raise ValueError('The method provided {} does not exists. The ones available are the following: {}'.format(method, _methods))

    smallmol_list = np.array([sm.copy() for sm in smallmol_list])

    if removeHs:
        for sm in smallmol_list:
            sm._removeAtoms(sm.get('element H', 'idx'))

    rdkitMols_list = [sm.toRdkitMol(includeConformer=True) for sm in smallmol_list]

    clustmethod = getattr(this_module, '_{}Clustering'.format(method))

    if method not in ['shape', 'mcs']:
        matrix = clustmethod(rdkitMols_list)

    else:
        aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
        matrix = aprun(total=len(rdkitMols_list), desc='{} Distance'.format(method)) \
                            (delayed(clustmethod)(mol1, rdkitMols_list) for mol1 in rdkitMols_list)

        matrix = np.array(matrix)


    db = DBSCAN(eps=distThresholds, min_samples=0, metric='precomputed' ).fit(matrix)

    labels = db.labels_

    populations = np.bincount(labels)
    n_clusters = np.max(labels)

    clusters_idx = np.empty((n_clusters,), dtype=object)
    clusters_smallmols = np.empty((n_clusters,), dtype=object)

    for n_cl in np.arange(n_clusters):
        idxs = np.where(labels == n_cl)[0]
        clusters_idx[n_cl] = idxs
        clusters_smallmols[n_cl] = smallmol_list[idxs]


    if returnDetails:
        details = {'numClusters':n_clusters,
                   'populations':populations,
                   'clusters':clusters_idx}
        return clusters_smallmols, details

    return clusters_smallmols




def _maccsClustering( rdkit_mols):
    from rdkit.Chem import MACCSkeys  # calcola MACCS keys

    fps = []
    for m in tqdm(rdkit_mols):
       fps.append(MACCSkeys.GenMACCSKeys(m))

    aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
    tanimoto_matrix = aprun(total=len(fps), desc='MACCS Distance') \
            (delayed(TanimotoDistances)(fp1, fps) for fp1 in fps)


    return np.array(tanimoto_matrix)

def _pathFingerprintsClustering(rdkit_mols):
    from rdkit.Chem.Fingerprints import FingerprintMols  # calcola path fingerprints

    fps = [ ]
    for m in tqdm(rdkit_mols):
        fps.append(FingerprintMols.FingerprintMol(m))

    aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
    tanimoto_matrix = aprun(total=len(fps), desc='PathFingerprints Distance') \
        (delayed(TanimotoDistances)(fp1, fps) for fp1 in fps)

    return np.array(tanimoto_matrix)

def _atomsFingerprintsClustering(rdkit_mols):
    from rdkit.Chem.AtomPairs import Pairs  # Atom pairs

    fps = []
    for m in tqdm(rdkit_mols):
        fps.append(Pairs.GetAtomPairFingerprint(m))

    aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
    dice_matrix = aprun(total=len(fps), desc='AtomsFingerprints Distance') \
        (delayed(DiceDistances)(fp1, fps) for fp1 in fps)

    return np.array(dice_matrix)

def _torsionsFingerprintsClustering(rdkit_mols):
    from rdkit.Chem.AtomPairs import Torsions  # Topological Torsions

    fps = [ ]
    for m in tqdm(rdkit_mols):
        fps.append(Torsions.GetHashedTopologicalTorsionFingerprint(m))

    aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
    dice_matrix = aprun(total=len(fps), desc='TorsionsFingerprints Distance') \
        (delayed(DiceDistances)(fp1, fps) for fp1 in fps)

    return np.array(dice_matrix)

def _circularFingerprintsClustering(rdkit_mols, radius=2):
    from rdkit.Chem import AllChem  # calcola circular fingerprints

    fps = []
    for m in rdkit_mols:
        fps.append(AllChem.GetMorganFingerprint(m, radius))

    aprun = ParallelExecutor(n_jobs=-1)  # _config['ncpus'])
    dice_matrix = aprun(total=len(fps), desc='CircularFingerprints Distance') \
        (delayed(DiceDistances)(fp1, fps) for fp1 in fps)

    return np.array(dice_matrix)

def _shapeClustering(mol1, rdkit_mols):
    from rdkit.Chem import rdMolAlign, rdShapeHelpers
    tanimoto_shape_row = []
    for mol2 in rdkit_mols:
        oa3 = rdMolAlign.GetO3A(mol1, mol2)
        oa3.Align()
        tani_shape = rdShapeHelpers.ShapeTanimotoDist(mol1, mol2)
        tanimoto_shape_row.append(tani_shape)
    return tanimoto_shape_row

def _mcsClustering(mol1, rdkit_mols):
    from rdkit.Chem import rdFMCS

    MCS_row = []
    for mol2 in rdkit_mols:
        sum_numHeavyAtoms = mol1.GetNumHeavyAtoms() + mol2.GetNumHeavyAtoms()
        mcsHeavyAtoms = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=True, completeRingsOnly=True, timeout=5)
        mcsNumHeavyAtoms = float(Chem.MolFromSmarts(mcsHeavyAtoms.smartsString).GetNumHeavyAtoms())
        distance = 1 - mcsNumHeavyAtoms * 2/sum_numHeavyAtoms
        MCS_row.append(distance)

    return MCS_row

def TanimotoDistances(fp1, fps):
    from rdkit import DataStructs  # fingerprint similarity

    tanimoto_row = []

    for fp2 in fps:
        tani = 1 - DataStructs.FingerprintSimilarity(fp1, fp2)
        tanimoto_row.append(tani)
    return tanimoto_row


def DiceDistances(fp1, fps):
    from rdkit import DataStructs  # fingerprint similarity

    dice_row = []

    for fp2 in fps:
        dice = 1 - DataStructs.DiceSimilarity(fp1, fp2)
        dice_row.append(dice)

    return dice_row
