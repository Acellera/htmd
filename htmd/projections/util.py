# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from scipy.sparse import lil_matrix
import numpy as np
from IPython.core.debugger import Tracer


def pp_calcDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
    distances = _distanceArray(mol, sel1, sel2, mol.numFrames, pbc)
    distances = _postProcessDistances(distances, sel1, sel2, truncate)

    if metric == 'contacts':
        #metric = lil_matrix(distances <= threshold)
        metric = distances <= threshold
    elif metric == 'distances':
        metric = distances.astype(dtype=np.float32)
    else:
        raise NameError('The metric you asked for is not supported. Check spelling and documentation')
    return metric


def pp_calcMinDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
    raise NameError('TODO: Not implemented yet')
    return


def _postProcessDistances(distances, sel1, sel2, truncate):
    # Setting upper triangle to -1 if same selections
    if np.array_equal(sel1, sel2):
        for i in range(len(distances)):
            distances[i][:, range(i+1)] = -1

    if np.ndim(distances[0]) > 1:  # 2D data
        distances = np.concatenate(distances, axis=1)
    else: # 1D data
        distances = np.vstack(distances).transpose()

    if np.array_equal(sel1, sel2):
        distances = distances[:, np.all(distances != -1, 0)]

    if truncate is not None:
        distances[distances > truncate] = truncate
    return np.atleast_1d(np.squeeze(distances))


def _distanceArray(mol, sel1, sel2, numframes, pbc):
    numsel1 = np.sum(sel1)
    numsel2 = np.sum(sel2)
    coords1 = mol.coords[sel1, :, :]
    coords2 = mol.coords[sel2, :, :]

    distances = []
    for j in range(numsel2):
        coo2 = coords2[j, :, :]  # 3 x numframes array
        dists = coords1 - coo2
        if pbc:
            if mol.box is None or np.sum(mol.box) == 0:
                raise NameError('No periodic box dimensions given in the molecule/trajectory. If you want to calculate distance without wrapping, set the pbc option to False')
            dists = _wrapDistances(mol.box, dists, _findDiffChain(mol, sel1, sel2, j, range(numsel1)))
        dists = np.transpose(np.sqrt(np.sum(dists * dists, 1)))
        distances.append(dists)
    return distances

'''
def _minDistanceArray(mol, sel1, sel2, group1, group2, numframes, pbc):
    distances = _distanceArray(mol, sel1, sel2, numframes, pbc)
    if np.ndim(distances[0]) > 1:  # 2D data
        distances = np.concatenate(distances, axis=1)
    else: # 1D data
        distances = np.vstack(distances).transpose()

    numsel1 = np.sum(sel1)
    numsel2 = np.sum(sel2)
    combo = []
    for s1 in group1:
        s1x = [s1 + (numsel1 * n2) for n2 in range(numsel2)]
        for s2 in group2:
            combo.append()

    newdistances = np.zeros((distances.shape[0], len(combo)))
    for c in combo:
        newdistances[:, c] = np.min(distances[:, c], axis=1)
'''


def _findDiffChain(mol, sel1, sel2, i, others):
    if np.array_equal(sel1, sel2) and len(mol.chain) > 0:
        chain = mol.get('chain', sel=sel1)
        diffchain = chain[others] != chain[i]
    else:
        diffchain = None
    return diffchain


def _wrapDistances(box, dist, diffchain):
    if diffchain is not None:
        numatoms = np.sum(diffchain)
        dist[diffchain, :, :] -= box * np.round(dist[diffchain, :, :] / box)
    else:
        numatoms = np.size(dist, 0)
        dist = dist - box * np.round(dist / box)
    return dist



def exportProjectionData(data, filename):
    """ Export results of a projection into an R-friendly data frame

    The format of the written data is:
      Trajectory Frame   CV1  CV2 ...
      <TrajName> <Frame> <V1> <V2> ...
      ...

    Parameters
    ----------
    data : htmd.metricdata.MetricData
        The results of a metric.project() operation
    filename : str
        The filename to be written.

    """

    out_file = open(filename, "w")

    nTrajs=len(data.simlist)

    if nTrajs==0:
        raise Exception("MetricData does not contain any trajectory")

    (junk,nVars)=data.dat[0].shape
    # Can we recover the mapping?
    # TODO check if combined trajectories work
    out_file.write("\t".join(["TrajName","Frame"]+
                             ["CV"+str(i) for i in range(nVars)]))
    out_file.write("\n")

    for tr in range(nTrajs):
        (nf, junk)=data.dat[tr].shape
        for fr in range(nf):
            fields = [data.simlist[tr].trajectory[0]]
            fields.append(fr)
            fields.extend(data.dat[tr][fr,])
            out_file.write("\t".join(str(el) for el in fields))
            out_file.write("\n")

    out_file.close()

