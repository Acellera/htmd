# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from scipy.sparse import lil_matrix
import numpy as np
from IPython.core.debugger import Tracer


def pp_calcDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
    distances = _distanceArray(mol, sel1, sel2, pbc)
    distances = _postProcessDistances(distances, sel1, sel2, truncate)

    if metric == 'contacts':
        #metric = lil_matrix(distances <= threshold)
        metric = distances <= threshold
    elif metric == 'distances':
        metric = distances.astype(dtype=np.float32)
    else:
        raise NameError('The metric you asked for is not supported. Check spelling and documentation')
    return metric


# def pp_calcMinDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
#     from scipy.spatial.distance import cdist
#     if pbc:
#         if mol.box is None or np.sum(mol.box) == 0:
#             raise NameError(
#                 'No periodic box dimensions given in the molecule/trajectory. If you want to calculate distance without wrapping, set the pbc option to False')
#         coords = _wrapCoords(mol.coords, mol.box)
#     else:
#         coords = mol.coords
#
#     mindist = np.zeros((mol.numFrames, len(sel1) * len(sel2)))
#     for i, s1 in enumerate(sel1):
#         for j, s2 in enumerate(sel2):
#             #print(i*len(sel1)+j)
#             for f in range(mol.numFrames):
#                 mindist[f, i*len(sel1)+j] = np.min(cdist(coords[s1, :, f], coords[s2, :, f])[:])
#
#     if metric == 'contacts':
#         mindist = mindist <= threshold
#     elif metric == 'distances':
#         mindist = mindist.astype(dtype=np.float32)
#     else:
#         raise NameError('The metric you asked for is not supported. Check spelling and documentation')
#     return mindist


# def pp_calcMinDistances_C(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
#     import os
#     import ctypes
#     from htmd.home import home
#     if pbc:
#         if mol.box is None or np.sum(mol.box) == 0:
#             raise NameError(
#                 'No periodic box dimensions given in the molecule/trajectory. If you want to calculate distance without wrapping, set the pbc option to False')
#         coords = _wrapCoords(mol.coords, mol.box)
#     else:
#         coords = mol.coords
#
#     # Converting from 2D boolean atomselect array to 2D int array where each row starts with the indexes of the boolean
#     groups1 = np.ones((sel1.shape[0], mol.numAtoms), dtype=np.int32) * -1
#     groups2 = np.ones((sel2.shape[0], mol.numAtoms), dtype=np.int32) * -1
#     for i in range(sel1.shape[0]):
#         idx = np.where(sel1[i, :])[0]
#         groups1[i, 0:len(idx)] = idx
#     for i in range(sel2.shape[0]):
#         idx = np.where(sel2[i, :])[0]
#         groups2[i, 0:len(idx)] = idx
#
#     # Running the actual calculations
#     lib = ctypes.cdll.LoadLibrary(os.path.join(home(), 'projections', 'mindist', 'mindist_ext.so'))
#     mindist = np.zeros((mol.numFrames, len(groups1) * len(groups2)))
#
#     for f in range(mol.numFrames):
#         # print('Frame {}'.format(f))
#         dist = np.zeros((len(groups1) * len(groups2),), dtype=np.float32)
#         coordsframecopy = coords[:, :, f].astype(np.float32)  # This is critically important! Otherwise trajectory slice is not copied and indexing fails in C
#         lib.mindist_single_frame(coordsframecopy.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
#                                  groups1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
#                                  groups2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
#                                  ctypes.c_int(len(groups1)),
#                                  ctypes.c_int(len(groups2)),
#                                  ctypes.c_int(mol.numAtoms),
#                                  dist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
#         # print(dist)
#         mindist[f, :] = dist
#
#     if metric == 'contacts':
#         mindist = mindist <= threshold
#     elif metric == 'distances':
#         mindist = mindist.astype(dtype=np.float32)
#     else:
#         raise NameError('The metric you asked for is not supported. Check spelling and documentation')
#     return mindist


def pp_calcMinDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
    import os
    import ctypes
    from htmd.home import home

    # Converting non-grouped boolean atomselection to group-style atomselections
    if np.ndim(sel1) != 2:
        sel1idx = tuple(np.where(sel1)[0])
        sel1 = np.zeros((len(sel1idx), len(sel1)), dtype=bool)
        sel1[range(sel1.shape[0]), sel1idx] = True
    if np.ndim(sel2) != 2:
        sel2idx = tuple(np.where(sel2)[0])
        sel2 = np.zeros((len(sel2idx), len(sel2)), dtype=bool)
        sel2[range(sel2.shape[0]), sel2idx] = True

    box = np.array([0, 0, 0], dtype=np.float32)
    if pbc:
        if mol.box is None or np.sum(mol.box) == 0:
            raise NameError('No periodic box dimensions given in the molecule/trajectory. '
                            'If you want to calculate distance without wrapping, set the pbc option to False')
        box = mol.box[:, 0]  # TODO: make it work for varying box size
        if np.max(mol.box.T - mol.box[:, 0]) != 0:
            raise NameError('Different box sizes per frame. Still unsupported by mindist. Contact Stefan Doerr.')

    coords = mol.coords

    # Converting from 2D boolean atomselect array to 2D int array where each row starts with the indexes of the boolean
    groups1 = np.ones((sel1.shape[0], mol.numAtoms), dtype=np.int32) * -1
    groups2 = np.ones((sel2.shape[0], mol.numAtoms), dtype=np.int32) * -1
    for i in range(sel1.shape[0]):
        idx = np.where(sel1[i, :])[0]
        groups1[i, 0:len(idx)] = idx
    for i in range(sel2.shape[0]):
        idx = np.where(sel2[i, :])[0]
        groups2[i, 0:len(idx)] = idx

    # Running the actual calculations
    lib = ctypes.cdll.LoadLibrary(os.path.join(home(libDir=True), 'mindist_ext.so'))
    mindist = np.zeros((mol.numFrames, len(groups1) * len(groups2)), dtype=np.float32)  # Preparing the return array
    lib.mindist_trajectory(coords.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                           box.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                           groups1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           groups2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           ctypes.c_int(len(groups1)),
                           ctypes.c_int(len(groups2)),
                           ctypes.c_int(mol.numAtoms),
                           ctypes.c_int(mol.numFrames),
                           ctypes.c_int(int(pbc)),
                           mindist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))

    if metric == 'contacts':
        mindist = mindist <= threshold
    elif metric == 'distances':
        mindist = mindist.astype(dtype=np.float32)
    else:
        raise NameError('The metric you asked for is not supported. Check spelling and documentation')
    return mindist


def _wrapCoords(coords, box):
    return coords - box * np.round(coords / box)


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


def _distanceArray(mol, sel1, sel2, pbc):
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


def _findDiffChain(mol, sel1, sel2, i, others):
    if np.array_equal(sel1, sel2) and len(mol.chain) > 0:
        chain = mol.get('chain', sel=sel1)
        diffchain = chain[others] != chain[i]
    else:
        diffchain = None
    return diffchain


def _wrapDistances(box, dist, diffchain):
    if diffchain is not None:
        dist[diffchain, :, :] -= box * np.round(dist[diffchain, :, :] / box)
    else:
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


def convertProjectionDataToPandas(md):
    """ Export results of a projection into a pandas data frame

    The format of the returned data is:
      Trajectory Frame   CV1  CV2 ...
      <TrajName> <Frame> <V1> <V2> ...
      ...

    Parameters
    ----------
    md : htmd.metricdata.MetricData
        The results of a metric.project() operation

    """

    import pandas as pd

    nTrajs=len(md.simlist)
    if nTrajs==0:
        raise Exception("MetricData does not contain any trajectory")

    dflist = []

    for tr in range(nTrajs):
        df0 = pd.DataFrame(md.dat[tr])
        df0.insert(0,'Trajectory', md.simlist[tr].trajectory[0])
        df0.insert(1,'Frame',range(len(df0)))
        dflist.append(df0)

    return pd.concat(dflist)



