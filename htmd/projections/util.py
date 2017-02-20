# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from scipy.sparse import lil_matrix
import numpy as np
from IPython.core.debugger import Tracer
import logging
import htmd.molecule.molecule
import htmd.progress.progress

logger = logging.getLogger(__name__)


def pp_calcDistances(mol, sel1, sel2, metric='distances', threshold=8, pbc=True, gap=1, truncate=None):
    distances = _distanceArray(mol, sel1, sel2, pbc)
    distances = _postProcessDistances(distances, sel1, sel2, truncate)

    if metric == 'contacts':
        # metric = lil_matrix(distances <= threshold)
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

    selfdist = np.array_equal(sel1, sel2)

    # Running the actual calculations
    lib = ctypes.cdll.LoadLibrary(os.path.join(home(libDir=True), 'mindist_ext.so'))
    mindist = np.zeros((mol.numFrames, len(groups1) * len(groups2)), dtype=np.float32)  # Preparing the return array
    if selfdist:
        mindist = np.zeros((mol.numFrames, int((len(groups1) * (len(groups2)-1))/2)), dtype=np.float32)

    #import time
    #t = time.time()
    lib.mindist_trajectory(coords.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                           box.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                           groups1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           groups2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           ctypes.c_int(len(groups1)),
                           ctypes.c_int(len(groups2)),
                           ctypes.c_int(mol.numAtoms),
                           ctypes.c_int(mol.numFrames),
                           ctypes.c_int(int(pbc)),
                           ctypes.c_int(int(selfdist)),
                           mindist.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
    #print(time.time() - t)

    if truncate is not None:
        mindist[mindist > truncate] = truncate

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
    # distances is a list of numpy arrays. Each numpy array is numFrames x numSel1. The list is length numSel2
    # Setting upper triangle to -1 if same selections
    if np.array_equal(sel1, sel2):
        for i in range(len(distances)):
                distances[i][:, range(i + 1)] = -1

    if np.ndim(distances[0]) > 1:  # 2D data
        distances = np.concatenate(distances, axis=1)
    else:  # 1D data
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
                raise NameError(
                    'No periodic box dimensions given in the molecule/trajectory. If you want to calculate distance without wrapping, set the pbc option to False')
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




def convertProjectionToDataFrame(md):
    """ Export results of a projection into a pandas data frame

    The format of the returned data contains:
      - TrajectoryID     (also index)
      - TrajectoryFrame  (also index)
      - Piece
      - PieceFile
      - PieceFrame
      - Columns containing the values

    Parameters
    ----------
    md : htmd.metricdata.MetricData
        The results of a metric.project() operation

    Returns
    -------
    df : pandas.DataFrame
        A DataFrame containing the results of the projection

    """

    import pandas as pd

    nTrajs = len(md.simlist)
    if nTrajs == 0:
        raise Exception("MetricData does not contain any trajectory")

    bar = htmd.progress.progress.ProgressBar(nTrajs, description="Converting {:d} trajectories".format(nTrajs))
    dflist = []

    curf=0
    for tr in range(nTrajs):
        df0 = pd.DataFrame(md.dat[tr])
        nf = len(df0)
        nfl = curf+np.array(range(nf))
        nfs = md.abs2sim(nfl)
        df0.insert(0, 'TrajectoryID', tr)
        # df0.insert(1, 'TrajectoryFile', md.simlist[tr].trajectory[0])
        df0.insert(1, 'TrajectoryFrame', range(len(df0)))
        df0.insert(2, 'Piece', [x.piece for x in nfs])
        df0.insert(3, 'PieceFile', [x.sim.trajectory[x.piece] for x in nfs])
        df0.insert(4, 'PieceFrame', [x.frame for x in nfs])
        dflist.append(df0)
        curf += nf
        bar.progress()

    df = pd.concat(dflist)

    df['PieceFile']=df['PieceFile'].astype('category')

    df.set_index(["TrajectoryID", "TrajectoryFrame"],
                 drop=False, inplace=True, verify_integrity=True)
    return df


def readSimlistIndices(prj, selector):
    """Convert a list of boolean values to a Molecule containing frames from the given simlist.

    Limitation: all trajectories must have the same (or compatible) topology.

    Parameters
    ----------
    prj : htmd.metricdata.MetricData
        The results of a metric.project() operation
    selector : list
        List of boolean values, i.e. whether to extract that frame.

    Returns
    -------
    mol : Molecule
        A molecule containing the frames selected.

    """
    idx = [i for i, x in enumerate(selector) if x]
    frs = prj.abs2sim(idx)
    nf = len(frs)

    molfile = prj.simlist[0].molfile
    mol = htmd.Molecule(molfile)
    mol.dropFrames([])

    tset = set()
    i = 0
    bar = htmd.progress.progress.ProgressBar(nf, description="Reading {:d} frames".format(nf))

    for i, f in enumerate(frs):
        tn = f.sim.trajectory[f.piece]
        tset.add(tn)
        bar.progress()
        # print("Read frame {:d} from trajectory {:s}, frame {:d}".format(i, tn, f.frame))
        mol.read(filename=tn,
                 frames=f.frame,
                 append=True)

    logger.info("Read {:d} frames from {:d} distinct trajectories".format(i + 1, len(tset)))
    return mol
