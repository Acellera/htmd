# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
# from numba import double
# from numba.decorators import jit


# @jit('void(f8[:,:], f8[:,:], f8[:,:], f8[:])', nopython=True, target='cpu')
def wrapping_dist_numba(X, Y, D, box):
    numX = X.shape[0]
    numY = Y.shape[0]
    numDim = X.shape[1]
    for i in range(numX):
        for j in range(numY):
            d = 0.0
            for k in range(numDim):
                dist = X[i, k] - Y[j, k]
                tmp = dist - box[k] * round(dist / box[k])
                d += tmp * tmp
            D[i, j] = np.sqrt(d)


def wrapping_dist_python(coor1, coor2, box):
    assert (coor1.ndim == 1) or (coor2.ndim == 1)
    dist = coor1 - coor2
    dist = dist - box * np.round(dist / box)
    return np.sqrt(np.sum(dist * dist, 1))


if __name__ == '__main__':
    pass
    # from wrappingdist import *
    # import time
    # from scipy.spatial.distance import cdist

    # allcoords = np.load('/tmp/data1.npy')
    # neighcoor = np.load('/tmp/data2.npy')
    # box = np.load('/tmp/box.npy')

    # #allcoords = np.random.random((10000, 3))
    # #neighcoor = np.random.random((10000, 3))

    # t = time.time()
    # for ac in allcoords:
    #     dists = wrapping_dist_python(ac, neighcoor, box)
    # print('python', time.time() - t)
    # t = time.time()
    # dists = cdist(allcoords, neighcoor)
    # print('cdist', time.time() - t)
    # t = time.time()
    # dists = np.zeros((allcoords.shape[0], neighcoor.shape[0]))
    # wrapping_dist_numba(allcoords, neighcoor, dists, box.astype(float))
    # print('numba', time.time() - t)