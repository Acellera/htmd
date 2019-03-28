# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from numba import jit
from math import sqrt, atan2


@jit(nopython=True)
def cross(vec1, vec2):
    """ Calculate the dot product of two 3d vectors. """
    a1, a2, a3 = vec1[0], vec1[1], vec1[2]
    b1, b2, b3 = vec2[0], vec2[1], vec2[2]
    result = np.zeros(3)
    result[0] = a2 * b3 - a3 * b2
    result[1] = a3 * b1 - a1 * b3
    result[2] = a1 * b2 - a2 * b1
    return result


@jit(nopython=True)
def dot(vec1, vec2):
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]


@jit(nopython=True)
def norm(vec):
    n = 0
    for i in range(3):
        n += vec[i] * vec[i]
    return sqrt(n)


@jit(nopython=True)
def pairwiseRMSD(coords):
    nframes = coords.shape[2]
    pairwise = np.zeros(int(nframes*(nframes-1)/2), dtype=np.float32)

    k = 0
    for i in range(nframes):
        for j in range(i+1, nframes):
            pairwise[k] = np.sqrt(np.mean((coords[:, :, i] - coords[:, :, j]) ** 2))
            k +=1
    return pairwise
