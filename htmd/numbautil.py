import numpy as np
from numba import jit
from math import sqrt, atan2


@jit('float32(float32, float32)', nopython=True)
def wrapDistance(d, box):
    if box == 0:
        return d
    return d - box * round(d / box)


@jit('float32(float32, float32)', nopython=True)
def wrapBondedDistance(d, box):
    if box == 0:
        return d
    # Assuming bonds can't cross multiple periodic boxes this is marginally faster
    hbox = box / 2
    if d < -hbox:
        return d + box
    elif d > hbox:
        return d - box
    return d


def dihedralAngle(pos, box=None):
    """ Calculates a dihedral angle.

    Parameters
    ----------
    pos: np.ndarray
        An array of 4x3 size where each row are the coordinates of an atom defining the dihedral angle
    box: np.ndarray
        The size of the periodic box

    Returns
    -------
    angle: float
        The angle in radians
    """
    if pos.ndim == 3 and pos.shape[2] > 1:
        if box is None:
            box = np.zeros((3, pos.shape[2]), dtype=pos.dtype)
        return dihedralAngleFrames(pos, box)
    else:
        if box is not None:
            box = box.squeeze()
        return dihedralAngleFull(pos.squeeze(), box)[0]

@jit(nopython=True)
def dihedralAngleFrames(pos, box):
    res = np.zeros(pos.shape[2], dtype=pos.dtype)
    for f in range(pos.shape[2]):
        res[f] = dihedralAngleFull(pos[:, :, f], box[:, f])[0]
    return res

@jit(nopython=True)
def dihedralAngleFull(pos, box=None):
    """ Calculates a dihedral angle.

    Parameters
    ----------
    pos: np.ndarray
        An array of 4x3 size where each row are the coordinates of an atom defining the dihedral angle
    box: np.ndarray
        The size of the periodic box
    """
    if pos.shape[0] != 4 or pos.shape[1] != 3:
        raise RuntimeError('dihedralAngles requires a 4x3 sized coordinate matrix as input.')
    if box is None:
        box = np.zeros(3, dtype=pos.dtype)

    r12 = np.zeros(3)
    r23 = np.zeros(3)
    r34 = np.zeros(3)

    r12[0] = wrapBondedDistance(pos[0, 0] - pos[1, 0], box[0])
    r12[1] = wrapBondedDistance(pos[0, 1] - pos[1, 1], box[1])
    r12[2] = wrapBondedDistance(pos[0, 2] - pos[1, 2], box[2])
    r23[0] = wrapBondedDistance(pos[1, 0] - pos[2, 0], box[0])
    r23[1] = wrapBondedDistance(pos[1, 1] - pos[2, 1], box[1])
    r23[2] = wrapBondedDistance(pos[1, 2] - pos[2, 2], box[2])
    r34[0] = wrapBondedDistance(pos[2, 0] - pos[3, 0], box[0])
    r34[1] = wrapBondedDistance(pos[2, 1] - pos[3, 1], box[1])
    r34[2] = wrapBondedDistance(pos[2, 2] - pos[3, 2], box[2])

    # A = cross(r12, r23)
    A = np.zeros(3)
    A[0] = r12[1] * r23[2] - r12[2] * r23[1]
    A[1] = r12[2] * r23[0] - r12[0] * r23[2]
    A[2] = r12[0] * r23[1] - r12[1] * r23[0]

    # B = cross(r23, r34)
    B = np.zeros(3)
    B[0] = r23[1] * r34[2] - r23[2] * r34[1]
    B[1] = r23[2] * r34[0] - r23[0] * r34[2]
    B[2] = r23[0] * r34[1] - r23[1] * r34[0]

    # C = cross(r23, A)
    C = np.zeros(3)
    C[0] = r23[1] * A[2] - r23[2] * A[1]
    C[1] = r23[2] * A[0] - r23[0] * A[2]
    C[2] = r23[0] * A[1] - r23[1] * A[0]

    rA = 1 / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2])
    rB = 1 / sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2])
    rC = 1 / sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2])

    B[0] *= rB
    B[1] *= rB
    B[2] *= rB

    cos_phi = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]) * rA
    sin_phi = (C[0] * B[0] + C[1] * B[1] + C[2] * B[2]) * rC

    phi = -atan2(sin_phi, cos_phi)

    return phi, r12, r23, r34, A, B, C, rA, rB, rC, sin_phi, cos_phi


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