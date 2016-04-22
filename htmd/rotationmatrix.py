# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from math import cos, sqrt, sin
import logging
logger = logging.getLogger(__name__)


def rotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta / 2)
    b, c, d = -axis * sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def rotationmatrix(axis, theta):
    logger.warning('The rotationmatrix method is deprecated. '
                   'It has been renamed to rotationMatrix to follow the naming style. Please change all uses.')
    return rotationMatrix(axis, theta)


if __name__ == '__main__':
    v = [3.0, 5., 0.]
    axis = [4.0, 4., 1.]
    theta = 1.2
    print(np.dot(rotationMatrix(axis, theta), v))
