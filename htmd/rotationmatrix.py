# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from math import cos, sqrt, sin
import logging
logger = logging.getLogger(__name__)


def rotationMatrix(axis, theta):
    """ Produces a rotation matrix given an axis and radians

    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    Parameters
    ----------
    axis: list
        The axis around which to rotate
    theta: float
        The rotation angle in radians

    Returns
    -------
    M: numpy.ndarray
        The rotation matrix.

    Examples
    --------
    >>> M = rotationMatrix([0, 0, 1], 1.5708)
    >>> M.round(4)
    array([[-0., -1.,  0.],
           [ 1., -0.,  0.],
           [ 0.,  0.,  1.]])

    >>> axis = [4.0, 4., 1.]
    >>> theta = 1.2
    >>> v = [3.0, 5., 0.]
    >>> np.dot(rotationMatrix(axis, theta), v).round(2)
    array([ 2.75,  4.77,  1.92])
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


if __name__ == '__main__':

    import doctest
    doctest.testmod()
