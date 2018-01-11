from math import cos, sqrt, sin
import numpy as np


def get_rotationMatrix(axis, theta):
    """ Generates a rotation matrix given an axis and radians
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


def rotate(coords, rotMat, center=(0,0,0)):
    """ Rotate a selection of atoms by a given rotation around a center
    Parameters
    ----------
    coords : np.ndarray
        Coordinates to rotate
    M : np.ndarray
        The rotation matrix
    center : list
        The rotation center
    sel :
        Atomselection for atoms to rotate
    """

    newcoords = coords - center
    return np.dot(newcoords, np.transpose(rotMat)) + center


def drawIsoSurface(values3d, resolution=1., plot_center=None, viewer=None):
    from htmd.vmdviewer import getCurrentViewer
    from htmd.molecule.util import writeVoxels
    # plot_center should be - molecule.get_center() + 12
    if len(values3d.shape) != 3:
        raise ValueError("Your provided a box of {} dimensions."
                         "\nThis only works with dimension of 3".format(len(values3d.shape)))
    from htmd.util import tempname
    if viewer is None:
        viewer = getCurrentViewer()
    mincoor = np.zeros(3, dtype=np.float64)
    maxcoor = np.array(values3d.shape, dtype=np.float64)
    rescoor = np.array([resolution] * 3)

    # Adjust the plotting center
    if plot_center is None:
        plot_center = maxcoor / 2.
    else:
        plot_center = np.array(plot_center)
    mincoor -= (plot_center + 0.5)  # TODO: Fix so it will work in case resolution != 1.
    maxcoor -= (plot_center + 0.5)

    outf = tempname(suffix='.cube')
    writeVoxels(values3d, outf, mincoor, maxcoor, rescoor)
    viewer.send('mol new {} type cube first 0 last -1 step 1 waitfor 1 volsets {{0 }}'.format(outf))
    viewer.send('mol modstyle 0 top Isosurface 0.75 0 2 0 1 1')