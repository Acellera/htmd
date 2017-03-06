# -> Notes: line 69 explanation missing

from htmd.molecule.molecule import Molecule
from htmd.vmdviewer import VMD
import numpy as np
import copy
import os
import logging
logger = logging.getLogger(__name__)


def _get_cellloc(center, angles, box_size):
    '''Computes the Unit Cell's location of the molecule object.
    Takes as argument the Euler angles of the Unit Cell's sides and performs
    a small correction to account for oblique sides.
    Allows to correct x position dependent on y and z.'''
    change = np.zeros(3)
    change[0] = center[1] / np.tan(angles[2]) + center[2] / np.tan(
        angles[1])  # get min_value x unit cell at given y,z
    return np.floor((center - change) / box_size)  # correct x axis boundaries according to angles


def _is_outside(center, cellloc):
    '''Checks if a molecule object is inside the defined unit cell.
    Takes as arguments the size vector [numpy array with 3 coordinates]
    and a vector [numpy array with three floats] with the Euler angles of
    the Unit Cell.
    Modifies the attribute:
           -> center: mean center of the molecule object.

    Creates the attributes:
           -> cellloc: location of the molecule object within the Unit Cell.
    Returns True if the molecule object is outside the Unit Cell, else
    returns False.'''
    return np.sum(np.abs(cellloc)) > 0, cellloc, center


def _move_inside(mol, axes, cellloc, center):
    '''Translates an HTMD molecule object inside another Unit Cell.
    Takes as argument the vectors defining the Unit Cell the molecule
    object is currently in and translates it inside another Unit Cell.
    '''
    neworigin = np.multiply(axes.transpose(), cellloc).transpose()  # computes translation to be applied to copy
    neworigin = np.sum(neworigin, axis=0)
    newcenter = center - neworigin  # computes difference between copy center and origin of its unit cell
    mol.center(loc=newcenter)  # places copy in the target unit cell (at an equivalent position to that of original unit cell)


def _place_crystal(mol, size, angles, axes):
    '''Places MoleculeCopy object inside crystal.
    Uses methods 'is_outside' and 'move_inside'.'''
    center = np.mean(mol.get('coords'), axis=0)
    cellloc = _get_cellloc(center, angles, size)  # integer from division for each axis = unit cell identifier
    if _is_outside(center, cellloc):
        _move_inside(mol, axes, cellloc, center)


def viewCrystalPacking(mol, hexagonal=False, style_display='NewCartoon'):
    """ Views the crystal packing of a protein

    Parameters
    ----------
    pdbfile : str
        Path to the pdb file which to read. Can also be a 4-letter PDB ID
    """
    if mol.crystalinfo is None or 'numcopies' not in mol.crystalinfo:
        raise RuntimeError('No crystallography data found in Molecule.')
    ci = mol.crystalinfo

    alpha, beta, gamma, a, b, c = ci['alpha'], ci['beta'], ci['gamma'], ci['a'], ci['b'], ci['c']
    alpha = np.deg2rad(float(alpha))
    beta = np.deg2rad(float(beta))
    gamma = np.deg2rad(float(gamma))

    caux = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    axes = np.array([[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0], [c * np.cos(beta), c * caux,
                           c * np.sqrt(1 - np.cos(beta) ** 2 - caux ** 2)]])
    size = np.array([axes[0][0], axes[1][1], axes[2][2]])

    molunit = Molecule()
    viewer = VMD()

    _draw_cell(axes, ci['sGroup'], viewer, hexagonal=hexagonal)

    # Creates copies of the molecule and places them correctly inside the complete Unit Cell
    hexagonal_molunit = None
    for i in range(ci['numcopies']):
        molecule = mol.copy()
        # apply SMTRY (Crystal Symmetry) operations
        molecule.rotateBy(ci['rotations'][i])
        molecule.moveBy(ci['translations'][i])
        # apply translation to inside of same Unit Cell.
        _place_crystal(molecule, size, [alpha, beta, gamma], axes)
        # pack copies to target Unit Cell
        molunit.append(molecule)
    if ci['sGroup'][0] == 'H' and hexagonal:
        hexagonal_molunit = Molecule()
        _build_hexagon(molunit, hexagonal_molunit)

    if hexagonal_molunit is not None:
        hexagonal_molunit.view(style=style_display, viewerhandle=viewer)
    else:
        molunit.view(style=style_display, viewerhandle=viewer)


def _draw_cell(axes, group, viewer, hexagonal=False):
    """ Draws lines in the viewer that represent Unit Cell.

    Parameters
    ----------
    draw_hexagon : bool
        If True, draws the Unit Cell with hexagonal boundaries if it is a hexagon. Else it draws a rectangle or
        trapezoid.
    """
    baseorigin = np.array([0, 0, 0])
    base_lo_l = axes[0]
    base_up_r = axes[1]
    base_lo_r = axes[0] + axes[1]
    top_origin = axes[2]
    top_up_r = axes[1] + axes[2]
    top_lo_l = axes[0] + axes[2]
    top_lo_r = axes.sum(axis=0)
    points = np.vstack((baseorigin, base_lo_l, base_up_r, base_lo_r))
    points = np.vstack((points, top_origin, top_lo_l, top_up_r, top_lo_r))
    if group[0] == 'H' and hexagonal:
        # Lattice is hexagonal, can be constructed by applying two 120º
        # rotations over the z axis
        alpha = np.deg2rad(120)
        rot_mat = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])
        rotpoints1 = np.dot(rot_mat, points.transpose())
        rotpoints2 = np.dot(rot_mat, rotpoints1)
        points = np.vstack((points, rotpoints1.transpose(),
                            rotpoints2.transpose()))
        lines_draw = np.array([[1, 3], [1, 5], [2, 3], [2, 6], [3, 7],
                               [5, 7], [6, 7]])
        lines_draw = np.vstack((lines_draw, lines_draw + 8, lines_draw + 16))
    else:
        lines_draw = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3],
                      [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]]
    for row in lines_draw:
        viewer.send('draw line {{{0} {1} {2}}} {{{3} {4} {5}}}'.format(
            *(tuple(points[row[0]]) + tuple(points[row[1]]))))


def _build_hexagon(molunit, hexagonal_molunit):
    """
    Draws the copies within the Unit Cell stored in molunit as a hexagon by performing three consecutive rotations
    on the Z-axis with an angle of 120 degrees and places them in a new Molecule Object.
    Creates a new attribute, hexagonal_molunit, which stores the created Unit Cell HTMD molecule object.
    """

    for i in range(1, 4):
        cellunit = copy.deepcopy(molunit)  # TODO: This cant possibly work
        alpha = np.deg2rad(120 * i)
        rot_mat = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])
        cellunit.rotateBy(rot_mat)
        hexagonal_molunit.append(cellunit)
