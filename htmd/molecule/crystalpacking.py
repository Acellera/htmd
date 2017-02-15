# -> Notes: line 69 explanation missing

from htmd.molecule.molecule import Molecule
from htmd.vmdviewer import VMD
import numpy as np
import copy


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


class UnitCell:
    """
    Class containing methods to read crystal information from PDB, generate outline of unit cell and draw crystal copies.
    """

    def __init__(self, pdbfile):
        """ Creates a UnitCell object.

        Parameters
        ----------
        pdbfile : str
            Path to the pdb file which to read.
        """
        '''
        Its attributes upon initialization are the following:
           -> mol: htmd molecule object of the original PDB.
           -> axes: axes of the Unit Cell.
           -> size: furthest coordinates from the origin of the Unit Cell.
           -> molunit: empty molecule object where the Unit Cell will be stored.
           -> viewer: HTMD viewer used.
        Furthermore, it calls the method 'readPDB' upon initializing, also
        acquiring the attributes defined there.'''
        self.mol = Molecule(pdbfile)
        self.readPDB(pdbfile)
        caux = (np.cos(self.alpha) - np.cos(self.beta) * np.cos(self.gamma)) / np.sin(self.gamma)
        self.axes = np.array([[self.a, 0, 0], [self.b * np.cos(self.gamma),
                                               self.b * np.sin(self.gamma), 0],
                              [self.c * np.cos(self.beta), self.c * caux,
                               self.c * np.sqrt(1 - np.cos(self.beta) ** 2 - caux ** 2)
                               ]])
        self.size = self.get_size()
        self.molunit = Molecule()
        self.viewer = VMD()

    def readPDB(self, pdbfile):
        """ Parses PDB file and stores useful data to build a Unit Cell within the instance's attributes.

        Parameters
        ----------
        pdbfile : str
            The path to the PDB file
        """
        '''
           Defines the following attributes:
               -> rot_mat: rotation matrix specified in the PDB file.
               -> trans_v: translation vector specified in the PDB file.
               -> num_copies: number of copies of the molecule per Unit Cell.
                              obtained from computing the length of trans_v.
               -> group: type of Unit Cell.
               -> alpha, beta and gamma: Euler angles of the Unit Cell.
               -> a, b and c: point distances (in Armstrongs) from origin to furthest point.'''
        f = open(pdbfile, "r")
        line = f.readline()
        while not line.startswith("EXPDTA"):
            line = f.readline()
        if (line.rfind("CRYSTAL") == -1) and (line.rfind("DIFFRACTION") == -1):
            raise RuntimeError(
                "The input pdb does not correspond to a crystallographic structure; it was obtained by {}".format(
                    " ".join(line.split()[1:])))
        while line.rfind("REMARK 290       1555") == -1:  # find start crystal structures
            line = f.readline()
        nstr = 0
        while len(line.split()) > 2:  # count structures
            nstr += 1
            line = f.readline()
        self.rot_mat = np.zeros((nstr, 3, 3))  # generate empty rotation matrices
        self.trans_v = np.zeros((nstr, 3))  # generate empty translation vectors
        while line.rfind("SMTRY1") == -1:
            line = f.readline()
        # fill rotation matrices and translation vectors
        n = 0
        dim = 0
        while len(line.split()) > 2:
            values = line.split()
            self.rot_mat[n][dim] = values[4:7]
            self.trans_v[n][dim] = values[7]
            line = f.readline()
            if dim == 2:
                n += 1
                dim = 0
            else:
                dim += 1
        self.num_copies = len(self.trans_v)
        # get axes and angles information
        while not line.startswith("CRYST1"):
            line = f.readline()
        (a, b, c, alpha, beta, gamma) = line.split()[1:7]
        self.group = line.split()[7:-2]
        # angles to radians
        self.alpha = np.deg2rad(float(alpha))
        self.beta = np.deg2rad(float(beta))
        self.gamma = np.deg2rad(float(gamma))
        # axes to float numbers
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)

    def get_size(self):
        """ Returns furthest point from origin to determine Unit Cell size.

        Returns
        -------
        point : np.ndarray
            Furthest point from origin.
        """
        return np.array([self.axes[0][0], self.axes[1][1], self.axes[2][2]])

    def draw_cell(self, draw_hexagon):
        """ Draws lines in the viewer that represent Unit Cell.

        Parameters
        ----------
        draw_hexagon : bool
            If True, draws the Unit Cell with hexagonal boundaries if it is a hexagon. Else it draws a rectangle or
            trapezoid.
        """
        baseorigin = np.array([0, 0, 0])
        base_lo_l = self.axes[0]
        base_up_r = self.axes[1]
        base_lo_r = self.axes[0] + self.axes[1]
        top_origin = self.axes[2]
        top_up_r = self.axes[1] + self.axes[2]
        top_lo_l = self.axes[0] + self.axes[2]
        top_lo_r = self.axes.sum(axis=0)
        points = np.vstack((baseorigin, base_lo_l, base_up_r, base_lo_r))
        points = np.vstack((points, top_origin, top_lo_l, top_up_r, top_lo_r))
        if self.group[0] == 'H' and draw_hexagon == True:
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
            self.viewer.send('draw line {{{0} {1} {2}}} {{{3} {4} {5}}}'.format(
                *(tuple(points[row[0]]) + tuple(points[row[1]]))))

    def build_hexagon(self):
        """
        Draws the copies within the Unit Cell stored in molunit as a hexagon by performing three consecutive rotations
        on the Z-axis with an angle of 120 degrees and places them in a new Molecule Object.
        Creates a new attribute, hexagonal_molunit, which stores the created Unit Cell HTMD molecule object.
        """

        for i in range(1, 4):
            cellunit = copy.deepcopy(self.molunit)  # TODO: This cant possibly work
            alpha = np.deg2rad(120 * i)
            rot_mat = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                                [np.sin(alpha), np.cos(alpha), 0],
                                [0, 0, 1]])
            cellunit.rotateBy(rot_mat)
            self.hexagonal_molunit.append(cellunit)

    def draw_copies(self, draw_all):
        """ Creates copies of the molecule and places them correctly inside the complete Unit Cell

        Parameters
        ----------
        draw_all : bool
            If True the Unit Cell is hexagonal
        """
        for i in range(self.num_copies):
            molecule = self.mol.copy()
            # apply SMTRY (Crystal Symmetry) operations
            molecule.rotateBy(self.rot_mat[i])
            molecule.moveBy(self.trans_v[i])
            # apply translation to inside of same Unit Cell.
            _place_crystal(molecule, self.size, [self.alpha, self.beta, self.gamma], self.axes)
            # pack copies to target Unit Cell
            self.molunit.append(molecule)
        if self.group[0] == 'H' and draw_all == True:
            self.hexagonal_molunit = Molecule()
            self.build_hexagon()

    def view_unit_cell(self, style_display):
        """ Displays the Unit Cell's molecule object in the viewer.

        Parameters
        ----------
        style_display : str
            Display style for the molecule un VMD
        """
        if hasattr(self, 'hexagonal_molunit'):
            self.hexagonal_molunit.view(style=style_display, viewerhandle=self.viewer)
        else:
            self.molunit.view(style=style_display, viewerhandle=self.viewer)


if __name__ == "__main__":
    # ###---VARIABLE DEFINITION---###
    # filename = "./3v03.pdb"
    # display_type = 'NewCartoon'
    # draw_crystal_type = False  # This variable determines whether an hexagonal crystal will be drawn as hexagonal or not.
    # full_build = False  # This variable determines whether, for an hexagonal crystal, the three different crystal pieces will be rotated.
    #
    # ###---MAIN PROCESS---###
    #
    # # First, a UnitCell object is built with the input PDB in filename is created.
    # # Next, the Unit Cell is created.
    # # Lastly, the copies of each molecule are drawn inside the Unit Cell and its visualization is called.
    # unitcell = UnitCell(filename)
    # unitcell.draw_cell(draw_crystal_type)
    # unitcell.draw_copies(full_build)
    # unitcell.view_unit_cell(display_type)
    pass
