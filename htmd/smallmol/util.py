# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import math
import numpy as np
from htmd.molecule.voxeldescriptors import _getGridCenters
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig


_highlight_colors = [(1.00,0.50,0.00), (0.00,0.50,1.00), (0.00,1.00,0.50),
                     (1.00,0.00,0.50), (0.50,0.00,1.00), (0.50,1.00,0.00),
                     (1.00,0.00,0.25), (0.00,0.25,1.00), (0.25,1.00,0.00)]

array_cache = {(16, 1.): _getGridCenters(np.array([-8] * 3),
                                               [16, 16, 16], 1.).reshape(16**3, 3),
                   (24, 1.): _getGridCenters(np.array([-12] * 3),
                                               [24, 24, 24], 1.).reshape(24**3, 3)}



fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


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
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def _getRotationMatrix(axis, theta, deg=False):
    """
    Rotational matrix used by Builder
    axis: np.array
        The axis of rotation
    theta: float
        The angle of the rotation
    deg: bool
        If True, the angle is considered in degree and will be converted into radians
        Default: False

    Returns
    -------
    rotmatrix: np.array
        The rotational matrix
    """
    # Development note: I will test the above rotationl matrix. If same behaviour this one will be deprecated

    if deg:
        theta = math.radians(theta)

    c, s = math.cos(theta), math.sin(theta)

    t = 1 - c
    X = axis[0]
    Y = axis[1]
    Z = axis[2]

    return np.array([ [t*X*X+c, t*X*Y-s*Z, t*X*Z+s*Y],
                      [t*X*Y+s*Z, t*Y*Y+c, t*Y*Z-s*X],
                      [t*X*Z-s*Y, t*Y*Z+s*X, t*Z*Z+c] ])

def _normalizeVector(v):
    """
    Returns the normalize vector

    Parameters
    ----------
    v: np.array
        The vector to normalize

    Returns
    -------
    V: np.array
        The normalized vector
    """
    from math import sqrt
    l = sqrt(v[0] ** 2 + v[1] ** 2. + v[2] ** 2)

    if l == 0:
        print('warning normalized vector goes to 0')
        return np.array([0,0,0])
    return v / l

def _getPerpendicular(v):
    """
    Returns a perpendicular vector

    Parameters
    ----------
    v: np.array
        The vector you want a perpendicular one

    Returns
    -------
    V:  np.array
        The perpendicular vector
    """
    V = [0,0,0]
    X, Y, Z = v[0], v[1], v[2]
    if X != 0:
        if Y != 0:
            V[0] = Y
            V[1] = -X
        elif Z != 0:
            V[2] = -X
            V[0] = Z
        else:
            V[1] = 1
    elif Y != 0:
        if Z != 0:
            V[2] = -Y
            V[1] = Z
        else:
            V[0] = 1
    elif Z !=0:
        V[0] = 1

    return np.array(V)

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
        Default: [0,0,0]
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


def alignMol(smallmol, refmol):
    """
    Return a new SmallMol object aligned to a refmol that can be a htmd.smallmol.smallmol.SmallMol
    or rdkit.Chem.rdchem.Mol. It removes all the conformers stored in the original object.


    Parameters
    ----------
    smallmol: htmd.smallmol.smallmol.SmallMol
     The SmallMol object to align
    refmol: htmd.smallmol.smallmol.SmallMol or rdkit.Chem.rdchem.Mol
        The molecule to align to

    Return
    ------
    newsmallmol: htmd.smallmol.smallmol.SmallMol
        a new SmallMol aligned to reference molecule
    """

    from htmd.smallmol.smallmol import SmallMol
    from rdkit.Chem.rdMolAlign import GetO3A

    if isinstance(refmol, SmallMol):
        refmol = refmol.toRdkitMol(includeConformer=True)

    sm_rdkit = smallmol.toRdkitMol(includeConformer=True)

    pyO3A = GetO3A(sm_rdkit, refmol)
    rmsd = pyO3A.Align()
    print('Alignment with a RMSD of {}'.format(rmsd))
    coords_new = sm_rdkit.GetConformer().GetPositions()

    sm_new = SmallMol(smallmol, fixHs=False)
    sm_new.removeConformers()
    sm_new.coords = coords_new[:, :, np.newaxis]

    return sm_new

def inputToOutput(input_file, input_format, output_format):
    """
    Converts the file from the input format to the output format specified. It uses the openbabel features

    Parameters
    ----------
    input_file: str
        The path of the input file to convert
    input_format: str
        The input file format
    output_format: str
        The output file format

    Returns
    -------
    outfile: str
        The output file generated
    """

    import openbabel
    import tempfile
    input_format = input_format[1:] if input_format.startswith('.') else input_format

    file = tempfile.NamedTemporaryFile(delete=True, suffix='.' + output_format)
    file.close()
    outfile = file.name

    _ = openbabel.OBConversion()
    _.SetInAndOutFormats(input_format, output_format)
    _mol = openbabel.OBMol()
    _.ReadFile(_mol, input_file)
    _.WriteFile(_mol, outfile)

    return outfile


def convertToString(arr):

    if isinstance(arr, list):
        arr_str = " ".join([str(i) for i in arr])
    elif isinstance(arr, tuple):
        arr_str = " ".join([str(i) for i in arr])
    else:
        arr_str = " ".join([ str(i) for i in arr[0] ])

    return arr_str

def _depictMol(mol, filename=None, ipython=False, atomlabels=None, highlightAtoms=None):
    """
    Returns the image or the ipython rendering.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol
        The rdkit molecule to depict
    filename: str
        The filename of the image
        Default: None
    ipython: bool
        If True, the SVG rendering for jupiter-nootebook are returned
        Default: False
    atomlabels: list
        List of the label to use for each atom
        Default: None
    highlightAtoms: list
        List of atom index to highlight. Can be also list of list for different selection-colors
        Default: None

    Returns
    -------
    svg: SVG
        If ipython set as True, the SVG rendering is returned

    """
    from os.path import splitext
    from rdkit.Chem.Draw import rdMolDraw2D
    from IPython.display import SVG

    if highlightAtoms is not None and not isinstance(highlightAtoms, list):
        raise ValueError('highlightAtoms should be a list of atom idx or a list of atom idx list ')

    # init the drawer object
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
    # get the drawer options
    opts = drawer.drawOptions()

    # add atomlabels
    if atomlabels is not None:
        for n, a in enumerate(atomlabels):
            opts.atomLabels[n] = a

    # draw molecule
    sel_atoms = []
    sel_colors = {}
    # highlight atoms
    if highlightAtoms is not None:
        if isinstance(highlightAtoms[0], list ):
            sel_atoms = [aIdx for subset in highlightAtoms for aIdx in subset]
            sel_colors = {aIdx: _highlight_colors[n%len(_highlight_colors)] for n, subset in enumerate(highlightAtoms) for aIdx in subset}
        else:
            sel_atoms = highlightAtoms
            sel_colors = { aIdx:_highlight_colors[0] for aIdx in sel_atoms }

    drawer.DrawMolecule(mol, highlightAtoms=sel_atoms, highlightBonds=[], highlightAtomColors=sel_colors)

    drawer.FinishDrawing()

    # svg object
    svg = drawer.GetDrawingText()

    # activate saving into a file
    if filename != None:
        ext = splitext(filename)[-1]
        filename = filename if ext != '' else filename + '.svg'
        f = open(filename, 'w')
        f.write(svg)
        f.close()

    # activate jupiter-notebook rendering
    if ipython:
        svg = svg.replace('svg:', '')
        return SVG(svg)
    else:
        return None

def depictMultipleMols(mols_list, filename=None, ipython=False, legends=None, highlightAtoms=None, mols_perrow=3):
    """
        Returns the image or the ipython rendering.

        Parameters
        ----------
        mols_list: list
            The list of the rdkit molecules to depict
        filename: str
            The filename of the image
            Default: None
        ipython: bool
            If True, the SVG rendering for jupiter-nootebook are returned
            Default: False
        legends: list
            List of titles subfigure for each molecule
            Default: None
        highlightAtoms: list
            List of list of atom index to highlight.
            Default: None
        mols_perrow: int
            The number of subfigures per row
            Default: 3

        Returns
        -------
        svg: SVG
            If ipython set as True, the SVG rendering is returned

        """
    from rdkit.Chem.Draw import MolsToGridImage
    from IPython.display import SVG
    from os.path import splitext

    sel_atoms = []
    sel_colors = []
    if highlightAtoms is not None:
        if isinstance(highlightAtoms[0][0], list):
            sel_atoms = [ [a for a in subset] for mol_set in highlightAtoms for subset in mol_set ]
            sel_colors = [ {aIdx:_highlight_colors[n%len(_highlight_colors)] for aIdx in subset } for mol_set in highlightAtoms for n, subset in enumerate(mol_set)  ]
        else:
            sel_atoms = highlightAtoms
            sel_colors = [ {aIdx: _highlight_colors[0] for aIdx in subset} for subset in highlightAtoms ]

    svg = MolsToGridImage(mols_list, highlightAtomLists=sel_atoms, highlightBondLists=[], highlightAtomColors=sel_colors,
                                                                legends=legends, molsPerRow=mols_perrow, useSVG=True)

    if filename:
        ext = splitext(filename)[-1]
        filename = filename if ext != '' else filename + '.svg'
        f = open(filename, 'w')
        f.write(svg)
        f.close()

    if ipython:
        return SVG(svg)
    else:
        return None
