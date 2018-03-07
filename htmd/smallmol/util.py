from math import cos, sqrt, sin
import numpy as np


_highlight_colors = [(1.00,0.50,0.00), (0.00,0.50,1.00), (0.00,1.00,0.50),
                     (1.00,0.00,0.50), (0.50,0.00,1.00), (0.50,1.00,0.00),
                     (1.00,0.00,0.25), (0.00,0.25,1.00), (0.25,1.00,0.00)]

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

def InputToOutput(input_file, input_format, output_format):
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

def _depictMol(mol, sketch=False, filename=None, ipython=False, optimize=False, optimizemode='std', removeHs=True, atomlabels=None, highlightAtoms=None):
    from rdkit.Chem import RemoveHs
    from rdkit.Chem.AllChem import  Compute2DCoords, EmbedMolecule, MMFFOptimizeMolecule, ETKDG
    from rdkit.Chem.Draw import rdMolDraw2D
    from IPython.display import SVG
    from copy import deepcopy
    from os.path import splitext

    if sketch and optimize:
        raise ValueError('Impossible to use optmization in  2D sketch representation')

    if optimizemode not in ['std', 'mmff']:
        raise ValueError('Optimization mode {} not understood. Can be "std" or  "ff"'.format(optimizemode))

    if highlightAtoms is not None and not isinstance(highlightAtoms, list):
        raise ValueError('highlightAtoms should be a list of atom idx or a list of atom idx list ')

    _mol = deepcopy(mol)

    # 2D representation. Set z coords to 0
    if sketch:
        Compute2DCoords(_mol)
    # Clean representation without hydrogens
    if removeHs:
        _mol = RemoveHs(_mol)

    # init the drawer object
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
    # get the drawer options
    opts = drawer.drawOptions()

    # add atomlabels
    if atomlabels is not None:
        for n, a in enumerate(_mol.GetAtoms()):
            _atomIdx = a.GetIdx()
            opts.atomLabels[n] = atomlabels[_atomIdx]

    # activate 3D coords optimization
    if optimize:
        if optimizemode == 'std':
            EmbedMolecule(_mol, ETKDG())
        elif optimizemode == 'mmff':
            MMFFOptimizeMolecule(_mol)

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

    drawer.DrawMolecule(_mol, highlightAtoms=sel_atoms, highlightBonds=[], highlightAtomColors=sel_colors)

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

def depictMultipleMols(mols_list, sketch=False, filename=None, ipython=False, optimize=False, optimizemode='std',
                       removeHs=True,  legends=None, highlightAtoms=None, mols_perrow=3):

    from rdkit.Chem import RemoveHs
    from rdkit.Chem.Draw import MolsToGridImage
    from rdkit.Chem.AllChem import Compute2DCoords, EmbedMolecule, MMFFOptimizeMolecule, ETKDG
    from copy import deepcopy
    from IPython.display import SVG
    from os.path import splitext

    if sketch and optimize:
        raise ValueError('Impossible to use optmization in  2D sketch representation')

    _mols = [ deepcopy(m) for m in mols_list ]

    if sketch:
        for _m in _mols: Compute2DCoords(_m)

    if removeHs:
        _mols = [ RemoveHs(_m) for _m in _mols ]

    # activate 3D coords optimization
    if optimize:
        if optimizemode == 'std':
            for _m in _mols: EmbedMolecule(_m)
        elif optimizemode == 'mmff':
            for _m in _mols: MMFFOptimizeMolecule(_m, ETKDG())


    sel_atoms = []
    sel_colors = []
    if highlightAtoms is not None:
        if isinstance(highlightAtoms[0][0], list):
            sel_atoms = [ [a for a in subset] for mol_set in highlightAtoms for subset in mol_set ]
            sel_colors = [ {aIdx:_highlight_colors[n%len(_highlight_colors)] for aIdx in subset } for mol_set in highlightAtoms for n, subset in enumerate(mol_set)  ]
        else:
            sel_atoms = highlightAtoms
            sel_colors = [ {aIdx: _highlight_colors[0] for aIdx in subset} for subset in highlightAtoms ]


    svg = MolsToGridImage(_mols, highlightAtomLists=sel_atoms, highlightBondLists=[], highlightAtomColors=sel_colors,
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