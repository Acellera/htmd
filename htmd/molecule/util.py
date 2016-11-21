# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial.distance import cdist
import ctypes as ct
import os
import htmd.home
import logging
logger = logging.getLogger(__name__)

libdir = htmd.home(libDir=True)
tmalignlib = ct.cdll.LoadLibrary(os.path.join(libdir, "tmalign.so"))


def molTMscore(mol, ref, selCAmol, selCAref):
    """ Calculates the TMscore between two Molecules

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A Molecule containing a single or multiple frames
    ref : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A reference Molecule containing a single frame. Will automatically keep only ref.frame.
    selCAmol : numpy.ndarray
        An atomselection array of booleans or indexes of the CA atoms for mol
    selCAref : numpy.ndarray
        An atomselection array of booleans or indexes of the CA atoms for ref

    Returns
    -------
    tmscoreRef : numpy.ndarray
        TMscore normalized by length of ref
    rmsd : numpy.ndarray
        RMSD for all frames
    """
    from htmd.builder.builder import sequenceID
    from htmd.molecule.molecule import _residueNameTable

    def calculateVariables(currmol):
        res = sequenceID((currmol.resid, currmol.insertion, currmol.segid, currmol.chain))
        caidx = currmol.name == 'CA'
        res = np.unique(res)
        reslen = len(res)
        res = res.astype(np.int32).ctypes.data_as(ct.POINTER(ct.c_int))

        # Calculate the protein sequence
        seq = ''.join([_residueNameTable[x] for x in currmol.resname[caidx]])
        seq = ct.c_char_p(seq.encode('utf-8'))

        # Keep only CA coordinates
        coords = currmol.coords[caidx, :, :].copy()
        return reslen, res, seq, coords

    mol = mol.copy()
    ref = ref.copy()
    mol.filter(selCAmol, _logger=False)
    ref.filter(selCAref, _logger=False)
    ref.dropFrames(keep=ref.frame)

    reslenMOL, residMOL, seqMOL, coordsMOL = calculateVariables(mol)
    reslenREF, residREF, seqREF, coordsREF = calculateVariables(ref)

    # void tmalign(int xlen, int ylen, int* xresno, int* yresno, char* seqx, char* seqy, float* xcoor, float* ycoor, int nframes, double &TM1, double &TM2, double &rmsd)
    resTM1 = (ct.c_double * mol.numFrames)()
    resTM2 = (ct.c_double * mol.numFrames)()
    resRMSD = (ct.c_double * mol.numFrames)()
    res = tmalignlib.tmalign(ct.c_int(reslenREF),
                       ct.c_int(reslenMOL),
                       residREF,
                       residMOL,
                       seqREF,
                       seqMOL,
                       coordsREF.ctypes.data_as(ct.POINTER(ct.c_float)),
                       coordsMOL.ctypes.data_as(ct.POINTER(ct.c_float)),
                       ct.c_int(mol.numFrames),
                       ct.byref(resTM1),
                       ct.byref(resTM2),
                       ct.byref(resRMSD))
    return np.ctypeslib.as_array(resTM1), np.ctypeslib.as_array(resRMSD)



def molRMSD(mol, refmol, rmsdsel1, rmsdsel2):
    """ Calculates the RMSD between two Molecules

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
    refmol
    rmsdsel1
    rmsdsel2

    Returns
    -------
    rmsd : float
        The RMSD between the two structures
    """
    dist = mol.coords[rmsdsel1, :, :] - refmol.coords[rmsdsel2, :, :]
    rmsd = np.sqrt(np.mean(np.sum(dist * dist, axis=1), axis=0))
    return np.squeeze(rmsd)


def sequenceID(field, prepend=None):
    """ Array of integers which increments at value change of another array

    Parameters
    ----------
    field : np.ndarray or tuple
        An array of values. Once a change in value happens, a new ID will be created in `seq`.
        If a tuple of ndarrays is passed, a change in any of them will cause an increase in `seq`.
    prepend : str
        A string to prepend to the incremental sequence

    Returns
    -------
    seq : np.ndarray
        An array of equal size to `field` containing integers which increment every time there is a change in `field`

    Examples
    --------
    >>> # A change in resid or segid will cause an increase in the sequence
    >>> sequenceID((mol.resid, mol.insertion, mol.segid))
    array([  1,   1,   1, ..., 285, 286, 287])
    """
    if isinstance(field, tuple):
        fieldlen = len(field[0])
    else:
        fieldlen = len(field)

    if prepend is None:
        seq = np.zeros(fieldlen, dtype=int)
    else:
        seq = np.empty(fieldlen, dtype=object)

    c = int(1)
    if prepend is None:
        seq[0] = c
    else:
        seq[0] = prepend + str(c)

    for i in range(1, fieldlen):
        if isinstance(field, tuple):  # Support tuples of multiple fields. Change in any of them will cause an increment
            for t in field:
                if t[i-1] != t[i]:
                    c += 1  # new sequence id
                    break
        elif field[i-1] != field[i]:
            c += 1  # new sequence id
        if prepend is None:
            seq[i] = c
        else:
            seq[i] = prepend + str(c)
    return seq


def _missingChain(mol):
    if mol.chain is None or np.size(mol.chain) == 0:
        raise NameError('Segid fields have to be set for all atoms in the Molecule object before building.')
    empty = [True if len(c) == 0 else False for c in mol.chain]
    if np.any(empty):
        idx = np.where(empty)[0]
        if len(idx) == 1:
            raise NameError('Atom ' + str(idx) + ' does not have a chain defined.')
        elif len(idx) <= 5:
            raise NameError('Atoms ' + str(idx) + ' do not have a chain defined.')
        else:
            raise NameError('Atoms [' + str(idx[0]) + ',' + str(idx[1]) + ',...,' + str(idx[-1]) + '] do not have chain defined.')


def _missingSegID(mol):
    if mol.segid is None or np.size(mol.segid) == 0:
        raise NameError('Segid fields have to be set for all atoms in the Molecule object before building.')
    empty = [True if len(s) == 0 else False for s in mol.segid]
    if np.any(empty):
        idx = np.where(empty)[0]
        if len(idx) == 1:
            raise NameError('Atom ' + str(idx) + ' does not have a segid defined.')
        elif len(idx) <= 5:
            raise NameError('Atoms ' + str(idx) + ' do not have a segid defined.')
        else:
            raise NameError('Atoms [' + str(idx[0]) + ',' + str(idx[1]) + ',...,' + str(idx[-1]) + '] do not have segid defined.')


def maxDistance(mol, sel='all', origin=[0, 0, 0]):
    """ Calculates the max distance of a set of atoms from an origin

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule containing the atoms
    sel : str
        Atomselection for atoms for which to calculate distances
    origin : list
        The origin x,y,z coordinates

    Returns
    -------
    maxd : float
        The maximum distance in Angstrom

    Example
    -------
    >>> y = maxDistance(mol, sel='protein', origin=[0, 0, 0])
    >>> print(round(y,2))
    48.39
    """
    coors = mol.get('coords', sel=sel)
    dists = cdist(np.atleast_2d(coors), np.atleast_2d(origin))
    return np.max(dists)


def boundingBox(mol, sel='all'):
    """ Calculates the bounding box of a selection of atoms.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule containing the atoms
    sel : str
        An Atomselection string of atoms.

    Returns
    -------
    bbox : np.ndarray
        The bounding box around the atoms selected in `sel`.

    Example
    -------
    >>> boundingBox(mol, sel='chain A')
    array([[-17.3390007 , -10.43700027,  -1.43900001],
           [ 25.40600014,  27.03800011,  46.46300125]], dtype=float32)

    """
    coords = mol.get('coords', sel=sel)
    maxc = np.squeeze(np.max(coords, axis=0))
    minc = np.squeeze(np.min(coords, axis=0))
    return np.vstack((minc, maxc))


def uniformRandomRotation():
    """ Return a uniformly distributed rotation matrix

    Returns
    -------
    M : np.ndarray
        A uniformly distributed rotation matrix
    """
    q, r = np.linalg.qr(np.random.normal(size=(3, 3)))
    return np.dot(q, np.diag(np.sign(np.diag(r))))


def writeVoxels(arr, filename, vecMin, vecMax, vecRes):
    """ Writes grid free energy to cube file

    Parameters
    ----------
    arr: np.ndarray
            array with volumetric data
    filename: str
            string with the name of the cubefile
    vecMin: np.ndarray
            3D vector denoting the minimal corner of the grid
    vecMax np.ndarray
            3D vector denoting the maximal corner of the grid
    vecRes: np.ndarray
            3D vector denoting the resolution of the grid in each dimension
    """

    outFile = open(filename, 'w')

    # conversion to gaussian units
    L = 1/0.52917725
    gauss_bin = vecRes*L
    #minCorner = 0.5*L*(vecMin - vecMax + vecRes)
    minCorner = L * (vecMin + 0.5 * vecRes)

    ngrid = np.array(np.floor((vecMax - vecMin) / vecRes), dtype=int)

    # write header
    outFile.write("CUBE FILE\n")
    outFile.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (1, minCorner[0], minCorner[1], minCorner[2]))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[0], gauss_bin[0], 0, 0))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[1], 0, gauss_bin[1], 0))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[2], 0, 0, gauss_bin[2]))
    outFile.write("%5d %12.6f %12.6f %12.6f %12.6f\n" % (1, 0, minCorner[0], minCorner[1], minCorner[2]))

    # main loop
    cont = 0

    for i in range(ngrid[0]):
        for j in range(ngrid[1]):
            for k in range(ngrid[2]):
                outFile.write("%13.5g" % arr[i][j][k])
                if np.mod(cont, 6) == 5:
                    outFile.write("\n")
                cont += 1

    outFile.close()


def drawCube(mi, ma, viewer=None):
    from htmd.vmdviewer import getCurrentViewer
    if viewer is None:
        viewer = getCurrentViewer()

    viewer.send('draw materials off')
    viewer.send('draw color red')
    c = ''
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], ma[0], mi[1], mi[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], ma[1], mi[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], mi[1], ma[2])

    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], ma[1], mi[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], mi[1], ma[2])

    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], ma[0], ma[1], mi[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], mi[0], ma[1], ma[2])

    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], ma[0], mi[1], ma[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], mi[0], ma[1], ma[2])

    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], ma[1], mi[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], mi[0], ma[1], ma[2])
    c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], mi[1], ma[2])
    viewer.send(c)

    """
    draw line "$minx $miny $minz" "$maxx $miny $minz"
    draw line "$minx $miny $minz" "$minx $maxy $minz"
    draw line "$minx $miny $minz" "$minx $miny $maxz"
    draw line "$maxx $miny $minz" "$maxx $maxy $minz"
    draw line "$maxx $miny $minz" "$maxx $miny $maxz"
    draw line "$minx $maxy $minz" "$maxx $maxy $minz"
    draw line "$minx $maxy $minz" "$minx $maxy $maxz"
    draw line "$minx $miny $maxz" "$maxx $miny $maxz"
    draw line "$minx $miny $maxz" "$minx $maxy $maxz"
    draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
    draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
    draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
    """


# A test method
if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import numpy as np
    from os import path
    import doctest
    #doctest.testmod(extraglobs={"mol" : Molecule("3PTB")})

    expectedTMscore = np.array([ 0.21418524,  0.2367377 ,  0.23433833,  0.21362964,  0.20935164,
        0.20279461,  0.27012895,  0.22675238,  0.21230793,  0.2372011 ])
    expectedRMSD = np.array([ 3.70322128,  3.43637027,  3.188193  ,  3.84455877,  3.53053882,
        3.46781854,  2.93777629,  2.97978692,  2.70792428,  2.63051318])

    mol = Molecule(path.join(home(), 'data', 'tmscore', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'tmscore', 'traj.xtc'))
    ref = Molecule(path.join(home(), 'data', 'tmscore', 'ntl9_2hbb.pdb'))
    tmscore, rmsd = molTMscore(mol, ref, mol.atomselect('protein'), ref.atomselect('protein'))

    assert np.allclose(tmscore, expectedTMscore)
    assert np.allclose(rmsd, expectedRMSD)
