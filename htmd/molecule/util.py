# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial.distance import cdist


def molRMSD(mol, refmol, rmsdsel1, rmsdsel2):
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
    >>> sequenceID((mol.resid, mol.segid))
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
    >>> bbox = boundingBox(mol, sel='chain A')
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
