# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

from htmd.molecule.util import sequenceID
import numpy as np
import logging
logger = logging.getLogger(__name__)


class DisulfideBridge(object):
    def __init__(self, segid1, resid1, segid2, resid2):
        if not isinstance(segid1, str) or not isinstance(segid2, str):
            raise NameError('segment1 and segment2 options need to be strings')
        if not isinstance(resid1, (int, np.int32, np.int64)) or not isinstance(resid2, (int, np.int32, np.int64)):
            raise NameError('residue1 and residue2 options need to be integers')
        self.segid1 = segid1
        self.resid1 = resid1
        self.segid2 = segid2
        self.resid2 = resid2


def embed(mol1, mol2, gap=1.3):
    '''Embeds one molecule into another removing overlaps.

    Will remove residues of mol1 which have collisions with atoms of mol2.

    Parameters
    ----------
    mol1 : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The first Molecule object
    mol2 : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The second Molecule object
    gap : float
        Minimum space in A between atoms of the two molecules

    Return
    ------
    newmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The resulting Molecule object

    Example
    -------
    >>> all = embed(memb, prot)
    '''
    mol1 = mol1.copy()
    mol2 = mol2.copy()
    #Set different occupancy to separate atoms of mol1 and mol2
    occ1 = mol1.get('occupancy')
    occ2 = mol2.get('occupancy')
    mol1.set('occupancy', 1)
    mol2.set('occupancy', 2)

    mol2.append(mol1)
    s1 = mol2.atomselect('occupancy 1')
    s2 = mol2.atomselect('occupancy 2')
    # Give unique "residue" beta number to all resids
    beta = mol2.get('beta')
    mol2.set('beta', sequenceID(mol2.resid))
    # Calculate overlapping atoms
    overlaps = mol2.atomselect('(occupancy 2) and same beta as exwithin '+ str(gap) + ' of (occupancy 1)')
    # Restore original beta and occupancy
    mol2.set('beta', beta)
    mol2.set('occupancy', occ1, s1)
    mol2.set('occupancy', occ2, s2)

    # Remove the overlaps
    mol2.remove(overlaps, _logger=False)
    return mol2


def detectDisulfideBonds(mol, thresh=3):
    """ Detect disulfide bonds in a molecule

    Automatically detects disulfide bonds and stores them in the Molecule.disubonds field. Returns the pairs of
    atom indexes.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule for which to detect disulfide bonds
    thresh : float
        The threshold under which two sulfurs are considered as bonded

    Returns
    -------
    pairs : np.ndarray
        The pairs of atom indexes forming a disulfide bond
    """
    disubonds = []

    idx = np.where(mol.atomselect('resname "CY.*" and name SG'))[0]
    if len(idx) == 0:
            return disubonds

    # Used only for displaying nice messages
    if mol.chain is None:
        chains = [''] * np.size(mol.serial, 0)
    else:
        chains = mol.chain
    if mol.segid is None:
        raise NameError('Cannot detect disulfide bonds without segment names defined.')
        #segids = [''] * np.size(mol.serial, 0)
    else:
        segids = mol.segid

    for sg in idx:
        resid = mol.resid[sg]
        sel = '(not resid "{0}") and resname "CY.*" and name SG and index > {1} and exwithin {2} of index {3}'.format(resid, sg, thresh, sg)
        idx2 = np.where(mol.atomselect(sel))[0]

        if len(idx2) == 0:
            continue

        for sg2 in idx2:
            disubonds.append(DisulfideBridge(mol.segid[sg], mol.resid[sg], mol.segid[sg2], mol.resid[sg2]))

        msg = 'Bond between A: [serial {0} resid {1} resname {2} chain {3} segid {4}]\n' \
              '             B: [serial {5} resid {6} resname {7} chain {8} segid {9}]\n'.format(
            mol.serial[sg], mol.resid[sg], mol.resname[sg], chains[sg], segids[sg],
            mol.serial[sg2], mol.resid[sg2], mol.resname[sg2], chains[sg2], segids[sg2]
        )
        # logger.info(msg)
        print(msg)
    if len(disubonds) == 1:
        logger.info('One disulfide bond was added')
    else:
        logger.info('{} disulfide bonds were added'.format(len(disubonds)))
    return disubonds


# TODO: Remove in upcoming versions
def segmentgaps(mol, sel='all', basename='P', spatial=True, spatialgap=4):
    logger.warning('segmentgaps will be deprecated in next versions of HTMD. '
                   'It is replaced by autoSegment to reflect the usage rather than the underlying method. '
                   'Please change all usages.')
    return autoSegment(mol, sel=sel, basename=basename, spatial=spatial, spatialgap=spatialgap)


def autoSegment(mol, sel='all', basename='P', spatial=True, spatialgap=4):
    """ Detects resid gaps in a selection and assigns incrementing segid to each fragment

    !!!WARNING!!! If you want to use atom selections like 'protein' or 'fragment',
    use this function on a Molecule containing only protein atoms, otherwise the protein selection can fail.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object
    sel : str
        Atom selection on which to check for gaps.
    basename : str
        The basename for segment ids. For example if given 'P' it will name the segments 'P1', 'P2', ...
    spatial : bool
        Only considers a discontinuity in resid as a gap of the CA atoms have distance more than `spatialgap` Angstrom
    spatialgap : float
        The size of a spatial gap which validates a discontinuity

    Returns
    -------
    newmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A new Molecule object with modified segids

    Example
    -------
    >>> newmol = autoSegment(mol,'chain B','P')
    """
    mol = mol.copy()

    idx = mol.atomselect(sel, indexes=True)
    rid = mol.get('resid', sel)  # TODO:maybe easier without sel, rid = mol.get('resid') and then no need of idx
    idxdiff = np.diff(idx)
    residiff = np.diff(rid)
    gappos = np.where((residiff != 1) & (residiff != 0))[0]  # Points to the index before the gap!

    idxstartseg = [idx[0]] + idx[gappos + 1].tolist()
    idxendseg = idx[gappos].tolist() + [idx[-1]]

    from IPython.core.debugger import Tracer
    Tracer()()

    mol.set('segid', basename, sel)

    if len(gappos) == 0:
        return mol

    if spatial:
        residbackup = mol.get('resid')
        mol.set('resid', sequenceID(mol.resid))  # Assigning unique resids to be able to do the distance selection

        todelete = []
        i = 0
        for s, e in zip(idxstartseg[1:], idxendseg[:-1]):
            coords = mol.get('coords', sel='resid "{}" "{}" and name CA'.format(mol.resid[e], mol.resid[s]))
            if np.shape(coords) == (2, 3):
                dist = np.sqrt(np.sum((coords[0, :] - coords[1, :]) ** 2))
                if dist < spatialgap:
                    todelete.append(i)
            i += 1
        # Join the non-real gaps into segments
        idxstartseg = np.delete(idxstartseg, np.array(todelete)+1)
        idxendseg = np.delete(idxendseg, todelete)

        mol.set('resid', residbackup)  # Restoring the original resids

    i = 0
    for s, e in zip(idxstartseg, idxendseg):
        newsegid = basename + str(i)
        if np.any(mol.segid == newsegid):
            raise RuntimeError('Segid {} already exists in the molecule. Please choose different prefix.'.format(newsegid))
        logger.info('Created segment {} between resid {} and {}.'.format(newsegid, mol.resid[s], mol.resid[e]))
        mol.segid[s:e+1] = newsegid
        i += 1

    return mol


def removeLipidsInProtein(prot, memb):
    """ Calculates the convex hull of the protein. If a lipid lies inside the hull it gets removed.

    This does not work well for lipids crossing out of the hull. If even one atom of the lipid is outside it will
    change the hull and will not get removed. I assume it will get removed by the clashes with the protein though.
    """
    # TODO: Do the same with Morphological Snakes
    from scipy.spatial import ConvexHull
    memb = memb.copy()

    # Convex hull of the protein
    cacoords = prot.get('coords', 'name CA')
    hull = ConvexHull(cacoords)

    sequence = sequenceID((memb.resid, memb.segid))
    uqres = np.unique(sequence)

    toremove = np.zeros(len(sequence), dtype=bool)
    numlipsrem = 0
    for res in uqres:  # For each lipid check if it's atoms lie within the convex hull
        atoms = np.where(sequence == res)[0]
        newhull = ConvexHull(np.append(cacoords, np.squeeze(memb.coords[atoms, :, :]), axis=0))

        # If the hull didn't change by adding the lipid, it lies within convex hull. Remove it.
        if list(hull.vertices) == list(newhull.vertices):
            toremove[atoms] = True
            numlipsrem += 1

    lipids = memb.atomselect('lipids')  # Only remove lipids, waters are ok
    memb.remove(toremove & lipids)
    return memb, numlipsrem


def removeHET(prot):
    prot = prot.copy()
    hetatoms = np.unique(prot.resname[prot.record == 'HETATM'])
    for het in hetatoms:
        logger.info('Found resname ''{}'' in structure. Removed assuming it is a bound ligand.'.format(het))
        prot.remove('resname {}'.format(het))
    return prot


def tileMembrane(memb, xmin, ymin, xmax, ymax):
    """ Tile the membrane in the X and Y dimensions to reach a specific size.
    Returns
    -------
    megamemb :
        A big membrane Molecule
    """
    from htmd.progress.progress import ProgressBar
    memb = memb.copy()
    memb.resid = sequenceID(memb.resid)

    minmemb = np.min(memb.get('coords', 'water'), axis=0).flatten()

    size = np.max(memb.get('coords', 'water'), axis=0) - np.min(memb.get('coords', 'water'), axis=0)
    size = size.flatten()
    xreps = int(np.ceil((xmax - xmin) / size[0]))
    yreps = int(np.ceil((ymax - ymin) / size[1]))

    logger.info('Replicating Membrane {}x{}'.format(xreps, yreps))

    from htmd.molecule.molecule import Molecule
    megamemb = Molecule()
    bar = ProgressBar(xreps * yreps, description='Replicating Membrane')
    k = 0
    for x in range(xreps):
        for y in range(yreps):
            tmpmemb = memb.copy()
            xpos = xmin + x * size[0]
            ypos = ymin + y * size[1]

            tmpmemb.moveBy([-float(minmemb[0]) + xpos, -float(minmemb[1]) + ypos, 0])
            sel = 'same resid as (x > {} or y > {})'.format(xmax, ymax)
            tmpmemb.remove(sel, _logger=False)
            tmpmemb.set('segid', 'M{}'.format(k))

            megamemb.append(tmpmemb)
            k += 1
            bar.progress()
    bar.stop()
    return megamemb


def minimalRotation(prot):
    """ Find the rotation around Z that minimizes the X and Y dimensions of the protein to best fit in a box.

    Essentially PCA in 2D
    """
    from numpy.linalg import eig
    from numpy import cov

    xycoords = prot.coords[:, 0:2]

    c = cov(np.transpose(np.squeeze(xycoords)))
    values, vectors = eig(c)
    idx = np.argsort(values)

    xa = vectors[0, idx[-1]]
    ya = vectors[1, idx[-1]]

    def cart2pol(x, y):
        # Cartesian to polar coordinates. Rho is the rotation angle
        rho = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)
        return rho, phi

    angle, _ = cart2pol(xa, ya)
    return angle + np.radians(45)


if __name__ == "__main__":
    from htmd import *
    from os import path

    p = Molecule(path.join(home(), 'data', 'building-protein-membrane', '4dkl.pdb'))
    p.filter('(chain B and protein) or water')
    p = autoSegment(p, 'protein', 'P')
    m = Molecule(path.join(home(), 'data', 'building-protein-membrane', 'membrane.pdb'))
    a = embed(p, m)
    print(np.unique(m.get('segid')))

    mol = Molecule('1ITG')
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '1ITG.pdb'))
    mol = autoSegment(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

    mol = Molecule('3PTB')
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '3PTB.pdb'))
    mol = autoSegment(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

