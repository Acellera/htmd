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
        sel = '(not resid {0}) and resname "CY.*" and name SG and index > {1} and exwithin {2} of index {3}'.format(resid, sg, thresh, sg)
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


def segmentgaps(mol, sel='all', basename='P', spatial=True, spatialgap=4):
    """ Detects resid gaps in a selection and assigns incrementing segid to each fragment

    !!!WARNING!!! If you want to use selections like protein or fragment,
    use this function on a pure protein Molecule, otherwise the protein
    selection will fail.

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
    >>> segmentgaps(mol,'chain B','P')
    """
    mol = mol.copy()

    idx = mol.atomselect(sel, indexes=True)
    rid = mol.get('resid', sel)  # TODO:maybe easier without sel, rid = mol.get('resid') and then no need of idx
    residiff = np.diff(rid)
    gappos = idx[np.where((residiff != 1) & (residiff != 0))[0]]

    if len(gappos) == 0:
        return mol

    gaps = np.array([idx[0]-1] + idx[gappos].tolist() + [idx[-1]])
    mol.set('segid', '', sel)

    if spatial:
        todelete = []
        for i in range(1, len(gaps)-1):
            gapidx = gaps[i]
            coords = mol.get('coords', sel='resid {} {} and name CA'.format(mol.resid[gapidx], mol.resid[gapidx+1]))
            if np.shape(coords) == (2, 3):
                dist = np.sqrt(np.sum((coords[0, :] - coords[1, :]) ** 2))
                if dist < spatialgap:
                    todelete.append(i)
        gaps = np.delete(gaps, todelete)

    start = gaps[0]+1
    for i in range(len(gaps) - 1):
        stop = gaps[i+1]
        newsegid = basename + str(i)
        logger.info('Created segment {} between resid {} and {}.'.format(newsegid, mol.resid[start], mol.resid[stop]))
        mol.segid[start:stop+1] = newsegid
        start = stop + 1

    return mol


if __name__ == "__main__":
    from htmd import *
    from os import path

    p = Molecule(path.join(home(), 'data', 'building-protein-membrane', '4dkl.pdb'))
    p.filter('(chain B and protein) or water')
    p = segmentgaps(p, 'protein', 'P')
    m = Molecule(path.join(home(), 'data', 'building-protein-membrane', 'membrane.pdb'))
    a = embed(p, m)
    print(np.unique(m.get('segid')))

    mol = Molecule('1ITG')
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '1ITG.pdb'))
    mol = segmentgaps(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

    mol = Molecule('3PTB')
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '3PTB.pdb'))
    mol = segmentgaps(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

