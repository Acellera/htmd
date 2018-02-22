# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

from htmd.molecule.util import sequenceID
import numpy as np
import logging
import string

logger = logging.getLogger(__name__)


class BuildError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        if isinstance(self.value, str):
            return repr(self.value)
        elif isinstance(self.value, list):
            return '\n'.join([v if isinstance(v, str) else repr(v) for v in self.value])


class MixedSegmentError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ResidueInsertionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class UnknownResidueError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MissingParameterError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MissingTorsionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MissingBondError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MissingAngleError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MissingAtomTypeError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


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

    def __str__(self):
        return 'Disulfide bridge between (segid, resid) ({}, {}) and ({}, {})'.format(self.segid1, self.resid1, self.segid2, self.resid2)

    def __repr__(self):
        return self.__str__()


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
    """ Automatically detects disulfide bonds in a molecule

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule for which to detect disulfide bonds
    thresh : float
        The threshold under which two sulfurs are considered as bonded

    Returns
    -------
    disubonds : np.ndarray
        A list of :class:`DisulfideBridge <htmd.builder.builder.DisulfideBridge>` objects
    """
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse.csgraph import connected_components
    import pandas as pd
    disubonds = []

    # Find all SG atoms belonging to resnames starting with CY
    idx = np.where([(rn[0:2] == 'CY') and (n == 'SG') for rn, n in zip(mol.resname, mol.name)])[0] # 'resname "CY.*" and name SG'
    if len(idx) == 0:
            return disubonds

    # Used only for displaying nice messages
    if mol.chain is None:
        chains = [''] * len(idx)
    else:
        chains = mol.chain[idx]
    segids = mol.segid[idx]
    resids = mol.resid[idx]
    serials = mol.serial[idx]
    resnames = mol.resname[idx]
    if np.any([len(s) == 0 for s in segids]):
        raise RuntimeError('Cannot detect disulfide bonds without segment names defined.')

    # Check that there is only one SG atom in the same resid/segid combo
    df = pd.DataFrame({'segids': segids, 'resids': resids, 'indexes': idx})
    groups = df.groupby(['resids','segids']).groups
    for k in groups:
        if len(df['indexes'][groups[k]].tolist()) != 1:
            raise RuntimeError('Multiple SG atoms detected in segment {} resid {}. Can\'t guess disulfide bridges.'.format(k[1], k[0]))

    sd = squareform(pdist(mol.coords[idx, :, mol.frame]))
    sd[np.diag_indices(sd.shape[0])] = thresh+1  # Set the diagonal over threshold
    close = sd < thresh
    rows, cols = np.where(close)

    numbonds = np.sum(close, axis=0)
    if np.any(numbonds > 1):
        pairs = [(s, r) for r, s in zip(resids[np.where(numbonds > 1)[0]], segids[np.where(numbonds > 1)[0]])]
        raise RuntimeError('SG atoms with (segid, resid) pairs {} have multiple possible bonds. Cannot guess disulfide bonds. Please specify them manually.'.format(pairs))

    uniquerowcols = list(set([tuple(sorted((r, c))) for r, c in zip(rows, cols)]))
    for rc in uniquerowcols:
        disubonds.append(DisulfideBridge(segids[rc[0]], resids[rc[0]], segids[rc[1]], resids[rc[1]]))
        msg = 'Bond between A: [serial {0} resid {1} resname {2} chain {3} segid {4}]\n' \
              '             B: [serial {5} resid {6} resname {7} chain {8} segid {9}]\n'.format(
            serials[rc[0]], resids[rc[0]], resnames[rc[0]], chains[rc[0]], segids[rc[0]],
            serials[rc[1]], resids[rc[1]], resnames[rc[1]], chains[rc[1]], segids[rc[1]]
        )
        print(msg)
        #logger.info(msg)

    if len(disubonds) == 1:
        logger.info('One disulfide bond was added')
    else:
        logger.info('{} disulfide bonds were added'.format(len(disubonds)))
    return sorted(disubonds, key=lambda x: x.resid1)


def autoSegment(mol, sel='all', basename='P', spatial=True, spatialgap=4.0, field="segid"):
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
        The size of a spatial gap which validates a discontinuity (A)
    field : str
        Field to fix. Can be "segid" (default), "chain", or "both"

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
    rid = mol.get('resid', sel)
    residiff = np.diff(rid)
    gappos = np.where((residiff != 1) & (residiff != 0))[0]  # Points to the index before the gap!

    # Letters to be used for chains, if free: 0123456789abcd...ABCD..., minus chain symbols already used
    used_chains = set(mol.chain)
    chain_alphabet = list(string.digits + string.ascii_letters)
    available_chains = [x for x in chain_alphabet if x not in used_chains]

    idxstartseg = [idx[0]] + idx[gappos + 1].tolist()
    idxendseg = idx[gappos].tolist() + [idx[-1]]

    mol.set('segid', basename, sel)

    if len(gappos) == 0:
        mol.set('segid', basename+'0', sel)
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
        # Fixup segid
        if field in ['segid', 'both']:
            newsegid = basename + str(i)
            if np.any(mol.segid == newsegid):
                raise RuntimeError('Segid {} already exists in the molecule. Please choose different prefix.'.format(newsegid))
            logger.info('Created segment {} between resid {} and {}.'.format(newsegid, mol.resid[s], mol.resid[e]))
            mol.segid[s:e+1] = newsegid
        # Fixup chain
        if field in ['chain', 'both']:
            newchainid = available_chains[i]
            logger.info('Set chain {} between resid {} and {}.'.format(newchainid, mol.resid[s], mol.resid[e]))
            mol.chain[s:e+1] = newchainid

        i += 1

    return mol


def autoSegment2(mol, sel='(protein or resname ACE NME)', basename='P', fields=('segid',), residgaps=False, residgaptol=1, chaingaps=True):
    """ Detects bonded segments in a selection and assigns incrementing segid to each segment

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object
    sel : str
        Atom selection on which to check for gaps.
    basename : str
        The basename for segment ids. For example if given 'P' it will name the segments 'P1', 'P2', ...
    fields : tuple of strings
        Field to fix. Can be "segid" (default) or any other Molecule field or combinations thereof.
    residgaps : bool
        Set to True to consider gaps in resids as structural gaps. Set to False to ignore resids
    residgaptol : int
        Above what resid difference is considered a gap. I.e. with residgaptol 1, 235-233 = 2 > 1 hence is a gap. We set
        default to 2 because in many PDBs single residues are missing in the proteins without any gaps.
    chaingaps : bool
        Set to True to consider changes in chains as structural gaps. Set to False to ignore chains

    Returns
    -------
    newmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A new Molecule object with modified segids

    Example
    -------
    >>> newmol = autoSegment2(mol)
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    if isinstance(fields, str):
        fields = (fields,)

    sel += ' and backbone or (resname NME ACE and name N C O CH3)'  # Looking for bonds only over the backbone of the protein
    idx = mol.atomselect(sel, indexes=True)  # Keep the original atom indexes to map from submol to mol
    submol = mol.copy()  # We filter out everything not on the backbone to calculate only those bonds
    submol.filter(sel, _logger=False)
    bonds = submol._getBonds()  # Calculate both file and guessed bonds

    if residgaps:
        # Remove bonds between residues without continuous resids
        bondresiddiff = np.abs(submol.resid[bonds[:, 0]] - submol.resid[bonds[:, 1]])
        bonds = bonds[bondresiddiff <= residgaptol, :]
    else:
        # Warning about bonds bonding non-continuous resids
        bondresiddiff = np.abs(submol.resid[bonds[:, 0]] - submol.resid[bonds[:, 1]])
        if np.any(bondresiddiff > 1):
            for i in np.where(bondresiddiff > residgaptol)[0]:
                logger.warning('Bonds found between resid gaps: resid {} and {}'.format(submol.resid[bonds[i, 0]],
                                                                                        submol.resid[bonds[i, 1]]))
    if chaingaps:
        # Remove bonds between residues without same chain
        bondsamechain = submol.chain[bonds[:, 0]] == submol.chain[bonds[:, 1]]
        bonds = bonds[bondsamechain, :]
    else:
        # Warning about bonds bonding different chains
        bondsamechain = submol.chain[bonds[:, 0]] == submol.chain[bonds[:, 1]]
        if np.any(bondsamechain == False):
            for i in np.where(bondsamechain == False)[0]:
                logger.warning('Bonds found between chain gaps: resid {}/{} and {}/{}'.format(submol.resid[bonds[i, 0]],
                                                                                              submol.chain[bonds[i, 0]],
                                                                                              submol.resid[bonds[i, 1]],
                                                                                              submol.chain[bonds[i, 1]]
                                                                                              ))

    # Calculate connected components using the bonds
    sparsemat = csr_matrix((np.ones(bonds.shape[0] * 2),  # Values
                            (np.hstack((bonds[:, 0], bonds[:, 1])),  # Rows
                             np.hstack((bonds[:, 1], bonds[:, 0])))), shape=[submol.numAtoms, submol.numAtoms])  # Columns
    numcomp, compidx = connected_components(sparsemat, directed=False)

    # Letters to be used for chains, if free: 0123456789abcd...ABCD..., minus chain symbols already used
    used_chains = set(mol.chain)
    chain_alphabet = list(string.digits + string.ascii_letters)
    available_chains = [x for x in chain_alphabet if x not in used_chains]

    mol = mol.copy()
    prevsegres = None
    for i in range(numcomp):  # For each connected component / segment
        segid = basename + str(i)
        backboneSegIdx = idx[compidx == i]  # The backbone atoms of the segment
        segres = mol.atomselect('same residue as index {}'.format(' '.join(map(str, backboneSegIdx)))) # Get whole residues

        # Warning about separating segments with continuous resids
        if i > 0 and (np.min(mol.resid[segres]) - np.max(mol.resid[prevsegres])) == 1:
            logger.warning('Separated segments {} and {}, despite continuous resids, due to lack of bonding.'.format(
                            basename + str(i-1), segid))

        # Add the new segment ID to all fields the user specified
        for f in fields:
            if f != 'chain':
                if np.any(mol.__dict__[f] == segid):
                    raise RuntimeError('Segid {} already exists in the molecule. Please choose different prefix.'.format(segid))
                mol.__dict__[f][segres] = segid  # Assign the segid to the correct atoms
            else:
                mol.__dict__[f][segres] = available_chains[i % len(available_chains)]
        logger.info('Created segment {} between resid {} and {}.'.format(segid, np.min(mol.resid[segres]),
                                                                         np.max(mol.resid[segres])))
        prevsegres = segres  # Store old segment atom indexes for the warning about continuous resids

    return mol


def _checkMixedSegment(mol):
    prot = mol.atomselect('protein')
    acenme = (mol.resname == 'ACE') | (mol.resname == 'NME')
    sel1 = prot | acenme  # 'protein or resname ACE NME'
    sel2 = ~prot & ~acenme  # 'not protein and not resname ACE NME'
    segsProt = np.unique(mol.segid[sel1])
    segsNonProt = np.unique(mol.segid[sel2])
    intersection = np.intersect1d(segsProt, segsNonProt)
    if len(intersection) != 0:
        logger.warning('Segments {} contain both protein and non-protein atoms. '
                       'Please assign separate segments to them or the build procedure might fail.'.format(intersection))


def _checkResidueInsertions(mol):
    ins = np.unique([x for x in mol.insertion if x != ''])
    if len(ins) > 0:
        raise ResidueInsertionError('Residue insertions "{}" were detected in input molecule and are not supported for'
                                    ' building. Please use mol.renumberResidues() on your protein '
                                    'Molecule.'.format(' '.join(ins)))
    pass


def removeLipidsInProtein(prot, memb, lipidsel='lipids'):
    """ Calculates the convex hull of the protein. If a lipid lies inside the hull it gets removed.

    This does not work well for lipids crossing out of the hull. If even one atom of the lipid is outside it will
    change the hull and will not get removed. I assume it will get removed by the clashes with the protein though.
    """
    return removeAtomsInHull(prot, memb, 'name CA', lipidsel)


def removeAtomsInHull(mol1, mol2, hullsel, removesel):
    """ Calculates the convex hull of an atom selection in mol1 and removes atoms within that hull in mol2.

    Parameters
    ----------
    mol1 : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        Molecule for which to calculate the convex hull
    mol2 : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        Molecule which contains the atoms which we check if they are within the hull
    hullsel : str
        Atomselection for atoms in mol1 from which to calculate the convex hull.
    removesel : str
        Atomselection for atoms in mol2 from which to remove the ones which are within the hull

    Returns
    -------
    newmol2 : Molecule
        mol2 but without any atoms located within the convex hull
    numrem : int
        Number of fragments removed
    """
    # TODO: Look into Morphological Snakes
    from scipy.spatial import ConvexHull

    mol2 = mol2.copy()
    # Convex hull of the protein
    hullcoords = mol1.get('coords', hullsel)
    hull = ConvexHull(hullcoords)

    sequence = sequenceID((mol2.resid, mol2.segid))
    uqres = np.unique(sequence)

    toremove = np.zeros(len(sequence), dtype=bool)
    numlipsrem = 0
    for res in uqres:  # For each fragment check if it's atoms lie within the convex hull
        atoms = np.where(sequence == res)[0]
        newhull = ConvexHull(np.vstack((hullcoords, mol2.get('coords', sel=atoms))))

        # If the hull didn't change by adding the fragment, it lies within convex hull. Remove it.
        if list(hull.vertices) == list(newhull.vertices):
            toremove[atoms] = True
            numlipsrem += 1

    rematoms = mol2.atomselect(removesel)
    mol2.remove(toremove & rematoms)
    return mol2, numlipsrem


def removeHET(prot):
    prot = prot.copy()
    hetatoms = np.unique(prot.resname[prot.record == 'HETATM'])
    for het in hetatoms:
        logger.info('Found resname ''{}'' in structure. Removed assuming it is a bound ligand.'.format(het))
        prot.remove('resname {}'.format(het))
    return prot


def tileMembrane(memb, xmin, ymin, xmax, ymax, buffer=1.5):
    """ Tile a membrane in the X and Y dimensions to reach a specific size.

    Parameters
    ----------
    memb : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The membrane to be tiled
    xmin : float
        Minimum x coordinate
    ymin : float
        Minimum y coordinate
    xmax : float
        Maximum x coordinate
    ymax : float
        Maximum y coordinate
    buffer : float
        Buffer distance between tiles

    Returns
    -------
    megamemb :
        A big membrane Molecule
    """
    from tqdm import tqdm
    memb = memb.copy()
    memb.resid = sequenceID((memb.resid, memb.insertion, memb.chain, memb.segid))

    minmemb = np.min(memb.get('coords', 'water'), axis=0).flatten()

    size = np.max(memb.get('coords', 'water'), axis=0) - np.min(memb.get('coords', 'water'), axis=0)
    size = size.flatten()
    xreps = int(np.ceil((xmax - xmin) / size[0]))
    yreps = int(np.ceil((ymax - ymin) / size[1]))

    logger.info('Replicating Membrane {}x{}'.format(xreps, yreps))

    from htmd.molecule.molecule import Molecule
    megamemb = Molecule()
    bar = tqdm(total=xreps * yreps, desc='Replicating Membrane')
    k = 0
    for x in range(xreps):
        for y in range(yreps):
            tmpmemb = memb.copy()
            xpos = xmin + x * (size[0] + buffer)
            ypos = ymin + y * (size[1] + buffer)

            tmpmemb.moveBy([-float(minmemb[0]) + xpos, -float(minmemb[1]) + ypos, 0])
            tmpmemb.remove('same resid as (x > {} or y > {})'.format(xmax, ymax), _logger=False)
            if tmpmemb.numAtoms == 0:
                continue

            tmpmemb.set('segid', 'M{}'.format(k), sel='not water')
            tmpmemb.set('segid', 'MW{}'.format(k), sel='water')

            megamemb.append(tmpmemb)
            k += 1
            bar.update(1)
    bar.close()

    # Membranes don't tile perfectly. Need to remove waters that clash with lipids of other tiles
    # Some clashes will still occur between periodic images however
    megamemb.remove('same resid as water and within 1.5 of not water', _logger=False)
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
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path

    p = Molecule(path.join(home(), 'data', 'building-protein-membrane', '4dkl.pdb'))
    p.filter('(chain B and protein) or water')
    p = autoSegment(p, 'protein', 'P')
    m = Molecule(path.join(home(), 'data', 'building-protein-membrane', 'membrane.pdb'))
    a = embed(p, m)
    print(np.unique(m.get('segid')))

    mol = Molecule(path.join(home(), 'data', 'building-protein-membrane', '1ITG_clean.pdb'))
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '1ITG.pdb'))
    mol = autoSegment(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

    mol = Molecule(path.join(home(), 'data', 'building-protein-membrane', '3PTB_clean.pdb'))
    ref = Molecule(path.join(home(), 'data', 'building-protein-membrane', '3PTB.pdb'))
    mol = autoSegment(mol, sel='protein')
    assert np.all(mol.segid == ref.segid)

