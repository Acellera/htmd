# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

from moleculekit.util import sequenceID
import numpy as np
import logging
import string
from htmd.decorators import _Deprecated

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


@_Deprecated('1.12.0', '<Read builder documentation on argument `disulfide`>')
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

    Will remove residues of mol2 which have collisions with atoms of mol1.

    Parameters
    ----------
    mol1 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The first Molecule object
    mol2 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The second Molecule object
    gap : float
        Minimum space in A between atoms of the two molecules

    Return
    ------
    newmol : :class:`Molecule <moleculekit.molecule.Molecule>` object
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


def convertDisulfide(mol, disu):
    from moleculekit.molecule import UniqueResidueID
    newdisu = []
    for d in disu:
        if not isinstance(d[0], str) or not isinstance(d[1], str):
            raise RuntimeError('All disulfide selections should be strings')
        newdisu.append([UniqueResidueID.fromMolecule(mol, d[0]), UniqueResidueID.fromMolecule(mol, d[1])])
    return newdisu


def detectDisulfideBonds(mol, thresh=3):
    """ Automatically detects disulfide bonds in a molecule

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The molecule for which to detect disulfide bonds
    thresh : float
        The threshold under which two sulfurs are considered as bonded

    Returns
    -------
    disubonds : np.ndarray
        A list of :class:`DisulfideBridge <htmd.builder.builder.DisulfideBridge>` objects
    """
    from scipy.spatial.distance import pdist, squareform
    from moleculekit.molecule import UniqueResidueID
    disubonds = []

    # Find all SG atoms belonging to resnames starting with CY
    idx = np.where([(rn[0:2] == 'CY') and (n == 'SG') for rn, n in zip(mol.resname, mol.name)])[0] # 'resname "CY.*" and name SG'
    if len(idx) == 0:
        return disubonds

    if np.any([len(s) == 0 for s in mol.segid[idx]]):
        raise RuntimeError('Cannot detect disulfide bonds without segment names defined.')

    residues = [UniqueResidueID.fromMolecule(mol, idx=i) for i in idx]
    for r1 in range(len(residues)):
        for r2 in range(r1+1, len(residues)):
            if residues[r1] == residues[r2]:
                raise RuntimeError('Multiple SG atoms detected in the same residue {}. '
                                   'Can\'t guess disulfide bridges.'.format(residues[r1]))

    sd = squareform(pdist(mol.coords[idx, :, mol.frame]))
    sd[np.diag_indices(sd.shape[0])] = thresh+1  # Set the diagonal over threshold
    close = sd < thresh
    rows, cols = np.where(close)

    numbonds = np.sum(close, axis=0)
    if np.any(numbonds > 1):
        multibonded_idx1 = np.where(numbonds > 1)[0]
        multibonded_indexes = np.where(close[multibonded_idx1])
        multibonded_idx1 = multibonded_idx1[multibonded_indexes[0]]
        multibonded_idx2 = multibonded_indexes[1]
        pairs = [(str(residues[r]), str(residues[c])) for r, c in zip(multibonded_idx1, multibonded_idx2)]
        raise RuntimeError('Sulphur atoms between pairs {} have multiple possible bonds. Cannot guess disulfide bonds. '
                           'Please specify them manually.'.format(pairs))

    uniquerowcols = list(set([tuple(sorted((r, c))) for r, c in zip(rows, cols)]))
    for rc in uniquerowcols:
        disubonds.append([residues[rc[0]], residues[rc[1]]])
        msg = 'Disulfide Bond between: {}\n' \
              '                   and: {}\n'.format(residues[rc[0]], residues[rc[1]])
        print(msg)

    if len(disubonds) == 1:
        logger.info('One disulfide bond was added')
    else:
        logger.info('{} disulfide bonds were added'.format(len(disubonds)))
    return sorted(disubonds, key=lambda x: x[0].resid)
    

def detectCisPeptideBonds(mol):
    from moleculekit.projections.metricdihedral import MetricDihedral, Dihedral

    protsel = mol.atomselect('protein and backbone and name C CA N')
    if np.sum(protsel) < 4: # Less atoms than dihedral
        return
        
    dih = Dihedral.proteinDihedrals(mol, sel=protsel, dih=('omega',))

    metr = MetricDihedral(dih=dih, sincos=False)
    data = metr.project(mol)
    mapping = metr.getMapping(mol)

    frames, idxs = np.where(np.abs(data) < 120)

    for ii in np.unique(idxs):
        currframes = frames[idxs == ii]
        nframes = len(currframes)
        description = mapping.loc[ii].description
        atomIndexes = mapping.loc[ii].atomIndexes

        currframes_str = "{}".format(currframes)
        if nframes> 5:
            currframes_str = "[{} ... {}]".format(currframes[0], currframes[-1])

        logger.warning("Found cis peptide bond in {} frames: {} in the omega diheral \"{}\" with indexes {}".format(nframes, currframes_str, description, atomIndexes))


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


def _checkLongResnames(mol, aliasresidues):
    for resname in np.unique(mol.resname):
        if len(resname) > 4 and resname not in aliasresidues:
            raise RuntimeError('Too long residue names in Molecule. Please give a 4-letter alias to these residues with the aliasresidues option.')


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
    mol1 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        Molecule for which to calculate the convex hull
    mol2 : :class:`Molecule <moleculekit.molecule.Molecule>` object
        Molecule which contains the atoms which we check if they are within the hull
    hullsel : str
        Atom selection string for atoms in mol1 from which to calculate the convex hull.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    removesel : str
        Atom selection string for atoms in mol2 from which to remove the ones which are within the hull.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__

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
    memb : :class:`Molecule <moleculekit.molecule.Molecule>` object
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

    from moleculekit.molecule import Molecule
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
    from moleculekit.molecule import Molecule
    from moleculekit.tools.autosegment import autoSegment
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

