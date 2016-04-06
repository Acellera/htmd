from htmd import *
import logging
logger = logging.getLogger(__name__)


def removeLipidsInProtein(prot, memb):
    """ Calculates the convex hull of the protein. If a lipid lies inside the hull it gets removed.

    This does not work well for lipids crossing out of the hull. If even one atom of the lipid is outside it will
    change the hull and will not get removed. I assume it will get removed by the clashes with the protein though.
    """
    # TODO: Do the same with Morphological Snakes
    from scipy.spatial import ConvexHull
    from htmd.builder.builder import sequenceID
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


def tilemembrane(memb, xmin, ymin, xmax, ymax):
    """ Tile the membrane in the X and Y dimensions to reach a specific size.
    Returns
    -------
    megamemb :
        A big membrane Molecule
    """
    from htmd.progress.progress import ProgressBar
    from htmd.builder.builder import sequenceID
    memb = memb.copy()
    memb.resid = sequenceID(memb.resid)

    minmemb = np.min(memb.get('coords', 'water'), axis=0).flatten()

    size = np.max(memb.get('coords', 'water'), axis=0) - np.min(memb.get('coords', 'water'), axis=0)
    size = size.flatten()
    xreps = int(np.ceil((xmax - xmin) / size[0]))
    yreps = int(np.ceil((ymax - ymin) / size[1]))

    logger.info('Replicating Membrane {}x{}'.format(xreps, yreps))

    megamemb = Molecule()
    bar = ProgressBar(xreps*yreps, description='Replicating Membrane')
    k = 0
    for x in range(xreps):
        for y in range(yreps):
            tmpmemb = memb.copy()
            xpos = xmin + x * size[0]
            ypos = ymin + y * size[1]

            tmpmemb.moveBy([-float(minmemb[0])+xpos, -float(minmemb[1])+ypos, 0])
            sel = 'same resid as (x > {} or y > {})'.format(xmax, ymax)
            tmpmemb.remove(sel, _logger=False)
            tmpmemb.set('segid', 'M{}'.format(k))

            megamemb.append(tmpmemb)
            k += 1
            bar.progress()
    bar.stop()
    return megamemb


def minimalrotation(prot):
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
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return rho, phi

    angle, _ = cart2pol(xa, ya)
    return angle + np.radians(45)



