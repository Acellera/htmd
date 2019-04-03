# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import logging
logger = logging.getLogger(__name__)


def least_square_fit_plane(coords):
    com = coords.mean(axis=0)
    numatm = coords.shape[0]

    axes = np.zeros((numatm, 3))
    for i in range(numatm):
        p1 = coords[i]
        if i == numatm-1:
            p2 = coords[0]
        else:
            p2 = coords[i + 1]
        a = np.cross(p1, p2)
        axes += a
    u, d, v = np.linalg.svd(axes)
    return v[0], com


def getAllNeighbours(lipids, idx):
    neigh = list(lipids[idx].neighbours)
    for j, l in enumerate(lipids):
        if idx in l.neighbours:
            neigh.append(j)
    return neigh


def wrap(coor, com2, box):
    com1 = coor.mean(axis=0)
    for i in range(2):
        if (com1[i] - com2[i]) > (box[i] / 2):
            coor[:, i] -= box[i]
        if (com1[i] - com2[i]) < -(box[i] / 2):
            coor[:, i] += box[i]
    return coor


def moveLipidToPos(mol, lip):
    from moleculekit.util import rotationMatrix
    mol = mol.copy()
    headpos = mol.coords[mol.name == lip.headname].flatten()[np.newaxis, :]
    mol.moveBy(-headpos)
    mol.rotateBy(rotationMatrix([0, 0, 1], np.deg2rad(lip.rot)))
    mol.moveBy(lip.xyz)
    return np.squeeze(mol.coords[:, :, 0])


def _getRingProperties(ringcoords):
    axis, com = least_square_fit_plane(ringcoords)

    # project atoms to the least square fit plane
    for i, atom in enumerate(ringcoords):
        w = np.dot(axis, atom - com) * axis + com
        ringcoords[i] = com + (atom - w)

    maxd = np.max(np.sqrt(np.sum(np.square(ringcoords - com), axis=1)))
    return axis, com, maxd, ringcoords


def _checkAtomsPenetrating(coords2, bonds2, ringcoords, axis, ringcom, maxd, box):
    # Taken from https://github.com/sunhwan/lipid-pentest
    # find two bonded atoms that are at the opposite side of the plane
    coords2 = wrap(coords2, ringcom, box)

    d = np.sqrt(np.sum(np.square(coords2 - ringcom), axis=1))
    closeidx = np.where(d < 3)[0]
    if len(closeidx) == 0:
        return

    flag = False
    penetrating = None
    for i in closeidx:
        bonded = np.hstack((bonds2[bonds2[:, 0] == i, 1], bonds2[bonds2[:, 1] == i, 0]))
        for j in bonded:
            crd1 = coords2[i, :]
            crd2 = coords2[j, :]
            v1 = np.dot(crd1 - ringcom, axis)
            v2 = np.dot(crd2 - ringcom, axis)
            if v1 * v2 > 0:
                continue

            # point of intersection of the least square fit plane
            s = -np.dot(axis, crd1 - ringcom) / np.dot(axis, crd2 - crd1)
            p = crd1 + s * (crd2 - crd1)

            d = np.sqrt(np.sum(np.square(p - ringcom)))
            if d > maxd:
                continue

            d = 0
            for k in range(0, ringcoords.shape[0]):
                p1 = ringcoords[k] - p
                try:
                    p2 = ringcoords[k + 1] - p
                except:
                    p2 = ringcoords[0] - p
                d += np.arccos(np.dot(p1, p2) / np.linalg.norm(p1) / np.linalg.norm(p2))

            wn = d / 2 / np.pi
            if 1.1 > wn > 0.9:
                from IPython.core.debugger import Tracer
                # Tracer()()
                penetrating = [i, j]
                flag = True
                break
        if flag:
            break
    return penetrating


def _detectRingPenetration(l1, lipids, box):
    lip1 = lipids[l1]
    neighbours = getAllNeighbours(lipids, l1)
    mol1 = lip1.mol
    coords1 = moveLipidToPos(mol1, lip1)
    bonds1 = mol1._getBonds()

    for l2, lip2 in zip(neighbours, lipids[neighbours]):
        if lip2.rings is None:
            continue
        mol2 = lip2.mol
        coords2 = moveLipidToPos(mol2, lip2)

        for r, ring in enumerate(lip2.rings):
            axis, ringcom, maxd, ringcoords = _getRingProperties(coords2[ring, :])

            pen = _checkAtomsPenetrating(coords1, bonds1, ringcoords, axis, ringcom, maxd, box)
            if pen is not None:
                logger.info('Lipid {} ring {} is being penetrated by lipid {} atoms {} {}'.format(l2, r, l1, pen[0], pen[1]))
                return True
    return False


def resolveRingPenetrations(lipids, box):
    lipids = np.array(lipids)
    while True:
        penetrators = []
        for l1 in range(len(lipids)):
            if _detectRingPenetration(l1, lipids, box):
                penetrators.append(l1)
        for p in penetrators:
            lipids[p].rot += 10
        if len(penetrators) == 0:
            logger.info('{} penetrating molecule(s) remaining'.format(len(penetrators)))
            break
        logger.info('{} penetrating molecule(s) remaining'.format(len(penetrators)))

