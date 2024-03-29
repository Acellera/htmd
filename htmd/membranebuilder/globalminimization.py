# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np


def _wrapping_dist(X, Y, box):
    vecdiff = np.swapaxes(X[..., None] - Y.T, 1, 2)
    wrappeddiff = vecdiff - box * np.round(vecdiff / box)
    return np.sqrt(np.sum(wrappeddiff * wrappeddiff, axis=2))


def _totalContacts(x, *args):
    from moleculekit.util import rotationMatrix

    lipids = args[0]
    headnames = args[1]
    zpos = args[2]
    thresh = args[3]
    neighbours = args[4]
    numlips = args[5]
    box = args[6]
    newxy = np.array(x[: numlips * 2]).reshape((-1, 2))
    rots = x[numlips * 2 :]

    allcoords = []
    for i in range(numlips):
        cc = np.squeeze(lipids[i].mol.coords)
        headpos = cc[lipids[i].mol.name == headnames[i]].flatten()[np.newaxis, :]
        cc = cc - headpos  # Center it on the head position
        newloc = np.hstack((newxy[i], zpos[i]))[
            np.newaxis, :
        ]  # Calculate new head position

        # Doing the rotation
        M = rotationMatrix([0, 0, 1], np.deg2rad(rots[i]))
        newcoords = np.dot(cc, np.transpose(M))

        allcoords.append(newcoords + newloc)

    allcoords = np.array(allcoords, dtype=object)
    numcontacts = 0
    for i in range(len(allcoords)):
        if len(neighbours[i]) == 0:
            continue
        neighcoor = np.vstack(allcoords[neighbours[i]])
        dists = _wrapping_dist(allcoords[i], neighcoor, np.array(box, dtype=float))
        numcontacts += np.count_nonzero(dists < thresh)

    return numcontacts


class _RandomDisplacementBounds(object):
    """random displacement with bounds"""

    def __init__(self, bounds, stepsizes):
        self.bounds = bounds
        self.stepsizes = stepsizes

    def __call__(self, x):
        """take a random step but ensure the new position is within the bounds"""
        xnew = np.clip(
            x + (np.random.uniform(size=x.shape) - 0.5) * self.stepsizes,
            self.bounds[:, 0],
            self.bounds[:, 1],
        )
        return xnew


def minimize(lipids, box, stepxy=0.5, steprot=50, contactthresh=2.6):
    from scipy.optimize import basinhopping

    # rotate in 10deg increments
    # translate in a 2x2 box in 0.25A increments
    # swap lipid conformer?
    headnames = [ll.headname for ll in lipids]
    zpos = [ll.xyz[2] for ll in lipids]
    neighbours = [ll.neighbours for ll in lipids]

    pos = np.vstack([ll.xyz[:2] for ll in lipids])
    x0 = pos.flatten().tolist()
    numlips = len(lipids)
    bounds = [(x - 1, x + 1) for x in x0]
    x0 += [180] * numlips  # Add the rotations
    bounds += [(0, 360)] * numlips  # Add the rotations
    stepsizes = np.ones(numlips * 2) * stepxy
    stepsizes = np.hstack((stepsizes, np.ones(numlips) * steprot))  # Add the rotations

    # define the new step taking routine and pass it to basinhopping
    take_step = _RandomDisplacementBounds(np.vstack(bounds), stepsizes=stepsizes)
    minimizer_kwargs = dict(
        method="L-BFGS-B",
        bounds=bounds,
        args=(lipids, headnames, zpos, contactthresh, neighbours, numlips, box),
    )
    res = basinhopping(
        _totalContacts,
        x0,
        minimizer_kwargs=minimizer_kwargs,
        disp=True,
        take_step=take_step,
        niter=10,
    )
    newpos = res.x[: numlips * 2].reshape((-1, 2))
    newrot = res.x[numlips * 2 :]
    return newpos, newrot
