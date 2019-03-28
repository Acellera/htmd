# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Rapidly exploring random tree to find reaction coordinate for umbrella sampling!
# Use VMD draw line to visualize! (dashed for failed, solid for connected)
import numpy as np
from moleculekit.util import maxDistance, uniformRandomRotation
from moleculekit.molecule import Molecule
from moleculekit.vmdviewer import getCurrentViewer
from scipy.spatial.distance import cdist


class Tree:
    def __init__(self, initialpoints):
        if isinstance(initialpoints, np.ndarray):
            initialpoints = initialpoints.tolist()
        self.points = initialpoints
        self.parent = [None] * len(self.points)
        self.costs = [0] * len(self.points)

    def addPoint(self, point, parent, dist):
        if isinstance(point, np.ndarray):
            point = point.tolist()
        if isinstance(dist, np.ndarray):
            dist = dist[0]
        self.points.append(point)
        self.parent.append(parent)
        self.costs.append(self.costs[parent] + dist)


def _getCoordinates(mol, ligandsel, ligcom):
    mol.center(sel='protein')
    protcoor = mol.get('coords', sel='protein and not {}'.format(ligandsel))
    if ligcom:
        ligcoor = np.mean(mol.get('coords', sel=ligandsel), axis=0)[np.newaxis, :]
    else:
        ligcoor = mol.get('coords', sel=ligandsel)
    return protcoor, ligcoor


def _newPoint(p_target, p_begin, step):
    v = p_target - p_begin
    norm_v = v / np.linalg.norm(v)
    if isinstance(step, np.ndarray):
        return p_begin + (step[:, np.newaxis] * norm_v)
    else:
        return p_begin + (step * norm_v)


def _dist(a, b):
    if isinstance(a, list):
        a = np.array(a)
    if isinstance(b, list):
        b = np.array(b)
    return cdist(np.atleast_2d(a), np.atleast_2d(b))


def _getNearest(coords, point):
    dist = _dist(coords, point)
    idx = np.argmin(dist)
    return idx, coords[idx]


def _collision(coords, point, buffer=2):
    dist = _dist(coords, point)
    if np.any(dist < buffer):
        return True
    else:
        return False


def _pointsOnSphere(radius, numsamples=1000):
    pointcoords = np.zeros((numsamples, 3))
    for i in range(numsamples):
        point = Molecule()
        point.empty(1)
        point.moveBy([0, 0, radius])
        point.rotateBy(uniformRandomRotation())
        pointcoords[i, :] = np.squeeze(point.coords)

    return pointcoords


def _collisionFreePath(obstacles, start, stop, colldist, step):
    linepoints = _newPoint(stop, start, np.arange(step, _dist(start, stop), step=step))
    distances = _dist(linepoints, obstacles)
    if np.any(distances <= colldist):
        return False
    else:
        return True


def _collisionFreeNeighbours(tree, point, radius, obstacles, colldist, step):
    dists = _dist(tree.points, point)
    near = np.where(dists < radius)[0]
    nn = []
    for n in near:
        if _collisionFreePath(obstacles, point, tree.points[n], colldist, step):
            nn.append(n)
    near = nn
    return near, dists[near]


def _rewire(tree, near, neardist, newidx):
    # Remove parent from near neighbours
    parentidx = np.where(near == tree.parent[newidx])[0]
    near = np.delete(near, parentidx)
    neardist = np.delete(neardist, parentidx)
    if len(near) == 0:
        return

    costs = np.array(tree.costs)
    newcosts = np.squeeze(neardist) + costs[newidx]
    shorter = np.where(newcosts < costs[near])[0]
    for s in shorter:
        #from IPython.core.debugger import Tracer
        #Tracer()()
        tree.costs[near[s]] = newcosts[s]
        tree.parent[near[s]] = newidx


def _chooseParent(tree, near, neardist):
    if len(near) == 0:
        raise AssertionError('Should include at least p_near.')

    costs = np.array(tree.costs)
    closecosts = costs[near] + np.squeeze(neardist)
    shortest = np.argmin(closecosts)
    return near[shortest], neardist[shortest]


def _endCondition(protcoor, target, curr, outdist, method='target'):
    if method == 'target':
        return _dist(target, curr) < 2
    if method == 'exited':
        return np.all(_dist(protcoor, curr) >= outdist)


def _pathOptimize(tree, currnode, obstacles, colldist, step):
    start = np.array(tree.points[currnode])
    endnode = tree.parent[currnode]
    prevend = endnode

    while True:
        if endnode is None:
            tree.parent[currnode] = prevend
            tree.costs[currnode] = tree.costs[prevend] + _dist(tree.points[currnode], tree.points[prevend])
            return

        stop = np.array(tree.points[endnode])
        if _collisionFreePath(obstacles, start, stop, colldist, step):
            prevend = endnode
            endnode = tree.parent[endnode]
        else:
            tree.parent[currnode] = prevend
            tree.costs[currnode] = tree.costs[prevend] + _dist(tree.points[currnode], tree.points[prevend])
            _pathOptimize(tree, prevend, obstacles, colldist, step)
            return


def _getBeacons(tree, endnode):
    beacons = []
    currnode = tree.parent[endnode]
    while tree.parent[currnode] is not None:
        beacons.append(currnode)
        currnode = tree.parent[currnode]
    return beacons


def _randomPoint(min, max):
    lens = max - min
    return np.random.rand(3) * lens + min


def rrtstarsmart(mol, ligandsel, step=1, maxiter=int(1E6), ligcom=False, colldist=2, outdist=8, radius=2.5):
    protcoor, ligcoor = _getCoordinates(mol, ligandsel, ligcom)
    viewer = _prepareViewer(mol, ligandsel)

    mincoor = np.squeeze(np.min(mol.coords, axis=0) - 10)
    maxcoor = np.squeeze(np.max(mol.coords, axis=0) + 10)
    print(mincoor, maxcoor)

    spherecoor = _pointsOnSphere(maxDistance(mol) + 5)

    tree = Tree(ligcoor)
    nocol = None

    for i in range(maxiter):
        print(i)
        #p_rand = spherecoor[np.random.randint(spherecoor.shape[0]), :]
        p_rand = _randomPoint(mincoor, maxcoor)
        parent_idx, p_near = _getNearest(tree.points, p_rand)
        p_new = _newPoint(p_rand, p_near, step)

        if _collision(protcoor, p_new, buffer=colldist):
            continue

        print('no collision')
        near, neardist = _collisionFreeNeighbours(tree, p_new, radius, protcoor, colldist, step)
        parent_idx, dist = _chooseParent(tree, near, neardist)
        tree.addPoint(p_new, parent_idx, dist)
        _rewire(tree, near, neardist, len(tree.points)-1)

        if _endCondition(protcoor, p_rand, p_new, outdist, method='exited'):
            break

    for i in range(len(tree.points)):
        if tree.parent[i] is not None:
            _drawline(viewer, tree.points[i], tree.points[tree.parent[i]])

    print('optimizing')
    _pathOptimize(tree, len(tree.points)-1, protcoor, colldist, step)
    print('finished optimization')
    # TODO: The algorithm is not finished. Normally you would get beacons and continue sampling around them to optimize
    # TODO: the path even more. However, we don't really need this in our case, unless I decide to make a real full impl
    # TODO: http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=6284384

    currnode = len(tree.points)-1
    viewer.send('draw color green')
    while tree.parent[currnode] is not None:
        _drawline(viewer, tree.points[currnode], tree.points[tree.parent[currnode]])
        currnode = tree.parent[currnode]

    return


def rrt(mol, ligandsel, step=1, maxiter=int(1E6), ligcom=False, colldist=2, outdist=8):
    protcoor, ligcoor = _getCoordinates(mol, ligandsel, ligcom)
    viewer = _prepareViewer(mol, ligandsel)

    spherecoor = _pointsOnSphere(maxDistance(mol) + 5)

    tree = Tree(ligcoor)
    nocol = None

    for i in range(maxiter):
        print(i)
        if nocol is None:
            p_rand = spherecoor[np.random.randint(spherecoor.shape[0]), :]  # TODO: Make it clever. Discard failed points
        else:  # Reusing the same point that caused no collision
            p_rand = nocol
        parent_idx, p_near = _getNearest(tree.points, p_rand)
        p_new = _newPoint(p_rand, p_near, step)

        if _collision(protcoor, p_new, buffer=colldist):
            nocol = None
            continue
        nocol = p_rand

        print('no collision')
        tree.addPoint(p_new, parent_idx, step)

        _drawline(viewer, p_near, p_new)
        if _dist(p_rand, p_new) < 2:
            break


def raytracing(mol, ligandsel, step=1, colldist=2, outdist=8, ligcom=False, numsamples=2000, ratioexposed=0, vmd=True):
    protcoor, ligcoor = _getCoordinates(mol, ligandsel, ligcom)

    spherecoor = _pointsOnSphere(maxDistance(mol) + 5, numsamples=numsamples)
    if vmd:
        _viewSphere(spherecoor)
        viewer = _prepareViewer(mol, ligandsel)

    distances = _dist(ligcoor, spherecoor)

    from joblib import Parallel, delayed
    results = Parallel(n_jobs=-2, verbose=11)(
        delayed(parallelfunc)(j, spherecoor[j, :], ligcoor, protcoor, step, colldist, outdist, distances[:, j])
        for j in range(spherecoor.shape[0]))

    points = []
    pointdist = []
    shortpoints = []
    for r in results:
        points += r[0]
        pointdist += r[1]
        shortpoints += r[2]

    if len(points) == 0:
        raise RuntimeError('No ligand atoms can exit the pocket without {}A clashes.'.format(colldist))

    points = np.array(points)
    numexposed = len(np.unique(points[:, 0]))
    percentexposed = (numexposed / ligcoor.shape[0]) * 100
    if numexposed < (ratioexposed * ligcoor.shape[0]):
        raise RuntimeError('Only {:.1f}% ligand atoms can exit the pocket without {}A clashes. '
                           'This collides with the user-defined required {:.1f}% exposed ligand atoms.'.format(
                            percentexposed, colldist, ratioexposed*100))
    print('{:.1f}% ligand atoms can exit the pocket without {}A clashes.'.format(percentexposed, colldist))


    from scipy.stats import mode
    modesphere = mode(points[:, 1]).mode[0]
    idx = points[:, 1] == modesphere
    meanlig = np.mean(ligcoor[points[idx, 0], :], axis=0)

    if vmd:
        viewer.send('draw color green')
        _drawline(viewer, meanlig, spherecoor[modesphere, :])

        for idx, p in enumerate(points.tolist()):
            if p[1] == modesphere:
                viewer.send('draw color red')
                _drawline(viewer, shortpoints[idx], spherecoor[p[1], :])
                viewer.send('draw color green')
                _drawline(viewer, ligcoor[p[0], :], shortpoints[idx])

    return spherecoor[modesphere, :] - meanlig


def parallelfunc(j, spherecoor, ligcoor, protcoor, step, colldist, outdist, distances):
    points = []
    pointdist = []
    shortpoints = []
    for i in range(ligcoor.shape[0]):
        linepoints = _newPoint(spherecoor, ligcoor[i, :], np.arange(step, distances[i], step=step))
        collisions = _dist(linepoints, protcoor)
        if not np.any(collisions <= colldist):
            linep_mincoll = np.min(collisions, axis=1)
            idx = np.min(np.where(linep_mincoll >= outdist)[0])
            points.append([i, j])
            pointdist.append(step * (idx+1))
            shortpoints.append(linepoints[idx, :])
    return points, pointdist, shortpoints


def _drawline(viewer, start, end, dashed=True):
    if dashed:
        append = 'style dashed'
    else:
        append = ''
    viewer.send('draw line {{ {} }} {{ {} }} {}'.format(' '.join(map(str, start)), ' '.join(map(str, end)), append))


def _prepareViewer(mol, ligandsel):
    mol.view(sel='protein', style='lines', hold=True)
    mol.view(sel=ligandsel, style='licorice', hold=True)
    mol.view()

    viewer = getCurrentViewer()
    viewer.send('draw color red')
    viewer.send('color Display Background black')
    viewer.send('draw materials off')
    return viewer


def _viewSphere(spherecoor):
    spheremol = Molecule()
    spheremol.empty(spherecoor.shape[0])
    spheremol.coords = np.atleast_3d(spherecoor)
    spheremol.view(guessBonds=False)
