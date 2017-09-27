from htmd.molecule.molecule import Molecule
from glob import glob
import numpy as np
from htmd.membranebuilder.wrappingdist import wrapping_dist_numba, wrapping_dist_python

# Based on: http://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0134(199601)24:1%3C92::AID-PROT7%3E3.0.CO;2-Q/epdf
# Structure, energetics, and dynamics of lipid–protein interactions: A molecular dynamics study
# of the gramicidin A channel in a DMPC bilayer


class _Lipid:
    def __init__(self, resname=None, headname=None, mol=None, xyz=None, rot=None, neighbours=None, rings=None,
                 area=None):
        self.resname = resname
        self.headname = headname
        self.mol = mol
        self.xyz = xyz
        self.rot = rot
        self.neighbours = neighbours
        self.rings = rings
        self.area = area

    def __repr__(self):
        return '<{}.{} object at {} {}>'.format(self.__class__.__module__, self.__class__.__name__, hex(id(self)),
                                                self.__str__())

    def __str__(self):
        s = ''
        if self.resname is not None:
            s += 'resname: {} '.format(self.resname)
        if self.headname is not None:
            s += 'headname: {} '.format(self.headname)
        if self.xyz is not None:
            s += 'xyz: {} '.format(self.xyz)
        if self.neighbours is not None:
            s += 'neigh: {} '.format(len(self.neighbours))
        if self.mol is not None:
            s += 'mol: {} '.format(id(self.mol))
        return s[:-1]


# def findNonClashingConformation(othermols, molname, resid, p, layer, positions, bilayerbuffer=0, thresh=2):
#     from scipy.spatial.distance import cdist
#
#     if len(othermols) == 0:
#         mol, _ = loadLipid(molname, np.random.choice(files[molname]), resid, positions[p], layer, bilayerbuffer)
#         return mol
#
#     reprcoords = [om.coords[0][np.newaxis, :] for om in othermols]  # Get first coordinate as representative
#     reprcoords = np.concatenate(reprcoords, axis=0).squeeze()
#     allcoords = np.array([om.coords for om in othermols], dtype=object)
#     # from IPython.core.debugger import Tracer
#     # Tracer()()
#     minclash = 999
#     minclashrot = None
#     for k in np.random.permutation(len(files[molname])):
#         mol, head = loadLipid(molname, files[molname][k], resid, positions[p], layer, bilayerbuffer)
#
#         for i in range(100):
#             headpos = mol.coords[mol.name == head, :, :].squeeze()
#             mol.rotateBy(rotationMatrix([0, 0, 1], 0.06), headpos)  # Rotate around the lipid head by 3.4 degrees at a time
#             notp = np.ones(len(positions), dtype=bool)
#             notp[p:] = False
#             repd = cdist(np.atleast_2d(positions[p]), np.atleast_2d(positions[notp]))
#
#             closelipids = np.where(repd.squeeze() < 20)[0]
#             if len(closelipids) == 0:
#                 return mol
#
#             # from IPython.core.debugger import Tracer
#             # Tracer()()
#
#             alld = cdist(mol.coords.squeeze(), np.concatenate(allcoords[closelipids], axis=0).squeeze())
#             if alld.min() > thresh:
#                 return mol
#             if alld.min() < minclash:
#                 minclash = alld.min()
#     print('min clash ' + str(minclash))
#     return mol


# def loadLipid(molname, fname, resid, pos, layer, bilayerbuffer):
#     def positionMolecule(molname, mol, head, pos, layer):
#         headpos = mol.coords[mol.name == head, :, :].squeeze()  # TODO: I will need a dictionary of the head atoms
#         zpos = thickness[molname] / 2 + bilayerbuffer
#         mol.moveBy([pos[0], pos[1], zpos] - headpos)
#         if layer == 'lower':
#             mol.coords[:, 2, :] *= -1
#
#     mol = Molecule(fname)
#     mol.remove('water', _logger=False)
#     mol.resid[:] = resid
#     head = headatoms[molname]
#     positionMolecule(molname, mol, head, pos, layer)
#     return mol, head


# def createMols(molpos, lipids, layer, resid, bilayerbuffer=0):
#     allmols = []
#     totalatoms = 0
#     for m, positions in enumerate(molpos):
#         molname = lipids[m][0]
#         for p in range(len(positions)):
#             resid += 1
#             mol, _ = loadLipid(molname, np.random.choice(files[molname]), resid, positions[p], layer, bilayerbuffer)
#             # mol = findNonClashingConformation(allmols, molname, resid, p, layer, positions, bilayerbuffer=bilayerbuffer)
#             allmols.append(mol)
#             totalatoms += mol.numAtoms
#     return allmols, totalatoms, resid


def _createLipids(lipidratio, area, lipiddb, files, leaflet=None):
    lipiddb = lipiddb.to_dict(orient='index')
    ratiosAPL = np.array([x[1] * lipiddb[x[0]]['APL'] for x in lipidratio])
    # Calculate the total areas per lipid type
    areaspl = area * (ratiosAPL / ratiosAPL.sum())
    # Calculate the counts from the total areas
    counts = np.round(areaspl / np.array([lipiddb[x[0]]['APL'] for x in lipidratio])).astype(int)

    lipids = []
    for i in range(len(lipidratio)):
        resname = lipidratio[i][0]
        rings = _detectRings(Molecule(files[resname][0]))
        for k in range(counts[i]):
            if leaflet == 'upper':
                xyz = np.array([np.nan, np.nan, lipiddb[resname]['Thickness'] / 2])
            elif leaflet == 'lower':
                xyz = np.array([np.nan, np.nan, - lipiddb[resname]['Thickness'] / 2])
            lipids.append(_Lipid(resname=resname, headname=lipiddb[resname]['Head'], rings=rings, area=lipiddb[resname]['APL'], xyz=xyz))
    return lipids


def _setPositionsLJSim(width, lipids):
    from htmd.membranebuilder.ljfluid import distributeLipids

    sigmas = np.array([2 * np.sqrt(l.area / np.pi) for l in lipids])
    resnames = [l.resname for l in lipids]

    cutoff = min(np.min(width) / 2, 3 * np.max(sigmas))
    pos = distributeLipids(width + [2 * cutoff], resnames, sigmas, cutoff=cutoff)
    for i in range(len(lipids)):
        lipids[i].xyz[:2] = pos[i, :2]


def _createMembraneMolecule(lipids):
    from htmd.rotationmatrix import rotationMatrix

    allmols = []
    for i, l in enumerate(lipids):
        mol = l.mol.copy()
        headpos = mol.coords[mol.name == l.headname].flatten()[np.newaxis, :]
        mol.moveBy(-headpos)
        mol.rotateBy(rotationMatrix([0, 0, 1], np.deg2rad(l.rot)))
        mol.moveBy(l.xyz)
        mol.resid[:] = i
        allmols.append(mol)
    
    def mergeMols(mollist):  # Divide and conquer approach for merging molecules
        while len(mollist) > 1:
            mollist[0].append(mollist[1])
            mollist = [mollist[0]] + mergeMols(mollist[2:])
        if len(mollist) == 1:
            return mollist
        return []

    return mergeMols(allmols)[0]


def _detectRings(mol):
    import networkx as nx
    bonds = mol._guessBonds()

    G = nx.Graph()
    G.add_edges_from(bonds)
    cycles = np.array(nx.cycle_basis(G))
    if len(cycles) == 0:
        return None

    cyclelen = np.array([len(c) for c in cycles])
    fivesix = cycles[np.where((cyclelen == 5) | (cyclelen == 6))[0]]

    return fivesix


def _findNeighbours(lipids, box):
    xypos = np.vstack([l.xyz[:2] for l in lipids])
    # from scipy.spatial.distance import pdist, squareform
    # dists = pdist(xypos)
    # dists = np.triu(squareform(dists))  # Only store each neighbour once by zeroing lower triangle
    # for i in range(len(lipids)):
    #     lipids[i].neighbours = np.where((dists[i] < 11) & (dists[i] != 0))[0]

    for i in range(len(lipids)):
        dist = wrapping_dist_python(xypos[i, :], xypos[i + 1:, :], box)
        lipids[i].neighbours = i + 1 + np.where(dist < 11)[0]


def _loadMolecules(lipids, files):
    from htmd.rotationmatrix import rotationMatrix
    # Create Molecules
    for l in lipids:
        randidx = np.random.randint(len(files[l.resname]))
        mol = Molecule(files[l.resname][randidx])
        mol.filter('not water', _logger=False)
        if l.xyz[2] < 0:
            mol.rotateBy(rotationMatrix([1, 0, 0], np.deg2rad(180)))  # Rotate the lower leaflet lipids upside down
        l.mol = mol
        l.rot = np.random.random() * 360 - 180  # Random starting rotation


def _locateLipidFiles(folder, lipidnames):
    import os
    files = {}
    for m in lipidnames:
        files[m] = glob(os.path.join(folder, m, '*', '*.crd'))
        if len(files[m]) == 0:
            raise RuntimeError('Could not locate crd files for lipid "{}" in folder {}'.format(m, folder))
    return files


def buildMembrane(xysize, ratioupper, ratiolower, waterbuff=20, equilibrate=True, outdir='./buildmemb/', lipidf=None):
    """ Construct a membrane containing arbitrary lipids and ratios of them.

    Parameters
    ----------
    xysize : list
        A list containing the size in x and y dimensions of the membrane in Angstroms
    ratioupper : list of lists
        A list containing sublists indicating the molecule name and the ratio of that molecule for the upper layer
    ratiolower : list of lists
        Same as ratioupper but for the lower layer
    waterbuff : float
        The z-dimension size of the water box above and below the membrane
    equilibrate : bool
        If True it equilibrates the membrane
    outdir : str
        A folder in which to store the psf and pdb files
    lipidf : str
        The path to the starting lipid conformations

    Returns
    -------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule`
        The resulting membrane including surrounding waters

    Examples
    --------
    >>> lipidratioupper = [['popc', 10], ['chl1', 1]]
    >>> lipidratiolower = [['popc', 8], ['chl1', 2]]
    >>> width = [50, 100]
    >>> res = buildMembrane(width, lipidratioupper, lipidratiolower)
    """
    from htmd.membranebuilder.ringpenetration import resolveRingPenetrations
    from htmd.builder.solvate import solvate
    from htmd.builder.charmm import build
    from htmd.util import tempname
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import os
    import pandas as pd

    if lipidf is None:
        lipidf = os.path.join(home(), 'membranebuilder', 'lipids')
    lipiddb = pd.read_csv(os.path.join(home(), 'membranebuilder', 'lipiddb.csv'), index_col='Name')

    uqlip = np.unique([x[0] for x in ratioupper] + [x[0] for x in ratiolower])
    files = _locateLipidFiles(lipidf, uqlip)

    area = np.prod(xysize)
    lipids = _createLipids(ratioupper, area, lipiddb, files, leaflet='upper')
    lipids += _createLipids(ratiolower, area, lipiddb, files, leaflet='lower')

    _setPositionsLJSim(xysize, [l for l in lipids if l.xyz[2] > 0])
    _setPositionsLJSim(xysize, [l for l in lipids if l.xyz[2] < 0])

    _findNeighbours(lipids, xysize)

    _loadMolecules(lipids, files)

    resolveRingPenetrations(lipids, xysize)
    memb = _createMembraneMolecule(lipids)

    # from globalminimization import minimize
    # newpos, newrot = minimize(lipids, xysize + [100], stepxy=0.5, steprot=50, contactthresh=1)
    # for i in range(len(lipids)):
    #     lipids[i].xyz[:2] = newpos[i]
    #     lipids[i].rot = newrot[i]
    #
    # resolveRingPenetrations(lipids, xysize)
    # endmemb = _createMembraneMolecule(lipids)

    minc = memb.get('coords', 'name P').min(axis=0) - 5
    maxc = memb.get('coords', 'name P').max(axis=0) + 5

    mm = [[0, 0, maxc[2] - 2], [xysize[0], xysize[1], maxc[2] + waterbuff]]
    smemb = solvate(memb, minmax=mm)
    mm = [[0, 0, minc[2] - waterbuff], [xysize[0], xysize[1], minc[2] + 2]]
    smemb = solvate(smemb, minmax=mm)

    smemb.moveBy([0, 0, -smemb.coords[:, 2, 0].min()])

    if outdir is None:
        outdir = tempname()
        print('Outdir ', outdir)
    res = build(smemb, ionize=False, stream=['str/lipid/toppar_all36_lipid_cholesterol_model_1.str'], outdir=outdir)

    if equilibrate:
        from htmd.membranebuilder.simulate_openmm import equilibrateSystem
        from shutil import copy, move
        outpdb = tempname(suffix='.pdb')
        charmmf = os.path.join(home(), 'membranebuilder', 'charmm-toppar')
        equilibrateSystem(os.path.join(outdir, 'structure.pdb'), os.path.join(outdir, 'structure.psf'), outpdb, charmmfolder=charmmf)
        res = Molecule(outpdb)
        res.center()
        move(os.path.join(outdir, 'structure.pdb'), os.path.join(outdir, 'starting_structure.pdb'))
        copy(outpdb, os.path.join(outdir, 'structure.pdb'))

    return res


def _findLeastAreaLipid(folder):
    """
    Use this to select the single least stretched conformation from a CHARMM-GUI lipid library
    """
    from glob import glob
    from scipy.spatial.distance import cdist
    ff = glob(folder + '/*/*.crd')
    maxdist = []
    for f in ff:
        m = Molecule(f)
        center = m.coords.mean(axis=0)
        dists = cdist(m.coords[:, :2, 0], center[:2].T)
        maxdist.append(dists.max())
    return ff[np.argmin(maxdist)], np.min(maxdist)


if __name__ == '__main__':
    pass



