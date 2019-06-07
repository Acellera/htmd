# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.molecule import Molecule
from glob import glob
import numpy as np
import logging

logger = logging.getLogger(__name__)

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


def listLipids():
    """ Lists all available lipids

    Examples
    --------
    >>> from htmd.membranebuilder.build_membrane import buildMembrane 
    >>> build_membrane.listLipids()
    ---- Lipids list: ...

    """

    from htmd.home import home
    import os
    from natsort import natsorted

    membranebuilderhome = os.path.join( home(shareDir=True), 'membranebuilder', 'lipids', '' )
    lipids = natsorted(glob( os.path.join(membranebuilderhome, '*', '')))
    print('---- Lipids list: ' + membranebuilderhome + ' ----')
    for l in lipids:
        print('- ', os.path.basename(os.path.abspath(l)))
    print('* Lipid DB file: ' + os.path.join(membranebuilderhome, 'lipiddb.csv'))



def _createLipids(lipidratio, area, lipiddb, files, leaflet=None):
    lipiddb = lipiddb.to_dict(orient='index')
    lipidnames = list(lipidratio.keys())
    ratiosAPL = np.array([lipidratio[lipn] * lipiddb[lipn]['APL'] for lipn in lipidnames])
    # Calculate the total areas per lipid type
    areaspl = area * (ratiosAPL / ratiosAPL.sum())
    # Calculate the counts from the total areas
    counts = np.round(areaspl / np.array([lipiddb[lipn]['APL'] for lipn in lipidnames])).astype(int)

    lipids = []
    for i in range(len(lipidnames)):
        resname = lipidnames[i]
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
    from moleculekit.util import rotationMatrix

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


def wrapping_dist_python(coor1, coor2, box):
    assert (coor1.ndim == 1) or (coor2.ndim == 1)
    dist = coor1 - coor2
    dist = dist - box * np.round(dist / box)
    return np.sqrt(np.sum(dist * dist, 1))


def _findNeighbours(lipids, box):
    xypos = np.vstack([l.xyz[:2] for l in lipids])

    for i in range(len(lipids)):
        dist = wrapping_dist_python(xypos[i, :], xypos[i + 1:, :], box)
        lipids[i].neighbours = i + 1 + np.where(dist < 11)[0]


def _loadMolecules(lipids, files):
    from moleculekit.util import rotationMatrix
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
        files[m] = glob(os.path.join(folder, m, '*.pdb'))
        if len(files[m]) == 0:
            raise RuntimeError('Could not locate pdb files for lipid "{}" in folder {}'.format(m, folder))
    return files


def buildMembrane(xysize, ratioupper, ratiolower, waterbuff=20, minimplatform='CPU', equilibrate=True, equilplatform='CUDA', outdir=None, lipidf=None):
    """ Construct a membrane containing arbitrary lipids and ratios of them.

    Parameters
    ----------
    xysize : list
        A list containing the size in x and y dimensions of the membrane in Angstroms
    ratioupper : dict
        A dict with keys the molecule names and the ratio of that molecule for the upper layer
    ratiolower : dict
        Same as ratioupper but for the lower layer
    waterbuff : float
        The z-dimension size of the water box above and below the membrane
    minimplatform : str
        The platform on which to run the minimization ('CUDA' or 'CPU')
    equilibrate : bool
        If True it equilibrates the membrane
    equilplatform : str
        The platform on which to run the equilibration ('CUDA' or 'CPU')
    outdir : str
        A folder in which to store the psf and pdb files
    lipidf : str
        The path to the folder containing the single-lipid PDB structures as well as the lipid DB file

    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule`
        The resulting membrane including surrounding waters

    Examples
    --------
    >>> lipidratioupper = {'popc': 10, 'chl1': 1}
    >>> lipidratiolower = {'popc': 8, 'chl1': 2}
    >>> width = [50, 100]
    >>> res = buildMembrane(width, lipidratioupper, lipidratiolower)
    """
    from htmd.membranebuilder.ringpenetration import resolveRingPenetrations
    from htmd.builder.solvate import solvate
    from htmd.builder.charmm import build
    from htmd.util import tempname
    from moleculekit.molecule import Molecule
    from htmd.home import home
    import os
    import pandas as pd

    if lipidf is None:
        lipidf = os.path.join(home(shareDir=True), 'membranebuilder', 'lipids')
    lipiddb = pd.read_csv(os.path.join(lipidf, 'lipiddb.csv'), index_col='Name')

    uqlip = np.unique(list(ratioupper.keys()) + list(ratiolower.keys()))
    files = _locateLipidFiles(lipidf, uqlip)

    area = np.prod(xysize)
    lipids = _createLipids(ratioupper, area, lipiddb, files, leaflet='upper')
    lipids += _createLipids(ratiolower, area, lipiddb, files, leaflet='lower')

    _setPositionsLJSim(xysize, [l for l in lipids if l.xyz[2] > 0])
    _setPositionsLJSim(xysize, [l for l in lipids if l.xyz[2] < 0])

    _findNeighbours(lipids, xysize)

    _loadMolecules(lipids, files)

    # from globalminimization import minimize
    # newpos, newrot = minimize(lipids, xysize + [100], stepxy=0.5, steprot=50, contactthresh=1)
    # for i in range(len(lipids)):
    #     lipids[i].xyz[:2] = newpos[i]
    #     lipids[i].rot = newrot[i]

    resolveRingPenetrations(lipids, xysize)
    memb = _createMembraneMolecule(lipids)

    minc = memb.get('coords', 'name P').min(axis=0) - 5
    maxc = memb.get('coords', 'name P').max(axis=0) + 5

    mm = [[0, 0, maxc[2] - 2], [xysize[0], xysize[1], maxc[2] + waterbuff]]
    smemb = solvate(memb, minmax=mm)
    mm = [[0, 0, minc[2] - waterbuff], [xysize[0], xysize[1], minc[2] + 2]]
    smemb = solvate(smemb, minmax=mm)

    smemb.moveBy([0, 0, -smemb.coords[:, 2, 0].min()])

    if outdir is None:
        outdir = tempname()
        logger.info('Outdir {}'.format(outdir))
    res = build(smemb, ionize=False, stream=['str/lipid/toppar_all36_lipid_cholesterol_model_1.str'], outdir=outdir)

    if equilibrate:
        from htmd.membranebuilder.simulate_openmm import equilibrateSystem
        from shutil import copy, move
        outpdb = tempname(suffix='.pdb')
        charmmf = os.path.join(home(shareDir=True), 'membranebuilder', 'charmm-toppar')
        equilibrateSystem(os.path.join(outdir, 'structure.pdb'), os.path.join(outdir, 'structure.psf'), outpdb,
                          charmmfolder=charmmf, equilplatform=equilplatform, minimplatform=minimplatform)
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



