# (c) 2015-2022 Acellera Ltd http://www.acellera.com
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
    def __init__(
        self,
        resname=None,
        headname=None,
        mol=None,
        xyz=None,
        rot=None,
        neighbours=None,
        rings=None,
        area=None,
    ):
        self.resname = resname
        self.headname = headname
        self.mol = mol
        self.xyz = xyz
        self.rot = rot
        self.neighbours = neighbours
        self.rings = rings
        self.area = area

    def __repr__(self):
        return "<{}.{} object at {} {}>".format(
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self)),
            self.__str__(),
        )

    def __str__(self):
        s = ""
        if self.resname is not None:
            s += "resname: {} ".format(self.resname)
        if self.headname is not None:
            s += "headname: {} ".format(self.headname)
        if self.xyz is not None:
            s += "xyz: {} ".format(self.xyz)
        if self.neighbours is not None:
            s += "neigh: {} ".format(len(self.neighbours))
        if self.mol is not None:
            s += "mol: {} ".format(id(self.mol))
        return s[:-1]


def listLipids():
    """Lists all available lipids

    Examples
    --------
    >>> from htmd.membranebuilder.build_membrane import buildMembrane
    >>> build_membrane.listLipids()
    ---- Lipids list: ...

    """

    from htmd.home import home
    import os
    from natsort import natsorted

    membranebuilderhome = os.path.join(
        home(shareDir=True), "membranebuilder", "lipids", ""
    )
    lipids = natsorted(glob(os.path.join(membranebuilderhome, "*", "")))
    print("---- Lipids list: " + membranebuilderhome + " ----")
    for ll in lipids:
        print("- ", os.path.basename(os.path.abspath(ll)))
    print("* Lipid DB file: " + os.path.join(membranebuilderhome, "lipiddb.csv"))


def _createLipids(lipidratio, area, lipiddb, files, leaflet=None):
    if leaflet not in ("upper", "lower"):
        raise ValueError(
            f"leaflet must be 'upper' or 'lower', got {leaflet!r}"
        )

    lipiddb = lipiddb.to_dict(orient="index")
    lipidnames = list(lipidratio.keys())
    ratiosAPL = np.array(
        [lipidratio[lipn] * lipiddb[lipn]["APL"] for lipn in lipidnames]
    )
    # Calculate the total areas per lipid type
    areaspl = area * (ratiosAPL / ratiosAPL.sum())
    # Calculate the counts from the total areas
    counts = np.round(
        areaspl / np.array([lipiddb[lipn]["APL"] for lipn in lipidnames])
    ).astype(int)

    if (counts == 0).any():
        zero_lipids = [lipidnames[i] for i in range(len(counts)) if counts[i] == 0]
        raise RuntimeError(
            f"Computed lipid count is 0 for {zero_lipids} in {leaflet} leaflet "
            f"given the requested xysize and ratios. Increase the membrane size "
            f"or adjust the ratios."
        )

    z_sign = 1 if leaflet == "upper" else -1
    lipids = []
    for i in range(len(lipidnames)):
        resname = lipidnames[i]
        rings = _detectRings(Molecule(files[resname][0]))
        for k in range(counts[i]):
            xyz = np.array(
                [np.nan, np.nan, z_sign * lipiddb[resname]["Thickness"] / 2]
            )
            lipids.append(
                _Lipid(
                    resname=resname,
                    headname=lipiddb[resname]["Head"],
                    rings=rings,
                    area=lipiddb[resname]["APL"],
                    xyz=xyz,
                )
            )
    return lipids


def _setPositionsLJSim(width, lipids):
    from htmd.membranebuilder.ljfluid import distributeLipids

    sigmas = np.array([2 * np.sqrt(ll.area / np.pi) for ll in lipids])
    resnames = [ll.resname for ll in lipids]

    cutoff = min(np.min(width) / 2, 3 * np.max(sigmas))
    pos = distributeLipids(width + [2 * cutoff], resnames, sigmas, cutoff=cutoff)
    for i in range(len(lipids)):
        lipids[i].xyz[:2] = pos[i, :2]


def _createMembraneMolecule(lipids):
    from moleculekit.util import rotationMatrix

    allmols = []
    numAtoms = 0
    numBonds = 0
    for i, l in enumerate(lipids):
        mol = l.mol.copy()
        headpos = mol.coords[mol.name == l.headname].flatten()[np.newaxis, :]
        mol.moveBy(-headpos)
        mol.rotateBy(rotationMatrix([0, 0, 1], np.deg2rad(l.rot)))
        mol.moveBy(l.xyz)
        mol.resid[:] = i
        numAtoms += mol.numAtoms
        numBonds += len(mol.bonds)
        allmols.append(mol)

    # Merge all the lipids into a single Molecule
    mol = Molecule().empty(numAtoms)
    mol.coords = np.zeros((numAtoms, 3, 1), dtype=np.float32)
    mol.bonds = np.zeros((numBonds, 2), dtype=np.uint32)
    mol.bondtype = np.empty(numBonds, dtype=object)
    start_idx = 0
    bond_idx = 0
    for mm in allmols:
        for prop in ["name", "resname", "resid", "segid", "chain", "element", "coords"]:
            mol.__dict__[prop][start_idx : start_idx + mm.numAtoms] = getattr(mm, prop)
        nb = len(mm.bonds)
        if nb:
            mol.bonds[bond_idx : bond_idx + nb] = mm.bonds + start_idx
            mol.bondtype[bond_idx : bond_idx + nb] = mm.bondtype
            bond_idx += nb
        start_idx += mm.numAtoms

    return mol


def _detectRings(mol):
    import networkx as nx

    bonds = mol._guessBonds()

    G = nx.Graph()
    G.add_edges_from(bonds)
    cycles = nx.cycle_basis(G)
    if len(cycles) == 0:
        return None

    fivesix = [c for c in cycles if len(c) in (5, 6)]
    return fivesix


def wrapping_dist_python(coor1, coor2, box):
    assert (coor1.ndim == 1) or (coor2.ndim == 1)
    dist = coor1 - coor2
    dist = dist - box * np.round(dist / box)
    return np.sqrt(np.sum(dist * dist, 1))


def _findNeighbours(lipids, box):
    xypos = np.vstack([ll.xyz[:2] for ll in lipids])
    leaflet = np.array([np.sign(ll.xyz[2]) for ll in lipids])

    for i in range(len(lipids)):
        dist = wrapping_dist_python(xypos[i, :], xypos[i + 1 :, :], box)
        same_leaflet = leaflet[i + 1 :] == leaflet[i]
        lipids[i].neighbours = i + 1 + np.where((dist < 11) & same_leaflet)[0]


def _loadMolecules(lipids, files):
    from moleculekit.util import rotationMatrix

    # Create Molecules
    for ll in lipids:
        randidx = np.random.randint(len(files[ll.resname]))
        mol = Molecule(files[ll.resname][randidx])
        mol.filter("not water", _logger=False)
        if ll.xyz[2] < 0:
            mol.rotateBy(
                rotationMatrix([1, 0, 0], np.deg2rad(180))
            )  # Rotate the lower leaflet lipids upside down
        ll.mol = mol
        ll.rot = np.random.random() * 360 - 180  # Random starting rotation


def _locateLipidFiles(folder, lipidnames):
    import os

    files = {}
    for mm in lipidnames:
        files[mm] = glob(os.path.join(folder, mm, "*.cif"))
        if len(files[mm]) == 0:
            raise RuntimeError(
                f'Could not locate cif files for lipid "{mm}" in folder {folder}'
            )
    return files


def _equilibrateOpenMM(
    smemb,
    minimize=0,
    equilibrate_ns=0,
    platform_name="CUDA",
    forcefield_files=None,
    temperature=300,
):
    """Minimize and/or equilibrate ``smemb`` with OpenMM.

    The htmd lipid library uses Lipid17-compatible atom names on merged
    single-residue lipids (POPC, POPE, CHL1, ...), which match the AMBER
    Lipid17 OpenMM XML directly. Water from htmd.solvate uses CHARMM names
    (TIP3 / OH2), so on a working copy we rename TIP3 -> HOH and OH2 -> O so
    the AMBER tip3p XML matches. Equilibrated coordinates and the final box
    are written back into the original Molecule.
    """
    import os
    from openmm import app, unit, Platform
    from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat
    from htmd.util import tempname

    if forcefield_files is None:
        forcefield_files = ["amber14/lipid17.xml", "amber14/tip3p.xml"]

    work = smemb.copy()
    is_water = work.resname == "TIP3"
    work.resname[is_water] = "HOH"
    work.name[is_water & (work.name == "OH2")] = "O"

    # Make sure the box is set so PDB CRYST1 gets written and OpenMM can do PME.
    coord_min = work.coords[:, :, 0].min(axis=0)
    coord_max = work.coords[:, :, 0].max(axis=0)
    box = (coord_max - coord_min).astype(np.float32)
    work.box = box.reshape(3, 1)
    work.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    pdb_path = tempname(suffix=".pdb")
    work.write(pdb_path)
    try:
        pdb = app.PDBFile(pdb_path)
    finally:
        os.remove(pdb_path)

    pdb.topology.setPeriodicBoxVectors(
        [
            (box[0], 0, 0) * unit.angstrom,
            (0, box[1], 0) * unit.angstrom,
            (0, 0, box[2]) * unit.angstrom,
        ]
    )

    ff = app.ForceField(*forcefield_files)
    system = ff.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=10 * unit.angstrom,
        constraints=app.HBonds,
    )

    if equilibrate_ns > 0:
        barostat = MonteCarloMembraneBarostat(
            1 * unit.bar,
            0 * unit.bar * unit.nanometer,
            temperature * unit.kelvin,
            MonteCarloMembraneBarostat.XYIsotropic,
            MonteCarloMembraneBarostat.ZFree,
        )
        system.addForce(barostat)

    integrator = LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        2 * unit.femtoseconds,
    )
    platform = Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    if minimize > 0:
        simulation.minimizeEnergy(maxIterations=int(minimize))

    if equilibrate_ns > 0:
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        nsteps = int(round(equilibrate_ns * 1_000_000 / 2))  # ns -> fs / 2 fs per step
        simulation.step(nsteps)

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    smemb.coords[:, :, 0] = positions.astype(np.float32)

    a, b, c = state.getPeriodicBoxVectors(asNumpy=True)
    smemb.box = np.array(
        [
            [a[0].value_in_unit(unit.angstrom)],
            [b[1].value_in_unit(unit.angstrom)],
            [c[2].value_in_unit(unit.angstrom)],
        ],
        dtype=np.float32,
    )


def buildMembrane(
    xysize,
    ratioupper,
    ratiolower,
    waterbuff=20,
    platform="CUDA",
    minimize=0,
    equilibrate=0,
    outdir=None,
    lipidf=None,
    forcefield_files=None,
    seed=None,
):
    """Construct a membrane containing arbitrary lipids and ratios of them.

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
    platform : str
        OpenMM platform on which to run the minimization/equilibration
        ('CUDA', 'OpenCL', 'CPU', or 'Reference')
    minimize : int
        If not 0 it minimizes the membrane for the given number of steps
    equilibrate : float
        If not 0 it equilibrates the membrane for the given number of nanoseconds
    outdir : str
        A folder in which to store the output PDB files
    lipidf : str
        The path to the folder containing the single-lipid PDB structures as well as the lipid DB file
    forcefield_files : list[str] or None
        OpenMM ForceField XML files used to parameterize the membrane during
        minimization/equilibration. Defaults to AMBER Lipid17 + TIP3P
        (``["amber14/lipid17.xml", "amber14/tip3p.xml"]``).
    seed : int or None
        Seed for the numpy global RNG. If provided, the build is reproducible
        (lipid conformer choice, initial rotations, and the LJ-fluid Halton
        shuffle). The OpenMM minimization/dynamics step is not seeded here.

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
    from htmd.util import tempname
    from htmd.home import home
    import os
    import pandas as pd

    if isinstance(equilibrate, bool):
        raise ValueError("equilibrate must be a float")

    if seed is not None:
        np.random.seed(seed)

    if lipidf is None:
        lipidf = os.path.join(home(shareDir=True), "membranebuilder", "lipids")
    lipiddb = pd.read_csv(os.path.join(lipidf, "lipiddb.csv"), index_col="Name")

    uqlip = np.unique(list(ratioupper.keys()) + list(ratiolower.keys()))
    files = _locateLipidFiles(lipidf, uqlip)

    area = np.prod(xysize)
    lipids = _createLipids(ratioupper, area, lipiddb, files, leaflet="upper")
    lipids += _createLipids(ratiolower, area, lipiddb, files, leaflet="lower")

    _setPositionsLJSim(xysize, [ll for ll in lipids if ll.xyz[2] > 0])
    _setPositionsLJSim(xysize, [ll for ll in lipids if ll.xyz[2] < 0])

    _findNeighbours(lipids, xysize)

    _loadMolecules(lipids, files)

    resolveRingPenetrations(lipids, xysize)
    memb = _createMembraneMolecule(lipids)

    head_mask = np.zeros(memb.numAtoms, dtype=bool)
    for resname, headname in {(ll.resname.upper(), ll.headname) for ll in lipids}:
        head_mask |= (memb.resname == resname) & (memb.name == headname)
    head_coords = memb.coords[head_mask, :, 0]
    minc = head_coords.min(axis=0) - 5
    maxc = head_coords.max(axis=0) + 5

    mm = [[0, 0, maxc[2] - 2], [xysize[0], xysize[1], maxc[2] + waterbuff]]
    smemb = solvate(memb, minmax=mm)
    mm = [[0, 0, minc[2] - waterbuff], [xysize[0], xysize[1], minc[2] + 2]]
    smemb = solvate(smemb, minmax=mm)

    smemb.moveBy([0, 0, -smemb.coords[:, 2, 0].min()])

    if outdir is None:
        outdir = tempname()
        logger.info(f"Outdir {outdir}")
    os.makedirs(outdir, exist_ok=True)

    if equilibrate > 0 or minimize > 0:
        smemb.write(os.path.join(outdir, "starting_structure.pdb"))
        _equilibrateOpenMM(
            smemb,
            minimize=minimize,
            equilibrate_ns=equilibrate,
            platform_name=platform,
            forcefield_files=forcefield_files,
        )

    smemb.write(os.path.join(outdir, "structure.pdb"))
    return smemb


def _findLeastAreaLipid(folder):
    """
    Use this to select the single least stretched conformation from a CHARMM-GUI lipid library
    """
    from glob import glob
    from scipy.spatial.distance import cdist

    ff = glob(folder + "/*/*.crd")
    maxdist = []
    for f in ff:
        m = Molecule(f)
        center = m.coords.mean(axis=0)
        dists = cdist(m.coords[:, :2, 0], center[:2].T)
        maxdist.append(dists.max())
    return ff[np.argmin(maxdist)], np.min(maxdist)
