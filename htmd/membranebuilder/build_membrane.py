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


def _solute_footprint(solute, head_z, slab=5.0, offset=4.0, buffer=5.0):
    """XY footprint of solute heavy atoms in a slab inside the leaflet.

    The slab is one-sided (extends toward the bilayer center), starting
    ``offset`` Angstrom inward from the head plane and continuing for
    ``slab`` Angstrom. The semantics is "TM structure that lipids must
    pack around inside this leaflet", not "anything near the head plane":

    - Atoms above the upper head plane (or below the lower head plane)
      are excluded, so extramembrane domains (cytoplasmic tails,
      extracellular loops, helices on top of the membrane, ...) do not
      contribute.
    - With ``offset`` ~ 4 A the slab also misses amphipathic helices
      that lie flat on the membrane surface. Lipid heads still sit
      under such helices, so they should not displace lipid placement.

    For the upper leaflet (``head_z > 0``): atoms in
    ``z in [head_z - offset - slab, head_z - offset]``. For the lower
    leaflet (``head_z < 0``): atoms in
    ``z in [head_z + offset, head_z + offset + slab]``.

    The returned per-atom radii are ``vdw_radius + buffer``. The buffer
    accounts for the finite size of the lipid head: with the Halton seed
    representing the head center, a seed must be at least ~one head-radius
    away from any solute atom to keep the head atoms outside protein pores
    or near surfaces.
    """
    from moleculekit.periodictable import periodictable

    coords = solute.coords[:, :, 0]
    is_heavy = solute.element != "H"
    z = coords[:, 2]
    if head_z >= 0:
        in_slab = (z >= head_z - offset - slab) & (z <= head_z - offset)
    else:
        in_slab = (z >= head_z + offset) & (z <= head_z + offset + slab)
    mask = is_heavy & in_slab
    if not mask.any():
        return None
    xy = coords[mask, :2]
    vdw = np.array(
        [
            periodictable[el].vdw_radius if el in periodictable else 1.7
            for el in solute.element[mask]
        ]
    )
    return xy, vdw + buffer


def _solute_area_fraction(footprint, xysize, n_samples=10000):
    """Monte Carlo estimate of the box-area fraction occupied by the footprint.

    Samples the centered XY box ``[-Lx/2, Lx/2] x [-Ly/2, Ly/2]`` uniformly and
    returns the fraction of points falling within any per-atom disk.
    """
    if footprint is None:
        return 0.0
    xy, radii = footprint
    samples = np.random.uniform(
        low=[-xysize[0] / 2, -xysize[1] / 2],
        high=[xysize[0] / 2, xysize[1] / 2],
        size=(n_samples, 2),
    )
    diffs = samples[:, None, :] - xy[None, :, :]
    dists2 = np.sum(diffs * diffs, axis=2)
    inside = (dists2 < (radii * radii)[None, :]).any(axis=1)
    return float(inside.mean())


def _createLipids(
    lipidratio, area, lipiddb, files, leaflet=None,
    area_fraction_used=0.0, head_z=15.0,
):
    if leaflet not in ("upper", "lower"):
        raise ValueError(
            f"leaflet must be 'upper' or 'lower', got {leaflet!r}"
        )
    if area_fraction_used >= 1.0:
        raise RuntimeError(
            f"Solute occupies the entire {leaflet} leaflet area "
            f"(fraction={area_fraction_used:.3f}); cannot place any lipids."
        )

    lipiddb = lipiddb.to_dict(orient="index")
    lipidnames = list(lipidratio.keys())
    ratiosAPL = np.array(
        [lipidratio[lipn] * lipiddb[lipn]["APL"] for lipn in lipidnames]
    )
    # Calculate the total areas per lipid type, scaled by the available area
    available_area = area * (1.0 - area_fraction_used)
    areaspl = available_area * (ratiosAPL / ratiosAPL.sum())
    # Calculate the counts from the total areas
    counts = np.round(
        areaspl / np.array([lipiddb[lipn]["APL"] for lipn in lipidnames])
    ).astype(int)

    if area_fraction_used > 0:
        logger.info(
            f"{leaflet} leaflet: solute occupies {area_fraction_used:.1%} of XY "
            f"area; placing {counts.sum()} lipids "
            f"({dict(zip(lipidnames, counts.tolist()))})."
        )

    if (counts == 0).any():
        zero_lipids = [lipidnames[i] for i in range(len(counts)) if counts[i] == 0]
        raise RuntimeError(
            f"Computed lipid count is 0 for {zero_lipids} in {leaflet} leaflet "
            f"given the requested xysize and ratios. Increase the membrane size "
            f"or adjust the ratios."
        )

    # All lipids in a leaflet sit at the same head plane (z = +-head_z).
    # The per-lipid Thickness column in lipiddb is the equilibrium
    # head-to-head distance for a *pure* bilayer of that lipid; using it
    # per-lipid in a mixed bilayer gives a stepped initial surface that's
    # worse than just picking a common plane and letting NPT settle each
    # species to its own depth.
    z_sign = 1 if leaflet == "upper" else -1
    lipids = []
    for i in range(len(lipidnames)):
        resname = lipidnames[i]
        rings = _detectRings(Molecule(files[resname][0]))
        for k in range(counts[i]):
            xyz = np.array([np.nan, np.nan, z_sign * head_z])
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


def _setPositionsLJSim(width, lipids, footprint=None):
    from htmd.membranebuilder.ljfluid import distributeLipids

    # Sigma is chosen so a hex-packed monolayer at the LJ minimum spacing
    # gives exactly the requested area-per-lipid (APL). For nearest-neighbor
    # distance a in a hex lattice the cell area is a^2 * sqrt(3)/2 = APL, so
    # a = sqrt(2*APL/sqrt(3)). The LJ minimum sits at a = 2^(1/6) * sigma,
    # giving sigma = sqrt(2*APL/sqrt(3)) / 2^(1/6). Setting sigma as the
    # disk diameter for a circle of area APL (the previous formula) is too
    # large by ~18% and makes the LJ sim try to space lipids ~40% further
    # apart than APL allows, which then squeezes them onto protein obstacles.
    sigmas = np.array(
        [np.sqrt(2 * ll.area / np.sqrt(3)) / (2 ** (1 / 6)) for ll in lipids]
    )
    resnames = [ll.resname for ll in lipids]

    cutoff = min(np.min(width) / 2, 3 * np.max(sigmas))
    forbidden_xy = footprint[0] if footprint is not None else None
    forbidden_radii = footprint[1] if footprint is not None else None
    pos, pos_initial = distributeLipids(
        width + [2 * cutoff],
        resnames,
        sigmas,
        cutoff=cutoff,
        forbidden_xy=forbidden_xy,
        forbidden_radii=forbidden_radii,
    )
    for i in range(len(lipids)):
        lipids[i].xyz[:2] = pos[i, :2]
        lipids[i].xyz_initial = pos_initial[i, :2].copy()


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


def _optimizeLipidRotations(lipids, solute, n_angles=24, search_radius=15.0):
    """Replace per-lipid random z-rotations with the rotation that maximizes
    the minimum distance to any solute heavy atom.

    Only lipids whose head is within ``search_radius`` of the solute (in 3D)
    are optimized; the rest keep their initial random rotation. The lipid is
    rotated around the z-axis in place (the same operation
    :func:`_createMembraneMolecule` later applies based on ``ll.rot``).
    """
    from scipy.spatial import cKDTree
    from moleculekit.util import rotationMatrix

    heavy = solute.element != "H"
    sxyz = solute.coords[heavy, :, 0]
    if len(sxyz) == 0:
        return
    tree = cKDTree(sxyz)

    angles = np.linspace(-180, 180, n_angles, endpoint=False)
    R = np.array(
        [rotationMatrix([0, 0, 1], np.deg2rad(a)) for a in angles]
    )  # (n_angles, 3, 3)

    n_optimized = 0
    for ll in lipids:
        head_d, _ = tree.query(ll.xyz[None, :], k=1)
        if head_d[0] > search_radius:
            continue

        mol = ll.mol
        lipid_heavy = mol.element != "H"
        head = mol.coords[mol.name == ll.headname].flatten()
        local = mol.coords[lipid_heavy, :, 0] - head  # (N, 3)

        # Rotate all-at-once: (n_angles, N, 3) = einsum(local @ R.T)
        rotated = np.einsum("nij,kj->nki", R, local)  # (n_angles, N, 3)
        placed = rotated + ll.xyz  # (n_angles, N, 3)

        # KDTree only queries 2D; flatten to (n_angles*N, 3) and reshape.
        flat = placed.reshape(-1, 3)
        d, _ = tree.query(flat, k=1, distance_upper_bound=search_radius)
        d = d.reshape(len(angles), -1)
        per_angle_min = d.min(axis=1)

        best = int(np.argmax(per_angle_min))
        ll.rot = float(angles[best])
        n_optimized += 1

    if n_optimized:
        logger.info(
            f"Optimized rotation for {n_optimized} lipids near solute "
            f"(within {search_radius:.1f} A)."
        )


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


def _writeLJDebugPDB(
    path, lipids, upper_fp, lower_fp, com_xy, upper_z, lower_z,
    use_initial=False,
):
    """Write the lipid head positions and per-leaflet obstacle positions to a
    single PDB for visual inspection.

    With ``use_initial=False`` the post-LJ-sim head positions
    (``ll.xyz``) are written; with ``use_initial=True`` the pre-sim
    Halton-filter positions (``ll.xyz_initial``) are written. Lipid heads
    use the lipid's own resname/headname; obstacle pseudo-atoms use ``Au``
    so they're easy to pick out in a viewer.
    """
    n_lipids = len(lipids)
    n_upper = 0 if upper_fp is None else len(upper_fp[0])
    n_lower = 0 if lower_fp is None else len(lower_fp[0])
    n_total = n_lipids + n_upper + n_lower

    mol = Molecule().empty(n_total)
    mol.coords = np.zeros((n_total, 3, 1), dtype=np.float32)

    for i, ll in enumerate(lipids):
        if use_initial:
            mol.coords[i, :2, 0] = ll.xyz_initial
            mol.coords[i, 2, 0] = ll.xyz[2]
        else:
            mol.coords[i, :, 0] = ll.xyz
        mol.name[i] = ll.headname
        mol.resname[i] = ll.resname.upper()
        mol.resid[i] = i
        mol.element[i] = ll.headname[0]
        mol.segid[i] = "M"

    idx = n_lipids
    for fp, z, segid in [(upper_fp, upper_z, "U"), (lower_fp, lower_z, "L")]:
        if fp is None:
            continue
        xy, _ = fp
        for j in range(len(xy)):
            mol.coords[idx, :2, 0] = xy[j] + com_xy
            mol.coords[idx, 2, 0] = z if z is not None else 0.0
            mol.name[idx] = "X"
            mol.resname[idx] = "OBS"
            mol.resid[idx] = idx
            mol.element[idx] = "Au"
            mol.segid[idx] = segid
            idx += 1

    mol.write(path)


def _equilibrateOpenMM(
    smemb,
    minimize=0,
    equilibrate_ns=0,
    platform_name="CUDA",
    forcefield_files=None,
    temperature=300,
    timestep_fs=2.0,
    solute=None,
    head_anchors=None,
    head_restraint_k=0.0,
):
    """Minimize and/or equilibrate ``smemb`` with OpenMM.

    The htmd lipid library uses Lipid17-compatible atom names on merged
    single-residue lipids (POPC, POPE, CHL1, ...), which match the AMBER
    Lipid17 OpenMM XML directly. Water from htmd.solvate uses CHARMM names
    (TIP3 / OH2), so on a working copy we rename TIP3 -> HOH and OH2 -> O so
    the AMBER tip3p XML matches. Equilibrated coordinates and the final box
    are written back into the original Molecule.

    When ``solute`` is provided, its heavy atoms are appended to the OpenMM
    System as frozen ``mass=0`` particles in the standard ``NonbondedForce``
    so the lipids relax around them. The ghost atoms are discarded before
    writing equilibrated coordinates back, so they never appear in the
    output Molecule.
    """
    import os
    import openmm
    from openmm import app, unit, Platform
    from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat
    from moleculekit.periodictable import periodictable

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
        flexibleConstraints=True,
    )

    n_lipid_atoms = system.getNumParticles()
    ghost_xyz = None
    if solute is not None:
        heavy = solute.element != "H"
        ghost_xyz = solute.coords[heavy, :, 0]
        ghost_elements = solute.element[heavy]
        # Find the standard NonbondedForce so we can add ghost particles.
        nb = next(
            f for f in (system.getForce(i) for i in range(system.getNumForces()))
            if isinstance(f, openmm.NonbondedForce)
        )
        for el in ghost_elements:
            # vdW radius is r_min/2; AMBER sigma = r_min / 2^(1/6).
            r_vdw = (
                periodictable[el].vdw_radius if el in periodictable else 1.7
            )
            sigma = 2.0 * r_vdw / (2 ** (1 / 6))
            system.addParticle(0.0)  # mass=0 -> frozen
            nb.addParticle(
                0.0 * unit.elementary_charge,
                sigma * 0.1 * unit.nanometer,
                0.4 * unit.kilojoule_per_mole,  # typical LJ strength
            )

    if head_anchors and head_restraint_k > 0:
        # Harmonic z-restraint on each lipid head atom toward its target
        # head plane (+head_z for upper, -head_z for lower). Lets tails
        # relax around the protein during minimization without dragging
        # heads off their chosen plane.
        anchor = openmm.CustomExternalForce("k_h * (z - z0)^2")
        anchor.addGlobalParameter("k_h", float(head_restraint_k) * 100.0)
        anchor.addPerParticleParameter("z0")
        for atom_idx, z_target in head_anchors:
            anchor.addParticle(int(atom_idx), [float(z_target) * 0.1])
        system.addForce(anchor)

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
        timestep_fs * unit.femtoseconds,
    )
    platform = Platform.getPlatformByName(platform_name)
    # Topology atom count must match System particle count; extend topology
    # with a dummy chain for the ghost atoms so OpenMM accepts the positions.
    if ghost_xyz is not None:
        ghost_chain = pdb.topology.addChain()
        carbon = app.Element.getBySymbol("C")
        for _ in range(len(ghost_xyz)):
            res = pdb.topology.addResidue("GHOST", ghost_chain)
            pdb.topology.addAtom("X", carbon, res)
        all_positions = np.vstack(
            [np.asarray(pdb.positions.value_in_unit(unit.angstrom)), ghost_xyz]
        ) * unit.angstrom
    else:
        all_positions = pdb.positions

    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(all_positions)

    if minimize > 0:
        from acemd.minimizer import minimize as _acemd_cg_minimize

        kcal = unit.kilocalorie_per_mole

        def _e():
            return simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal)

        e_initial = _e()
        # CG handles catastrophic starting energies (severe clashes that make
        # L-BFGS plateau); then L-BFGS to convergence cleans up the residual
        # large forces that would otherwise blow up MD on the first timestep.
        _acemd_cg_minimize(system, simulation.context, int(minimize))
        e_after_cg = _e()
        simulation.minimizeEnergy()
        e_final = _e()
        logger.info(
            f"Minimization: {e_initial:.4g} -> {e_after_cg:.4g} "
            f"(CG, {int(minimize)} steps) -> {e_final:.4g} "
            f"(L-BFGS to convergence) [kcal/mol]"
        )

    if equilibrate_ns > 0:
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        nsteps = int(round(equilibrate_ns * 1_000_000 / timestep_fs))
        simulation.step(nsteps)

    # enforcePeriodicBox=False keeps molecules intact across the PBC; we
    # wrap with moleculekit afterwards which knows about the bonds set on
    # smemb and won't split lipid molecules across the box edge.
    state = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=False
    )
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    smemb.coords[:, :, 0] = positions[:n_lipid_atoms].astype(np.float32)

    a, b, c = state.getPeriodicBoxVectors(asNumpy=True)
    smemb.box = np.array(
        [
            [a[0].value_in_unit(unit.angstrom)],
            [b[1].value_in_unit(unit.angstrom)],
            [c[2].value_in_unit(unit.angstrom)],
        ],
        dtype=np.float32,
    )
    smemb.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    # Wrap. With a solute, center wrapping on the (equilibrated) solute COM
    # so lipids/water cluster around the protein in the periodic image
    # rather than the box origin.
    if solute is not None:
        solute_com = positions[n_lipid_atoms:].mean(axis=0)
        smemb.wrap(wrapcenter=solute_com)
    else:
        smemb.wrap()


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
    solute=None,
    timestep_fs=2.0,
    head_z=15.0,
    head_restraint_k=0.0,
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
    solute : :class:`Molecule <moleculekit.molecule.Molecule>` or None
        Optional pre-positioned solute (typically a protein) around which the
        membrane is built. Coordinates must already be in the membrane frame:
        XY centered on the origin (i.e. spanning ``[-Lx/2, Lx/2] x [-Ly/2,
        Ly/2]``) and bilayer center at z=0. The XY footprint of the solute in
        each leaflet's head plane is used to reduce per-leaflet lipid counts
        proportionally so the resulting membrane has the correct area-per-lipid.
    timestep_fs : float
        Integrator timestep in femtoseconds for the OpenMM equilibration.
        Default 2.0 (compatible with ``constraints=HBonds``).

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
    upper_fraction = 0.0
    lower_fraction = 0.0
    upper_fp = None
    lower_fp = None
    com_xy = np.zeros(2, dtype=np.float32)
    if solute is not None:
        # The solute defines where the membrane should sit in XY: take the
        # COM of its membrane-embedded heavy atoms (|z| < mean_thickness/2),
        # falling back to the full COM for peripheral solutes lying on top
        # of a bilayer. The membrane below is built in the centered frame
        # [-Lx/2, Lx/2] and only the footprint xy is translated; the user's
        # Molecule is never modified.
        mean_thickness = float(
            np.mean([lipiddb.loc[name, "Thickness"] for name in uqlip])
        )
        z_solute = solute.coords[:, 2, 0]
        embedded = np.abs(z_solute) < mean_thickness / 2
        anchor_mask = embedded if embedded.any() else np.ones(solute.numAtoms, bool)
        com_xy = solute.coords[anchor_mask, :2, 0].mean(axis=0).astype(np.float32)

        upper_fp = _solute_footprint(solute, head_z)
        lower_fp = _solute_footprint(solute, -head_z)
        # Translate footprint xy from the user's frame to the centered LJ
        # frame so Halton/obstacles see the protein at the box origin.
        if upper_fp is not None:
            upper_fp = (upper_fp[0] - com_xy, upper_fp[1])
        if lower_fp is not None:
            lower_fp = (lower_fp[0] - com_xy, lower_fp[1])
        upper_fraction = _solute_area_fraction(upper_fp, xysize)
        lower_fraction = _solute_area_fraction(lower_fp, xysize)

    lipids = _createLipids(
        ratioupper, area, lipiddb, files, leaflet="upper",
        area_fraction_used=upper_fraction, head_z=head_z,
    )
    lipids += _createLipids(
        ratiolower, area, lipiddb, files, leaflet="lower",
        area_fraction_used=lower_fraction, head_z=head_z,
    )

    _setPositionsLJSim(xysize, [ll for ll in lipids if ll.xyz[2] > 0], footprint=upper_fp)
    _setPositionsLJSim(xysize, [ll for ll in lipids if ll.xyz[2] < 0], footprint=lower_fp)

    # Translate lipid xy from the centered LJ frame to the user's frame
    # (shifting by the solute's COM) so everything downstream - rotation
    # optimization, ring penetration, membrane assembly, solvation - already
    # lives in the user's frame and the user's solute can be appended as is.
    if solute is not None:
        for ll in lipids:
            ll.xyz[:2] += com_xy
            ll.xyz_initial += com_xy

    _findNeighbours(lipids, xysize)

    _loadMolecules(lipids, files)

    if solute is not None:
        _optimizeLipidRotations(lipids, solute)

    resolveRingPenetrations(lipids, xysize)
    memb = _createMembraneMolecule(lipids)

    head_mask = np.zeros(memb.numAtoms, dtype=bool)
    for resname, headname in {(ll.resname.upper(), ll.headname) for ll in lipids}:
        head_mask |= (memb.resname == resname) & (memb.name == headname)
    head_coords = memb.coords[head_mask, :, 0]
    minc = head_coords.min(axis=0) - 5
    maxc = head_coords.max(axis=0) + 5

    mm = [
        [minc[0] - 5, minc[1] - 5, maxc[2] - 2],
        [maxc[0] + 5, maxc[1] + 5, maxc[2] + waterbuff],
    ]
    smemb = solvate(memb, minmax=mm)
    mm = [
        [minc[0] - 5, minc[1] - 5, minc[2] - waterbuff],
        [maxc[0] + 5, maxc[1] + 5, minc[2] + 2],
    ]
    smemb = solvate(smemb, minmax=mm)

    if outdir is None:
        outdir = tempname()
        logger.info(f"Outdir {outdir}")
    os.makedirs(outdir, exist_ok=True)

    _writeLJDebugPDB(
        os.path.join(outdir, "lj_packing_initial.pdb"),
        lipids,
        upper_fp,
        lower_fp,
        com_xy,
        head_z if solute is not None else None,
        -head_z if solute is not None else None,
        use_initial=True,
    )
    _writeLJDebugPDB(
        os.path.join(outdir, "lj_packing.pdb"),
        lipids,
        upper_fp,
        lower_fp,
        com_xy,
        head_z if solute is not None else None,
        -head_z if solute is not None else None,
    )

    if equilibrate > 0 or minimize > 0:
        smemb.write(os.path.join(outdir, "starting_structure.pdb"))

        head_anchors = None
        if head_restraint_k > 0:
            head_anchors = []
            for i, ll in enumerate(lipids):
                m = (smemb.resid == i) & (smemb.name == ll.headname)
                idx = np.where(m)[0]
                if len(idx) != 1:
                    continue
                z_target = head_z if ll.xyz[2] > 0 else -head_z
                head_anchors.append((int(idx[0]), float(z_target)))

        _equilibrateOpenMM(
            smemb,
            minimize=minimize,
            equilibrate_ns=equilibrate,
            platform_name=platform,
            forcefield_files=forcefield_files,
            timestep_fs=timestep_fs,
            solute=solute,
            head_anchors=head_anchors,
            head_restraint_k=head_restraint_k,
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
