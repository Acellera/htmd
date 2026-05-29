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
        return (
            f"<{self.__class__.__module__}.{self.__class__.__name__} "
            f"object at {hex(id(self))} {self.__str__()}>"
        )

    def __str__(self):
        s = ""
        if self.resname is not None:
            s += f"resname: {self.resname} "
        if self.headname is not None:
            s += f"headname: {self.headname} "
        if self.xyz is not None:
            s += f"xyz: {self.xyz} "
        if self.neighbours is not None:
            s += f"neigh: {len(self.neighbours)} "
        if self.mol is not None:
            s += f"mol: {id(self.mol)} "
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
    ``slab`` Angstrom. The slab samples the membrane-embedded part of the
    solute that lipids must pack around, not extramembrane structure:

    - Atoms beyond the head plane (cytoplasmic tails, extracellular
      loops, helices on top of the membrane) are excluded.
    - With ``offset`` ~ 4 A the slab is set inside the hydrophobic core
      and misses amphipathic helices that lie flat on the membrane
      surface; lipid heads still sit under such helices.

    For the upper leaflet (``head_z > 0``): atoms in
    ``z in [head_z - offset - slab, head_z - offset]``. For the lower
    leaflet (``head_z < 0``): atoms in
    ``z in [head_z + offset, head_z + offset + slab]``.

    The returned per-atom radii are ``vdw_radius + buffer``. The buffer
    accounts for the finite size of the lipid head: with the Halton seed
    representing the head center, a seed must be at least one head-radius
    away from any solute atom to keep the head atoms outside protein
    pores or near surfaces.
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

    # All lipids in a leaflet share the same head plane (z = +-head_z) so
    # mixed bilayers start with a flat surface. Per-species depths emerge
    # naturally during NPT equilibration.
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


def _setPositionsLJSim(width, lipids, footprint=None, platform_name=None):
    from htmd.membranebuilder.ljfluid import distributeLipids

    # Sigma is chosen so a hex-packed monolayer at the LJ minimum spacing
    # gives exactly the requested area-per-lipid (APL). For nearest-neighbor
    # distance a in a hex lattice the cell area is a^2 * sqrt(3)/2 = APL,
    # so a = sqrt(2*APL/sqrt(3)). The LJ minimum sits at a = 2^(1/6) * sigma,
    # giving sigma = sqrt(2*APL/sqrt(3)) / 2^(1/6).
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
        platform_name=platform_name,
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


# Common water residue and oxygen atom names across the major force fields.
# Used by _detectWaterNaming to identify water in arbitrary input Molecules
# (e.g. CHARMM-GUI output uses TIP3/OH2, AMBER/PDB uses HOH/O, GROMACS uses
# SOL/OW).
_WATER_RESNAMES = (
    "TIP3", "TIP3P", "TIP4", "TIP4P", "TIP5", "TIP5P",
    "HOH", "WAT", "SOL", "H2O", "T3P", "T4P",
)
_WATER_OXYGEN_NAMES = ("OH2", "O", "OW")


def _detectWaterNaming(mol):
    """Detect the (resname, oxygen_atom_name) used for water in ``mol``.

    Returns ``(resname, oxygen_name)`` or ``(None, None)`` if no water is
    found. When multiple water resnames are present, the most common one
    wins; the oxygen name is the most common atom name within that resname
    whose element is O.
    """
    is_water = np.isin(mol.resname, _WATER_RESNAMES)
    if not is_water.any():
        return None, None
    resnames, counts = np.unique(mol.resname[is_water], return_counts=True)
    water_resname = str(resnames[counts.argmax()])
    sel = (mol.resname == water_resname) & np.isin(mol.name, _WATER_OXYGEN_NAMES)
    if not sel.any():
        sel = (mol.resname == water_resname) & (mol.element == "O")
    if not sel.any():
        return water_resname, None
    names, counts = np.unique(mol.name[sel], return_counts=True)
    return water_resname, str(names[counts.argmax()])


def _defaultHeadAtoms():
    """Default ``{resname: head_atom_name}`` map from the shipped lipiddb.csv."""
    from htmd.home import home
    import os
    import pandas as pd

    csv = os.path.join(
        home(shareDir=True), "membranebuilder", "lipids", "lipiddb.csv"
    )
    db = pd.read_csv(csv, index_col="Name")
    return {name.upper(): str(db.loc[name, "Head"]) for name in db.index}


def _headAnchorsFromMolecule(mol, head_atoms=None, head_z=None, midplane_z=None):
    """Internal: build head-atom z-restraint anchors for a membrane Molecule.

    Returns ``[(atom_idx, z_target_A), ...]``. Upper/lower leaflet is split
    at ``midplane_z`` (median head-atom z by default). Each anchor targets
    either the per-leaflet mean head z (when ``head_z`` is None) or
    ``midplane_z +- |head_z|``.
    """
    if head_atoms is None:
        head_atoms = _defaultHeadAtoms()

    head_mask = np.zeros(mol.numAtoms, dtype=bool)
    for resname, atomname in head_atoms.items():
        head_mask |= (mol.resname == resname.upper()) & (mol.name == atomname)

    if not head_mask.any():
        return []

    head_idx = np.where(head_mask)[0]
    head_z_coords = mol.coords[head_idx, 2, 0]

    if midplane_z is None:
        midplane_z = float(np.median(head_z_coords))
    upper = head_z_coords >= midplane_z

    if head_z is None:
        z_upper = float(head_z_coords[upper].mean()) if upper.any() else midplane_z
        z_lower = float(head_z_coords[~upper].mean()) if (~upper).any() else midplane_z
    else:
        z_upper = midplane_z + abs(float(head_z))
        z_lower = midplane_z - abs(float(head_z))

    anchors = []
    for idx, z in zip(head_idx, head_z_coords):
        target = z_upper if z >= midplane_z else z_lower
        anchors.append((int(idx), float(target)))
    return anchors


def equilibrateMembrane(
    mol,
    minimize=0,
    equilibrate_ns=0,
    platform_name="CUDA",
    forcefield_files=None,
    temperature=300,
    timestep_fs=2.0,
    solute=None,
    head_anchors=None,
    head_restraint_k=0.0,
    head_atoms=None,
    water_resname="HOH",
    water_oxygen_name="O",
):
    """Minimize and/or equilibrate a solvated membrane ``mol`` with OpenMM.

    Works both for membranes built by :func:`buildMembrane` and for
    user-supplied membranes (e.g. from CHARMM-GUI or another tool). The
    function detects the water naming used by ``mol`` and renames it to
    match the chosen force field XMLs. When ``head_restraint_k > 0`` and
    no explicit ``head_anchors`` are given, anchors are auto-derived via
    :func:`_headAnchorsFromMolecule`.

    Equilibrated coordinates and the final box are written back into
    ``mol`` in place.

    Parameters
    ----------
    mol : moleculekit.molecule.Molecule
        Solvated membrane. Modified in place.
    minimize : int
        Conjugate-gradient minimization steps before L-BFGS. 0 to skip.
    equilibrate_ns : float
        Length of NPT equilibration in nanoseconds. 0 to skip.
    platform_name : str
        OpenMM platform name (``CUDA``, ``OpenCL``, ``CPU``, ``Reference``).
    forcefield_files : list of str, optional
        Force field XMLs. Defaults to ``["amber14/lipid17.xml",
        "amber14/tip3p.xml"]``.
    temperature : float
        Temperature in kelvin.
    timestep_fs : float
        Integrator timestep in fs.
    solute : moleculekit.molecule.Molecule, optional
        If given, heavy atoms of its membrane-spanning slab are added to
        the OpenMM System as frozen ``mass=0`` particles so the lipids
        relax around the solute. The solute Molecule itself is not
        modified and ghost atoms are discarded before writing back.
    head_anchors : list of (int, float), optional
        Explicit ``(atom_idx, z_target_A)`` anchors. If ``None`` and
        ``head_restraint_k > 0``, anchors are auto-derived from ``mol``.
    head_restraint_k : float
        Harmonic z-restraint constant (kcal/mol/A^2) on head atoms.
        0 disables the restraint.
    head_atoms : dict, optional
        ``{resname: head_atom_name}`` used when auto-deriving anchors.
        Defaults to the htmd lipid library map.
    water_resname, water_oxygen_name : str
        Target water resname and oxygen atom name expected by the chosen
        force field XML. Defaults match AMBER ``tip3p.xml`` (``HOH``/``O``).
        Whatever water naming is detected in ``mol`` is renamed to these
        on a working copy; ``mol`` itself is not renamed.
    """
    import os
    import openmm
    from openmm import app, unit, Platform
    from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat
    from moleculekit.periodictable import periodictable

    from htmd.util import tempname

    if forcefield_files is None:
        forcefield_files = ["amber14/lipid17.xml", "amber14/tip3p.xml"]

    work = mol.copy()
    detected_resname, detected_oxygen = _detectWaterNaming(work)
    if detected_resname is not None and detected_resname != water_resname:
        is_water = work.resname == detected_resname
        work.resname[is_water] = water_resname
        if detected_oxygen is not None and detected_oxygen != water_oxygen_name:
            work.name[is_water & (work.name == detected_oxygen)] = water_oxygen_name

    # buildMembrane leaves mol.box as all zeros (solvate does not set it),
    # so we fall back to coord extents in that case. For user-supplied
    # membranes that already have a valid box, respect it.
    if (
        work.box is not None
        and work.box.size == 3
        and np.all(work.box > 0)
    ):
        box = work.box[:, 0].astype(np.float32)
    else:
        coord_min = work.coords[:, :, 0].min(axis=0)
        coord_max = work.coords[:, :, 0].max(axis=0)
        box = (coord_max - coord_min).astype(np.float32)
        work.box = box.reshape(3, 1)
    work.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    if head_anchors is None and head_restraint_k > 0:
        head_anchors = _headAnchorsFromMolecule(work, head_atoms=head_atoms)

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
        # Restrict ghost atoms to the membrane-spanning slab so we only
        # equilibrate lipids around the protein's TM region. Extramembrane
        # domains add atoms with no nearby lipids to influence.
        known_lipid_resnames = [
            r.upper() for r in (head_atoms or _defaultHeadAtoms()).keys()
        ]
        is_lipid = np.isin(work.resname, known_lipid_resnames)
        if not is_lipid.any():
            raise ValueError(
                "No lipid residues found in mol; pass head_atoms with the "
                "lipid resnames present in your membrane."
            )
        lipid_z = work.coords[is_lipid, 2, 0]
        # Buffer of a few angstrom so that if the membrane expands a bit
        # during equilibration the lipids still feel the solute and don't
        # spill past its TM edge.
        z_buffer = 5.0
        z_lo = float(lipid_z.min()) - z_buffer
        z_hi = float(lipid_z.max()) + z_buffer

        heavy = solute.element != "H"
        in_slab = (solute.coords[:, 2, 0] >= z_lo) & (solute.coords[:, 2, 0] <= z_hi)
        mask = heavy & in_slab
        ghost_xyz = solute.coords[mask, :, 0]
        ghost_elements = solute.element[mask]
        logger.info(
            f"Ghost atoms: {mask.sum()}/{heavy.sum()} solute heavy atoms in "
            f"lipid Z range [{z_lo:.1f}, {z_hi:.1f}]"
        )
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
            system.addParticle(0.0)  # frozen
            nb.addParticle(
                0.0 * unit.elementary_charge,
                sigma * 0.1 * unit.nanometer,
                0.4 * unit.kilojoule_per_mole,
            )

    if head_anchors and head_restraint_k > 0:
        # Harmonic z-restraint on each lipid head atom toward its target
        # head plane (+head_z for upper, -head_z for lower).
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
        # CG copes with catastrophic starting energies; L-BFGS then drives
        # forces low enough that the first MD timestep is stable.
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

    # enforcePeriodicBox=False keeps molecules intact across the PBC;
    # mol.wrap() afterwards uses the bonds on mol to wrap whole lipids.
    state = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=False
    )
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    mol.coords[:, :, 0] = positions[:n_lipid_atoms].astype(np.float32)

    a, b, c = state.getPeriodicBoxVectors(asNumpy=True)
    mol.box = np.array(
        [
            [a[0].value_in_unit(unit.angstrom)],
            [b[1].value_in_unit(unit.angstrom)],
            [c[2].value_in_unit(unit.angstrom)],
        ],
        dtype=np.float32,
    )
    mol.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    # With a solute, center the wrap on the solute COM so lipids/water
    # cluster around the protein rather than the box origin.
    if solute is not None:
        solute_com = positions[n_lipid_atoms:].mean(axis=0)
        mol.wrap(wrapcenter=solute_com)
    else:
        mol.wrap()


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
        OpenMM platform on which to run the LJ-fluid lipid placement and the
        minimization/equilibration ('CUDA', 'OpenCL', 'CPU', or 'Reference')
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
        membrane is built. The solute must be aligned with bilayer center at
        z=0; the membrane is shifted in XY to follow the solute's
        membrane-embedded COM. The user's Molecule is not modified.
    timestep_fs : float
        Integrator timestep in femtoseconds for the OpenMM equilibration.
        Default 2.0 (compatible with ``constraints=HBonds``).
    head_z : float
        Half-bilayer head-plane Z (Angstrom). Every lipid head is placed at
        ``+head_z`` (upper leaflet) or ``-head_z`` (lower) regardless of
        species, so a mixed bilayer starts with a flat head plane that NPT
        relaxes to per-species depths. Default 15.0.
    head_restraint_k : float
        If > 0, apply a harmonic z-restraint of strength ``head_restraint_k``
        (kJ/mol/A^2) to each lipid head atom toward its initial head plane
        during minimization and equilibration. Default 0.0 (no restraint).

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
        # XY anchor for the membrane: the COM of solute heavy atoms inside
        # the bilayer (|z| < mean_thickness/2). Falls back to the full COM
        # for peripheral solutes that don't span the bilayer.
        mean_thickness = float(
            np.mean([lipiddb.loc[name, "Thickness"] for name in uqlip])
        )
        z_solute = solute.coords[:, 2, 0]
        embedded = np.abs(z_solute) < mean_thickness / 2
        anchor_mask = embedded if embedded.any() else np.ones(solute.numAtoms, bool)
        com_xy = solute.coords[anchor_mask, :2, 0].mean(axis=0).astype(np.float32)

        upper_fp = _solute_footprint(solute, head_z)
        lower_fp = _solute_footprint(solute, -head_z)
        # Lipid placement (Halton/obstacles) happens in the box-centered
        # frame, so translate the footprint by -com_xy.
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

    _setPositionsLJSim(
        xysize,
        [ll for ll in lipids if ll.xyz[2] > 0],
        footprint=upper_fp,
        platform_name=platform,
    )
    _setPositionsLJSim(
        xysize,
        [ll for ll in lipids if ll.xyz[2] < 0],
        footprint=lower_fp,
        platform_name=platform,
    )

    # Move lipids into the solute's frame so everything downstream
    # (rotation optimization, ring penetration, assembly, solvation) lives
    # in the user's coordinate frame.
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

    # Per-lipid templates ship with empty segid, which makes
    # htmd.builder.amber.build reject the structure ("Atoms ... do not have
    # segid defined"). Assign a default lipid segid so downstream builders
    # get a fully-segmented molecule. The waters added by solvate() below
    # already get their own segids (W0, W1, ...).
    memb.segid[memb.segid == ""] = "MEMB"

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

    # Set the PBC cell to the lipid head bbox (minc/maxc already include
    # a +-5 A tail-buffer) in xy and the water-padded lipid extent in z.
    # solvate fills [minc-5, maxc+5] in xy, so the outer ~5 A of water
    # wraps to the opposite edge at startup; the MonteCarloMembraneBarostat
    # (XYIsotropic, ZFree) then converges the cell to the natural APL.
    smemb.box = np.array(
        [
            [maxc[0] - minc[0]],
            [maxc[1] - minc[1]],
            [maxc[2] - minc[2] + 2.0 * waterbuff],
        ],
        dtype=np.float32,
    )
    smemb.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

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

        # Use the lipid placement's design intent (xyz[2] sign) to assign
        # leaflets, not the current head z. This is more robust than the
        # generic helper for membranes fresh out of LJ packing, where a
        # head could occasionally cross the midplane.
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

        equilibrateMembrane(
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
    """Select the single least-stretched conformation from a folder of
    per-conformer ``.crd`` files (one subfolder per lipid)."""
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
