# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

from moleculekit.util import sequenceID
import numpy as np
import logging

if TYPE_CHECKING:
    from moleculekit.molecule import Molecule

logger = logging.getLogger(__name__)


class BuildError(Exception):
    def __init__(self, text, errors=()):
        self.text = text
        self.errors = errors

    def __str__(self):
        if isinstance(self.text, str):
            return repr(self.text)
        elif isinstance(self.text, list):
            return "\n".join([v if isinstance(v, str) else repr(v) for v in self.text])


class _MissingErrorType(Exception):
    def __init__(self, text, values=()):
        self.text = text
        self.values = values
        if isinstance(self.values, np.ndarray):
            self.values = self.values.tolist()

    def __str__(self):
        return self.text


class MixedSegmentError(Exception):
    pass


class ResidueInsertionError(Exception):
    pass


class MissingResidueError(_MissingErrorType):
    pass


class MissingParameterError(_MissingErrorType):
    pass


class MissingTorsionError(_MissingErrorType):
    pass


class MissingBondError(_MissingErrorType):
    pass


class MissingAngleError(_MissingErrorType):
    pass


class MissingAtomTypeError(_MissingErrorType):
    pass


def _has_cell(mol: "Molecule") -> bool:
    """Whether a Molecule carries a usable periodic cell.

    A fresh Molecule has ``box`` of shape ``(3, 0)``, and a Molecule read from
    a file with no cell can have an all-zero box, so both size and content are
    checked.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The Molecule to inspect.

    Returns
    -------
    has_cell : bool
        True if `mol` has three positive box lengths.
    """
    return (
        mol.box is not None and mol.box.size >= 3 and bool(np.all(mol.box[:3, 0] > 0))
    )


def _cell_angles(mol: "Molecule") -> list:
    """A Molecule's cell angles, defaulting to 90 degrees when it has none.

    ``boxangles`` cannot be tested for presence by size alone. A reader that
    supplies a box but no angles leaves it as ``np.zeros((3, 1))``
    (``moleculekit/readers.py:470``), whose size is 3, and emitting those zeros
    as angles would describe a degenerate zero-volume cell.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The Molecule to inspect.

    Returns
    -------
    angles : list of float
        ``[alpha, beta, gamma]`` in degrees, or ``[90.0, 90.0, 90.0]`` when
        `mol` carries no usable angles.
    """
    if (
        mol.boxangles is not None
        and mol.boxangles.size >= 3
        and bool(np.any(mol.boxangles[:3, 0]))
    ):
        return [float(v) for v in mol.boxangles[:3, 0]]
    return [90.0, 90.0, 90.0]


def _write_topology_input(mol: "Molecule", path: str, **kwargs) -> None:
    """Write `mol` to `path` with no cell, for tleap/psfgen input.

    Both tleap's ``loadpdb`` and psfgen discard CRYST1 on read anyway, and a
    non-rectangular cell in the intermediate PDB is more likely to be
    misread than ignored, so it is dropped before writing.

    `mol`'s ``box``/``boxangles`` are saved and restored around the write
    rather than working on a ``mol.copy()``: a full Molecule copy duplicates
    every atom field just to protect two small arrays, and `mol` here is
    sometimes the very same object called once per segment in a loop whose
    caller reads it again right after (see ``charmm._write_segments``).

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The Molecule to write. Its box/boxangles are unchanged on return.
    path : str
        The output file path.
    **kwargs
        Forwarded to :meth:`Molecule.write <moleculekit.molecule.Molecule.write>`.
    """
    box, boxangles = mol.box, mol.boxangles
    mol.box = np.zeros((3, 1), dtype=np.float32)
    mol.boxangles = np.zeros((3, 1), dtype=np.float32)
    try:
        mol.write(path, **kwargs)
    finally:
        mol.box, mol.boxangles = box, boxangles


def embed(mol1: "Molecule", mol2: "Molecule", gap: float = 1.3) -> "Molecule":
    """Embed one molecule into another, removing overlapping residues.

    Removes residues of mol2 that have collisions with atoms of mol1, then
    appends mol1 into mol2.

    Parameters
    ----------
    mol1 : :class:`Molecule <moleculekit.molecule.Molecule>`
        The first Molecule object (embedded into mol2).
    mol2 : :class:`Molecule <moleculekit.molecule.Molecule>`
        The second Molecule object (residues clashing with mol1 are removed).
    gap : float, optional
        Minimum distance in Angstroms between atoms of the two molecules below
        which a residue is considered to clash.

    Returns
    -------
    newmol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The resulting Molecule object with mol1 embedded into mol2.

    Examples
    --------
    >>> merged = embed(memb, prot)
    """
    from scipy.spatial.distance import cdist

    dists = cdist(mol2.coords[:, :, 0], mol1.coords[:, :, 0])
    s2close = np.where(dists < gap)[0]

    mol2res = sequenceID(mol2.resid)
    s2closeres = np.isin(mol2res, np.unique(mol2res[s2close]))

    mol2 = mol2.copy()
    mol2.remove(s2closeres, _logger=False)
    mol2.append(mol1)
    return mol2


def convertDisulfide(mol: "Molecule", disu: list) -> list:
    """Convert disulfide selection-string pairs to UniqueResidueID pairs.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule the selections refer to.
    disu : list
        A list of pairs of atom selection strings, each pair identifying the two
        residues forming a disulfide bond.

    Returns
    -------
    newdisu : list
        The same pairs with each selection string resolved to a
        :class:`UniqueResidueID <moleculekit.molecule.UniqueResidueID>`.
    """
    from moleculekit.molecule import UniqueResidueID

    newdisu = []
    for d in disu:
        if not isinstance(d[0], str) or not isinstance(d[1], str):
            raise RuntimeError("All disulfide selections should be strings")
        newdisu.append(
            [
                UniqueResidueID.fromMolecule(mol, d[0]),
                UniqueResidueID.fromMolecule(mol, d[1]),
            ]
        )
    return newdisu


def detectDisulfideBonds(mol: "Molecule", thresh: float = 3) -> list:
    """Automatically detect disulfide bonds in a molecule.

    Finds all SG atoms in cysteine-like residues (resnames starting with "CY")
    and returns pairs whose inter-sulfur distance is below `thresh`.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule for which to detect disulfide bonds.
    thresh : float, optional
        Distance threshold in Angstroms below which two sulfur atoms are
        considered bonded.

    Returns
    -------
    disubonds : list
        A list of pairs of
        :class:`UniqueResidueID <moleculekit.molecule.UniqueResidueID>` objects
        representing the detected disulfide bonds, sorted by residue ID.

    Raises
    ------
    RuntimeError
        If segment names are not defined, if multiple SG atoms are found in the
        same residue, or if a sulfur atom has more than one possible bond partner.
    """
    from scipy.spatial.distance import pdist, squareform
    from moleculekit.molecule import UniqueResidueID

    disubonds = []

    # Find all SG atoms belonging to resnames starting with CY
    # 'resname "CY.*" and name SG'
    idx = np.where(
        [(rn[0:2] == "CY") and (n == "SG") for rn, n in zip(mol.resname, mol.name)]
    )[0]
    if len(idx) == 0:
        return disubonds

    if np.any([len(s) == 0 for s in mol.segid[idx]]):
        raise RuntimeError(
            "Cannot detect disulfide bonds without segment names defined."
        )

    residues = [UniqueResidueID.fromMolecule(mol, idx=i) for i in idx]
    for r1 in range(len(residues)):
        for r2 in range(r1 + 1, len(residues)):
            if residues[r1] == residues[r2]:
                raise RuntimeError(
                    f"Multiple SG atoms detected in the same residue {residues[r1]}. "
                    "Can't guess disulfide bridges."
                )

    sd = squareform(pdist(mol.coords[idx, :, mol.frame]))
    sd[np.diag_indices(sd.shape[0])] = thresh + 1  # Set the diagonal over threshold
    close = sd < thresh
    rows, cols = np.where(close)

    numbonds = np.sum(close, axis=0)
    if np.any(numbonds > 1):
        multibonded_idx1 = np.where(numbonds > 1)[0]
        multibonded_indexes = np.where(close[multibonded_idx1])
        multibonded_idx1 = multibonded_idx1[multibonded_indexes[0]]
        multibonded_idx2 = multibonded_indexes[1]
        pairs = [
            (str(residues[r]), str(residues[c]))
            for r, c in zip(multibonded_idx1, multibonded_idx2)
        ]
        raise RuntimeError(
            f"Sulphur atoms between pairs {pairs} have multiple possible bonds. Cannot guess disulfide bonds. "
            "Please specify them manually."
        )

    uniquerowcols = list(set([tuple(sorted((r, c))) for r, c in zip(rows, cols)]))
    for rc in uniquerowcols:
        disubonds.append([residues[rc[0]], residues[rc[1]]])
        msg = (
            f"Disulfide Bond between: {residues[rc[0]]}\n"
            f"                   and: {residues[rc[1]]}\n"
        )
        print(msg)

    if len(disubonds) == 1:
        logger.info("One disulfide bond was added")
    else:
        logger.info(f"{len(disubonds)} disulfide bonds were added")
    return sorted(disubonds, key=lambda x: x[0].resid)


def detectCisPeptideBonds(mol: "Molecule", respect_bonds: bool = False) -> None:
    """Detect and warn about cis peptide bonds in a protein.

    Projects the protein omega dihedrals and logs a warning for every frame and
    residue whose omega angle indicates a cis peptide bond.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule to check for cis peptide bonds.
    respect_bonds : bool, optional
        If True, only report cis bonds whose backbone atoms form a connected
        component in the molecule's bond graph.
    """
    from moleculekit.projections.metricdihedral import MetricDihedral, Dihedral
    import networkx as nx

    protsel = mol.atomselect("protein and backbone and name C CA N")
    if np.sum(protsel) < 4:  # Less atoms than dihedral
        return

    dih = Dihedral.proteinDihedrals(mol, sel=protsel, dih=("omega",))

    metr = MetricDihedral(dih=dih, sincos=False)
    data = metr.project(mol)
    mapping = metr.getMapping(mol)

    frames, idxs = np.where(np.abs(data) < 120)
    if respect_bonds:
        bond_graph = nx.Graph()
        bond_graph.add_edges_from(mol.bonds)

    for ii in np.unique(idxs):
        currframes = frames[idxs == ii]
        nframes = len(currframes)
        description = mapping.loc[ii].description
        atomIndexes = [int(x) for x in mapping.loc[ii].atomIndexes]

        if respect_bonds:
            sub_bond_graph = bond_graph.subgraph(atomIndexes)
            if not nx.is_connected(sub_bond_graph):
                continue

        currframes_str = "{}".format(currframes)
        if nframes > 5:
            currframes_str = f"[{currframes[0]} ... {currframes[-1]}]"

        logger.warning(
            f'Found cis peptide bond in {nframes} frames: {currframes_str} in the omega diheral "{description}" with indexes {atomIndexes}'
        )


def _checkMixedSegment(mol):
    prot = mol.atomselect("protein")
    acenme = (mol.resname == "ACE") | (mol.resname == "NME")
    sel1 = prot | acenme  # 'protein or resname ACE NME'
    sel2 = ~prot & ~acenme  # 'not protein and not resname ACE NME'
    segsProt = np.unique(mol.segid[sel1])
    segsNonProt = np.unique(mol.segid[sel2])
    intersection = np.intersect1d(segsProt, segsNonProt)
    if len(intersection) != 0:
        logger.warning(
            f"Segments {intersection} contain both protein and non-protein atoms. "
            "Please assign separate segments to them or the build procedure might fail."
        )


def _checkLongResnames(mol, aliasresidues):
    for resname in np.unique(mol.resname):
        if len(resname) > 4 and resname not in aliasresidues:
            raise RuntimeError(
                "Too long residue names in Molecule. Please give a 4-letter alias to these residues with the aliasresidues option."
            )


def removeLipidsInProtein(prot, memb, lipidsel="lipids"):
    """Calculates the convex hull of the protein. If a lipid lies inside the hull it gets removed.

    This does not work well for lipids crossing out of the hull. If even one atom of the lipid is outside it will
    change the hull and will not get removed. I assume it will get removed by the clashes with the protein though.
    """
    return removeAtomsInHull(prot, memb, "name CA", lipidsel)


def removeAtomsInHull(
    mol1: "Molecule",
    mol2: "Molecule",
    hullsel: str | np.ndarray,
    removesel: str | np.ndarray,
) -> tuple:
    """Calculate the convex hull of an atom selection in mol1 and remove atoms within that hull in mol2.

    Parameters
    ----------
    mol1 : :class:`Molecule <moleculekit.molecule.Molecule>`
        Molecule for which to calculate the convex hull.
    mol2 : :class:`Molecule <moleculekit.molecule.Molecule>`
        Molecule containing atoms to check for hull membership.
    hullsel : str or np.ndarray
        An atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`) for
        atoms in mol1 from which to calculate the convex hull.
    removesel : str or np.ndarray
        An atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`) for
        atoms in mol2 from which to remove those located within the hull.

    Returns
    -------
    newmol2 : :class:`Molecule <moleculekit.molecule.Molecule>`
        mol2 without atoms located within the convex hull.
    numrem : int
        Number of fragments removed.
    """
    # TODO: Look into Morphological Snakes
    from scipy.spatial import ConvexHull

    mol2 = mol2.copy()
    # Convex hull of the protein
    hullcoords = mol1.get("coords", hullsel)
    hull = ConvexHull(hullcoords)

    sequence = sequenceID((mol2.resid, mol2.segid))
    uqres = np.unique(sequence)

    toremove = np.zeros(len(sequence), dtype=bool)
    numlipsrem = 0
    for (
        res
    ) in uqres:  # For each fragment check if it's atoms lie within the convex hull
        atoms = np.where(sequence == res)[0]
        newhull = ConvexHull(np.vstack((hullcoords, mol2.get("coords", sel=atoms))))

        # If the hull didn't change by adding the fragment, it lies within convex hull. Remove it.
        if list(hull.vertices) == list(newhull.vertices):
            toremove[atoms] = True
            numlipsrem += 1

    rematoms = mol2.atomselect(removesel)
    mol2.remove(toremove & rematoms)
    return mol2, numlipsrem


# Distance in Angstrom from a lipid's head atom to where its hydrocarbon
# starts. htmd's lipiddb records head-to-head thickness (phosphate to
# phosphate for the phospholipids: 38.5 for POPC), while the tails that pack
# around a solute occupy only the core between the two headgroup regions.
# 6 A is the offset that turns POPC's 38.5 into its ~26.5 A hydrocarbon core.
HEAD_TO_CORE_OFFSET = 6.0


def lipid_inaccessible_points(
    solute: "Molecule",
    thickness: float,
    probe: float = 2.5,
    resolution: float = 1.0,
    min_volume: float = 1200.0,
) -> np.ndarray:
    """Points inside the bilayer slab that no lipid can reach.

    The free space in the slab is flood-filled inward from its lateral edges
    with a probe the size of an acyl methylene. Whatever the fill cannot reach
    is walled off from the bulk lipid by the solute, so no lipid could ever
    diffuse into it and none should be placed there. Space that is merely
    concave, such as a groove between protomers or a lateral fenestration,
    stays reachable and is left alone. That reachability is what separates a
    sealed pore from a homodimer interface lipids legitimately fill; the
    shape of the space does not.

    The solute must be pre-aligned with the bilayer centered at z=0, which is
    how :func:`get_opm_pdb <moleculekit.opm.get_opm_pdb>` returns it.

    Parameters
    ----------
    solute : :class:`Molecule <moleculekit.molecule.Molecule>`
        The membrane-embedded solute, centered on the bilayer midplane.
    thickness : float
        Full bilayer thickness in Angstrom. Only ``|z| < thickness / 2`` is
        searched, since that is the span lipids occupy. This is the
        hydrocarbon core, not a head-to-head thickness, and it must not reach
        past the solute's hydrophobic belt: a taller slab runs into the
        vestibules above and below the belt, where the fill enters from the
        side and comes back down into the region being tested. Erring thin is
        safe, erring tall is not, and the error is one of omission either way
        (nothing reported sealed rather than something sealed that is not).
    probe : float
        Radius in Angstrom of the probe that must reach a point for a lipid to
        occupy it, i.e. one acyl methylene.
    resolution : float
        Grid spacing in Angstrom.
    min_volume : float
        Sealed regions below this volume in cubic Angstrom are ignored, since
        nothing large enough to matter fits inside them. One POPC displaces
        roughly 1200 A^3.

    Returns
    -------
    points : np.ndarray
        An ``(N, 3)`` array of grid points, one per sealed voxel. Empty when
        the solute seals nothing off, which is the common case.

    Examples
    --------
    >>> from moleculekit.opm import get_opm_pdb
    >>> mol, thickness = get_opm_pdb("5YIL")             # doctest: +SKIP
    >>> pts = lipid_inaccessible_points(mol, thickness)  # doctest: +SKIP
    """
    from moleculekit.periodictable import periodictable
    from scipy.spatial import cKDTree
    from scipy import ndimage

    half = thickness / 2
    coords = solute.coords[:, :, 0]
    # Reach past the slab so atoms just outside it still wall the grid in.
    sel = (solute.element != "H") & (np.abs(coords[:, 2]) < half + probe + 2)
    xyz = coords[sel]
    radii = np.array([periodictable[el].vdw_radius for el in solute.element[sel]])

    # Far enough out that every boundary point clears the largest atom, so
    # the faces the fill seeds from are bulk lipid by construction.
    margin = radii.max() + probe
    # np.arange stops short of its endpoint, so pad it by one step to keep
    # the far faces a full margin out rather than wherever they happen to land.
    axes = [
        np.arange(
            xyz[:, 0].min() - margin, xyz[:, 0].max() + margin + resolution, resolution
        ),
        np.arange(
            xyz[:, 1].min() - margin, xyz[:, 1].max() + margin + resolution, resolution
        ),
        np.arange(-half, half, resolution),
    ]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)

    dist, nearest = cKDTree(xyz).query(grid.reshape(-1, 3), k=1, workers=-1)
    free = (dist >= radii[nearest] + probe).reshape(grid.shape[:3])
    labels = ndimage.label(free)[0]

    # Bulk lipid is whatever the slab's lateral faces open onto. A region that
    # is free but touches none of them is enclosed by the solute.
    bulk = np.unique(
        np.concatenate(
            [
                labels[0].ravel(),
                labels[-1].ravel(),
                labels[:, 0].ravel(),
                labels[:, -1].ravel(),
            ]
        )
    )
    sealed = np.bincount(labels.ravel()) * resolution**3 >= min_volume
    sealed[0] = False  # the occupied label
    sealed[bulk] = False
    return grid[sealed[labels]].astype(np.float32)


def removeHET(prot: "Molecule") -> "Molecule":
    """Remove all HETATM residues from a structure.

    Each unique HETATM residue name is removed, assuming it is a bound ligand or
    other heteroatom group.

    Parameters
    ----------
    prot : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule to clean.

    Returns
    -------
    prot : :class:`Molecule <moleculekit.molecule.Molecule>`
        A copy of the input molecule with all HETATM residues removed.
    """
    prot = prot.copy()
    hetatoms = np.unique(prot.resname[prot.record == "HETATM"])
    for het in hetatoms:
        logger.info(
            "Found resname "
            "{}"
            " in structure. Removed assuming it is a bound ligand.".format(het)
        )
        prot.remove("resname {}".format(het))
    return prot


def tileMembrane(
    memb: "Molecule",
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    buffer: float = 1.5,
) -> "Molecule":
    """Tile a membrane in the X and Y dimensions to reach a specific size.

    Parameters
    ----------
    memb : :class:`Molecule <moleculekit.molecule.Molecule>`
        The membrane to be tiled.
    xmin : float
        Minimum x coordinate.
    ymin : float
        Minimum y coordinate.
    xmax : float
        Maximum x coordinate.
    ymax : float
        Maximum y coordinate.
    buffer : float, optional
        Buffer distance in Angstroms between tiles.

    Returns
    -------
    megamemb : :class:`Molecule <moleculekit.molecule.Molecule>`
        A tiled membrane Molecule covering the specified dimensions.
    """
    from tqdm import tqdm

    memb = memb.copy()
    memb.resid = sequenceID((memb.resid, memb.insertion, memb.chain, memb.segid))

    minmemb = np.min(memb.get("coords", "water"), axis=0).flatten()

    size = np.max(memb.get("coords", "water"), axis=0) - np.min(
        memb.get("coords", "water"), axis=0
    )
    size = size.flatten()
    xreps = int(np.ceil((xmax - xmin) / size[0]))
    yreps = int(np.ceil((ymax - ymin) / size[1]))

    logger.info("Replicating Membrane {}x{}".format(xreps, yreps))

    from moleculekit.molecule import Molecule

    megamemb = Molecule()
    bar = tqdm(total=xreps * yreps, desc="Replicating Membrane")
    k = 0
    for x in range(xreps):
        for y in range(yreps):
            tmpmemb = memb.copy()
            xpos = xmin + x * (size[0] + buffer)
            ypos = ymin + y * (size[1] + buffer)

            tmpmemb.moveBy([-float(minmemb[0]) + xpos, -float(minmemb[1]) + ypos, 0])
            tmpmemb.remove(
                "same resid as (x > {} or y > {})".format(xmax, ymax), _logger=False
            )
            if tmpmemb.numAtoms == 0:
                continue

            tmpmemb.set("segid", "M{}".format(k), sel="not water")
            tmpmemb.set("segid", "MW{}".format(k), sel="water")

            megamemb.append(tmpmemb)
            k += 1
            bar.update(1)
    bar.close()

    # Membranes don't tile perfectly. Need to remove waters that clash with lipids of other tiles
    # Some clashes will still occur between periodic images however
    megamemb.remove("same resid as water and within 1.5 of not water", _logger=False)
    return megamemb


def minimalRotation(prot):
    """Find the rotation around Z that minimizes the X and Y dimensions of the protein to best fit in a box.

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
