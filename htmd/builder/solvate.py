# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from moleculekit.molecule import Molecule
from moleculekit.unitcell import box_vectors_to_lengths_and_angles
import numpy as np
import logging

logger = logging.getLogger(__name__)


def _segid_gen(prefix, mol, mode="decimal"):
    import string

    segids = np.unique(mol.segid)
    for prefix in [prefix] + list(string.ascii_uppercase + string.digits):
        i = 0
        while True:
            if mode == "decimal":
                segid = f"{prefix}{i:d}"
            elif mode == "hex":
                segid = f"{prefix}{i:X}"
            elif mode == "alphanum":
                segid = "{0}{1:c}{2:c}{3:c}".format(
                    prefix,
                    int(np.floor(np.floor(i / 26) / 26) + 65),
                    int(np.mod(np.floor(i / 26), 26) + 65),
                    int(np.mod(i, 26) + 65),
                )
            if len(segid) > 4:
                break
            i += 1

            if segid not in segids:
                yield segid


_CELL_SHAPES = ("rectangular", "cube", "octahedron", "dodecahedron")


def _cell_vectors(shape: str, width: float) -> np.ndarray:
    """Box vectors for the equilateral representative of a cell shape.

    Each shape is returned with a=b=c and alpha=beta=gamma. That form is what
    the AMBER prmtop's single-angle ``BOX_DIMENSIONS`` field can store without
    loss, and it is accepted by OpenMM as a reduced cell.

    The vectors are built directly rather than through a lengths-and-angles
    conversion. That conversion returns 1.5000000000000004 for the 60 degree
    cell's ``b[0]``, which trips OpenMM's ``2*|b0| <= a0`` reduction check.

    Parameters
    ----------
    shape : str
        One of ``"cube"``, ``"octahedron"`` or ``"dodecahedron"``.
        ``"rectangular"`` is not accepted; it uses the per-axis region path.
    width : float
        Cell edge length in Angstroms. All three edges have this length.

    Returns
    -------
    vectors : np.ndarray
        A ``(3, 3)`` array whose rows are the lattice vectors and columns the
        Cartesian components.

    Raises
    ------
    ValueError
        If `shape` is ``"rectangular"`` or is not a known shape name.
    """
    d = float(width)
    s2, s3, s6 = np.sqrt(2.0), np.sqrt(3.0), np.sqrt(6.0)
    if shape == "cube":
        # 90/90/90, volume d**3
        return np.array([[d, 0.0, 0.0], [0.0, d, 0.0], [0.0, 0.0, d]])
    if shape == "octahedron":
        # Truncated octahedron: 109.4712206 x3, BCC lattice, volume 0.7698 d**3
        return np.array(
            [
                [d, 0.0, 0.0],
                [-d / 3.0, 2.0 * s2 * d / 3.0, 0.0],
                [-d / 3.0, -s2 * d / 3.0, s6 * d / 3.0],
            ]
        )
    if shape == "dodecahedron":
        # Rhombic dodecahedron: 60 x3, FCC lattice, volume 0.7071 d**3
        return np.array(
            [
                [d, 0.0, 0.0],
                [d / 2.0, s3 * d / 2.0, 0.0],
                [d / 2.0, s3 * d / 6.0, s6 * d / 3.0],
            ]
        )
    if shape == "rectangular":
        raise ValueError(
            "shape 'rectangular' has no single cell width. It uses the "
            "per-axis region path and never reaches _cell_vectors."
        )
    raise ValueError(f"Unknown cell shape '{shape}'. Valid: {_CELL_SHAPES}")


def _cell_lengths_and_angles(vectors: np.ndarray) -> tuple[list[float], list[float]]:
    """Convert box vectors to edge lengths and angles.

    Parameters
    ----------
    vectors : np.ndarray
        A ``(3, 3)`` array whose rows are the lattice vectors.

    Returns
    -------
    lengths : list of float
        Edge lengths ``[a, b, c]`` in Angstroms.
    angles : list of float
        Cell angles ``[alpha, beta, gamma]`` in degrees.
    """
    la = box_vectors_to_lengths_and_angles(*vectors)
    return [float(v) for v in la[:3]], [float(v) for v in la[3:]]


def _cell_width(coords: np.ndarray, center: np.ndarray, pad: float) -> float:
    """Equilateral cell edge length that pads a solute by `pad` per side.

    `pad` is per-side, so the image distance is ``2 * pad``, matching GROMACS
    ``-d``. The ``4 * pad`` floor is OpenMM's rule in that convention.

    Parameters
    ----------
    coords : np.ndarray
        Atom coordinates, shape ``(N, 3)``.
    center : np.ndarray
        Cell center, shape ``(3,)``.
    pad : float
        Padding in Angstroms, applied per side.

    Returns
    -------
    width : float
        Cell edge length in Angstroms.
    """
    radius = float(np.linalg.norm(coords - center, axis=1).max())
    return max(2.0 * radius + 2.0 * pad, 4.0 * pad)


def _ws_halfspaces(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Half-spaces bounding the Wigner-Seitz cell of a lattice.

    A point ``d``, measured from the cell center, is inside the cell when
    ``d @ L.T <= h`` holds for every row of `L`. That is the standard
    nearest-lattice-point condition written as half-spaces.

    The 26 non-zero offsets in ``{-1, 0, 1}**3`` are sufficient for the
    lattices this module builds. The BCC cell's faces come from its 8 nearest
    neighbors (+-a, +-b, +-c, +-(a+b+c)) and 6 second neighbors (+-(a+b),
    +-(b+c), +-(a+c)); the FCC cell's from its 12 nearest (+-a, +-b, +-c,
    +-(a-b), +-(b-c), +-(a-c)) and 6 second. All have coefficients in
    ``{-1, 0, 1}``.

    Parameters
    ----------
    vectors : np.ndarray
        A ``(3, 3)`` array whose rows are the lattice vectors.

    Returns
    -------
    normals : np.ndarray
        A ``(26, 3)`` array of lattice vectors to the neighboring cells.
    offsets : np.ndarray
        A ``(26,)`` array of half-space offsets, ``0.5 * |L|**2``.
    """
    import itertools

    offs = np.array(
        [n for n in itertools.product((-1, 0, 1), repeat=3) if any(n)], dtype=float
    )
    normals = offs @ vectors
    return normals, 0.5 * (normals**2).sum(axis=1)


def _ws_bbox(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bounding box of the Wigner-Seitz cell, relative to its center.

    The water tiling loop needs this rather than the primitive
    parallelepiped's bounding box, which is smaller along some axes. At width
    60 the parallelepiped spans z in +-24.49 while the WS cell spans +-36.74,
    so tiling the parallelepiped's box and culling to the WS cell would leave a
    vacuum slab at each z extreme.

    The box is computed exactly, by enumerating the WS cell's vertices as the
    solutions of every triple of the 26 bounding half-spaces that satisfies all
    26. Coincident duplicate vertices are returned for the cube and
    dodecahedron, which does not affect the box.

    Parameters
    ----------
    vectors : np.ndarray
        A ``(3, 3)`` array whose rows are the lattice vectors.

    Returns
    -------
    lo : np.ndarray
        A ``(3,)`` array of minimum offsets from the cell center.
    hi : np.ndarray
        A ``(3,)`` array of maximum offsets from the cell center.
    """
    import itertools

    normals, offsets = _ws_halfspaces(vectors)
    verts = []
    for i, j, k in itertools.combinations(range(len(normals)), 3):
        rows = normals[[i, j, k]]
        if abs(np.linalg.det(rows)) < 1e-9:
            continue
        point = np.linalg.solve(rows, offsets[[i, j, k]])
        if np.all(normals @ point <= offsets + 1e-6):
            verts.append(point)
    verts = np.array(verts)
    return verts.min(axis=0), verts.max(axis=0)


def solvate(
    mol: Molecule,
    pad: float | None = None,
    minmax: list | np.ndarray | None = None,
    centersel: str | np.ndarray | None = None,
    boxsize: float | list | np.ndarray | None = None,
    shape: str = "rectangular",
    exclude_z: tuple | list | None = None,
    negx: float = 0,
    posx: float = 0,
    negy: float = 0,
    posy: float = 0,
    negz: float = 0,
    posz: float = 0,
    buffer: float = 2.4,
    watsize: float = 65.4195,
    prefix: str = "W",
    rotate: bool = False,
    spdb: str | None = None,
) -> Molecule:
    """Solvate a molecular system in a water box.

    Places water molecules around the input molecule by tiling a pre-built
    water box and removing waters that clash with existing atoms or fall outside
    the specified box boundaries.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule to solvate.
    pad : float, optional
        Uniform padding in Angstroms to add around the molecule in all six
        directions. Overrides `negx`, `posx`, `negy`, `posy`, `negz`, `posz`.
    minmax : list or np.ndarray, optional
        Explicit box boundaries as a 2D array of the form
        ``[[minx, miny, minz], [maxx, maxy, maxz]]``. If None, derived from
        the molecule's own coordinates.
    centersel : str or np.ndarray, optional
        An atom selection string, a boolean mask, or an integer index array (see :meth:`Molecule.atomselect <moleculekit.molecule.Molecule.atomselect>`)
        defining the center of the solvation box. The geometric center of the
        selected atoms is used. With `shape` left as ``"rectangular"`` it must
        be combined with `boxsize`; with any other shape it may be combined
        with `pad` instead, and it then only moves the cell center, because
        the width is still measured from that center to the furthest atom of
        the whole molecule, so an off-center selection enlarges the cell
        rather than clipping the molecule.
    boxsize : float or list or np.ndarray, optional
        Dimensions of the solvation box. A single float creates a cubic box;
        a 3-element list ``[sx, sy, sz]`` creates an axis-aligned box. Must be
        combined with `centersel`.
    shape : str, optional
        Unit cell shape. ``"rectangular"`` (the default) uses the per-axis
        region defined by `pad`, `minmax`, `boxsize` or the `negx`-`posz`
        arguments, and reproduces the historical behavior. ``"cube"``,
        ``"octahedron"`` (truncated octahedron) and ``"dodecahedron"``
        (rhombic dodecahedron) are equilateral cells that need a single width,
        taken from `pad` or from a scalar `boxsize`. A non-rectangular cell
        reaches the same minimum image distance with less water: 77.0% of a
        cube's volume for the truncated octahedron, 70.7% for the rhombic
        dodecahedron.
    exclude_z : tuple or list, optional
        A ``(zlo, zhi)`` pair, strictly interpreted: water whose z lies
        strictly between the two values is removed, and water exactly on
        either boundary is kept. Water is not placed between these two z
        values.
        Use it to keep water out of a membrane's hydrophobic slab in a single
        `solvate` call: without it, a call spanning the full z range of a
        bilayer places water inside the tail region, where there is enough
        free volume for a water molecule to sit further than `buffer` from any
        lipid atom.
    negx : float, optional
        Padding in Angstroms in the -x direction.
    posx : float, optional
        Padding in Angstroms in the +x direction.
    negy : float, optional
        Padding in Angstroms in the -y direction.
    posy : float, optional
        Padding in Angstroms in the +y direction.
    negz : float, optional
        Padding in Angstroms in the -z direction.
    posz : float, optional
        Padding in Angstroms in the +z direction.
    buffer : float, optional
        Minimum distance in Angstroms between water molecules and other atoms.
    watsize : float, optional
        Edge length in Angstroms of the pre-built water box tile.
    prefix : str, optional
        Prefix string used for water segment names.
    rotate : bool, optional
        If True, rotate the molecule to minimize box volume (not yet implemented).
    spdb : str, optional
        Path to a custom solvent box file, in any format Molecule can read.
        If None, uses the built-in water box. Bonds are read from the file if
        present and guessed from the coordinates otherwise.

    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        A copy of the input molecule with water molecules added.

    Raises
    ------
    ValueError
        If `shape` is not a known shape name; if a non-rectangular `shape` is
        combined with `minmax`, a 3-element `boxsize`, or any of the per-axis
        padding arguments; if a non-rectangular `shape` is given neither `pad`
        nor `boxsize`; if the resulting cell width is not positive; if
        `centersel` matches no atoms; or if `exclude_z` is not an increasing
        ``(zlo, zhi)`` pair of finite numbers.

    Examples
    --------
    >>> smol = solvate(mol, pad=10)
    >>> smol = solvate(mol, minmax=[[-20, -20, -20],[20, 20, 20]])
    >>> smol = solvate(mol, centersel="protein", boxsize=100)
    >>> smol = solvate(mol, centersel="protein", boxsize=[80, 80, 120])
    >>> smol = solvate(mol, pad=12, shape="octahedron")
    """
    from tqdm import tqdm
    from htmd.home import home

    mol = mol.copy()
    if mol.numFrames > 1:
        logger.warning(
            "Multiple frames in Molecule. Solvate keeps only frame 0 and discards the rest."
        )
        mol.coords = np.atleast_3d(mol.coords[:, :, 0])

    if shape not in _CELL_SHAPES:
        raise ValueError(f"Unknown cell shape '{shape}'. Valid: {_CELL_SHAPES}")

    if exclude_z is not None:
        exclude_z = [float(v) for v in exclude_z]
        # `not a < b` also rejects NaN, which would silently disable filtering.
        if len(exclude_z) != 2 or not exclude_z[0] < exclude_z[1]:
            raise ValueError(
                f"exclude_z must be an increasing (zlo, zhi) pair of finite "
                f"numbers, got {exclude_z}"
            )

    cell = None
    if shape != "rectangular":
        per_axis = any(v != 0 for v in (negx, posx, negy, posy, negz, posz))
        if minmax is not None or per_axis:
            raise ValueError(
                f"shape '{shape}' needs a single cell width and cannot be "
                "combined with minmax or the per-axis padding arguments. "
                "Use pad or a scalar boxsize."
            )
        if boxsize is not None and np.atleast_1d(np.array(boxsize)).size != 1:
            raise ValueError(
                f"shape '{shape}' needs a single cell width, so boxsize must "
                "be a scalar. A 3-element boxsize only applies to "
                "shape='rectangular'."
            )

        if mol.numAtoms > 0:
            coords = mol.coords[:, :, 0]
            if centersel is not None:
                selatoms = mol.atomselect(centersel)
                if not np.any(selatoms):
                    raise ValueError(f"Atom selection '{centersel}' matched no atoms.")
                center = coords[selatoms].mean(axis=0)
            else:
                center = 0.5 * (coords.min(axis=0) + coords.max(axis=0))
        else:
            coords = np.zeros((1, 3))
            center = np.zeros(3)

        if boxsize is not None:
            width = float(np.atleast_1d(np.array(boxsize, dtype=float))[0])
        elif pad is not None:
            width = _cell_width(coords, center, pad)
        else:
            raise ValueError(f"shape '{shape}' needs either pad or boxsize.")

        if not width > 0:
            # Zero width would give nan angles from a 0/0 in the conversion.
            raise ValueError(
                f"shape '{shape}' needs a positive cell width, got {width}. "
                "Check pad and boxsize."
            )

        cell = _cell_vectors(shape, width)
        lengths, angles = _cell_lengths_and_angles(cell)
        logger.info(
            f"Cell shape '{shape}': width {width:.2f} A, "
            f"lengths [{lengths[0]:.2f}, {lengths[1]:.2f}, {lengths[2]:.2f}], "
            f"angles [{angles[0]:.2f}, {angles[1]:.2f}, {angles[2]:.2f}]"
        )
        # Tile the WS bounding box, which is larger than the parallelepiped's
        # on some axes; _outOfBoundaries carves the cell out of it.
        ws_lo, ws_hi = _ws_bbox(cell)
        minmax = np.array([center + ws_lo, center + ws_hi])
        pad = None
        centersel = None
        boxsize = None

    if (centersel is None) != (boxsize is None):
        raise ValueError("centersel and boxsize must both be specified together.")

    if centersel is not None:
        if minmax is not None or pad is not None:
            raise ValueError("centersel/boxsize cannot be combined with minmax or pad.")
        selatoms = mol.atomselect(centersel)
        if not np.any(selatoms):
            raise ValueError(f"Atom selection '{centersel}' matched no atoms.")
        center = mol.get("coords", sel=selatoms).mean(axis=0)
        boxsize = np.atleast_1d(np.array(boxsize, dtype=float))
        if boxsize.shape == (1,):
            boxsize = np.repeat(boxsize, 3)
        elif boxsize.shape != (3,):
            raise ValueError("boxsize must be a scalar or a 3-element iterable.")
        half = boxsize / 2.0
        minmax = np.array([center - half, center + half])
        logger.info(
            f"Box center from selection '{centersel}': [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}], "
            f"box size: [{boxsize[0]:.2f}, {boxsize[1]:.2f}, {boxsize[2]:.2f}]"
        )

    if spdb is None:
        spdb = os.path.join(home(shareDir=True), "solvate", "wat.bcif.gz")

    if os.path.isfile(spdb):
        logger.info("Using water file at: " + spdb)
        water = Molecule(spdb)
    else:
        raise NameError("No solvent file found in " + spdb)

    # Without bonds, Molecule.wrap treats each water atom as its own molecule
    # and splits waters across periodic boundaries.
    if water.bonds.shape[0] == 0:
        logger.info(f"No bonds found in {spdb}. Guessing them from the coordinates.")
        water.guessBonds()

    if pad is not None:
        negx = pad
        posx = pad
        negy = pad
        posy = pad
        negz = pad
        posz = pad

    if rotate:
        raise NameError("Rotation not implemented yet")

    # Calculate min max coordinates from molecule
    if mol.numAtoms > 0:
        minmol = np.min(mol.get("coords"), axis=0)
        maxmol = np.max(mol.get("coords"), axis=0)
    else:
        minmol = [np.inf, np.inf, np.inf]
        maxmol = [-np.inf, -np.inf, -np.inf]

    if minmax is None:
        minc = minmol
        maxc = maxmol
    else:
        if isinstance(minmax, list):
            minmax = np.array(minmax)
        minc = minmax[0, :]
        maxc = minmax[1, :]

    xmin = float(minc[0] - negx)
    xmax = float(maxc[0] + posx)
    ymin = float(minc[1] - negy)
    ymax = float(maxc[1] + posy)
    zmin = float(minc[2] - negz)
    zmax = float(maxc[2] + posz)

    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    nx = int(np.ceil((dx + 2 * buffer) / watsize))
    ny = int(np.ceil((dy + 2 * buffer) / watsize))
    nz = int(np.ceil((dz + 2 * buffer) / watsize))

    # Calculate number of preexisting water segments with given prefix
    if mol.numAtoms > 0:
        preexist = len(np.unique(mol.get("segid", sel=f'segid "{prefix}.*"')))
    else:
        preexist = 0

    numsegs = nx * ny * nz
    logger.info(f"Replicating {numsegs} water segments, {nx} by {ny} by {nz}")

    # Check that we won't run out of segment name characters, and switch to
    # using hexadecimal or alphanumeric naming schemes in cases where decimal
    # numbered segnames won't fit into the field width.
    testsegname = f"{prefix}{numsegs + preexist:d}"
    testsegnamehex = f"{prefix}{numsegs + preexist:X}"
    writemode = "decimal"
    if len(testsegname) > 4 and len(testsegnamehex) <= 4:
        writemode = "hex"
        logger.warning(
            "Warning: decimal naming would overrun segname field. Using hexadecimal segnames instead..."
        )
    elif len(testsegnamehex) > 4:
        writemode = "alphanum"
        logger.warning(
            "Warning: decimal or hex naming would overrun segname field. Using alphanumeric segnames instead..."
        )

    segid_gen = _segid_gen(prefix, mol, writemode)

    minx = minmol[0] - buffer
    miny = minmol[1] - buffer
    minz = minmol[2] - buffer
    maxx = maxmol[0] + buffer
    maxy = maxmol[1] + buffer
    maxz = maxmol[2] + buffer

    bar = tqdm(total=nx * ny * nz, desc="Solvating")
    waterboxes = np.empty(numsegs, dtype=object)
    n = preexist
    w = 0
    for i in range(nx):
        movex = xmin + i * watsize
        movexmax = movex + watsize
        xoverlap = True
        if movex > maxx or movexmax < minx:
            xoverlap = False

        for j in range(ny):
            movey = ymin + j * watsize
            moveymax = movey + watsize
            yoverlap = True
            if movey > maxy or moveymax < miny:
                yoverlap = False

            for k in range(nz):
                movez = zmin + k * watsize
                movezmax = movez + watsize
                zoverlap = True
                if movez > maxz or movezmax < minz:
                    zoverlap = False

                segname = next(segid_gen)

                waterboxes[w] = water.copy()
                waterboxes[w].moveBy([movex, movey, movez])
                waterboxes[w].set("segid", segname)

                mol.append(waterboxes[w])
                watsel = mol.segid == segname

                selover = np.zeros(len(watsel), dtype=bool)
                if (
                    xoverlap and yoverlap and zoverlap
                ):  # Remove water overlapping with other segids
                    selover = _overlapWithOther(mol, segname, buffer)
                # Remove water outside the boundaries
                selout = _outOfBoundaries(
                    mol,
                    segname,
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    zmin,
                    zmax,
                    cell=((cell, center) if cell is not None else None),
                    exclude_z=exclude_z,
                )
                sel = selover | selout

                # mol.write('temp.pdb')
                mol.filter(mol.segid != segname, _logger=False)
                waterboxes[w].filter(np.invert(sel[watsel]), _logger=False)
                # waterboxes[w].write('wat' + str(w) + '.pdb')
                n += 1
                w += 1
                bar.update(1)
    bar.close()

    waters = 0
    for i in range(numsegs):
        waters += waterboxes[i].numAtoms
        if waterboxes[i].numAtoms != 0:
            mol.append(waterboxes[i])

    logger.info(f"{int(waters / 3)} water molecules were added to the system.")

    if cell is not None:
        lengths, angles = _cell_lengths_and_angles(cell)
    else:
        lengths = [xmax - xmin, ymax - ymin, zmax - zmin]
        angles = [90.0, 90.0, 90.0]
    mol.box = np.array(lengths, dtype=np.float32).reshape(3, 1)
    mol.boxangles = np.array(angles, dtype=np.float32).reshape(3, 1)
    mol.crystalinfo = dict(
        zip(
            ("a", "b", "c", "alpha", "beta", "gamma"),
            [float(v) for v in lengths + angles],
        )
    )
    return mol


def _overlapWithOther(mol, segname, buffer):
    # Optimized version of this atomselection:
    # segid {segname} and same resid as (segid {segname} and within {buffer} of not segid {segname})
    from moleculekit.atomselect_utils import within_distance

    segmask = mol.segid == segname
    segidx = np.where(segmask)[0].astype(np.uint32)
    notsegidx = np.where(~segmask)[0].astype(np.uint32)

    contacts = np.zeros(len(segidx), dtype=bool)
    within_distance(
        mol.coords[:, :, 0],
        buffer,
        sel1=segidx,
        sel2=notsegidx,
        sel2_min_coords=mol.coords[notsegidx, :, 0].min(axis=0),
        sel2_max_coords=mol.coords[notsegidx, :, 0].max(axis=0),
        results=contacts,
    )

    close_resid = np.unique(mol.resid[segidx[contacts]])
    res = segmask & np.isin(mol.resid, close_resid)
    return res


def _outOfBoundaries(
    mol,
    segname,
    xmin,
    xmax,
    ymin,
    ymax,
    zmin,
    zmax,
    cell=None,
    exclude_z=None,
):
    # Implementing the following atomselection
    # segid {segname} and same resid as (segid {segname} and (x < {xmin} or x > {xmax} or y < {ymin} or y > {ymax} or z < {zmin} or z > {zmax}))

    segnamesel = mol.segid == segname
    if cell is None:
        oob = (
            (mol.x < xmin)
            | (mol.x > xmax)
            | (mol.y < ymin)
            | (mol.y > ymax)
            | (mol.z < zmin)
            | (mol.z > zmax)
        )
    else:
        # Cull the Wigner-Seitz cell, not the parallelepiped: equal volume but a
        # smaller inscribed diameter, which would leave solute atoms unwatered.
        vectors, center = cell
        normals, offsets = _ws_halfspaces(vectors)
        dist = (mol.coords[:, :, 0] - center) @ normals.T
        oob = np.any(dist > offsets + 1e-6, axis=1)

    if exclude_z is not None:
        oob = oob | ((mol.z > exclude_z[0]) & (mol.z < exclude_z[1]))

    residsel = np.isin(mol.resid, mol.resid[segnamesel & oob])
    return segnamesel & residsel
