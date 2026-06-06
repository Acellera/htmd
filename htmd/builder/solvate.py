# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from moleculekit.molecule import Molecule
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


def solvate(
    mol: Molecule,
    pad: float | None = None,
    minmax: list | np.ndarray | None = None,
    centersel: str | np.ndarray | None = None,
    boxsize: float | list | np.ndarray | None = None,
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
        selected atoms is used. Must be combined with `boxsize`.
    boxsize : float or list or np.ndarray, optional
        Dimensions of the solvation box. A single float creates a cubic box;
        a 3-element list ``[sx, sy, sz]`` creates an axis-aligned box. Must be
        combined with `centersel`.
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
        Path to a custom water PDB file. If None, uses the built-in water box.

    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        A copy of the input molecule with water molecules added.

    Examples
    --------
    >>> smol = solvate(mol, pad=10)
    >>> smol = solvate(mol, minmax=[[-20, -20, -20],[20, 20, 20]])
    >>> smol = solvate(mol, centersel="protein", boxsize=100)
    >>> smol = solvate(mol, centersel="protein", boxsize=[80, 80, 120])
    """
    from tqdm import tqdm
    from htmd.home import home

    mol = mol.copy()
    if mol.numFrames > 1:
        logger.warning(
            "Multiple frames in Molecule. Solvate keeps only frame 0 and discards the rest."
        )
        mol.coords = np.atleast_3d(mol.coords[:, :, 0])

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
        spdb = os.path.join(home(shareDir=True), "solvate", "wat.pdb")

    if os.path.isfile(spdb):
        logger.info("Using water pdb file at: " + spdb)
        water = Molecule(spdb)
    else:
        raise NameError("No solvent pdb file found in " + spdb)

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
                    mol, segname, xmin, xmax, ymin, ymax, zmin, zmax
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


def _outOfBoundaries(mol, segname, xmin, xmax, ymin, ymax, zmin, zmax):
    # Implementing the following atomselection
    # segid {segname} and same resid as (segid {segname} and (x < {xmin} or x > {xmax} or y < {ymin} or y > {ymax} or z < {zmin} or z > {zmax}))

    segnamesel = mol.segid == segname
    oob = (
        (mol.x < xmin)
        | (mol.x > xmax)
        | (mol.y < ymin)
        | (mol.y > ymax)
        | (mol.z < zmin)
        | (mol.z > zmax)
    )
    residsel = np.isin(mol.resid, mol.resid[segnamesel & oob])
    return segnamesel & residsel
