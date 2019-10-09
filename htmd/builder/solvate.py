# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from moleculekit.molecule import Molecule
import numpy as np
import sys
import logging
logger = logging.getLogger(__name__)


def solvate(mol, pad=None, minmax=None, negx=0, posx=0, negy=0, posy=0, negz=0, posz=0, buffer=2.4, watsize=65.4195,
            prefix='WT', keysel='name OH2', rotate=False, rotsel='all', rotinc=36, spdb=None,
            spsf=None, stop=None):
    """ Solvates the system in a water box


    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The molecule object we want to solvate
    pad : float
        The padding to add to the minmax in all dimensions. You can specify different padding in each dimension using
        the negx, negy, negz, posx, posy, posz options. This option will override any values in the neg and pos options.
    minmax : list
        Min and max dimensions. Should be a 2D matrix of the form [[minx, miny, minz], [maxx, maxy, maxz]]. If none is
        given, it is calculated from the minimum and maximum coordinates in the mol.
    negx : float
        The padding in the -x dimension
    posx : float
        The padding in the +x dimension
    negy : float
        The padding in the -y dimension
    posy : float
        The padding in the +y dimension
    negz : float
        The padding in the -z dimension
    posz : float
        The padding in the +z dimension
    buffer : float
        How much buffer space to leave empty between waters and other molecules
    watsize : float
        The size of the water box
    prefix : str
        The prefix used for water segments
    keysel : str
        The key selection for water atoms
    rotate : bool
        Enable automated rotation of molecule to fit best in box
    rotsel : str
        The selection of atoms to rotate
    rotinc : float
        The increment in degrees to rotate
    spdb : str
        The path to the water pdb file
    spsf : str
        The path to the water psf file
    stop : str
        The path to the water topology file

    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        A solvated molecule

    Examples
    --------
    >>> smol = solvate(mol, pad=10)
    >>> smol = solvate(mol, minmax=[[-20, -20, -20],[20, 20, 20]])
    """
    from tqdm import tqdm
    from htmd.home import home

    mol = mol.copy()
    if mol.numFrames > 1:
        logger.warning('Multiple frames in Molecule. Solvate keeps only frame 0 and discards the rest.')
        mol.coords = np.atleast_3d(mol.coords[:, :, 0])

    if spdb is None:
        spdb = os.path.join(home(shareDir=True), 'solvate', 'wat.pdb')

    if os.path.isfile(spdb):
        logger.info('Using water pdb file at: ' + spdb)
        water = Molecule(spdb)
    else:
        raise NameError('No solvent pdb file found in ' + spdb)

    if pad is not None:
        negx = pad; posx = pad; negy = pad; posy = pad; negz = pad; posz = pad

    if rotate:
        raise NameError('Rotation not implemented yet')

    # Calculate min max coordinates from molecule
    if mol.numAtoms > 0:
        minmol = np.min(mol.get('coords'), axis=0)
        maxmol = np.max(mol.get('coords'), axis=0)
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
        preexist = len(np.unique(mol.get('segid', sel='segid "{}.*"'.format(prefix))))
    else:
        preexist = 0

    numsegs = nx * ny * nz
    logger.info('Replicating ' + str(numsegs) + ' water segments, ' + str(nx) + ' by ' + str(ny) + ' by ' + str(nz))

    # Check that we won't run out of segment name characters, and switch to
    # using hexadecimal or alphanumeric naming schemes in cases where decimal
    # numbered segnames won't fit into the field width.
    testsegname    = '{0}{1:d}'.format(prefix, numsegs + preexist)
    testsegnamehex = '{0}{1:X}'.format(prefix, numsegs + preexist)
    writemode = 'decimal'
    if len(testsegname) > 4 and len(testsegnamehex) <= 4:
        writemode = 'hex'
        logger.warning('Warning: decimal naming would overrun segname field. Using hexadecimal segnames instead...')
    elif len(testsegnamehex) > 4:
        writemode = 'alphanum'
        logger.warning('Warning: decimal or hex naming would overrun segname field. Using alphanumeric segnames instead...')

    minx = minmol[0] - buffer; miny = minmol[1] - buffer; minz = minmol[2] - buffer
    maxx = maxmol[0] + buffer; maxy = maxmol[1] + buffer; maxz = maxmol[2] + buffer

    bar = tqdm(total=nx*ny*nz, desc='Solvating')
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

                if writemode == 'decimal':
                    segname = '{0}{1:d}'.format(prefix, n)
                elif writemode == 'hex':
                    segname = '{0}{1:x}'.format(prefix, n)
                elif writemode == 'alphanum':
                    segname = '{0}{1:c}{2:c}{3:c}'.format(prefix, int(np.floor(np.floor(n/26)/26) + 65), int(np.mod(np.floor(n/26), 26) + 65), int(np.mod(n, 26) + 65))

                waterboxes[w] = water.copy()
                waterboxes[w].moveBy([movex, movey, movez])
                waterboxes[w].set('segid', segname)

                mol.append(waterboxes[w])
                watsel = mol.segid == segname

                selover = np.zeros(len(watsel), dtype=bool)
                if xoverlap and yoverlap and zoverlap:  # Remove water overlapping with other segids
                    selover = _overlapWithOther(mol, segname, buffer)
                # Remove water outside the boundaries
                selout = _outOfBoundaries(mol, segname, xmin, xmax, ymin, ymax, zmin, zmax)
                sel = selover | selout

                #mol.write('temp.pdb')
                mol.filter(mol.segid != segname, _logger=False)
                waterboxes[w].filter(np.invert(sel[watsel]), _logger=False)
                #waterboxes[w].write('wat' + str(w) + '.pdb')
                n += 1
                w += 1
                bar.update(1)
    bar.close()

    waters = 0
    for i in range(numsegs):
        waters += waterboxes[i].numAtoms
        if waterboxes[i].numAtoms != 0:
            mol.append(waterboxes[i])

    logger.info('{} water molecules were added to the system.'.format(int(waters/3)))
    return mol


def _overlapWithOther(mol, segname, buffer):
    sel = 'segid {} and same resid as (segid {} and within {} of not segid {})'.format(segname, segname, buffer, segname)
    return mol.atomselect(sel)


def _outOfBoundaries(mol, segname, xmin, xmax, ymin, ymax, zmin, zmax):
    sel = 'segid {} and same resid as (segid {} and (x < {} or x > {} or y < {} or y > {} or z < {} or z > {}))'.format(
        segname, segname, xmin, xmax, ymin, ymax, zmin, zmax
    )
    return mol.atomselect(sel)
