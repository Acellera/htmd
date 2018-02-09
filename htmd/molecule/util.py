# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial.distance import cdist
import ctypes as ct
import os
from htmd.home import home
import platform
import logging
logger = logging.getLogger(__name__)

libdir = home(libDir=True)
if platform.system() == "Windows":
    tmalignlib = ct.cdll.LoadLibrary(os.path.join(libdir, "tmalign.dll"))
else:
    tmalignlib = ct.cdll.LoadLibrary(os.path.join(libdir, "tmalign.so"))


def molTMscore(mol, ref, selCAmol, selCAref):
    """ Calculates the TMscore between two Molecules

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A Molecule containing a single or multiple frames
    ref : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A reference Molecule containing a single frame. Will automatically keep only ref.frame.
    selCAmol : numpy.ndarray
        An atomselection array of booleans or indexes of the CA atoms for mol
    selCAref : numpy.ndarray
        An atomselection array of booleans or indexes of the CA atoms for ref

    Returns
    -------
    tmscoreRef : numpy.ndarray
        TMscore normalized by length of ref
    rmsd : numpy.ndarray
        RMSD only OF COMMON RESIDUES for all frames. This is not the same as a full protein RMSD!!!

    Examples
    --------
    tmscore, rmsd = molTMscore(mol, ref, mol.atomselect('protein'), ref.atomselect('protein'))
    """
    from htmd.builder.builder import sequenceID
    from htmd.molecule.molecule import _residueNameTable

    def calculateVariables(currmol):
        res = sequenceID((currmol.resid, currmol.insertion, currmol.segid, currmol.chain))
        caidx = currmol.name == 'CA'
        res = np.unique(res)
        reslen = len(res)
        # Calculate the protein sequence
        seq = ''.join([_residueNameTable[x] for x in currmol.resname[caidx]])
        seq = ct.c_char_p(seq.encode('utf-8'))

        # Keep only CA coordinates
        coords = currmol.coords[caidx, :, :].copy()
        return reslen, res.astype(np.int32), seq, coords

    mol = mol.copy()
    ref = ref.copy()
    mol.filter(selCAmol, _logger=False)
    ref.filter(selCAref, _logger=False)
    ref.dropFrames(keep=ref.frame)

    reslenMOL, residMOL, seqMOL, coordsMOL = calculateVariables(mol)
    reslenREF, residREF, seqREF, coordsREF = calculateVariables(ref)

    # DLLEXPORT void tmalign(int xlen, int ylen, int* xresno, int* yresno, char* seqx, char* seqy,
    # float* xcoor, float* ycoor, int nframes,
    # double *TM1, double *TM2, double *rmsd)
    # tmalignlib.tmalign.argtypes = [ct.c_int, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p,
    #                                ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.c_int,
    #                                ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
    resTM1 = (ct.c_double * mol.numFrames)()
    resTM2 = (ct.c_double * mol.numFrames)()
    resRMSD = (ct.c_double * mol.numFrames)()
    tmalignlib.tmalign(ct.c_int(reslenREF),
                       ct.c_int(reslenMOL),
                       residREF.ctypes.data_as(ct.POINTER(ct.c_int32)),
                       residMOL.ctypes.data_as(ct.POINTER(ct.c_int32)),
                       seqREF,
                       seqMOL,
                       coordsREF.ctypes.data_as(ct.POINTER(ct.c_float)),
                       coordsMOL.ctypes.data_as(ct.POINTER(ct.c_float)),
                       ct.c_int(mol.numFrames),
                       ct.byref(resTM1),
                       ct.byref(resTM2),
                       ct.byref(resRMSD))
    resTM1 = np.ctypeslib.as_array(resTM1)
    resRMSD = np.ctypeslib.as_array(resRMSD)
    return resTM1.astype(np.float32), resRMSD.astype(np.float32)



def molRMSD(mol, refmol, rmsdsel1, rmsdsel2):
    """ Calculates the RMSD between two Molecules

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
    refmol
    rmsdsel1
    rmsdsel2

    Returns
    -------
    rmsd : float
        The RMSD between the two structures
    """
    dist = mol.coords[rmsdsel1, :, :] - refmol.coords[rmsdsel2, :, :]
    rmsd = np.sqrt(np.mean(np.sum(dist * dist, axis=1), axis=0))
    return np.squeeze(rmsd)


def orient(mol, sel="all"):
    """Rotate a molecule so that its main axes are oriented along XYZ.

    The calculation is based on the axes of inertia of the given
    selection, but masses will be ignored. After the operation, the
    main axis will be parallel to the Z axis, followed by Y and X (the
    shortest axis). Only the first frame is oriented.  The reoriented
    molecule is returned.

    Parameters
    ----------
    mol :
        The Molecule to be rotated
    sel : 
        Atom selection on which the rotation is computed

    Examples
    --------
    >>> mol = Molecule("1kdx")
    >>> mol = orient(mol,"chain B")

    """
    if mol.numFrames != 1:
        logger.warning("Only the first frame is considered for the orientation")
    m = mol.copy()
    s = m.atomselect(sel)
    x = m.coords[s,:,0]
    c = np.cov(x.T)
    ei = np.linalg.eigh(c)
    logger.info("Moments of intertia: "+str(ei[0]))
    ev=ei[1].T
    if np.linalg.det(ev)<0: ev=-ev # avoid inversions
    m.rotateBy(ev)
    return(m)


def sequenceID(field, prepend=None):
    """ Array of integers which increments at value change of another array

    Parameters
    ----------
    field : np.ndarray or tuple
        An array of values. Once a change in value happens, a new ID will be created in `seq`.
        If a tuple of ndarrays is passed, a change in any of them will cause an increase in `seq`.
    prepend : str
        A string to prepend to the incremental sequence

    Returns
    -------
    seq : np.ndarray
        An array of equal size to `field` containing integers which increment every time there is a change in `field`

    Examples
    --------
    >>> # A change in resid, insertion, chain or segid will cause an increase in the sequence
    >>> sequenceID((mol.resid, mol.insertion, mol.chain, mol.segid))
    array([  1,   1,   1, ..., 285, 286, 287])
    >>> # it is typically used to renumber resids as follows
    >>> mol.set('resid', sequenceID((mol.resid, mol.insertion, mol.chain, mol.segid)))
    """
    if isinstance(field, tuple):
        fieldlen = len(field[0])
    else:
        fieldlen = len(field)

    if prepend is None:
        seq = np.zeros(fieldlen, dtype=int)
    else:
        seq = np.empty(fieldlen, dtype=object)

    c = int(0)
    if prepend is None:
        seq[0] = c
    else:
        seq[0] = prepend + str(c)

    for i in range(1, fieldlen):
        if isinstance(field, tuple):  # Support tuples of multiple fields. Change in any of them will cause an increment
            for t in field:
                if t[i-1] != t[i]:
                    c += 1  # new sequence id
                    break
        elif field[i-1] != field[i]:
            c += 1  # new sequence id
        if prepend is None:
            seq[i] = c
        else:
            seq[i] = prepend + str(c)
    return seq


def _missingChain(mol):
    if mol.chain is None or np.size(mol.chain) == 0:
        raise NameError('Segid fields have to be set for all atoms in the Molecule object before building.')
    empty = [True if len(c) == 0 else False for c in mol.chain]
    if np.any(empty):
        idx = np.where(empty)[0]
        if len(idx) == 1:
            raise NameError('Atom ' + str(idx) + ' does not have a chain defined.')
        elif len(idx) <= 5:
            raise NameError('Atoms ' + str(idx) + ' do not have a chain defined.')
        else:
            raise NameError('Atoms [' + str(idx[0]) + ',' + str(idx[1]) + ',...,' + str(idx[-1]) + '] do not have chain defined.')


def _missingSegID(mol):
    if mol.segid is None or np.size(mol.segid) == 0:
        raise NameError('Segid fields have to be set for all atoms in the Molecule object before building.')
    empty = [True if len(s) == 0 else False for s in mol.segid]
    if np.any(empty):
        idx = np.where(empty)[0]
        if len(idx) == 1:
            raise NameError('Atom ' + str(idx) + ' does not have a segid defined.')
        elif len(idx) <= 5:
            raise NameError('Atoms ' + str(idx) + ' do not have a segid defined.')
        else:
            raise NameError('Atoms [' + str(idx[0]) + ',' + str(idx[1]) + ',...,' + str(idx[-1]) + '] do not have segid defined.')


def maxDistance(mol, sel='all', origin=[0, 0, 0]):
    """ Calculates the max distance of a set of atoms from an origin

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule containing the atoms
    sel : str
        Atomselection for atoms for which to calculate distances
    origin : list
        The origin x,y,z coordinates

    Returns
    -------
    maxd : float
        The maximum distance in Angstrom

    Example
    -------
    >>> y = maxDistance(mol, sel='protein', origin=[0, 0, 0])
    >>> print(round(y,2))
    48.39
    """
    coors = mol.get('coords', sel=sel)
    dists = cdist(np.atleast_2d(coors), np.atleast_2d(origin))
    return np.max(dists)


def boundingBox(mol, sel='all'):
    """ Calculates the bounding box of a selection of atoms.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The molecule containing the atoms
    sel : str
        An Atomselection string of atoms.

    Returns
    -------
    bbox : np.ndarray
        The bounding box around the atoms selected in `sel`.

    Example
    -------
    >>> boundingBox(mol, sel='chain A')
    array([[-17.3390007 , -10.43700027,  -1.43900001],
           [ 25.40600014,  27.03800011,  46.46300125]], dtype=float32)

    """
    coords = mol.get('coords', sel=sel)
    maxc = np.squeeze(np.max(coords, axis=0))
    minc = np.squeeze(np.min(coords, axis=0))
    return np.vstack((minc, maxc))


def uniformRandomRotation():
    """
    Return a uniformly distributed rotation 3 x 3 matrix

    The initial description of the calculation can be found in the section 5 of "How to generate random matrices from
    the classical compact groups" of Mezzadri (PDF: https://arxiv.org/pdf/math-ph/0609050.pdf; arXiv:math-ph/0609050;
    and NOTICES of the AMS, Vol. 54 (2007), 592-604). Sample code is provided in that section as the ``haar_measure``
    function.

    Apparently this code can randomly provide flipped molecules (chirality-wise), so a fix found in
    https://github.com/tmadl/sklearn-random-rotation-ensembles/blob/5346f29855eb87241e616f6599f360eba12437dc/randomrotation.py
    was applied.

    Returns
    -------
    M : np.ndarray
        A uniformly distributed rotation 3 x 3 matrix
    """
    q, r = np.linalg.qr(np.random.normal(size=(3, 3)))
    M = np.dot(q, np.diag(np.sign(np.diag(r))))
    if np.linalg.det(M) < 0:  # Fixing the flipping
        M[:, 0] = -M[:, 0]  # det(M)=1
    return M


def writeVoxels(arr, filename, vecMin, vecMax, vecRes):
    """ Writes grid free energy to cube file

    Parameters
    ----------
    arr: np.ndarray
            array with volumetric data
    filename: str
            string with the name of the cubefile
    vecMin: np.ndarray
            3D vector denoting the minimal corner of the grid
    vecMax np.ndarray
            3D vector denoting the maximal corner of the grid
    vecRes: np.ndarray
            3D vector denoting the resolution of the grid in each dimension
    """

    outFile = open(filename, 'w')

    # conversion to gaussian units
    L = 1/0.52917725
    gauss_bin = vecRes*L
    #minCorner = 0.5*L*(vecMin - vecMax + vecRes)
    minCorner = L * (vecMin + 0.5 * vecRes)

    ngrid = np.array(np.floor((vecMax - vecMin) / vecRes), dtype=int)

    # write header
    outFile.write("CUBE FILE\n")
    outFile.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (1, minCorner[0], minCorner[1], minCorner[2]))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[0], gauss_bin[0], 0, 0))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[1], 0, gauss_bin[1], 0))
    outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[2], 0, 0, gauss_bin[2]))
    outFile.write("%5d %12.6f %12.6f %12.6f %12.6f\n" % (1, 0, minCorner[0], minCorner[1], minCorner[2]))

    # main loop
    cont = 0

    for i in range(ngrid[0]):
        for j in range(ngrid[1]):
            for k in range(ngrid[2]):
                outFile.write("%13.5g" % arr[i][j][k])
                if np.mod(cont, 6) == 5:
                    outFile.write("\n")
                cont += 1

    outFile.close()


def sequenceStructureAlignment(mol, ref, molseg=None, refseg=None, maxalignments=10, nalignfragment=1):
    """ Aligns two structures by their longests sequences alignment

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule we want to align
    ref : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The reference Molecule to which we want to align
    molseg : str
        The segment of `mol` we want to align
    refseg : str
        The segment of `ref` we want to align to
    maxalignments : int
        The maximum number of alignments we want to produce
    nalignfragment : int
        The number of fragments used for the alignment.

    Returns
    -------
    mols : list
        A list of Molecules each containing a different alignment.
    """
    from htmd.util import ensurelist
    try:
        from Bio import pairwise2
    except ImportError as e:
        raise ImportError('You need to install the biopython package to use this function. Try using `conda install biopython`.')
    from Bio.SubsMat import MatrixInfo as matlist

    if len([x for x in np.unique(mol.altloc) if len(x)]) > 1:
        raise RuntimeError('Alternative atom locations detected in `mol`. Please remove these before calling this function.')
    if len([x for x in np.unique(ref.altloc) if len(x)]) > 1:
        raise RuntimeError('Alternative atom locations detected in `ref`. Please remove these before calling this function.')

    seqmol = mol.sequence()
    seqref = ref.sequence()

    if len(seqmol) > 1:
        logger.info('Multiple segments ({}) detected in `mol`. Alignment will be done on all. Otherwise please specify which segment to align.'.format(list(seqmol.keys())))
        seqmol = mol.sequence(noseg=True)
    if len(seqref) > 1:
        logger.info('Multiple segments ({}) detected in `ref`. Alignment will be done on all. Otherwise please specify which segment to align.'.format(list(seqref.keys())))
        seqref = ref.sequence(noseg=True)

    if molseg is None:
        molseg = list(seqmol.keys())[0]
    if refseg is None:
        refseg = list(seqref.keys())[0]

    def getSegIdx(m, mseg):
        # Calculate the atoms which belong to the selected segments
        if isinstance(mseg, str) and mseg == 'protein':
            msegidx = m.atomselect('protein and name CA')
        else:
            msegidx = np.zeros(m.numAtoms, dtype=bool)
            for seg in ensurelist(mseg):
                msegidx |= (m.segid == seg) & (m.name == 'CA')
        return np.where(msegidx)[0]
    molsegidx = getSegIdx(mol, molseg)
    refsegidx = getSegIdx(ref, refseg)

    # Create fake residue numbers for the selected segment
    molfakeresid = sequenceID((mol.resid[molsegidx], mol.insertion[molsegidx], mol.chain[molsegidx]))
    reffakeresid = sequenceID((ref.resid[refsegidx], ref.insertion[refsegidx], ref.chain[refsegidx]))

    # TODO: Use BLOSUM62?
    alignments = pairwise2.align.globaldx(seqref[refseg], seqmol[molseg], matlist.blosum62)
    numaln = len(alignments)

    if numaln > maxalignments:
        logger.warning('{} alignments found. Limiting to {} as specified in the `maxalignments` argument.'.format(numaln, maxalignments))

    alignedstructs = []
    for i in range(min(maxalignments, numaln)):
        refaln = np.array(list(alignments[i][0]))
        molaln = np.array(list(alignments[i][1]))
 
        # By doing cumsum we calculate how many letters were before the current letter (i.e. residues before current)
        residref = np.cumsum(refaln != '-') - 1  # Start them from 0
        residmol = np.cumsum(molaln != '-') - 1  # Start them from 0

        # Find the region of maximum alignment between the molecules
        dsig = np.hstack(([False], (refaln != '-') & (molaln != '-'), [False])).astype(int)
        dsigdiff = np.diff(dsig)
        startIndex = np.where(dsigdiff > 0)[0]
        endIndex = np.where(dsigdiff < 0)[0]
        duration = endIndex - startIndex
        duration_sorted = np.sort(duration)[::-1]

        _list_starts = []
        _list_finish = []
        for n in range(nalignfragment):
            if n == len(duration):
                break
            idx = np.where(duration == duration_sorted[n])[0]
            start = startIndex[idx][0]
            finish = endIndex[idx][0]
            _list_starts.append(start)
            _list_finish.append(finish)

        # Get the "resids" of the aligned residues only
        refalnresid = np.concatenate([ residref[start:finish] for start, finish in zip(_list_starts,_list_finish)])
        molalnresid = np.concatenate([ residmol[start:finish] for start, finish in zip(_list_starts, _list_finish) ])
        refidx = []
        for r in refalnresid:
            refidx += list(refsegidx[reffakeresid == r])
        molidx = []
        for r in molalnresid:
            molidx += list(molsegidx[molfakeresid == r])        

        molboolidx = np.zeros(mol.numAtoms, dtype=bool)
        molboolidx[molidx] = True
        refboolidx = np.zeros(ref.numAtoms, dtype=bool)
        refboolidx[refidx] = True

        start_residues = np.concatenate([ mol.resid[molsegidx[molfakeresid == residmol[r]]] for r in _list_starts])
        finish_residues = np.concatenate([ mol.resid[molsegidx[molfakeresid == residmol[r-1]]] for r in _list_finish])
        logger.info('Alignment #{} was done on {} residues: mol segid {} resid {}'.format(
            i, len(refalnresid), np.unique(mol.segid[molidx])[0], ', '.join(['{}-{}'.format(s,f) for s, f in zip(start_residues,finish_residues)])  ))

        alignedmol = mol.copy()
        alignedmol.align(molboolidx, ref, refboolidx)
        alignedstructs.append(alignedmol)

    return alignedstructs


def rcsbFindMutatedResidues(pdbid):
    import requests
    try:
        from bs4 import BeautifulSoup
        import lxml
    except ImportError:
        raise ImportError('You need to install the \'beautifulsoup4\' and \'lxml\' packages to use this function.')
    tomutate = {}

    connected = False
    while not connected:
        try:
            res = requests.get('http://www.rcsb.org/pdb/explore.do?structureId={}'.format(pdbid))
        except requests.ConnectionError as coer:
            import time
            time.sleep(5)
            continue
        connected = True

    soup = BeautifulSoup(res.text, 'lxml')
    table = soup.find(id='ModifiedResidueTable')

    if table:
        trs = table.find_all('tr')

        for tr in trs:
            td = tr.find_all('td')
            if td:
                mutname = td[0].find_all('a')[0].text.strip()
                orgname = td[5].text.strip()
                tomutate[mutname] = orgname
    return tomutate


def rcsbFindLigands(pdbid):
    import requests
    from bs4 import BeautifulSoup
    ligands = []

    connected = False
    while not connected:
        try:
            res = requests.get('http://www.rcsb.org/pdb/explore.do?structureId={}'.format(pdbid))
        except requests.ConnectionError as coer:
            import time
            time.sleep(5)
            continue
        connected = True

    soup = BeautifulSoup(res.text, 'lxml')
    table = soup.find(id='LigandsTable')
    if table:
        trs = table.find_all('tr')

        for tr in trs:
            td = tr.find_all('td')
            if td:
                name = td[0].find_all('a')[0].text.strip()
                ligands.append(name)
    return ligands


# def drawCube(mi, ma, viewer=None):
#     from htmd.vmdviewer import getCurrentViewer
#     if viewer is None:
#         viewer = getCurrentViewer()
#
#     viewer.send('draw materials off')
#     viewer.send('draw color red')
#     c = ''
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], ma[0], mi[1], mi[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], ma[1], mi[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], mi[1], ma[2])
#
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], ma[1], mi[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], mi[1], ma[2])
#
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], ma[0], ma[1], mi[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], mi[0], ma[1], ma[2])
#
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], ma[0], mi[1], ma[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], mi[0], ma[1], ma[2])
#
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], ma[1], mi[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], mi[0], ma[1], ma[2])
#     c += 'draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], mi[1], ma[2])
#     viewer.send(c)
#
#     """
#     draw line "$minx $miny $minz" "$maxx $miny $minz"
#     draw line "$minx $miny $minz" "$minx $maxy $minz"
#     draw line "$minx $miny $minz" "$minx $miny $maxz"
#     draw line "$maxx $miny $minz" "$maxx $maxy $minz"
#     draw line "$maxx $miny $minz" "$maxx $miny $maxz"
#     draw line "$minx $maxy $minz" "$maxx $maxy $minz"
#     draw line "$minx $maxy $minz" "$minx $maxy $maxz"
#     draw line "$minx $miny $maxz" "$maxx $miny $maxz"
#     draw line "$minx $miny $maxz" "$minx $maxy $maxz"
#     draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
#     draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
#     draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
#     """


def guessAnglesAndDihedrals(bonds, cyclicdih=False):
    """
    Generate a guess of angle and dihedral N-body terms based on a list of bond index pairs.
    """

    import networkx as nx

    g = nx.Graph()
    g.add_nodes_from(np.unique(bonds))
    g.add_edges_from([tuple(b) for b in bonds])

    angles = []
    for n in g.nodes():
        neighbors = list(g.neighbors(n))
        for e1 in range(len(neighbors)):
            for e2 in range(e1+1, len(neighbors)):
                angles.append((neighbors[e1], n, neighbors[e2]))

    angles = sorted([sorted([angle, angle[::-1]])[0] for angle in angles])
    angles = np.array(angles)

    dihedrals = []
    for a1 in range(len(angles)):
        for a2 in range(a1+1, len(angles)):
            a1a = angles[a1]
            a2a = angles[a2]
            a2f = a2a[::-1]  # Flipped a2a. We don't need flipped a1a as it produces the flipped versions of these 4
            if np.all(a1a[1:] == a2a[:2]) and (cyclicdih or (a1a[0] != a2a[2])):
                dihedrals.append(list(a1a) + [a2a[2]])
            if np.all(a1a[1:] == a2f[:2]) and (cyclicdih or (a1a[0] != a2f[2])):
                dihedrals.append(list(a1a) + [a2f[2]])
            if np.all(a2a[1:] == a1a[:2]) and (cyclicdih or (a2a[0] != a1a[2])):
                dihedrals.append(list(a2a) + [a1a[2]])
            if np.all(a2f[1:] == a1a[:2]) and (cyclicdih or (a2f[0] != a1a[2])):
                dihedrals.append(list(a2f) + [a1a[2]])

    dihedrals = sorted([sorted([dihedral, dihedral[::-1]])[0] for dihedral in dihedrals])
    dihedrals = np.array(dihedrals)

    return angles, dihedrals


# A test method
if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import numpy as np
    from os import path
    import doctest
    #doctest.testmod(extraglobs={"mol" : Molecule("3PTB")})

    expectedTMscore = np.array([ 0.21418524,  0.2367377 ,  0.23433833,  0.21362964,  0.20935164,
        0.20279461,  0.27012895,  0.22675238,  0.21230793,  0.2372011 ])
    expectedRMSD = np.array([ 3.70322128,  3.43637027,  3.188193  ,  3.84455877,  3.53053882,
        3.46781854,  2.93777629,  2.97978692,  2.70792428,  2.63051318])

    mol = Molecule(path.join(home(), 'data', 'tmscore', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'tmscore', 'traj.xtc'))
    ref = Molecule(path.join(home(), 'data', 'tmscore', 'ntl9_2hbb.pdb'))
    tmscore, rmsd = molTMscore(mol, ref, mol.atomselect('protein'), ref.atomselect('protein'))

    assert np.allclose(tmscore, expectedTMscore)
    assert np.allclose(rmsd, expectedRMSD)

    #
    # rhodopsin = Molecule('1F88')
    # d3r = Molecule('3PBL')
    # alnmol = sequenceStructureAlignment(rhodopsin, d3r)
