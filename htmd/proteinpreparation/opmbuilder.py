
from htmd import *
from .util import removeHET, removeLipidsInProtein, minimalrotation, tilemembrane
from sysprep.systempreparation import systemPreparation
import logging
logger = logging.getLogger(__name__)


def opmbuilder(pdbid, protfile, membfile, outdir, dbcursor=None):
    try:
        prot = Molecule(protfile)
    except:
        raise RuntimeError('ERROR: File {} failed to even load.'.format(protfile))
    logger.info('Loaded protein {}'.format(protfile))

    # Get coordinates of membrane dummy atoms
    uqz = np.unique(prot.get('coords', 'resname DUM')[:, 2])
    numlayers = len(uqz)
    if numlayers == 1:
        if dbcursor:
            dbcursor.execute('UPDATE report SET numlayers=? WHERE pdbid=?', (1, pdbid))
    #        return
    #    else:
    #        raise RuntimeError('Found only one Dummy layer. Skipping protein...')

    #if len(uqz) != 2:
    #    raise RuntimeError('Found more than 2 z-coordinate values for dummy atoms!')
    if numlayers == 2:
        dumwidth = abs(uqz[0]) + abs(uqz[1])
        if dumwidth < 30 or dumwidth > 35:
            if dbcursor:
                dbcursor.execute('UPDATE report SET membsize=? WHERE pdbid=?', (str(dumwidth), pdbid))
        #        return
        #    else:
        #        raise RuntimeError('Dummy membrane width ({} A) indicates bigger or smaller membrane than the one used!'.format(dumwidth))
        if dbcursor:
            dbcursor.execute('UPDATE report SET membsize=? WHERE pdbid=?', (str(dumwidth), pdbid))

    # Remove Hetero atoms
    #tmp = prot.remove('resname DUM')
    #idx = prot.remove('hetero')  # hetero thinks some ligands are nucleic acids and keeps them
    #prot = removeHET(prot)  # My removeHET is too conservative. There are cases of broken PDBs where the bond is too large to a atom and it should be removed
    idx = prot.remove('not protein')

    if dbcursor:
        dbcursor.execute('UPDATE report SET numhet=? WHERE pdbid=?', (str(len(idx)), pdbid))

    # Detect engineered residues and allow renaming

    # Automatically detect protein segments (should be chains?)
    prot.set('segid', 'P')  # TODO: I should join these two into a single autosegments function
    prot = segmentgaps(prot)
    if dbcursor:
        dbcursor.execute('UPDATE report SET numseg=? WHERE pdbid=?', (str(len(np.unique(prot.segid))), pdbid))

    # Mutations

    # Build proteins for charmm to add caps (should there be any?)
    try:
        prot = systemPreparation(prot)
        if dbcursor:
            dbcursor.execute('UPDATE report SET syspreppass=? WHERE pdbid=?', (str(1), pdbid))
    except:
        if dbcursor:
            dbcursor.execute('UPDATE report SET syspreppass=? WHERE pdbid=?', (str(0), pdbid))
            return
        else:
            raise RuntimeError('pdb2pqr failed due to missing backbone atoms')

    # TODO: DANGEROUS: Needed the remove 'not protein' for case 1m6b which has some MET atoms not bonded (why?) after sysprep.
    tmp = prot.remove('not protein')
    if dbcursor:
        dbcursor.execute('UPDATE report SET numnotprot=?, numprotatm=? WHERE pdbid=?', (str(len(tmp)), str(prot.numAtoms), pdbid))

    # Build protein solo
    try:
        prot = charmm.build(prot, ionize=False, outdir='/tmp/build/')
        if dbcursor:
            dbcursor.execute('UPDATE report SET protbuildpass=? WHERE pdbid=?', (str(1), pdbid))
    except:
        if dbcursor:
            dbcursor.execute('UPDATE report SET protbuildpass=? WHERE pdbid=?', (str(0), pdbid))
            return
        else:
            raise RuntimeError('protein only build failed')

    # Rotate proteins around Z to have the maximum variance in the box XY diagonal
    r = minimalrotation(prot)
    prot.rotate([0, 0, 1], r)
    logger.info('Rotated the protein by {} degrees around the Z axis.'.format(np.degrees(r)))

    # Move protein to [0,0]. Don't use .center as it would ruin the Z coordinate
    meanpos = np.mean(prot.get('coords'), axis=0)
    prot.moveBy([-meanpos[0], -meanpos[1], 0])

    # Replicate the membrane
    minc = np.min(prot.coords, axis=0).flatten()
    maxc = np.max(prot.coords, axis=0).flatten()
    membrane = Molecule(membfile)
    buffer = 20  # Add 20 A of membrane around the protein
    memb = tilemembrane(membrane, minc[0]-buffer, minc[1]-buffer, maxc[0]+buffer, maxc[1]+buffer)

    # remove any lipids inside the protein hull
    memb, num = removeLipidsInProtein(prot, memb)
    if dbcursor:
        dbcursor.execute('UPDATE report SET numlipidsinhull=? WHERE pdbid=?', (str(num), pdbid))

    # Append the membrane removing collisions
    system = prot.copy()
    system.append(memb, collisions=True, coldist=1.3)

    # Calculate box dimensions
    minz = np.min(system.coords, axis=0)[2] - 5
    maxz = np.max(system.coords, axis=0)[2] + 5
    minxy = np.min(memb.get('coords', 'water'), axis=0)
    maxxy = np.max(memb.get('coords', 'water'), axis=0)
    logger.info([[minxy[0], minxy[1], minz[0]], [maxxy[0], maxxy[1], maxz[0]]])

    if dbcursor:
        dbcursor.execute('UPDATE report SET boxx=?, boxy=?, boxz=? WHERE pdbid=?', (str(maxxy[0]-minxy[0]), str(maxxy[1]-minxy[1]), str(maxz[0]-minz[0]), pdbid))

    # Solvate the system
    system = solvate(system, minmax=[[minxy[0], minxy[1], minz[0]], [maxxy[0], maxxy[1], maxz[0]]])

    # Build the system
    try:
        system = charmm.build(system, outdir=outdir)
        if dbcursor:
            dbcursor.execute('UPDATE report SET fullbuildpass=?, numsysatm=? WHERE pdbid=?', (str(1), str(system.numAtoms), pdbid))
    except:
        if dbcursor:
            dbcursor.execute('UPDATE report SET fullbuildpass=? WHERE pdbid=?', (str(0), pdbid))
            return
        else:
            raise RuntimeError('full build failed')
    return system
