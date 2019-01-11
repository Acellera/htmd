import os
import string
from tempfile import NamedTemporaryFile

import numpy as np
import pybel
from htmd.builder.builder import autoSegment2
from htmd.builder.preparation import proteinPrepare
from htmd.molecule.molecule import Molecule
from htmd.molecule.writers import _deduce_PDB_atom_name, checkTruncations
from htmd.util import ensurelist
from joblib import Parallel, delayed
from tqdm import tqdm


def getPDBQTAtomType(atype, aidx, mol, aromaticNitrogen=False):
    tmptype = ''
    # carbons
    if atype == 'Car':
        tmptype = 'A'
    elif atype.startswith('C'):
        tmptype = 'C'
    # nitrogens
    if atype.startswith('N'):
        tmptype = 'N'
        if atype in ['Nam', 'Npl', 'Ng+']:
            bs, bo = np.where(mol.bonds == aidx)
            if len(bs) == 2:
                tmptype += 'A'
        elif atype == 'Nar':
            # if mol.resname[aidx] == 'HIE':
            bs, bo = np.where(mol.bonds == aidx)
            if len(bs) == 2:
                tmptype += 'a' if aromaticNitrogen else 'A'
            else:
                tmptype += 'n' if aromaticNitrogen else ''
        elif atype[-1] != '+':
            tmptype += 'A'
    # elif atype.startswith('N'):
    #   print(atype, aidx)
    #  tmptype = 'NA'
    # oxygens
    if atype.startswith('O'):
        tmptype = 'OA'
    # sulfurs
    if atype.startswith('S'):
        tmptype = 'S'
        if atype not in ['Sox', 'Sac']:
            tmptype += 'A'

    # hydrogens
    if atype.startswith('H'):
        tmptype = 'H'
        # print(aidx)
        # print(np.where(mol.bonds == aidx))
        bond = np.where(mol.bonds == aidx)[0][0]
        oidx = [a for a in mol.bonds[bond] if a != aidx][0]
        if mol.element[oidx] not in ['C', 'A']:
            tmptype += 'D'
    if tmptype == '':
        tmptype = atype[0]

    return tmptype


def getProperties(mol):
    name = NamedTemporaryFile(suffix='.pdb').name
    mol.write(name)
    mpybel = next(pybel.readfile('pdb', name))
    # print(name)
    residues = pybel.ob.OBResidueIter(mpybel.OBMol)
    atoms = [[r.GetName(), r.GetNum(), r.GetAtomID(at), at.GetType(), round(at.GetPartialCharge(), 3)]
             for r in residues
             for at in pybel.ob.OBResidueAtomIter(r)]
    return atoms


def molPDBQT(mol, aromaticNitrogen=False, inplace=False):
    if mol.bonds.shape[0] < mol.numAtoms:
        raise ValueError('The protein does not have bonds. Assign them with htmd molecule _guessBonds or use the '
                         'prepProtForFeats')

    atomsProp = getProperties(mol)
    for n, a in enumerate(atomsProp):
        if a[0] == 'HIP':
            if a[2].strip().startswith('C') and a[2].strip() not in ['CA', 'C', 'CB']:
                a[3] = 'Car'

            #print(n, a)
    charges = ['{0:.3f}'.format(a[-1]) for a in atomsProp]
    pdbqtATypes = [getPDBQTAtomType(a[3], n, mol, aromaticNitrogen)
                   for n, a in enumerate(atomsProp)]
    cmol = mol
    if not inplace:
        cmol = mol.copy()
    cmol.element = np.array(pdbqtATypes, dtype='O')
    cmol.charge = np.array(charges, dtype='float32')

    return cmol


def writePDBQT(mol, fname):
    m = mol.copy()
    # m.segid = np.array([' '+str(c) if c > 0 else str(c) for c in mol.charge], dtype='O')
    m.segid = mol.charge  # np.array([], dtype='float32')
    # print(m.segid)
    PDBQTwrite(m, fname)
    # m.write(fname + '.pdbqt', type='pdb')


def PDBQTwrite(mol, filename, frames=None, writebonds=True):
    if frames is None and mol.numFrames != 0:
        frames = mol.frame
    else:
        frames = 0
    frames = ensurelist(frames)

    checkTruncations(mol)
    if mol.numFrames != 0:
        coords = np.atleast_3d(mol.coords[:, :, frames])
        box = mol.box[:, frames[0]]
    else:  # If Molecule only contains topology, PDB requires some coordinates so give it zeros
        coords = np.zeros((mol.numAtoms, 3, 1), dtype=np.float32)
        box = None

    numFrames = coords.shape[2]
    nAtoms = coords.shape[0]

    serial = np.arange(1, np.size(coords, 0) + 1).astype(object)
    serial[serial > 99999] = '*****'
    serial = serial.astype('U5')

    if nAtoms > 0:
        if coords.max() >= 1E8 or coords.min() <= -1E7:
            raise RuntimeError(
                'Cannot write PDB coordinates with values smaller than -1E7 or larger than 1E8')
        if mol.occupancy.max() >= 1E6 or mol.occupancy.min() <= -1E5:
            raise RuntimeError(
                'Cannot write PDB occupancy with values smaller than -1E5 or larger than 1E6')
        if mol.beta.max() >= 1E6 or mol.beta.min() <= -1E5:
            raise RuntimeError(
                'Cannot write PDB beta/temperature with values smaller than -1E5 or larger than 1E6')

    fh = open(filename, 'w')

    if box is not None and not np.all(mol.box == 0):
        fh.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 \n" %
                 (box[0], box[1], box[2], 90, 90, 90))

    for f in range(numFrames):
        fh.write("MODEL    %5d\n" % (frames[f] + 1))
        for i in range(0, len(mol.record)):
            name = _deduce_PDB_atom_name(mol.name[i], mol.resname[i])

            fh.write(
                "{!s:6.6}{!s:>5.5} {}{!s:>1.1}{!s:4.4}{!s:>1.1}{!s:>4.4}{!s:>1.1}   {}{}{}{}{}    {:6.3f} {!s:<2.2}  \n".format(
                    mol.record[i],
                    serial[i], name, mol.altloc[i],
                    mol.resname[i], mol.chain[i],
                    mol.resid[i],
                    mol.insertion[i],
                    '{:8.3f}'.format(coords[i, 0, f])[:8],
                    '{:8.3f}'.format(coords[i, 1, f])[:8],
                    '{:8.3f}'.format(coords[i, 2, f])[:8],
                    '{:6.2f}'.format(mol.occupancy[i])[:6],
                    '{:6.2f}'.format(mol.beta[i])[:6],
                    mol.segid[i],
                    mol.element[i]
                )
            )
            # TODO : convert charges to ints if we ever write them
            # if i < len(mol.record) - 1 and mol.segid[i] != mol.segid[i + 1]:
            #    fh.write("TER\n")

        if writebonds and mol.bonds is not None and len(mol.bonds) != 0:
            bondedatoms = np.unique(mol.bonds)
            # Don't print bonds over 99999 as it overflows the field
            bondedatoms = bondedatoms[bondedatoms < 99998]

            for a in bondedatoms:
                partners = mol.bonds[mol.bonds[:, 0] == a, 1]
                partners = np.unique(
                    np.append(partners, mol.bonds[mol.bonds[:, 1] == a, 0]))
                # Don't print bonds over 99999 as it overflows the field
                partners = partners[partners < 99998] + 1
                # I need to support multi-line printing of atoms with more than 4 bonds
                while len(partners) >= 3:  # Write bonds as long as they are more than 3 in fast more
                    fh.write("CONECT%5d%5d%5d%5d\n" %
                             (a + 1, partners[0], partners[1], partners[2]))
                    partners = partners[3:]
                if len(partners) > 0:  # Write the rest of the bonds
                    line = "CONECT%5d" % (a + 1)
                    for p in partners:
                        line = "%s%5d" % (line, p)
                    fh.write(line)
                    fh.write('\n')

        fh.write("ENDMDL\n")
    fh.write("END\n")

    fh.close()


def prepProtForFeats(mol, pH=7, thickness=None, autosegment=False, prepare=True):
    list_chains = list(string.ascii_uppercase)
    mol.filter('protein')

    if autosegment:
        mol = autoSegment2(mol)
        segids = np.unique(mol.segid)

        for n, s in enumerate(segids):
            idxs = np.where(mol.segid == s)
            mol.chain[idxs] = list_chains[n]

    if prepare:
        pmol = proteinPrepare(
            mol, pH=pH, hydrophobicThickness=thickness) if prepare else mol
        bonds = pmol._guessBonds()
        pmol.bonds = bonds
    else:
        pmol = mol

    return pmol


def _getHydrophonic(atypes):
    return atypes == 'C'


def _getAromatic(atypes):
    return (atypes == 'A') | (atypes == 'Na') | (atypes == 'Nn')


def _getAcceptor(atypes):
    return (atypes == 'OA') | (atypes == 'NA') | (atypes == 'SA') | (atypes == 'Na')


def _getDonors(atypes, bonds):
    hidxs = np.where((atypes == 'HD') | (atypes == 'HS'))[0]

    donors = np.zeros(atypes.shape[0], dtype=bool)
    for n, hidx in enumerate(hidxs):
        bidx, baidx = np.where(bonds == hidx)
        bidx = bidx[0]
        baidx = baidx[0]

        baidx = 0 if baidx == 1 else 1
        oaidx = bonds[bidx][baidx]
        oel = atypes[oaidx]
        isdonor = 1 if oel[0] in ['N', 'O', 'S'] else 0
        donors[n] = isdonor
    return donors


def _getPosIonizable(mol):
    # arginine, lysine and histidine
    posIonizables = np.zeros(mol.numAtoms, dtype=bool)

    # ARG
    n_idxs = np.where(((mol.resname == 'ARG') | (mol.resname == 'AR0')) & (
        mol.element == 'N') & (mol.name != 'N'))
    allc_idxs = np.where((mol.resname == 'ARG') & (
        mol.element == 'C') & (mol.name != 'C'))[0]
    c_idxs = []
    for c in allc_idxs:
        bs = np.where(mol.bonds == c)[0]
        if len(bs) == 3:
            c_idxs.append(c)

    aidxs = n_idxs[0].tolist() + c_idxs

    # LYS
    n_idxs = np.where(((mol.resname == 'LYS') | (mol.resname == 'LYN')) & (
        mol.element == 'N') & (mol.name != 'N'))
    aidxs += n_idxs[0].tolist()

    # HIS, HID, HIE, HIP, HSD, HSE
    n_idxs = np.where(((mol.resname == 'HIS') | (mol.resname == 'HID') |
                       (mol.resname == 'HIE') | (mol.resname == 'HIP') |
                       (mol.resname == 'HSE') | (mol.resname == 'HSD') |
                       (mol.resname == 'HSP')) &
                      ((mol.element == 'N') | (mol.element == 'NA') | (mol.element == 'Nn') | (mol.element == 'Na')) &
                      (mol.name != 'N'))

    c_idxs = np.where(((mol.resname == 'HIS') | (mol.resname == 'HID') |
                       (mol.resname == 'HIE') | (mol.resname == 'HIP') |
                       (mol.resname == 'HSE') | (mol.resname == 'HSD') |
                       (mol.resname == 'HSP')) &
                      (mol.element == 'A'))

    aidxs += n_idxs[0].tolist() + c_idxs[0].tolist()

    posIonizables[aidxs] = 1

    return posIonizables


def _getNegIonizable(mol):
    # aspartic and glutamate
    negIonizables = np.zeros(mol.numAtoms, dtype=bool)

    # ASP
    o_idxs = np.where(((mol.resname == 'ASP') | (mol.resname == 'ASH')) &
                      (mol.element == 'OA') & (mol.name != 'O'))
    allc_idxs = np.where(((mol.resname == 'ASP') | (mol.resname == 'ASH')) &
                         (mol.element == 'C') & (mol.name != 'C'))[0]
    c_idxs = []
    for c in allc_idxs:
        bs = np.where(mol.bonds == c)[0]
        if len(bs) == 3:
            c_idxs.append(c)
    aidxs = o_idxs[0].tolist() + c_idxs

    # Glutamate
    o_idxs = np.where(((mol.resname == 'GLU') | (mol.resname == 'GLH')) &
                      (mol.element == 'OA') & (mol.name != 'O'))

    allc_idxs = np.where(((mol.resname == 'GLU') | (mol.resname == 'GLH')) &
                         (mol.element == 'C') & (mol.name != 'C'))[0]
    c_idxs = []
    for c in allc_idxs:
        bs = np.where(mol.bonds == c)[0]
        if len(bs) == 3:
            c_idxs.append(c)
    aidxs += o_idxs[0].tolist() + c_idxs

    negIonizables[aidxs] = 1

    return negIonizables


def _getOccupancy(elements):
    return np.array(elements) != 'H'


def _getMetals(elements):
    return (elements == 'MG') | (elements == 'ZN') | (elements == 'MN') | \
           (elements == 'CA') | (elements == 'FE') | (elements == 'HG') | \
           (elements == 'CD') | (elements == 'NI') | (elements == 'CO') | \
           (elements == 'CU') | (elements == 'K') | (elements == 'LI') | \
           (elements == 'Mg') | (elements == 'Zn') | (elements == 'Mn') | \
           (elements == 'Ca') | (elements == 'Fe') | (elements == 'Hg') | \
           (elements == 'Cd') | (elements == 'Ni') | (elements == 'Co') | \
           (elements == 'Cu') | (elements == 'Li')


def getFeatures(mol):
    atypes = mol.element
    elements = [el[0] for el in atypes]

    hydr = _getHydrophonic(atypes)
    arom = _getAromatic(atypes)
    acc = _getAcceptor(atypes)
    don = _getDonors(atypes, mol.bonds)
    pos = _getPosIonizable(mol)
    neg = _getNegIonizable(mol)
    metals = _getMetals(atypes)
    occ = _getOccupancy(elements)

    return np.vstack((hydr, arom, acc, don, pos, neg, metals, occ)).T.copy()


def myPdbqt(pdb, outfolder):
    oname = pdb.split('/')[-2]
    outname = './{}/{}.pdbqt'.format(outfolder, oname)
    if os.path.isfile(outname):
        return
    try:
        m = Molecule(pdb)
    except:
        return (pdb, 'read')
    try:
        pmol = prepProtForFeats(m)  # , prepare=False)
        pmolPDBQT = molPDBQT(pmol)
        writePDBQT(pmolPDBQT, outname)
    except:
        return (pdb, 'process')


def parallel(func, listobj, n_cpus=-1, *args):
    results = Parallel(n_jobs=n_cpus)(delayed(func)(ob, *args)
                                      for ob in tqdm(listobj))
    return results
