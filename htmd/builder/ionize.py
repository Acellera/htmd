# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import scipy.spatial.distance as distance
from htmd.molecule.molecule import Molecule
from htmd.molecule.util import sequenceID
import logging
logger = logging.getLogger(__name__)


def ionize(netcharge, nwater, neutralize=True, saltconc=None, cation=None, anion=None, ff='charmm'):
    if cation is None:
        if ff == 'amber':
            cation = 'Na+'
        elif ff == 'charmm':
            cation = 'SOD'
        else:
            raise NameError('Invalid forcefield. Choose between ''charmm'' or ''amber''')
    if anion is None:
        if ff == 'amber':
            anion = 'Cl-'
        elif ff == 'charmm':
            anion = 'CLA'
        else:
            raise NameError('Invalid forcefield. Choose between ''charmm'' or ''amber''')

    if saltconc is not None:
        neutralize = False

    cationcharge = _ionGetCharge(cation, ff)
    anioncharge = _ionGetCharge(anion, ff)

    roundnetcharge = int(np.round(netcharge))
    done = False
    ncation = 0
    nanion = 0

    # First we neutralize the system
    if netcharge > 0:
        if anioncharge == -1:
            nanion = roundnetcharge
            ncation = 0
        elif anioncharge == -2:
            if np.mod(roundnetcharge, 2) == 0:  # Total charge is even
                nanion = int(roundnetcharge / 2)
                ncation = 0
            else:  # The charge is odd
                if cationcharge == 1:
                    nanion = 1 + int(np.floor(roundnetcharge / 2))
                    ncation = 1
                else:
                    raise NameError('Failed to neutralize because the cation and anion charges cannot reach 0')
    elif netcharge < 0:
        roundnetchargeabs = abs(roundnetcharge)
        if cationcharge == 1:
            nanion = 0
            ncation = roundnetchargeabs
        elif cationcharge == 2:
            if np.mod(roundnetchargeabs, 2) == 0:
                nanion = 0
                ncation = int(roundnetchargeabs / 2)
            else:
                if anioncharge == -1:
                    nanion = 1
                    ncation = 1 + int(np.floor(roundnetchargeabs / 2))
                else:
                    raise NameError('Failed to neutralize because the cation and anion charges cannot reach 0')
    else:
        logger.info('The system is already neutral.')
        done = True

    if not done:
        newtotalcharge = nanion * anioncharge + ncation * cationcharge + roundnetcharge
        if newtotalcharge != 0:
            raise NameError('Could not neutralize the system')
    if neutralize:
        return anion, cation, nanion, ncation

    # If we want a specific salt concentration, guess chemical formula
    if cationcharge == 1 and anioncharge == -1:  # e.g., NaCl, Kcl, ...
        cationstoich = 1
        anionstoich = 1
    elif cationcharge == 2 and anioncharge == -1:  # e.g., MgCl2
        cationstoich = 1
        anionstoich = 2
    elif cationcharge == 1 and anioncharge == -2:
        cationstoich = 2
        anionstoich = 1
    elif cationcharge == 2 and anioncharge == -2:
        cationstoich = 1
        anionstoich = 1
    else:
        raise NameError('Unsupported ion charge; cannot guess chamical formula.')

    num = int(np.floor(0.5 + 0.0187 * saltconc * nwater))
    logger.info('Adding {} anions + {} cations for neutralizing and {} ions for the given salt concentration.'.format(nanion, ncation, (cationstoich + anionstoich) * num))
    ncation += cationstoich * num
    nanion += anionstoich * num
    return anion, cation, _ionGetAtomname(anion), _ionGetAtomname(cation), nanion, ncation


def _ionGetAtomname(ion):
    if ion == 'ZN2':
        return 'ZN'
    if ion == 'CD2':
        return 'CD'
    return ion


def _ionGetCharge(ion, ff):
    charmmions = dict()
    charmmions['SOD'] = (1, 'sodium', 'Na+')
    charmmions['MG'] = (2, 'magnesium', 'Mg2+')
    charmmions['POT'] = (1, 'potassium', 'K+')
    charmmions['CES'] = (1, 'cesium', 'Cs+')
    charmmions['CAL'] = (2, 'calcium', 'Ca2+')
    charmmions['ZN2'] = (2, 'zinc', 'Zn2+')
    charmmions['CLA'] = (-1, 'chloride', 'Cl-')

    amberions = dict()
    amberions['Na+'] = (1, 'sodium', 'Na+')
    amberions['Cl-'] = (-1, 'chloride', 'Cl-')
    amberions['K+'] = (1, 'potassium', 'K+')
    amberions['Cs+'] = (1, 'cesium', 'Cs+')

    if ff == 'charmm':
        if ion not in charmmions:
            raise NameError('Ion ' + ion + ' is not supported in CHARMM. Supported ions are: ' + str(charmmions.keys()))
        return charmmions[ion][0]
    elif ff == 'amber':
        if ion not in amberions:
            raise NameError('Ion ' + ion + ' is not supported in AMBER. Supported ions are: ' + str(amberions.keys()))
        return amberions[ion][0]


def ionizePlace(mol, anion, cation, anionatom, cationatom, nanion, ncation, dfrom=5, dbetween=5, segname=None):
    newmol = mol.copy()

    logger.info('Min distance of ions from molecule: ' + str(dfrom) + 'A')
    logger.info('Min distance between ions: ' + str(dbetween) + 'A')
    logger.info('Placing ' + str(nanion+ncation) + ' ions.')

    if (nanion + ncation) == 0:
        return newmol

    segname = _getSegname(newmol, segname)
    nions = nanion + ncation

    betabackup = newmol.beta
    newmol.set('beta', sequenceID((newmol.resid, newmol.segid)))

    # Find water oxygens to replace with ions
    ntries = 0
    maxtries = 10
    while True:
        ionlist = np.empty(0, dtype=int)
        watindex = newmol.atomselect('noh and water and not (within ' + str(dfrom) + ' of not water)', indexes=True)
        watsize = len(watindex)

        if watsize == 0:
            raise NameError('No waters could be found further than ' + str(dfrom) + ' from other molecules to be replaced by ions. You might need to solvate with a bigger box.')

        while len(ionlist) < nions:
            if len(watindex) == 0:
                break
            randwat = np.random.randint(len(watindex))
            thision = watindex[randwat]
            addit = True
            if len(ionlist) != 0:  # Check for distance from precious ions
                ionspos = newmol.get('coords', sel=ionlist)
                thispos = newmol.get('coords', sel=thision)
                dists = distance.cdist(np.atleast_2d(ionspos), np.atleast_2d(thispos), metric='euclidean')

                if np.any(dists < dbetween):
                    addit = False
            if addit:
                ionlist = np.append(ionlist, thision)
                watindex = np.delete(watindex, randwat)
        if len(ionlist) == nions:
            break

        ntries += 1
        if ntries == maxtries:
            raise NameError('Failed to add ions after ' + str(maxtries) + ' attempts. Try decreasing the ''from'' and ''between'' parameters, decreasing ion concentration or making a larger water box.')

    # Delete waters but keep their coordinates
    waterpos = np.atleast_2d(newmol.get('coords', ionlist))
    stringsel = 'beta'
    for x in newmol.beta[ionlist]:
        stringsel += ' ' + str(int(x))
    atmrem = np.sum(newmol.atomselect(stringsel))
    atmput = 3 * len(ionlist)
    # assert atmrem == atmput, 'Removing {} atoms instead of {}. Report this bug.'.format(atmrem, atmput)
    sel = newmol.atomselect(stringsel, indexes=True)
    newmol.remove(sel, _logger=False)
    # assert np.size(sel) == atmput, 'Removed {} atoms instead of {}. Report this bug.'.format(np.size(sel), atmput)
    betabackup = np.delete(betabackup, sel)

    # Add the ions
    randidx = np.random.permutation(np.size(waterpos, 0))
    atom = Molecule()
    atom.empty(1)
    atom.set('chain', 'I')
    atom.set('segid', 'I')

    for i in range(nanion):
        atom.name = anionatom
        atom.resname = anion
        atom.resid = newmol.resid[-1] + 1
        atom.coords = waterpos[randidx[i], :]
        newmol.insert(atom, len(newmol.name))
    for i in range(ncation):
        atom.name = cationatom
        atom.resname = cation
        atom.resid = newmol.resid[-1] + 1
        atom.coords = waterpos[randidx[i+nanion], :]
        newmol.insert(atom, len(newmol.name))

    # Restoring the original betas
    sel = np.ones(len(betabackup) + nions, dtype=bool)
    sel[len(betabackup)::] = False
    newmol.set('beta', betabackup, sel=sel)
    return newmol


def _getSegname(mol, segname):
    defsegnames = ['ION','ION1','ION3','ION4','ION5','ION6','ION7','ION8','ION9','IN1','IN2','IN3','IN4','IN5','IN6','IN7','IN8','IN9']

    # Make sure segname is not taken
    if segname is None:
        for i in range(len(defsegnames)):
            sel = mol.atomselect('segname ' + defsegnames[i])
            if not np.any(sel):
                segname = defsegnames[i]
                break
        if segname is None:
            raise NameError('All default segnames taken! Provide your own.')
    else:
        sel = mol.atomselect('segname ' + segname)
        if np.any(sel):
            raise NameError('Segname already taken. Provide a different one.')
