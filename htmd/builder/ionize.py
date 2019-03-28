# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import scipy.spatial.distance as distance
from moleculekit.molecule import Molecule
from moleculekit.util import sequenceID
import logging
logger = logging.getLogger(__name__)

# Info about all ions, with charge, name, amber and charmm names (may
# need residues too). Only the first column is currently used.
_ions = {
    'NA': (1, 'sodium', 'Na+', 'SOD'),
    'MG': (2, 'magnesium', 'Mg2+', 'MG'),
    'ZN': (2, 'zinc', 'Zn2+', 'ZN'),
    'K': (1, 'potassium', 'K+', 'POT'),
    'CS': (1, 'cesium', 'Cs+', 'CES'),
    'CA': (2, 'calcium', 'Ca2+', 'CAL'),
    'CL': (-1, 'chloride', 'Cl-', 'CLA')
}


def ionize(netcharge, nwater, neutralize=True, saltconc=None, cation=None, anion=None, ff=None):
    if ff:
        logger.warning('The ff option is no longer needed and is deprecated. '
                       'Please use standard names (e.g. NA, CL, K) for `cation` and `anion` arguments.')
    
    if cation is None:
            cation = 'NA'
    if anion is None:
            anion = 'CL'

    if saltconc is not None:
        neutralize = False

    cationcharge = _ionGetCharge(cation)
    anioncharge = _ionGetCharge(anion)

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
    return anion, cation, anion, cation, nanion, ncation


def _ionGetCharge(ion):
    if ion not in _ions:
        raise NameError('Ion {:s} not in the database'.format(ion))
    return _ions[ion][0]


def ionizePlace(mol, anion_resname, cation_resname, anion_name, cation_name, nanion, ncation, dfrom=5, dbetween=5, segname=None):
    """Place a given number of negative and positive ions in the solvent.

    Replaces water molecules al long as they respect the given distance criteria.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object
    anion_resname : str
        Resname of the added anions
    cation_resname : str
        Resname of the added cations
    anion_name : str
        Name of the added anions
    cation_name : str
        Name of the added cations
    nanion : int
        Number of anions to add
    ncation : int
        Number of cations to add
    dfrom : float
        Min distance of ions from molecule
    dbetween : float
        Min distance between ions
    segname : str
        Segment name to add
        
    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The molecule with the ions added
    """

    newmol = mol.copy()

    logger.debug('Min distance of ions from molecule: ' + str(dfrom) + 'A')
    logger.debug('Min distance between ions: ' + str(dbetween) + 'A')
    logger.debug('Placing {:d} anions and {:d} cations.'.format(nanion,ncation))

    if (nanion + ncation) == 0:
        return newmol

    nions = nanion + ncation

    betabackup = newmol.beta.copy()
    newmol.set('beta', sequenceID((newmol.resid, newmol.insertion, newmol.segid)))

    # Find water oxygens to replace with ions
    ntries = 0
    maxtries = 10
    while True:
        ionlist = []
        watindex = newmol.atomselect('noh and water and not (within ' + str(dfrom) + ' of not water)', indexes=True)
        watsize = len(watindex)

        if watsize == 0:
            raise NameError('No waters could be found further than ' + str(dfrom) + ' from other molecules to be replaced by ions. You might need to solvate with a bigger box or disable the ionize property when building.')

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
                ionlist.append(thision)
                watindex = np.delete(watindex, randwat)
        if len(ionlist) == nions:
            break

        ntries += 1
        if ntries == maxtries:
            raise NameError('Failed to add ions after ' + str(maxtries) + ' attempts. Try decreasing the ''from'' and ''between'' parameters, decreasing ion concentration or making a larger water box.')

    # Delete waters but keep their coordinates
    waterpos = np.atleast_2d(newmol.get('coords', ionlist))
    betasel = np.zeros(newmol.numAtoms, dtype=bool)
    for b in newmol.beta[ionlist]:
        betasel |= newmol.beta == b
    atmrem = np.sum(betasel)
    atmput = 3 * len(ionlist)
    # assert atmrem == atmput, 'Removing {} atoms instead of {}. Report this bug.'.format(atmrem, atmput)
    sel = np.where(betasel)[0]
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
        atom.set('name', anion_name)
        atom.set('resname', anion_resname)
        atom.set('resid', newmol.resid[-1] + 1)
        atom.coords = waterpos[randidx[i], :]
        newmol.insert(atom, len(newmol.name))
    for i in range(ncation):
        atom.set('name', cation_name)
        atom.set('resname', cation_resname)
        atom.set('resid', newmol.resid[-1] + 1)
        atom.coords = waterpos[randidx[i+nanion], :]
        newmol.insert(atom, len(newmol.name))

    # Restoring the original betas
    newmol.beta[:len(betabackup)] = betabackup
    return newmol


def _getSegname(mol, segname):
    defsegnames = ['ION','ION1','ION3','ION4','ION5','ION6','ION7','ION8','ION9','IN1','IN2','IN3','IN4','IN5','IN6','IN7','IN8','IN9']

    # Make sure segname is not taken
    if segname is None:
        for i in range(len(defsegnames)):
            sel = mol.segid == defsegnames[i]
            if not np.any(sel):
                segname = defsegnames[i]
                break
        if segname is None:
            raise NameError('All default segnames taken! Provide your own.')
    else:
        sel = mol.segid == segname
        if np.any(sel):
            raise NameError('Segname already taken. Provide a different one.')
