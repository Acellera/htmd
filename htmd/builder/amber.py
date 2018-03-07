# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function


from htmd.home import home
import numpy as np
import os
from glob import glob
from htmd.molecule.util import _missingSegID, sequenceID
import shutil
from htmd.builder.builder import detectDisulfideBonds, BuildError, UnknownResidueError, MissingAngleError, MissingBondError, MissingParameterError, MissingTorsionError, MissingAtomTypeError
from htmd.builder.builder import _checkMixedSegment, _checkResidueInsertions
from subprocess import call, check_output, DEVNULL
from htmd.molecule.molecule import Molecule
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.util import ensurelist
from natsort import natsorted
import logging
logger = logging.getLogger(__name__)


def _findTleap():
    tleap = shutil.which("tleap", mode=os.X_OK)
    if not tleap:
        raise FileNotFoundError('tleap not found. You should either have AmberTools or ambermini installed '
                                '(to install ambermini do: conda install ambermini -c acellera)')
    if os.path.islink(tleap):
        if os.path.isabs(os.readlink(tleap)):
            tleap = os.readlink(tleap)
        else:
            tleap = os.path.join(os.path.dirname(tleap), os.readlink(tleap))
    return tleap

def listFiles():
    """ Lists all available AMBER forcefield files

    Example
    -------
    >>> from htmd.builder import amber
    >>> amber.listFiles()             # doctest: +ELLIPSIS
    ---- Forcefield files list: ... ----
    leaprc.constph
    leaprc.DNA.bsc1
    leaprc.DNA.OL15
    leaprc.ff14SB.redq
    ...
    """
    from os.path import join
    tleap = _findTleap()

    amberhome = os.path.normpath(os.path.join(os.path.dirname(tleap), '../'))

    # Original AMBER FFs
    ffdir = join(amberhome, 'dat', 'leap', 'cmd')
    ffs = glob(join(ffdir, '*'))
    print('---- Forcefield files list: ' + join(ffdir, '') + ' ----')
    for f in sorted(ffs, key=str.lower):
        if os.path.isdir(f):
            continue
        print(f.replace(join(ffdir, ''), ''))

    oldffdir = join(amberhome, 'dat', 'leap', 'cmd', 'oldff')
    ffs = glob(join(oldffdir, '*'))
    print('---- OLD Forcefield files list: ' + join(ffdir, '') + ' ----')
    for f in sorted(ffs, key=str.lower):
        print(f.replace(join(ffdir, ''), ''))

    topodir = os.path.join(amberhome, 'dat', 'leap', 'prep')
    topos = glob(join(topodir, '*'))
    print('---- Topology files list: ' + join(topodir, '') + ' ----')
    for f in sorted(topos, key=str.lower):
        if os.path.isdir(f):
            continue
        print(f.replace(join(topodir, ''), ''))

    # FRCMOD files
    frcmoddir = os.path.join(amberhome, 'dat', 'leap', 'parm')
    ffs = glob(join(frcmoddir, 'frcmod.*'))
    print('---- Parameter files list: ' + join(frcmoddir, '') + ' ----')
    for f in sorted(ffs, key=str.lower):
        print(f.replace(join(frcmoddir, ''), ''))

    # Extra AMBER FFs on HTMD
    htmdamberdir = os.path.abspath(os.path.join(home(), 'builder', 'amberfiles', ''))
    extraffs = glob(join(htmdamberdir, '*', 'leaprc.*'))
    print('---- Extra forcefield files list: ' + join(htmdamberdir, '') + ' ----')
    for f in sorted(extraffs, key=str.lower):
        print(f.replace(join(htmdamberdir, ''), ''))

    # Extra AMBER FFs on HTMD (*.frcmod, *.in) @cuzzo87
    extratopos = glob(join(htmdamberdir, '*', '*.in'))
    print('---- Extra topology files list: ' + join(htmdamberdir, '') + ' ----')
    for f in sorted(extratopos, key=str.lower):
        print(f.replace(join(htmdamberdir, ''), ''))

    extraparams = glob(join(htmdamberdir, '*', '*.frcmod'))
    print('---- Extra parameter files list: ' + join(htmdamberdir, '') + ' ----')
    for f in sorted(extraparams, key=str.lower):
        print(f.replace(join(htmdamberdir, ''), ''))


def _locateFile(fname, type, tleap):
    amberhome = os.path.normpath(os.path.join(os.path.dirname(tleap), '../'))
    htmdamberdir = os.path.abspath(os.path.join(home(), 'builder', 'amberfiles', ''))
    if type == 'topo':
        topodir = os.path.join(amberhome, 'dat', 'leap', 'prep')
        foundfile = glob(os.path.join(topodir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        foundfile = glob(os.path.join(htmdamberdir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        logger.warning('Was not able to find topology file {}'.format(fname))
    elif type == 'param':
        parmdir = os.path.join(amberhome, 'dat', 'leap', 'parm')
        foundfile = glob(os.path.join(parmdir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        foundfile = glob(os.path.join(htmdamberdir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        logger.warning('Was not able to find parameter file {}'.format(fname))
    elif type == 'ff':
        ffdir = os.path.join(amberhome, 'dat', 'leap', 'cmd')
        foundfile = glob(os.path.join(ffdir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        foundfile = glob(os.path.join(htmdamberdir, fname))
        if len(foundfile) != 0:
            return foundfile[0]
        logger.warning('Was not able to find forcefield file {}'.format(fname))


def defaultFf():
    """ Returns the default leaprc forcefield files used by amber.build """
    return ['leaprc.lipid14', os.path.join('oldff', 'leaprc.ff14SB'), 'leaprc.gaff']


def defaultTopo():
    """ Returns the default topology `prepi` files used by amber.build """
    return []


def defaultParam():
    """ Returns the default parameter `frcmod` files used by amber.build """

    # Common choices:
    #  * frcmod.ionsjc_tip3p:  Monovalent ion parameters for Ewald and TIP3P water from Joung & Cheatham JPCB (2008)
    #  * frcmod.ions1lm_126_tip3p + frcmod.ions234lm_126_tip3p :
    #       Li/Merz ion parameters of monovalent ions for TIP3P water model (12-6 normal usage set)
    #       Li/Merz ion parameters of divalent to tetravalent ions for TIP3P water model (12-6 normal usage set)
    #
    # See page 50 of Amber16 manual.

    return ['frcmod.ionsjc_tip3p', 'frcmod.ions234lm_126_tip3p']


def build(mol, ff=None, topo=None, param=None, prefix='structure', outdir='./build', caps=None, ionize=True, saltconc=0,
          saltanion=None, saltcation=None, disulfide=None, tleap=None, execute=True, atomtypes=None,
          offlibraries=None, gbsa=False, igb=2):
    """ Builds a system for AMBER

    Uses tleap to build a system for AMBER. Additionally it allows the user to ionize and add disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object containing the system
    ff : list of str
        A list of leaprc forcefield files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available forcefield files.
        Default: :func:`amber.defaultFf <htmd.builder.amber.defaultFf>`
    topo : list of str
        A list of topology `prepi/prep/in` files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available topology files.
        Default: :func:`amber.defaultTopo <htmd.builder.amber.defaultTopo>`
    param : list of str
        A list of parameter `frcmod` files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available parameter files.
        Default: :func:`amber.defaultParam <htmd.builder.amber.defaultParam>`
    prefix : str
        The prefix for the generated pdb and psf files
    outdir : str
        The path to the output directory
        Default: './build'
    caps : dict
        A dictionary with keys segids and values lists of strings describing the caps for a particular protein segment.
        e.g. caps['P'] = ['ACE', 'NME'] or caps['P'] = ['none', 'none']. Default: will apply ACE and NME caps to every
        protein segment.
    ionize : bool
        Enable or disable ionization
    saltconc : float
        Salt concentration to add to the system after neutralization.
    saltanion : {'Cl-'}
        The anion type. Please use only AMBER ion atom names.
    saltcation : {'Na+', 'K+', 'Cs+'}
        The cation type. Please use only AMBER ion atom names.
    disulfide : list of :class:`DisulfideBridge <htmd.builder.builder.DisulfideBridge>` objects
        If None it will guess disulfide bonds. Otherwise provide a list of `DisulfideBridge` objects.
    tleap : str
        Path to tleap executable used to build the system for AMBER
    execute : bool
        Disable building. Will only write out the input script needed by tleap. Does not include ionization.
    atomtypes : list of triplets
        Custom atom types defined by the user as ('type', 'element', 'hybrid') triplets
        e.g. (('C1', 'C', 'sp2'), ('CI', 'C', 'sp3')). Check `addAtomTypes` in AmberTools docs.
    offlibraries : str or list
        A path or a list of paths to OFF library files. Check `loadOFF` in AmberTools docs.
    gbsa : bool
        Modify radii for GBSA implicit water model
    igb : int
        GB model. Select: 1 for mbondi, 2 and 5 for mbondi2, 7 for bondi and 8 for mbondi3.
        Check section 4. The Generalized Born/Surface Area Model of the AMBER manual.

    Returns
    -------
    molbuilt : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The built system in a Molecule object

    Example
    -------
    >>> from htmd.ui import *
    >>> mol = Molecule("3PTB")
    >>> molbuilt = amber.build(mol, outdir='/tmp/build')  # doctest: +SKIP
    >>> # More complex example
    >>> disu = [DisulfideBridge('P', 157, 'P', 13), DisulfideBridge('K', 1, 'K', 25)]
    >>> molbuilt = amber.build(mol, outdir='/tmp/build', saltconc=0.15, disulfide=disu)  # doctest: +SKIP
    """
    # Remove pdb protein bonds as they can be regenerated by tleap. Keep non-protein bonds i.e. for ligands
    mol = mol.copy()
    _removeProteinBonds(mol)

    if tleap is None:
        tleap = _findTleap()
    else:
        if shutil.which(tleap) is None:
            raise NameError('Could not find executable: `{}` in the PATH. Cannot build for AMBER.'.format(tleap))

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    _cleanOutDir(outdir)
    if ff is None:
        ff = defaultFf()
    if topo is None:
        topo = defaultTopo()
    if param is None:
        param = defaultParam()
    if caps is None:
        caps = _defaultProteinCaps(mol)

    _missingSegID(mol)
    _checkMixedSegment(mol)
    _checkResidueInsertions(mol)

    mol = _charmmLipid2Amber(mol)

    _applyProteinCaps(mol, caps)

    f = open(os.path.join(outdir, 'tleap.in'), 'w')
    f.write('# tleap file generated by amber.build\n')

    # Printing out the forcefields
    if isinstance(ff, str):
        ff = [ff]
    for i, force in enumerate(ff):
        if os.path.isfile(force):
            f.write('source {}\n'.format(force))
        else: # Search amber directories for the file
            p = _locateFile(force, 'ff', tleap)
            if p is not None:
                newname = 'ff{}_{}'.format(i, os.path.basename(force))
                shutil.copy(p, os.path.join(outdir, newname))
                f.write('source {}\n'.format(newname))
    f.write('\n')

    if gbsa:
        gbmodels = {1: 'mbondi', 2: 'mbondi2', 5: 'mbondi2', 7: 'bondi', 8: 'mbondi3'}
        f.write('set default PBradii {}\n\n'.format(gbmodels[igb]))

    # Adding custom atom types
    if atomtypes is not None:
        atomtypes = ensurelist(tocheck=atomtypes[0], tomod=atomtypes)
        f.write('addAtomTypes {\n')
        for at in atomtypes:
            if len(at) != 3:
                raise RuntimeError('Atom type definitions have to be triplets. Check the AMBER documentation.')
            f.write('    {{ "{}" "{}" "{}" }}\n'.format(at[0], at[1], at[2]))
        f.write('}\n\n')

    # Loading OFF libraries
    if offlibraries is not None:
        if not isinstance(offlibraries, list) and not isinstance(offlibraries, tuple):
            offlibraries = [offlibraries, ]
        for off in offlibraries:
            f.write('loadoff {}\n\n'.format(off))

    # Loading frcmod parameters
    f.write('# Loading parameter files\n')
    for i, p in enumerate(param):
        if os.path.isfile(p):
            f.write('loadamberparams {}\n'.format(p))
        else: # Search amber directories for the file
            p = _locateFile(p, 'param', tleap)
            if p is not None:
                newname = 'param{}_{}'.format(i, os.path.basename(p))
                shutil.copy(p, os.path.join(outdir, newname))
                f.write('loadamberparams {}\n'.format(newname))
    f.write('\n')

    # Loading prepi topologies
    f.write('# Loading prepi topologies\n')
    for i, t in enumerate(topo):
        if os.path.isfile(t):
            f.write('loadamberprep {}\n'.format(t))
        else: # Search amber directories for the file
            t = _locateFile(t, 'topo', tleap)
            if t is not None:
                newname = 'topo{}_{}'.format(i, os.path.basename(t))
                shutil.copy(t, os.path.join(outdir, newname))
                f.write('loadamberprep {}\n'.format(newname))
    f.write('\n')

    # Detect disulfide bridges if not defined by user
    if disulfide is None and not ionize:
        logger.info('Detecting disulfide bonds.')
        disulfide = detectDisulfideBonds(mol)

    # Fix structure to match the disulfide patching
    if not ionize and len(disulfide) != 0:
        for d in disulfide:
            # Rename the residues to CYX if there is a disulfide bond
            atoms1 = (mol.segid == d.segid1) & (mol.resid == d.resid1)
            atoms2 = (mol.segid == d.segid2) & (mol.resid == d.resid2)
            mol.resname[atoms1] = 'CYX'
            mol.resname[atoms2] = 'CYX'
            # Remove (eventual) HG hydrogens on these CYS (from proteinPrepare)
            mol.remove((atoms1 & (mol.name == 'HG')) | (atoms2 & (mol.name == 'HG')), _logger=False)

    # Printing and loading the PDB file. AMBER can work with a single PDB file if the segments are separate by TER
    logger.debug('Writing PDB file for input to tleap.')
    pdbname = os.path.join(outdir, 'input.pdb')

    # mol2 files have atomtype, here we only write parts not coming from mol2
    mol.write(pdbname, mol.atomtype == '')
    if not os.path.isfile(pdbname):
        raise NameError('Could not write a PDB file out of the given Molecule.')
    f.write('# Loading the system\n')
    f.write('mol = loadpdb input.pdb\n\n')

    if np.sum(mol.atomtype != '') != 0:
        logger.debug('Writing mol2 files for input to tleap.')
        segs = np.unique(mol.segid[mol.atomtype != ''])
        combstr = 'mol = combine {mol'
        for s in segs:
            name = 'segment{}'.format(s)
            mol2name = os.path.join(outdir, '{}.mol2'.format(name))
            mol.write(mol2name, (mol.atomtype != '') & (mol.segid == s))
            if not os.path.isfile(mol2name):
                raise NameError('Could not write a mol2 file out of the given Molecule.')
            f.write('# Loading the rest of the system\n')
            f.write('{} = loadmol2 {}.mol2\n\n'.format(name, name))
            combstr += ' {}'.format(name)
        combstr += '}\n\n'
        f.write(combstr)

    # Write patches for disulfide bonds (only after ionizing)
    if not ionize and len(disulfide) != 0:
        f.write('# Adding disulfide bonds\n')
        for d in disulfide:
            # Convert to stupid amber residue numbering
            uqseqid = sequenceID((mol.resid, mol.insertion, mol.segid)) + mol.resid[0]
            uqres1 = int(np.unique(uqseqid[(mol.segid == d.segid1) & (mol.resid == d.resid1)]))
            uqres2 = int(np.unique(uqseqid[(mol.segid == d.segid2) & (mol.resid == d.resid2)]))
            f.write('bond mol.{}.SG mol.{}.SG\n'.format(uqres1, uqres2))
        f.write('\n')

    f.write('# Writing out the results\n')
    f.write('saveamberparm mol ' + prefix + '.prmtop ' + prefix + '.crd\n')
    f.write('quit')
    f.close()

    molbuilt = None
    if execute:
        # Source paths of extra dirs (our dirs, not amber default)
        htmdamberdir = os.path.abspath(os.path.join(home(), 'builder', 'amberfiles'))
        sourcepaths = [htmdamberdir]
        sourcepaths += [os.path.join(htmdamberdir, os.path.dirname(f))
                        for f in ff if os.path.isfile(os.path.join(htmdamberdir, f))]
        extrasource = []
        for p in sourcepaths:
            extrasource.append('-I')
            extrasource.append('{}'.format(p))
        logpath = os.path.abspath(os.path.join(outdir, 'log.txt'))
        logger.info('Starting the build.')
        currdir = os.getcwd()
        os.chdir(outdir)
        f = open(logpath, 'w')
        try:
            cmd = [tleap, '-f', './tleap.in']
            cmd[1:1] = extrasource
            call(cmd, stdout=f)
        except:
            raise NameError('tleap failed at execution')
        f.close()
        errors = _logParser(logpath)
        os.chdir(currdir)
        if errors:
            raise BuildError(errors + ['Check {} for further information on errors in building.'.format(logpath)])
        logger.info('Finished building.')

        if os.path.exists(os.path.join(outdir, 'structure.crd')) and \
                        os.path.getsize(os.path.join(outdir, 'structure.crd')) != 0 and \
                        os.path.getsize(os.path.join(outdir, 'structure.prmtop')) != 0:
            molbuilt = Molecule(os.path.join(outdir, 'structure.prmtop'))
            molbuilt.read(os.path.join(outdir, 'structure.crd'))
        else:
            raise BuildError('No structure pdb/prmtop file was generated. Check {} for errors in building.'.format(logpath))

        if ionize:
            shutil.move(os.path.join(outdir, 'structure.crd'), os.path.join(outdir, 'structure.noions.crd'))
            shutil.move(os.path.join(outdir, 'structure.prmtop'), os.path.join(outdir, 'structure.noions.prmtop'))
            totalcharge = np.sum(molbuilt.charge)
            nwater = np.sum(molbuilt.atomselect('water and noh'))
            anion, cation, anionatom, cationatom, nanion, ncation = ionizef(totalcharge, nwater, saltconc=saltconc, ff='amber', anion=saltanion, cation=saltcation)
            newmol = ionizePlace(mol, anion, cation, anionatom, cationatom, nanion, ncation)
            # Redo the whole build but now with ions included
            return build(newmol, ff=ff, topo=topo, param=param, prefix=prefix, outdir=outdir, caps={}, ionize=False,
                         execute=execute, saltconc=saltconc, disulfide=disulfide, tleap=tleap, atomtypes=atomtypes,
                         offlibraries=offlibraries)
    tmpbonds = molbuilt.bonds
    molbuilt.bonds = []  # Removing the bonds to speed up writing
    molbuilt.write(os.path.join(outdir, 'structure.pdb'))
    molbuilt.bonds = tmpbonds  # Restoring the bonds
    return molbuilt


def _applyProteinCaps(mol, caps):

    # AMBER capping
    # =============
    # This is the (horrible) way of adding caps in tleap:
    # For now, this is hardwired for ACE and NME
    # 1. Change one of the hydrogens of the N terminal (H[T]?[123]) to the ACE C atom, giving it a new resid
    # 1a. If no hydrogen present, create the ACE C atom.
    # 2. Change one of the oxygens of the C terminal ({O,OT1,OXT}) to the NME N atom, giving it a new resid
    # 2a. If no oxygen present, create the NME N atom.
    # 3. Reorder to put the new atoms first and last
    # 4. Remove the lingering hydrogens of the N terminal and oxygens of the C terminal.

    # Define the atoms to be replaced (0 and 1 corresponds to N- and C-terminal caps)
    terminalatoms = {'ACE': ['H1', 'H2', 'H3', 'HT1', 'HT2', 'HT3'], 'NME': ['OXT', 'OT1', 'O']}  # XPLOR names for H[123] and OXT are HT[123]
                                                                         # and OT1, respectively.
    capresname = ['ACE', 'NME']
    capatomtype = ['C', 'N']

    # For each caps definition
    for seg in caps:
        prot = mol.atomselect('protein')  # Can't move this out since we remove atoms in this loop
        # Get the segment
        segment = np.where(mol.segid == seg)[0]
        # Test segment
        if len(segment) == 0:
            raise RuntimeError('There is no segment {} in the molecule.'.format(seg))
        if not np.any(prot & (mol.segid == seg)):
            raise RuntimeError('Segment {} is not protein. Capping for non-protein segments is not supported.'.format(seg))
        # For each cap
        passed = False
        for i, cap in enumerate(caps[seg]):
            if cap is None or (isinstance(cap, str) and cap == 'none'):
                continue
            # Get info on segment and its terminals
            segidm = mol.segid == seg  # Mask for segid
            segididx = np.where(segidm)[0]
            resids = mol.resid[segididx]
            terminalids = [segididx[0], segididx[-1]]
            terminalresids = [resids[0], resids[-1]]
            residm = mol.resid == terminalresids[i]  # Mask for resid

            if not passed:
                orig_terminalresids = terminalresids
                passed = True

            if cap is None or cap == '':  # In case there is no cap defined
                logger.warning('No cap provided for resid {} on segment {}. Did not apply it.'.format(terminalresids[i], seg))
                continue
            elif cap not in capresname:  # If it is defined, test if supported
                raise RuntimeError('In segment {}, the {} cap is not supported. Try using {} instead.'.format(seg, cap, capresname))

            # Test if cap is already applied
            testcap = np.where(segidm & residm & (mol.resname == cap))[0]
            if len(testcap) != 0:
                logger.warning('Cap {} already exists on segment {}. Did not re-apply it.'.format(cap, seg))
                continue

            # Test if the atom to change exists
            termatomsids = np.zeros(residm.shape, dtype=bool)
            for atm in terminalatoms[cap]:
                termatomsids |= mol.name == atm
            termatomsids = np.where(termatomsids & segidm & residm)[0]

            if len(termatomsids) == 0:
                # Create new atom
                termcaid = np.where(segidm & residm & (mol.name == 'CA'))[0]
                termcenterid = np.where(segidm & residm & (mol.name == capatomtype[1-i]))[0]
                atom = Molecule()
                atom.empty(1)
                atom.record = np.array(['ATOM'], dtype=Molecule._dtypes['record'])
                atom.name = np.array([capatomtype[i]], dtype=Molecule._dtypes['name'])
                atom.resid = np.array([terminalresids[i]-1+2*i], dtype=Molecule._dtypes['resid'])
                atom.resname = np.array([cap], dtype=Molecule._dtypes['resname'])
                atom.segid = np.array([seg], dtype=Molecule._dtypes['segid'])
                atom.element = np.array([capatomtype[i]], dtype=Molecule._dtypes['element'])
                atom.chain = np.array([np.unique(mol.chain[segidm])], dtype=Molecule._dtypes['chain'])  # TODO: Assumption of single chain in a segment might be wrong
                atom.coords = mol.coords[termcenterid] + 0.33 * np.subtract(mol.coords[termcenterid],
                                                                            mol.coords[termcaid])
                mol.insert(atom, terminalids[i])
                # logger.info('In segment {}, resid {} had none of these atoms: {}. Capping was performed by creating '
                #             'a new atom for cap construction by tleap.'.format(seg, terminalresids[i],
                #                                                                ' '.join(terminalatoms[cap])))

            else:
                # Select atom to change, do changes to cap, and change resid
                newatom = np.max(termatomsids)
                mol.set('resname', cap, sel=newatom)
                mol.set('name', capatomtype[i], sel=newatom)
                mol.set('element', capatomtype[i], sel=newatom)
                mol.set('resid', terminalresids[i]-1+2*i, sel=newatom)  # if i=0 => resid-1; i=1 => resid+1

                # Reorder
                neworder = np.arange(mol.numAtoms)
                neworder[newatom] = terminalids[i]
                neworder[terminalids[i]] = newatom
                _reorderMol(mol, neworder)

        # For each cap
        for i, cap in enumerate(caps[seg]):
            if cap is None or (isinstance(cap, str) and cap == 'none'):
                continue
            # Remove lingering hydrogens or oxygens in terminals
            mol.remove('segid {} and resid "{}" and name {}'.format(seg, orig_terminalresids[i], ' '.join(terminalatoms[cap])),
                       _logger=False)


def _removeProteinBonds(mol):
    segs = np.unique(mol.segid[mol.atomtype != ''])  # Keeping bonds related to mol2 files
    mol.deleteBonds(~np.in1d(mol.segid, segs))


def _defaultProteinCaps(mol):
    # Defines ACE and NME (neutral terminals) as default for protein segments
    # Of course, this might not be ideal for proteins that require charged terminals

    segsProt = np.unique(mol.get('segid', sel='protein'))
    caps = dict()
    for s in segsProt:
        if len(np.unique(mol.resid[mol.segid == s])) < 10:
            logger.warning('Segment {} consists of a peptide with less than 10 residues. It will not be capped by '
                           'default. If you want to cap it use the caps argument of amber.build to manually define caps'
                           'for all segments.'.format(s))
            continue
        caps[s] = ['ACE', 'NME']
    return caps


def _cleanOutDir(outdir):
    from glob import glob
    files = glob(os.path.join(outdir, 'structure.*'))
    files += glob(os.path.join(outdir, 'log.*'))
    files += glob(os.path.join(outdir, '*.log'))
    for f in files:
        os.remove(f)


def _charmmLipid2Amber(mol):
    """ Convert a CHARMM lipid membrane to AMBER format

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object containing the membrane

    Returns
    -------
    newmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A new Molecule object with the membrane converted to AMBER
    """

    resdict = _readcsvdict(os.path.join(home(), 'builder', 'charmmlipid2amber.csv'))

    natoms = mol.numAtoms
    neworder = np.array(list(range(natoms)))  # After renaming the atoms and residues I have to reorder them

    begs = np.zeros(natoms, dtype=bool)
    fins = np.zeros(natoms, dtype=bool)
    begters = np.zeros(natoms, dtype=bool)
    finters = np.zeros(natoms, dtype=bool)

    # Iterate over the translation dictionary
    mol = mol.copy()
    incrresids = sequenceID((mol.resid, mol.insertion, mol.segid))

    for res in resdict.keys():
        molresidx = mol.resname == res
        if not np.any(molresidx):
            continue
        names = mol.name.copy()  # Need to make a copy or I accidentally double-modify atoms

        atommap = resdict[res]
        for atom in atommap.keys():
            rule = atommap[atom]

            molatomidx = np.zeros(len(names), dtype=bool)
            molatomidx[molresidx] = names[molresidx] == atom

            mol.set('resname', rule.replaceresname, sel=molatomidx)
            mol.set('name', rule.replaceatom, sel=molatomidx)
            neworder[molatomidx] = rule.order

            if rule.order == 0:  # First atom (with or without ters)
                begs[molatomidx] = True
            if rule.order == rule.natoms - 1:  # Last atom (with or without ters)
                fins[molatomidx] = True
            if rule.order == 0 and rule.ter:  # First atom with ter
                begters[molatomidx] = True
            if rule.order == rule.natoms - 1 and rule.ter:  # Last atom with ter
                finters[molatomidx] = True

    uqresids = np.unique(incrresids[begs])
    residuebegs = np.ones(len(uqresids), dtype=int) * -1
    residuefins = np.ones(len(uqresids), dtype=int) * -1
    for i in range(len(uqresids)):
        residuebegs[i] = np.where(incrresids == uqresids[i])[0][0]
        residuefins[i] = np.where(incrresids == uqresids[i])[0][-1]
    for i in range(len(residuebegs)):
        beg = residuebegs[i]
        fin = residuefins[i] + 1
        neworder[beg:fin] = neworder[beg:fin] + beg
    idx = np.argsort(neworder)

    _reorderMol(mol, idx)

    begters = np.where(begters[idx])[0]  # Sort the begs and ters
    finters = np.where(finters[idx])[0]

    #if len(begters) > 999:
    #    raise NameError('More than 999 lipids. Cannot define separate segments for all of them.')

    for i in range(len(begters)):
        map = np.zeros(len(mol.resid), dtype=bool)
        map[begters[i]:finters[i]+1] = True
        mol.set('resid', sequenceID(mol.get('resname', sel=map)), sel=map)
        mol.set('segid', 'L{}'.format(i % 2), sel=map)

    return mol


def _reorderMol(mol, order):
    for k in mol._atom_fields:
        if mol.__dict__[k] is not None and np.size(mol.__dict__[k]) != 0:
            if k == 'coords':
                mol.__dict__[k] = mol.__dict__[k][order, :, :]
            else:
                mol.__dict__[k] = mol.__dict__[k][order]


def _readcsvdict(filename):
    import csv
    from collections import namedtuple
    if os.path.isfile(filename):
        csvfile = open(filename, 'r')
    else:
        raise NameError('File ' + filename + ' does not exist')

    resdict = dict()

    Rule = namedtuple('Rule', ['replaceresname', 'replaceatom', 'order', 'natoms', 'ter'])

    # Skip header line of csv file. Line 2 contains dictionary keys:
    csvfile.readline()
    csvreader = csv.DictReader(csvfile)
    for line in csvreader:
        searchres = line['search'].split()[1]
        searchatm = line['search'].split()[0]
        if searchres not in resdict:
            resdict[searchres] = dict()
        resdict[searchres][searchatm] = Rule(line['replace'].split()[1], line['replace'].split()[0], int(line['order']), int(line['num_atom']), line['TER'] == 'True')
    csvfile.close()

    return resdict


def _logParser(fname):
    import re
    unknownres_regex = re.compile('Unknown residue:\s+(\w+)')
    missingparam_regex = re.compile('For atom: (.*) Could not find vdW \(or other\) parameters for type:\s+(\w+)')
    missingtorsion_regex = re.compile('No torsion terms for\s+(.*)$')
    missingbond_regex = re.compile('Could not find bond parameter for:\s+(.*)$')
    missingangle_regex = re.compile('Could not find angle parameter:\s+(.*)$')
    missingatomtype_regex = re.compile('FATAL:\s+Atom (.*) does not have a type')

    unknownres = []
    missingparam = []
    missingtorsion = []
    missingbond = []
    missingangle = []
    missingatomtype = []
    with open(fname, 'r') as f:
        for line in f:
            if unknownres_regex.search(line):
                unknownres.append(unknownres_regex.findall(line)[0])
            if missingparam_regex.search(line):
                missingparam.append(missingparam_regex.findall(line)[0])
            if missingtorsion_regex.search(line):
                missingtorsion.append(missingtorsion_regex.findall(line)[0])
            if missingbond_regex.search(line):
                missingbond.append(missingbond_regex.findall(line)[0])
            if missingangle_regex.search(line):
                missingangle.append(missingangle_regex.findall(line)[0])
            if missingatomtype_regex.search(line):
                missingatomtype.append(missingatomtype_regex.findall(line)[0])

    errors = []
    if len(unknownres):
        errors.append(UnknownResidueError('Unknown residue(s) {} found in the input structure. '
                                          'You are either missing a topology definition for the residue or you need to '
                                          'rename it to the correct residue name'.format(np.unique(unknownres))))
    if len(missingparam):
        errors.append(MissingParameterError('Missing parameters for atom {} and type {}'.format(missingparam, missingparam)))
    if len(missingtorsion):
        errors.append(MissingTorsionError('Missing torsion terms for {}'.format(missingtorsion)))
    if len(missingbond):
        errors.append(MissingBondError('Missing bond parameters for {}'.format(missingbond)))
    if len(missingangle):
        errors.append(MissingAngleError('Missing angle parameters for {}'.format(missingangle)))
    if len(missingatomtype):
        errors.append(MissingAtomTypeError('Missing atom type for {}'.format(missingatomtype)))

    return errors

if __name__ == '__main__':
    from htmd.molecule.molecule import Molecule, mol_equal
    from htmd.builder.solvate import solvate
    from htmd.builder.preparation import proteinPrepare
    from htmd.builder.builder import autoSegment2
    from htmd.home import home
    from htmd.util import tempname
    import os
    from os.path import join
    from glob import glob
    import numpy as np
    import filecmp
    import doctest

    failure_count, _ = doctest.testmod()
    # if failure_count != 0:
    #     raise Exception('Doctests failed')

    def _cutfirstline(infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(outfile, 'w') as fout:
            fout.writelines(data[1:])


    def _compareResultFolders(compare, tmpdir, pid):
        ignore_ftypes = ('.log', '.txt')
        files = []
        deletefiles = []
        for f in glob(join(compare, '*')):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith('prmtop'):
                _cutfirstline(f, join(compare, fname + '.mod'))
                _cutfirstline(join(tmpdir, fname), os.path.join(tmpdir, fname + '.mod'))
                files.append(os.path.basename(f) + '.mod')
                deletefiles.append(join(compare, fname + '.mod'))
            else:
                files.append(os.path.basename(f))

        match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError(
                'Different results produced by amber.build for test {} between {} and {} in files {}.'.format(pid, compare, tmpdir, mismatch))

        for f in deletefiles:
            os.remove(f)

    # Test with proteinPrepare
    pdbids = ['3PTB']
    # pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)
    for pid in pdbids:
        np.random.seed(1)
        mol = Molecule(pid)
        mol.filter('protein')
        if mol._checkInsertions():
            mol.renumberResidues()
        mol = proteinPrepare(mol)
        mol.filter('protein')  # Fix for bad proteinPrepare hydrogen placing
        smol = solvate(mol)
        ffs = defaultFf()
        tmpdir = tempname()
        bmol = build(smol, ff=ffs, outdir=tmpdir)

        refdir = home(dataDir=join('test-amber-build', pid))
        _compareResultFolders(refdir, tmpdir, pid)
        shutil.rmtree(tmpdir)

    # Test without proteinPrepare
    pdbids = ['3PTB']
    # pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
    for pid in pdbids:
        np.random.seed(1)
        mol = Molecule(pid)
        mol.filter('protein')
        if mol._checkInsertions():
            mol.renumberResidues()
        smol = solvate(mol)
        ffs = defaultFf()
        tmpdir = tempname()
        bmol = build(smol, ff=ffs, outdir=tmpdir)

        refdir = home(dataDir=join('test-amber-build-nopp', pid))
        _compareResultFolders(refdir, tmpdir, pid)
        shutil.rmtree(tmpdir)

    # Test protein ligand building with parametrized ligand
    mol = Molecule('3ptb')
    mol.filter('protein')
    mol.renumberResidues()
    lig = Molecule(join(home(dataDir='test-param'), 'h2o2_gaff2', 'parameters', 'GAFF2', 'B3LYP-cc-pVDZ-vacuum', 'mol.mol2'))
    lig.segid[:] = 'L'
    newmol = Molecule()
    newmol.append(lig)
    newmol.append(mol)
    smol = solvate(newmol)
    tmpdir = tempname()
    bmol = build(newmol, outdir=tmpdir, ionize=False)
    shutil.rmtree(tmpdir)

    # # Test protein-ligand building
    # folder = home(dataDir='building-protein-ligand')
    # prot = Molecule(os.path.join(folder, 'trypsin.pdb'))
    # prot.filter('protein')
    # prot.renumberResidues()
    # prot = autoSegment2(prot)
    # prot = proteinPrepare(prot)
    # prot1 = prot
    # prot2 = prot.copy()
    # lig1 = Molecule(os.path.join(folder, 'benzamidine.mol2'))
    # lig1.set('segid', 'L')
    # lig2 = Molecule(os.path.join(folder, 'benzamidine.pdb'))
    # lig2.set('segid', 'L')
    # prot1.append(lig1)
    # prot2.append(lig2)
    # smol1 = solvate(prot1)
    # smol2 = solvate(prot2)
    # tmpdir1 = tempname()
    # tmpdir2 = tempname()
    # np.random.seed(1)
    # params = amber.defaultParam() + [os.path.join(folder, 'benzamidine.frcmod')]
    # bmol1 = amber.build(smol1, param=params, outdir=tmpdir1)
    # np.random.seed(1)
    # params = amber.defaultParam() + [os.path.join(folder, 'benzamidineprepi.frcmod')]
    # topos = amber.defaultTopo() + [os.path.join(folder, 'benzamidine.prepi')]
    # bmol2 = amber.build(smol2, topo=topos, param=params, outdir=tmpdir2)
    # _compareResultFolders(tmpdir1, tmpdir2, 'ben-tryp')
    # shutil.rmtree(tmpdir1)
    # shutil.rmtree(tmpdir2)

