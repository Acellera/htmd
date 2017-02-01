# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function


from htmd.home import home
import numpy as np
import os
import os.path as path
from glob import glob
from htmd.molecule.util import _missingSegID, sequenceID
import shutil
from htmd.builder.builder import detectDisulfideBonds
from htmd.builder.builder import _checkMixedSegment
from subprocess import call, check_output, DEVNULL
from htmd.molecule.molecule import Molecule
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.util import ensurelist
import logging
logger = logging.getLogger(__name__)


def listFiles():
    """ Lists all available AMBER forcefield files
    """
    tleap = shutil.which("tleap")
    if not tleap:
        raise NameError('tleap not found. You should either have AmberTools or ambermini installed '
                        '(to install ambermini do: conda install ambermini -c acellera)')
    if os.path.islink(tleap):
        if path.isabs(os.readlink(tleap)):
            tleap = os.readlink(tleap)
        else:
            tleap = os.path.join(os.path.dirname(tleap), os.readlink(tleap))

    amberhome = path.normpath(path.join(path.dirname(tleap), '../'))

    # Original AMBER FFs
    amberdir = path.join(amberhome, 'dat', 'leap', 'cmd')
    ffs = [f for f in os.listdir(amberdir) if path.isfile(path.join(amberdir, f))]
    print('---- Forcefield files list: ' + path.join(amberdir, '') + ' ----')
    for f in ffs:
        print(f)

    amberdir = path.join(amberhome, 'dat', 'leap', 'cmd', 'oldff')
    ffs = [f for f in os.listdir(amberdir) if path.isfile(path.join(amberdir, f))]
    print('---- OLD Forcefield files list: ' + path.join(amberdir, '') + ' ----')
    for f in ffs:
        print(f)

    # FRCMOD files
    frcmoddir = path.join(amberhome, 'dat', 'leap', 'parm')
    ffs = glob(frcmoddir+"/frcmod.*")
    print('---- Parameter files list: ' + path.join(frcmoddir, '') + ' ----')
    for f in ffs:
        print(path.basename(f))

    # Extra AMBER FFs on HTMD
    htmdamberdir = path.abspath(path.join(home(), 'builder', 'amberfiles', ''))
    extraffs = [f + '/' + path.basename(glob(os.path.join(htmdamberdir, f) + '/leaprc.*')[0])
                for f in os.listdir(htmdamberdir) if os.path.isdir(os.path.join(htmdamberdir, f))
                and len(glob(os.path.join(htmdamberdir, f) + '/leaprc.*')) == 1]
    print('---- Extra forcefield files list: ' + path.join(htmdamberdir, '') + ' ----')
    for f in extraffs:
        print(f)



def build(mol, ff=None, topo=None, param=None, prefix='structure', outdir='./build', caps=None, ionize=True, saltconc=0,
          saltanion=None, saltcation=None, disulfide=None, tleap='tleap', execute=True, atomtypes=None,
          offlibraries=None):
    """ Builds a system for AMBER

    Uses tleap to build a system for AMBER. Additionally it allows the user to ionize and add disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object containing the system
    ff : list of str
        A list of leaprc forcefield files. Default: ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']
    topo : list of str
        A list of topology `prepi` files.
    param : list of str
        A list of parameter `frcmod` files. Default: ['frcmod.ionsjc_tip3p',]
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

    Returns
    -------
    molbuilt : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The built system in a Molecule object

    Example
    -------
    >>> from htmd import *
    >>> mol = Molecule("3PTB")
    >>> amber.listFiles()             # doctest: +ELLIPSIS
    ---- Topologies files list...
    top/top_all36_prot.rtf
    top/top_water_ions.rtf
    ...
    >>> molbuilt = amber.build(mol, outdir='/tmp/build')  # doctest: +SKIP
    ...
    >>> # More complex example
    >>> ffs = ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']
    >>> params = ['frcmod.ionsjc_tip3p',]
    >>> disu = [DisulfideBridge('P', 157, 'P', 13), DisulfideBridge('K', 1, 'K', 25)]
    >>> molbuilt = amber.build(mol, ff=ffs, param=params, outdir='/tmp/build', saltconc=0.15, disulfide=disu)  # doctest: +SKIP
    """
    # Remove pdb bonds!
    mol = mol.copy()
    mol.bonds = np.empty((0, 2), dtype=np.uint32)
    if shutil.which(tleap) is None:
        raise NameError('Could not find executable: `' + tleap + '` in the PATH. Cannot build for AMBER.')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    _cleanOutDir(outdir)
    if ff is None:
        ff = ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']
    if topo is None:
        topo = []
    if param is None:
        param = ['frcmod.ionsjc_tip3p',]
    if caps is None:
        caps = _defaultProteinCaps(mol)

    _missingSegID(mol)
    _checkMixedSegment(mol)

    logger.info('Converting CHARMM membranes to AMBER.')
    mol = _charmmLipid2Amber(mol)

    #_checkProteinGaps(mol)
    _applyProteinCaps(mol, caps)

    f = open(path.join(outdir, 'tleap.in'), 'w')
    f.write('# tleap file generated by amber.build\n')

    # Printing out the forcefields
    if isinstance(ff, str):
        ff = [ff]
    for force in ff:
        f.write('source ' + force + '\n')
    f.write('\n')

    # Adding custom atom types
    if atomtypes is not None:
        atomtypes = ensurelist(tocheck=atomtypes[0], tomod=atomtypes)
        f.write('addAtomTypes {\n')
        for at in atomtypes:
            if len(at) != 3:
                raise RuntimeError('Atom type definitions have to be triplets. Check the AMBER documentations.')
            f.write('    {{ "{}" "{}" "{}" }}\n'.format(at[0], at[1], at[2]))
        f.write('}\n\n')

    # Loading OFF libraries
    if offlibraries is not None:
        if not isinstance(offlibraries, list) and not isinstance(offlibraries, tuple):
            offlibraries = [offlibraries, ]
        for off in offlibraries:
            f.write('loadoff {}\n\n'.format(off))

    # Loading user parameters
    f.write('# Loading parameter files\n')
    for p in param:
        try:
            shutil.copy(p, outdir)
            f.write('loadamberparams ' + path.basename(p) + '\n')
        except:
            f.write('loadamberparams ' + p + '\n')
            logger.info("Path {:s} not found, assuming a standard AmberTools file.".
                        format(p))
    f.write('\n')

    # Printing out topologies
    f.write('# Loading prepi topologies\n')
    for t in topo:
        shutil.copy(t, outdir)
        f.write('loadamberprep ' + path.basename(t) + '\n')
    f.write('\n')

    # Detect disulfide bonds
    if disulfide is None and not ionize:
        logger.info('Detecting disulfide bonds.')
        disulfide = detectDisulfideBonds(mol)
        if len(disulfide) != 0:
            for d in disulfide:
                # Convert to stupid amber residue numbering
                uqseqid = sequenceID((mol.resid, mol.insertion, mol.segid)) + mol.resid[0]
                uqres1 = int(np.unique(uqseqid[mol.atomselect('segid {} and resid {}'.format(d.segid1, d.resid1))]))
                uqres2 = int(np.unique(uqseqid[mol.atomselect('segid {} and resid {}'.format(d.segid2, d.resid2))]))
                # Rename the CYS to CYX if there is a disulfide bond
                mol.set('resname', 'CYX', sel='segid {} and resid {}'.format(d.segid1, d.resid1))
                mol.set('resname', 'CYX', sel='segid {} and resid {}'.format(d.segid2, d.resid2))
                # Remove (eventual) HG hydrogens on these CYS (from proteinPrepare)
                mol.remove('name HG and segid {} and resid {}'.format(d.segid1, d.resid1), _logger=False)
                mol.remove('name HG and segid {} and resid {}'.format(d.segid2, d.resid2), _logger=False)

    # Printing and loading the PDB file. AMBER can work with a single PDB file if the segments are separate by TER
    logger.info('Writing PDB file for input to tleap.')
    pdbname = path.join(outdir, 'input.pdb')

    # mol2 files have atomtype, here we only write parts not coming from mol2
    mol.write(pdbname, mol.atomtype == '')
    if not os.path.isfile(pdbname):
        raise NameError('Could not write a PDB file out of the given Molecule.')
    f.write('# Loading the system\n')
    f.write('mol = loadpdb input.pdb\n\n')

    if np.sum(mol.atomtype != '') != 0:
        logger.info('Writing mol2 files for input to tleap.')
        segs = np.unique(mol.segid[mol.atomtype != ''])
        combstr = 'mol = combine {mol'
        for s in segs:
            name = 'segment{}'.format(s)
            mol2name = path.join(outdir, '{}.mol2'.format(name))
            mol.write(mol2name, (mol.atomtype != '') & (mol.segid == s))
            if not os.path.isfile(mol2name):
                raise NameError('Could not write a mol2 file out of the given Molecule.')
            f.write('# Loading the rest of the system\n')
            f.write('{} = loadmol2 {}.mol2\n\n'.format(name, name))
            combstr += ' {}'.format(name)
        combstr += '}\n\n'
        f.write(combstr)

    # Printing out patches for the disulfide bridges
    if disulfide is None and not ionize:
        logger.info('Detecting disulfide bonds.')
        disulfide = detectDisulfideBonds(mol)

    # Write patches for disulfide bonds (only after ionizing)
    if not ionize and len(disulfide) != 0:
        f.write('# Adding disulfide bonds\n')
        for d in disulfide:
            # Convert to stupid amber residue numbering
            uqseqid = sequenceID((mol.resid, mol.insertion, mol.segid)) + mol.resid[0]
            uqres1 = int(np.unique(uqseqid[mol.atomselect('segid {} and resid {}'.format(d.segid1, d.resid1))]))
            uqres2 = int(np.unique(uqseqid[mol.atomselect('segid {} and resid {}'.format(d.segid2, d.resid2))]))
            f.write('bond mol.{}.SG mol.{}.SG\n'.format(uqres1, uqres2))
        f.write('\n')

    f.write('# Writing out the results\n')
    f.write('saveamberparm mol ' + prefix + '.prmtop ' + prefix + '.crd\n')
    f.write('quit')
    f.close()

    molbuilt = None
    if execute:
        # Source paths of extra dirs (our dirs, not amber default)
        htmdamberdir = path.abspath(path.join(home(), 'builder', 'amberfiles'))
        sourcepaths = [htmdamberdir]
        sourcepaths += [path.join(htmdamberdir, path.dirname(f))
                        for f in ff if path.isfile(path.join(htmdamberdir, f))]
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
        os.chdir(currdir)
        logger.info('Finished building.')

        if path.getsize(path.join(outdir, 'structure.crd')) != 0 and path.getsize(path.join(outdir, 'structure.prmtop')) != 0:
            molbuilt = Molecule(path.join(outdir, 'structure.prmtop'))
            molbuilt.read(path.join(outdir, 'structure.crd'))
        else:
            raise NameError('No structure pdb/prmtop file was generated. Check {} for errors in building.'.format(logpath))

        if ionize:
            shutil.move(path.join(outdir, 'structure.crd'), path.join(outdir, 'structure.noions.crd'))
            shutil.move(path.join(outdir, 'structure.prmtop'), path.join(outdir, 'structure.noions.prmtop'))
            totalcharge = np.sum(molbuilt.charge)
            nwater = np.sum(molbuilt.atomselect('water and noh'))
            anion, cation, anionatom, cationatom, nanion, ncation = ionizef(totalcharge, nwater, saltconc=saltconc, ff='amber', anion=saltanion, cation=saltcation)
            newmol = ionizePlace(mol, anion, cation, anionatom, cationatom, nanion, ncation)
            # Redo the whole build but now with ions included
            return build(newmol, ff=ff, topo=topo, param=param, prefix=prefix, outdir=outdir, caps={}, ionize=False,
                         execute=execute, saltconc=saltconc, disulfide=disulfide, tleap=tleap)
    tmpbonds = molbuilt.bonds
    molbuilt.bonds = []  # Removing the bonds to speed up writing
    molbuilt.write(path.join(outdir, 'structure.pdb'))
    molbuilt.bonds = tmpbonds # Restoring the bonds
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
    terminalatoms = {'ACE': 'H1 H2 H3 HT1 HT2 HT3', 'NME': 'OXT OT1 O'}  # XPLOR names for H[123] and OXT are HT[123]
                                                                         # and OT1, respectively.
    capresname = ['ACE', 'NME']
    capatomtype = ['C', 'N']

    # For each caps definition
    for seg in caps:
        # Get the segment
        segment = mol.atomselect('segid {}'.format(seg), indexes=True)
        # Test segment
        if len(segment) == 0:
            raise RuntimeError('There is no segment {} in the molecule.'.format(seg))
        if len(mol.atomselect('protein and segid {}'.format(seg), indexes=True)) == 0:
            raise RuntimeError(
                'Segment {} is not protein. Capping for non-protein segments is not supported.'.format(seg))
        # For each cap
        for i, cap in enumerate(caps[seg]):
            if cap is None or (isinstance(cap, str) and cap == 'none'):
                continue
            # Get info on segment and its terminals
            segment = mol.atomselect('segid {}'.format(seg), indexes=True)
            resids = np.unique(mol.get('resid', sel=segment))
            terminalids = [segment[0], segment[-1]]
            terminalresids = [np.min(resids), np.max(resids)]
            if i == 0:
                orig_terminalresids = [np.min(resids), np.max(resids)]
            # In case there is no cap defined
            if cap is None or cap == '':
                logger.warning(
                    'No cap provided for resid {} on segment {}. Did not apply it.'.format(terminalresids[i], seg))
                continue
            # If it is defined, test if supported
            elif cap not in capresname:
                raise RuntimeError(
                    'In segment {}, the {} cap is not supported. Try using {} instead.'.format(seg, cap, capresname))
            # Test if cap is already applied
            testcap = mol.atomselect('segid {} and resid "{}" and resname {}'.format(seg, terminalresids[i], cap),
                                     indexes=True)
            if len(testcap) != 0:
                logger.warning('Cap {} already exists on segment {}. Did not re-apply it.'.format(cap, seg))
                continue
            # Test if the atom to change exists
            termatomsids = mol.atomselect('segid {} and resid "{}" and name {}'.format(seg,
                                                                                       terminalresids[i],
                                                                                       terminalatoms[cap]),
                                          indexes=True)
            if len(termatomsids) == 0:
                # Create new atom
                termcaid = mol.atomselect('segid {} and resid "{}" and name CA'.format(seg, terminalresids[i]),
                                        indexes=True)
                termcenterid = mol.atomselect('segid {} and resid "{}" and name {}'.format(seg, terminalresids[i],
                                                                                           capatomtype[-i+1]),
                                        indexes=True)  # if i=0 => capatomtype[1]; i=1 => capatomtype[0]
                atom = Molecule()
                atom.empty(1)
                atom.record = 'ATOM'
                atom.name = capatomtype[i]
                atom.resid = terminalresids[i]-1+2*i
                atom.resname = cap
                atom.segid = seg
                atom.element = capatomtype[i]
                atom.chain = np.unique(mol.get('chain', sel='segid {}'.format(seg)))
                atom.coords = mol.coords[termcenterid] + 0.33 * np.subtract(mol.coords[termcenterid],
                                                                               mol.coords[termcaid])
                mol.insert(atom, terminalids[i])
                # newatom = mol.numAtoms - 1
                logger.info('In segment {}, resid {} had none of these atoms: {}. Capping was performed by creating '
                            'a new atom for cap construction by tleap.'.format(seg, terminalresids[i],
                                                                              terminalatoms[cap]))
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
            mol.remove('segid {} and resid "{}" and name {}'.format(seg, orig_terminalresids[i], terminalatoms[cap]),
                       _logger=False)

def _defaultProteinCaps(mol):
    # Defines ACE and NME (neutral terminals) as default for protein segments
    # Of course, this might not be ideal for proteins that require charged terminals

    segsProt = np.unique(mol.get('segid', sel='protein'))
    caps = dict()
    for s in segsProt:
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
    resdict = _readcsvdict(path.join(home(), 'builder', 'charmmlipid2amber.csv'))

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

    if len(begters) > 999:
        raise NameError('More than 999 lipids. Cannot define separate segments for all of them.')

    for i in range(len(begters)):
        map = np.zeros(len(mol.resid), dtype=bool)
        map[begters[i]:finters[i]+1] = True
        mol.set('resid', sequenceID(mol.get('resname', sel=map)), sel=map)
        mol.set('segid', 'L' + str(i+1), sel=map)

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


if __name__ == '__main__':
    from htmd.molecule.molecule import Molecule, mol_equal
    from htmd.builder.solvate import solvate
    from htmd.builder.preparation import proteinPrepare
    from htmd.builder.builder import autoSegment2
    from htmd.home import home
    from htmd.util import tempname
    import os
    from glob import glob
    import numpy as np
    import filecmp

    def cutfirstline(infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(outfile, 'w') as fout:
            fout.writelines(data[1:])


    def _compareResultFolders(compare, tmpdir, pid):
        ignore_ftypes = ('.log', '.txt')
        files = []
        deletefiles = []
        for f in glob(os.path.join(compare, '*')):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith('prmtop'):
                cutfirstline(f, os.path.join(compare, fname + '.mod'))
                cutfirstline(os.path.join(tmpdir, fname), os.path.join(tmpdir, fname + '.mod'))
                files.append(os.path.basename(f) + '.mod')
                deletefiles.append(os.path.join(compare, fname + '.mod'))
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
    #pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)
    for pid in pdbids:
        np.random.seed(1)
        mol = Molecule(pid)
        mol.filter('protein')
        mol = proteinPrepare(mol)
        mol.filter('protein')  # Fix for bad proteinPrepare hydrogen placing
        smol = solvate(mol)
        ffs = ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']
        tmpdir = tempname()
        bmol = build(smol, ff=ffs, outdir=tmpdir)

        refdir = home(dataDir=os.path.join('test-amber-build', pid))
        _compareResultFolders(refdir, tmpdir, pid)
        shutil.rmtree(tmpdir)

    # Test without proteinPrepare
    pdbids = ['3PTB']
    #pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
    for pid in pdbids:
        np.random.seed(1)
        mol = Molecule(pid)
        mol.filter('protein')
        smol = solvate(mol)
        ffs = ['leaprc.lipid14', 'leaprc.ff14SB', 'leaprc.gaff']
        tmpdir = tempname()
        bmol = build(smol, ff=ffs, outdir=tmpdir)

        refdir = home(dataDir=os.path.join('test-amber-build-nopp', pid))
        _compareResultFolders(refdir, tmpdir, pid)
        shutil.rmtree(tmpdir)

    # # Test protein-ligand building
    # folder = home(dataDir='building-protein-ligand')
    # prot = Molecule(os.path.join(folder, 'trypsin.pdb'))
    # prot.filter('protein')
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
    # bmol1 = build(smol1, param=[os.path.join(folder, 'benzamidine.frcmod')], outdir=tmpdir1)
    # np.random.seed(1)
    # bmol2 = build(smol2, topo=[os.path.join(folder, 'benzamidine.prepi')], param=[os.path.join(folder, 'benzamidineprepi.frcmod')], outdir=tmpdir2)
    # _compareResultFolders(tmpdir1, tmpdir2, 'ben-tryp')
    # shutil.rmtree(tmpdir1)
    # shutil.rmtree(tmpdir2)

