# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import numpy as np
import os.path as path
import os
import re
import shutil
import textwrap

from subprocess import call
from htmd.home import home
from htmd.molecule.molecule import Molecule
from htmd.molecule.util import _missingChain, _missingSegID
from htmd.builder.builder import detectDisulfideBonds
from htmd.builder.builder import _checkMixedSegment
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.vmdviewer import getVMDpath
from glob import glob
from natsort import natsorted

import logging
logger = logging.getLogger(__name__)

__test__ = {'build-opm-1u5u': """
    >>> from htmd.proteinpreparation.proteinpreparation import prepareProtein
    >>> from htmd.util import diffMolecules

    >>> pdb = os.path.join(home(dataDir="test-charmm-build"), '1u5u_opm.pdb')
    >>> mol = Molecule(pdb)
    >>> mol.filter('protein')
    >>> mol.set('segid', 'P')

    >>> pmol = prepareProtein(mol)
    >>> bmol = build(pmol, outdir='/tmp/build/', ionize=False)

    >>> refpdb = os.path.join(home(dataDir="test-charmm-build"), '1u5u_built_protonated.pdb')
    >>> ref = Molecule(refpdb)

    >>> difflist = diffMolecules(bmol, ref, sel="name CA")
    >>> print(difflist)
    []
""" }

class BuildError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def listFiles():
    """ Lists all available Charmm topologies and parameter files
    """
    charmmdir = path.join(home(), 'builder', 'charmmfiles', '')  # maybe just lookup current module?
    topos = natsorted(glob(path.join(charmmdir, 'top', '*.rtf')))
    params = natsorted(glob(path.join(charmmdir, 'par', '*.prm')))
    streams = natsorted(glob(path.join(charmmdir, 'str', '*', '*.str')))
    print('---- Topologies files list: ' + path.join(charmmdir, 'top', '') + ' ----')
    for t in topos:
        t = t.replace(charmmdir, '')
        print(t)
    print('---- Parameters files list: ' + path.join(charmmdir, 'par', '') + ' ----')
    for p in params:
        p = p.replace(charmmdir, '')
        print(p)
    print('---- Stream files list: ' + path.join(charmmdir, 'str', '') + ' ----')
    for s in streams:
        s = s.replace(charmmdir, '')
        print(s)


def build(mol, topo=None, param=None, stream=None, prefix='structure', outdir='./', caps=None, ionize=True, saltconc=0,
          saltanion=None, saltcation=None, disulfide=None, patches=None, psfgen=None, execute=True):
    """ Builds a system for CHARMM

    Uses VMD and psfgen to build a system for CHARMM. Additionally it allows for ionization and adding of disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The Molecule object containing the system
    topo : list of str
        A list of topology `rtf` files.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available topology files.
        Default: ['top/top_all36_prot.rtf', 'top/top_all36_lipid.rtf', 'top/top_water_ions.rtf']
    param : list of str
        A list of parameter `prm` files.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available parameter files.
        Default: ['par/par_all36_prot_mod.prm', 'par/par_all36_lipid.prm', 'par/par_water_ions.prm']
    stream : list of str
        A list of stream `str` files containing topologies and parameters.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available stream files.
        Default: ['str/prot/toppar_all36_prot_arg0.str']
    prefix : str
        The prefix for the generated pdb and psf files
    outdir : str
        The path to the output directory
    caps : dict
        A dictionary with keys segids and values lists of strings describing the caps of that segment.
        e.g. caps['P'] = ['first ACE', 'last CT3']. Default: will apply ACE and CT3 caps to proteins and none caps
        to the rest
    ionize : bool
        Enable or disable ionization
    saltconc : float
        Salt concentration (in Molar) to add to the system after neutralization.
    saltanion : {'CLA'}
        The anion type. Please use only CHARMM ion atom names.
    saltcation : {'SOD', 'MG', 'POT', 'CES', 'CAL', 'ZN2'}
        The cation type. Please use only CHARMM ion atom names.
    disulfide : list of :class:`DisulfideBridge <htmd.builder.builder.DisulfideBridge>` objects
        If None it will guess disulfide bonds. Otherwise provide a list of `DisulfideBridge` objects.
    patches : list of str
        Any further patches the user wants to apply
    psfgen : str
        Path to psfgen executable used to build for CHARMM
    execute : bool
        Disable building. Will only write out the input script needed by psfgen. Does not include ionization.

    Returns
    -------
    molbuilt : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The built system in a Molecule object

    Example
    -------
    >>> mol = Molecule("3PTB")
    >>> charmm.listFiles()             # doctest: +ELLIPSIS
    ---- Topologies files list...
    top/top_all36_prot.rtf
    top/top_water_ions.rtf
    ...
    >>> topos  = ['top/top_all36_prot.rtf', './benzamidine.rtf']
    >>> params = ['par/par_all36_prot_mod.prm', './benzamidine.prm']
    >>> molbuilt = charmm.build(mol, topo=topos, param=params, outdir='/tmp/build', saltconc=0.15)  # doctest: +SKIP
    """

    mol = mol.copy()
    _missingSegID(mol)
    _checkMixedSegment(mol)
    if psfgen is None:
        try:
            psfgen = shutil.which('psfgen', mode=os.X_OK)
        except:
            raise FileNotFoundError('Could not find psfgen executable, or no execute permissions are given. '
                                    'Run `conda install psfgen`.')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    _cleanOutDir(outdir)
    if topo is None:
        topo = ['top/top_all36_prot.rtf', 'top/top_all36_lipid.rtf', 'top/top_water_ions.rtf']
    if param is None:
        param = ['par/par_all36_prot_mod.prm', 'par/par_all36_lipid.prm', 'par/par_water_ions.prm']
    if stream is None:
        stream = ['str/prot/toppar_all36_prot_arg0.str']
    if caps is None:
        caps = _defaultCaps(mol)

    # Splitting the stream files and adding them to the list of parameter and topology files
    charmmdir = path.join(home(), 'builder', 'charmmfiles')
    for s in stream:
        if s[0] != '.' and path.isfile(path.join(charmmdir, s)):
            s = path.join(charmmdir, s)
        outrtf, outprm = _prepareStream(s)
        topo.append(outrtf)
        param.append(outprm)

    #_missingChain(mol)
    #_checkProteinGaps(mol)
    if patches is None:
        patches = []
    if isinstance(patches, str):
        patches = [patches]
    # Find protonated residues and add patches for them
    patches += _protonationPatches(mol)

    f = open(path.join(outdir, 'build.vmd'), 'w')
    f.write('# psfgen file generated by charmm.build\n')
    f.write('package require psfgen;\n')
    f.write('psfcontext reset;\n\n')

    # Copying and printing out the topologies
    for i in range(len(topo)):
        if topo[i][0] != '.' and path.isfile(path.join(charmmdir, topo[i])):
            topo[i] = path.join(charmmdir, topo[i])
        localname = '{}.'.format(i) + path.basename(topo[i])
        shutil.copy(topo[i], path.join(outdir, localname))
        f.write('topology ' + localname + '\n')
    f.write('\n')

    _printAliases(f)

    # Printing out segments
    logger.info('Writing out segments.')
    segments = _getSegments(mol)
    for seg in segments:
        pdbname = 'segment' + seg + '.pdb'
        mol.write(path.join(outdir, pdbname), sel='segid '+seg)

        segatoms = mol.atomselect('segid {}'.format(seg))
        segwater = mol.atomselect('segid {} and water'.format(seg))

        f.write('segment ' + seg + ' {\n')
        if np.all(segatoms == segwater):  # If segment only contains waters, set: auto none
            f.write('\tauto none\n')
        f.write('\tpdb ' + pdbname + '\n')
        if caps is not None and seg in caps:
            for c in caps[seg]:
                f.write('\t' + c + '\n')
        f.write('}\n')
        f.write('coordpdb ' + pdbname + ' ' + seg + '\n\n')

    # Printing out patches for the disulfide bridges
    if disulfide is None:
        disulfide = detectDisulfideBonds(mol)

    if len(disulfide) != 0:
        for d in disulfide:
            f.write('patch DISU {}:{} {}:{}\n'.format(d.segid1, d.resid1, d.segid2, d.resid2))
        f.write('\n')

    # Printing out extra patches
    if len(patches) != 0:
        for p in patches:
            f.write(p + '\n')
        f.write('\n')

    f.write('guesscoord\n')
    f.write('writepsf ' + prefix + '.psf\n')
    f.write('writepdb ' + prefix + '.pdb\n')
    #f.write('quit\n')
    f.close()

    if param is not None:
        _charmmCombine(param, path.join(outdir, 'parameters'))

    molbuilt = None
    if execute:
        logpath = os.path.abspath('{}/log.txt'.format(outdir))
        logger.info('Starting the build.')
        currdir = os.getcwd()
        os.chdir(outdir)
        f = open(logpath, 'w')
        #call([vmd, '-dispdev', 'text', '-e', './build.vmd'], stdout=f)
        call([psfgen, './build.vmd'], stdout=f)
        f.close()
        _logParser(logpath)
        os.chdir(currdir)
        logger.info('Finished building.')

        if path.isfile(path.join(outdir, 'structure.pdb')) and path.isfile(path.join(outdir, 'structure.psf')):
            molbuilt = Molecule(path.join(outdir, 'structure.pdb'))
            molbuilt.read(path.join(outdir, 'structure.psf'))
        else:
            raise BuildError('No structure pdb/psf file was generated. Check {} for errors in building.'.format(logpath))

        if ionize:
            shutil.move(path.join(outdir, 'structure.pdb'), path.join(outdir, 'structure.noions.pdb'))
            shutil.move(path.join(outdir, 'structure.psf'), path.join(outdir, 'structure.noions.psf'))
            totalcharge = np.sum(molbuilt.charge)
            nwater = np.sum(molbuilt.atomselect('water and noh'))
            anion, cation, anionatom, cationatom, nanion, ncation = ionizef(totalcharge, nwater, saltconc=saltconc, ff='charmm', anion=saltanion, cation=saltcation)
            newmol = ionizePlace(mol, anion, cation, anionatom, cationatom, nanion, ncation)
            # Redo the whole build but now with ions included
            return build(newmol, topo=topo, param=param, prefix=prefix, outdir=outdir, ionize=False, caps=caps,
                         execute=execute, saltconc=saltconc, disulfide=disulfide, patches=patches, psfgen=psfgen)
    _checkFailedAtoms(molbuilt)
    return molbuilt


def _cleanOutDir(outdir):
    from glob import glob
    files = glob(os.path.join(outdir, 'structure.*'))
    files += glob(os.path.join(outdir, 'log.*'))
    files += glob(os.path.join(outdir, '*.log'))
    files += glob(os.path.join(outdir, 'segment*'))
    for f in files:
        os.remove(f)


def _getSegments(mol):
    # Calculate unique segments but keep sorting
    indexes = np.unique(mol.segid, return_index=True)[1]
    uqseg = [mol.segid[index] for index in sorted(indexes)]
    return uqseg


def _printAliases(f):
    lines = '''
        # Aliases
        pdbalias residue G GUA
        pdbalias residue C CYT
        pdbalias residue A ADE
        pdbalias residue T THY
        pdbalias residue U URA

        foreach bp { GUA CYT ADE THY URA } {
            pdbalias atom $bp "O5\*" O5'
            pdbalias atom $bp "C5\*" C5'
            pdbalias atom $bp "O4\*" O4'
            pdbalias atom $bp "C4\*" C4'
            pdbalias atom $bp "C3\*" C3'
            pdbalias atom $bp "O3\*" O3'
            pdbalias atom $bp "C2\*" C2'
            pdbalias atom $bp "O2\*" O2'
            pdbalias atom $bp "C1\*" C1'
        }

        pdbalias atom ILE CD1 CD
        pdbalias atom SER HG HG1
        pdbalias residue HIS HSD

        # Heme aliases
        pdbalias residue HEM HEME
        pdbalias atom HEME "N A" NA
        pdbalias atom HEME "N B" NB
        pdbalias atom HEME "N C" NC
        pdbalias atom HEME "N D" ND

        # Water aliases
        pdbalias residue HOH TIP3
        pdbalias atom TIP3 O OH2

        # Ion aliases
        pdbalias residue K POT
        pdbalias atom K K POT
        pdbalias residue ICL CLA
        pdbalias atom ICL CL CLA
        pdbalias residue INA SOD
        pdbalias atom INA NA SOD
        pdbalias residue CA CAL
        pdbalias atom CA CA CAL
        pdbalias residue ZN ZN2

        # Other aliases
        pdbalias atom LYS 1HZ HZ1
        pdbalias atom LYS 2HZ HZ2
        pdbalias atom LYS 3HZ HZ3

        pdbalias atom ARG 1HH1 HH11
        pdbalias atom ARG 2HH1 HH12
        pdbalias atom ARG 1HH2 HH21
        pdbalias atom ARG 2HH2 HH22

        pdbalias atom ASN 1HD2 HD21
        pdbalias atom ASN 2HD2 HD22

        # Aliases for PDB ions / elements
        pdbalias residue LI LIT
        pdbalias residue NA SOD
        pdbalias residue K POT
        pdbalias residue CA CAL
        pdbalias residue ZN ZN2
        pdbalias residue CL CLA
        pdbalias residue RB RUB
        pdbalias residue CD CD2
        pdbalias residue CS CES
        pdbalias residue BA BAR

        # Aliases for Maestro residues
        pdbalias residue AR0 ARG
        pdbalias residue GLH GLU
        pdbalias residue ASH ASP
        pdbalias residue LYN LYS
        pdbalias residue HIE HSE
        pdbalias residue HID HSD
        pdbalias residue HIP HSP
        pdbalias residue CYX CYS
        pdbalias residue WAT TIP3

        # Generated by gen_psfaliases.py (Toni) on 2016-03-11 10:21
        pdbalias atom ALA H HN
        pdbalias atom ARG H HN
        pdbalias atom ARG HB3 HB1
        pdbalias atom ARG HG3 HG1
        pdbalias atom ARG HD3 HD1
        pdbalias atom ASP H HN
        pdbalias atom ASP HB3 HB1
        pdbalias atom ASN H HN
        pdbalias atom ASN HB3 HB1
        pdbalias atom CYS H HN
        pdbalias atom CYS HB3 HB1
        pdbalias atom GLU H HN
        pdbalias atom GLU HB3 HB1
        pdbalias atom GLU HG3 HG1
        pdbalias atom GLN H HN
        pdbalias atom GLN HB3 HB1
        pdbalias atom GLN HG3 HG1
        pdbalias atom GLY H HN
        pdbalias atom GLY HA3 HA1
        pdbalias atom HIS H HN
        pdbalias atom HIS HB3 HB1
        pdbalias atom ILE H HN
        pdbalias atom ILE HG13 HG11
        pdbalias atom ILE HD11 HD1
        pdbalias atom ILE HD12 HD2
        pdbalias atom ILE HD13 HD3
        pdbalias atom LEU H HN
        pdbalias atom LEU HB3 HB1
        pdbalias atom LYS H HN
        pdbalias atom LYS HB3 HB1
        pdbalias atom LYS HG3 HG1
        pdbalias atom LYS HD3 HD1
        pdbalias atom LYS HE3 HE1
        pdbalias atom MET H HN
        pdbalias atom MET HB3 HB1
        pdbalias atom MET HG3 HG1
        pdbalias atom PHE H HN
        pdbalias atom PHE HB3 HB1
        pdbalias atom PRO H2 HT2
        pdbalias atom PRO H3 HT1
        pdbalias atom PRO HB3 HB1
        pdbalias atom PRO HG3 HG1
        pdbalias atom PRO HD3 HD1
        pdbalias atom SER H HN
        pdbalias atom SER HB3 HB1
        pdbalias atom THR H HN
        pdbalias atom TRP H HN
        pdbalias atom TRP HB3 HB1
        pdbalias atom TYR H HN
        pdbalias atom TYR HB3 HB1
        pdbalias atom VAL H HN
    '''
    f.write(textwrap.dedent(lines))
    f.write('\n\n')


def _defaultCaps(mol):
    # neutral for protein, nothing for any other segment
    # of course this might not be ideal for protein which require charged terminals

    segsProt = np.unique(mol.get('segid', sel='protein'))
    segsNonProt = np.unique(mol.get('segid', sel='not protein'))
    caps = dict()
    for s in segsProt:
        nter, cter = _removeCappedResidues(mol, s)
        caps[s] = ['first {}'.format(nter), 'last {}'.format(cter)]
    for s in segsNonProt:
        caps[s] = ['first none', 'last none']
    return caps


def _removeCappedResidues(mol, seg):
    # Default caps for charmm
    nter = 'ACE'
    cter = 'CT3'

    # Mapping from various residue caps to charmm patches
    '''
    CHARMM patches:
    # N-terminus patches
    NTER         1.00 ! standard N-terminus
    GLYP         1.00 ! Glycine N-terminus
    PROP         1.00 ! Proline N-Terminal
    ACE          0.00 ! acetylated N-terminus
    ACED         0.00 ! acetylated N-terminus (to create dipeptide)
    ACP          0.00 ! acetylated N-terminus for proline
    ACPD         0.00 ! acetylated N-terminus for proline (to create dipeptide)
    NNEU         0.00 ! neutral N-terminus; charges from LSN
    # C-Terminus patches
    CTER        -1.00 ! standard C-terminus
    CNEU         0.00 ! protonated (neutral) C-terminu, charges from ASPP
    CTP          0.00 ! protonated C-terminus
    CT1          0.00 ! methylated C-terminus from methyl acetate
    CT2          0.00 ! amidated C-terminus
    CT3          0.00 ! N-Methylamide C-terminus
    '''

    ntercaps = {'ACE': 'ACE'}
    ctercaps = {'NMA': 'CT3'}

    # If caps already exist, remove them and convert them to charmm caps
    for n in ntercaps:
        sel = 'segid {} and resname {}'.format(seg, n)
        asel = mol.atomselect(sel)
        if np.sum(asel) != 0:
            mol.remove(asel, _logger=False)
            nter = ntercaps[n]
            break  # I expect only one capped n-terminal residue in one segment!
    for c in ctercaps:
        sel = 'segid {} and resname {}'.format(seg, c)
        asel = mol.atomselect(sel)
        if np.sum(asel) != 0:
            mol.remove(asel, _logger=False)
            cter = ctercaps[c]
            break  # I expect only one capped c-terminal residue in one segment!
    return nter, cter


# Mapping Maestro protonated residue names to CHARMM patches
def _protonationPatches(mol):
    protonations = {'GLH': 'GLUP', 'ASH': 'ASPP', 'LYN': 'LSN', 'AR0': 'RN1'}
    aliases = {}  # Some protonations don't exist in CHARMM
    # TODO: Remember to alias all of them before applying patches
    patches = []

    for pro in protonations:
        pseg = mol.get('segid', sel='resname {} and name CA'.format(pro))
        pres = mol.get('resid', sel='resname {} and name CA'.format(pro))
        if len(pseg) == 0:
            continue
        for r in range(len(pseg)):
            #from IPython.core.debugger import Tracer
            #Tracer()()
            patch = 'patch {} {}:{}'.format(protonations[pro], pseg[r], pres[r])
            patches.append(patch)

    '''for pro in aliases:
        sel = mol.atomselect('resname {}'.format(pro))
        if np.sum(sel) != 0:
            logger.warning('Found resname {}. This protonation state does not exist in CHARMM '
                           'and will be reverted to {}.'.format(pro, aliases[pro]))
            mol.set('resname', aliases[pro], sel=sel)'''
    for pro in aliases:
        sel = mol.atomselect('resname {}'.format(pro))
        if np.sum(sel) != 0:
            raise RuntimeError('Found resname {}. This protonation state does not exist in CHARMM. Cannot build.')

    return patches


def _charmmCombine(prmlist, outfile):
    """ Combines CHARMM parameter files.

    Parameters
    ----------
    prmlist
    outfile

    Returns
    -------

    """
    # Process parameter files
    prm_list = ["!COMMENTS\n", "!ATOMS\n", "BONDS\n", "ANGLES\n", "DIHEDRALS\n", "IMPROPER\n", "CMAP\n", "NONBONDED\n", "NBFIX\n", "HBOND\n"]

    charmmdir = path.join(home(), 'builder', 'charmmfiles')
    for myfile in prmlist:
        if myfile[0] != '.' and path.isfile(path.join(charmmdir, myfile)):
            myfile = path.join(charmmdir, myfile)
        if not path.isfile(myfile):
            raise FileNotFoundError(myfile + ' file does not exist. Cannot create combined parameter file.')
        fn = os.path.basename(myfile)
        fh = open(myfile, "r")
        context = 0
        for line in fh:
            if re.search(r'^ATOMS', line):
                context = 1
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^BOND', line):
                context = 2
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^ANGL', line):
                context = 3
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^DIHE', line) or re.search(r'^THET', line):
                context = 4
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^IMPR', line) or re.search(r'^IMPH', line):
                context = 5
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^CMAP', line) or re.search(r'^NBON', line):
                context = 6
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^NONB', line):
                context = 7
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^NBFI', line):
                context = 8
                prm_list[context] += _sec_name(fn)
            elif re.search(r'^HBON', line):
                context = 9
                prm_list[context] += _sec_name(fn)
            else:
                if context == 0: # COMMENTS
                    if re.search(r'^\s*\!+',line,re.I) or re.search(r'^\s*\*+',line,re.I):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 1: # ATOMS
                    if re.search(r'^\s*\!+',line,re.I):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 2: # BONDS
                    if re.search(r'^\s*\!+',line) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 3: # ANGLES
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 4: # DIHEDRALS
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 5: # IMPROPER
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 6: # CMAP
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+',line) or \
                       re.search(r'^\s*[0-9\.-]+\s+[0-9\.-]+\s+[0-9\.-]+\s+[0-9\.-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 7: # NONBONDED
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[_A-Z0-9a-z\.]+\s+[0-9\.\-]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 8: # NBFIX
                    if re.search(r'^\s*\!+',line,re.I) or \
                       re.search(r'^\s*$',line) or \
                       re.search(r'^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+',line):
                        prm_list[context] += line
                    else:
                        prm_list[context] += "!"+line
                elif context == 9: # HBOND
                    if not re.search(r'^END',line,re.I):
                        prm_list[context] += line
                else:
                    continue

    prm = ''.join(map(str, prm_list))+"END"
    prmfh = open(outfile, "w")
    prmfh.write(prm)
    prmfh.close()


# Create string to indicate source of section
def _sec_name(filename):
    return "!Following lines added from %s\n" % (filename)


def _charmmSplit(filename, outdir):
    """ Splits a stream file into an rtf and prm file.

    Parameters
    ----------
    filename : str
        Stream file name
    """
    regex = re.compile('^toppar_(\w+)\.str$')
    base = os.path.basename(os.path.normpath(filename))
    base = regex.findall(base)[0]
    outrtf = os.path.join(outdir, 'top_{}.rtf'.format(base))
    outprm = os.path.join(outdir, 'par_{}.prm'.format(base))

    startrtf = re.compile('^read rtf card', flags=re.IGNORECASE)
    startprm = re.compile('^read param? card', flags=re.IGNORECASE)
    endsection = re.compile('^end', flags=re.IGNORECASE)

    rtfsection = 0
    prmsection = 0
    section = 'junk'

    rtfstr = ''
    prmstr = ''

    f = open(filename, 'r')
    for line in f:
        if startrtf.match(line):
            rtfsection += 1
            if rtfsection > 1:
                rtfstr += '! WARNING -- ANOTHER rtf SECTION FOUND\n'
            section = 'rtf'
        elif startprm.match(line):
            prmsection += 1
            if prmsection > 1:
                prmstr += '! WARNING -- ANOTHER para SECTION FOUND\n'
            section = 'prm'
        elif endsection.match(line):
            section = 'junk'
        elif section == 'rtf':
            rtfstr += line
        elif section == 'prm':
            prmstr += line
    f.close()

    if rtfsection > 1:
        logger.warning('Multiple ({}) rtf topology sections found in {} stream file.'.format(rtfsection, filename))
    if prmsection > 1:
        logger.warning('Multiple ({}) prm parameter sections found in {} stream file.'.format(prmsection, filename))

    f = open(outrtf, 'w')
    f.write(rtfstr + 'END\n')
    f.close()
    f = open(outprm, 'w')
    f.write(prmstr + 'END\n')
    f.close()
    return outrtf, outprm


def _prepareStream(filename):
    from htmd.util import tempname
    tmpout = tempname()
    os.makedirs(tmpout)
    return _charmmSplit(filename, tmpout)


def _logParser(fname):
    import re
    failedcoor = re.compile('Warning: failed to set coordinate for atom')
    failedangle = re.compile('Warning: failed to guess coordinate due to bad angle')
    poorlycoor = re.compile('Warning: poorly guessed coordinate(s?)')
    otherwarn = re.compile('Warning')

    failedcoorcount = 0
    failedanglecount = 0
    poorlycoorcount = -1  # Discount the summary report message in the log
    otherwarncount = 0
    f = open(fname, 'r')
    for line in f:
        if failedcoor.search(line):
            failedcoorcount += 1
        elif failedangle.search(line):
            failedanglecount += 1
        elif poorlycoor.search(line):
            poorlycoorcount += 1
        elif otherwarn.search(line):
            otherwarncount += 1

    f.close()
    warnings = False
    if failedcoorcount > 0:
        warnings = True
        logger.warning('Failed to set coordinates for {} atoms.'.format(failedcoorcount))
    if failedanglecount > 0:
        warnings = True
        logger.warning('Failed to guess coordinates for {} atoms due to bad angles.'.format(failedanglecount))
    if poorlycoorcount > 0:
        warnings = True
        logger.warning('Poorly guessed coordinates for {} atoms.'.format(poorlycoorcount))
    if otherwarncount > 0:
        warnings = True
        logger.warning('{} undefined warnings were produced during building.'.format(otherwarncount))
    if warnings:
        logger.warning('Please check {} for further information.'.format(fname))


def _checkFailedAtoms(mol):
    if mol is None:
        return
    idx = np.where(np.sum(mol.coords == 0, axis=1) == 3)[0]
    if len(idx) != 0:
        logger.critical('Atoms with indexes {} are positioned at [0,0,0]. This can cause simulations to crash. '
                        'Check log file for more details.'.format(idx))


if __name__ == '__main__':
    from htmd import *
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import numpy as np
    import os

    import doctest
    doctest.testmod()


