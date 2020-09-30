# (c) 2015-2018 Acellera Ltd http://www.acellera.com
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
from moleculekit.molecule import Molecule
from moleculekit.util import _missingChain, _missingSegID
from htmd.builder.builder import (
    detectDisulfideBonds,
    convertDisulfide,
    detectCisPeptideBonds,
    _checkMixedSegment,
    _checkLongResnames,
    UnknownResidueError,
    BuildError,
)
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from moleculekit.vmdviewer import getVMDpath
from glob import glob
from unittest import TestCase


import logging

logger = logging.getLogger(__name__)


def htmdCharmmHome():
    """ Returns the location of the CHARMM files distributed with HTMD"""
    return os.path.abspath(
        os.path.join(home(shareDir=True), "builder", "charmmfiles", "")
    )


def listFiles():
    """ Lists all available Charmm topologies and parameter files

    Examples
    --------
    >>> from htmd.builder import charmm
    >>> charmm.listFiles()             # doctest: +ELLIPSIS
    ---- Topologies files list...

    """
    from natsort import natsorted

    charmmdir = htmdCharmmHome()
    topos = natsorted(glob(path.join(charmmdir, "top", "*.rtf")))
    params = natsorted(glob(path.join(charmmdir, "par", "*.prm")))
    streams = natsorted(glob(path.join(charmmdir, "str", "*", "*.str")))
    print("---- Topologies files list: " + path.join(charmmdir, "top", "") + " ----")
    for t in topos:
        t = t.replace(charmmdir, "")
        print(t)
    print("---- Parameters files list: " + path.join(charmmdir, "par", "") + " ----")
    for p in params:
        p = p.replace(charmmdir, "")
        print(p)
    print("---- Stream files list: " + path.join(charmmdir, "str", "") + " ----")
    for s in streams:
        s = s.replace(charmmdir, "")
        print(s)


def search(key, name):
    """ Searches for CHARMM files containing a given definition.

    Parameters
    ----------
    key : str
        A key
    name : str
        The corresponding name
        
    Examples
    --------
    >>> charmm.search(key='RESI', name = 'CHL1')  # doctest: +SKIP
    """
    os.system(
        'find {} -type f -exec grep -n "{} {}" {{}} +'.format(
            htmdCharmmHome(), key, name
        )
    )


def defaultTopo():
    """ Returns the default topologies used by charmm.build """
    return [
        "top/top_all36_prot.rtf",
        "top/top_all36_lipid.rtf",
        "top/top_water_ions.rtf",
        "top/top_all36_cgenff.rtf",
    ]


def defaultParam():
    """ Returns the default parameters used by charmm.build """
    return [
        "par/par_all36m_prot.prm",
        "par/par_all36_lipid.prm",
        "par/par_water_ions.prm",
        "par/par_all36_cgenff.prm",
    ]


def defaultStream():
    """ Returns the default stream files used by charmm.build """
    return ["str/prot/toppar_all36_prot_arg0.str", "str/misc/toppar_ions_won.str"]


def build(
    mol,
    topo=None,
    param=None,
    stream=None,
    prefix="structure",
    outdir="./build",
    caps=None,
    ionize=True,
    saltconc=0,
    saltanion=None,
    saltcation=None,
    disulfide=None,
    regenerate=["angles", "dihedrals"],
    patches=None,
    noregen=None,
    aliasresidues=None,
    psfgen=None,
    execute=True,
    _clean=True,
):
    """ Builds a system for CHARMM

    Uses VMD and psfgen to build a system for CHARMM. Additionally it allows for ionization and adding of disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object containing the system
    topo : list of str
        A list of topology `rtf` files.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available topology files.
        Default: ['top/top_all36_prot.rtf', 'top/top_all36_lipid.rtf', 'top/top_water_ions.rtf']
    param : list of str
        A list of parameter `prm` files.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available parameter files.
        Default: ['par/par_all36_prot.prm', 'par/par_all36_lipid.prm', 'par/par_water_ions.prm']
    stream : list of str
        A list of stream `str` files containing topologies and parameters.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available stream files.
        Default: ['str/prot/toppar_all36_prot_arg0.str']
    prefix : str
        The prefix for the generated pdb and psf files
    outdir : str
        The path to the output directory
        Default: './build'
    caps : dict
        A dictionary with keys segids and values lists of strings describing the caps of that segment.
        e.g. caps['P'] = ['first ACE', 'last CT3'] or caps['P'] = ['first none', 'last none']. 
        Default: will apply ACE and CT3 caps to proteins and none caps to the rest.
    ionize : bool
        Enable or disable ionization
    saltconc : float
        Salt concentration (in Molar) to add to the system after neutralization.
    saltanion : {'CLA'}
        The anion type. Please use only CHARMM ion atom names.
    saltcation : {'SOD', 'MG', 'POT', 'CES', 'CAL', 'ZN2'}
        The cation type. Please use only CHARMM ion atom names.
    disulfide : list of pairs of atomselection strings
        If None it will guess disulfide bonds. Otherwise provide a list pairs of atomselection strings for each pair of
        residues forming the disulfide bridge.
    regenerate : None or list of strings of: ['angles', 'dihedrals']
        Disable angle/dihedral regeneration with `regenerate=None`, or enable it with `regenerate=['angles', 'diheldrals']` 
        or just one of the two options with `regenerate=['angles']` or `regenerate=['diheldrals']`.
    patches : list of str
        Any further patches the user wants to apply
    noregen : list of str
        A list of patches that must not be regenerated (angles and dihedrals)
        Default: ['FHEM', 'PHEM', 'PLOH', 'PLO2', 'PLIG', 'PSUL']
    aliasresidues : dict of aliases
        A dictionary of key: value pairs of residue names we want to alias
    psfgen : str
        Path to psfgen executable used to build for CHARMM
    execute : bool
        Disable building. Will only write out the input script needed by psfgen. Does not include ionization.

    Returns
    -------
    molbuilt : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The built system in a Molecule object

    Example
    -------
    >>> from htmd.ui import *
    >>> mol = Molecule("3PTB")
    >>> mol.filter("not resname BEN")
    >>> molbuilt = charmm.build(mol, outdir='/tmp/build', ionize=False)  # doctest: +ELLIPSIS
    Bond between A: [serial 185 resid 42 resname CYS chain A segid 0]
                 B: [serial 298 resid 58 resname CYS chain A segid 0]...
    >>> # More complex example
    >>> topos  = ['top/top_all36_prot.rtf', './benzamidine.rtf', 'top/top_water_ions.rtf']
    >>> params = ['par/par_all36_prot.prm', './benzamidine.prm', 'par/par_water_ions.prm']
    >>> disu = [['segid P and resid 157', 'segid P and resid 13'], ['segid K and resid 1', 'segid K and resid 25']]
    >>> ar = {'SAPI24': 'SP24'}  # Alias large resnames to a short-hand version
    >>> molbuilt = charmm.build(mol, topo=topos, param=params, outdir='/tmp/build', saltconc=0.15, disulfide=disu, aliasresidues=ar)  # doctest: +SKIP
    """

    mol = mol.copy()
    _missingSegID(mol)
    _checkMixedSegment(mol)
    _checkLongResnames(mol, aliasresidues)
    if psfgen is None:
        psfgen = shutil.which("psfgen", mode=os.X_OK)
        if not psfgen:
            raise FileNotFoundError(
                "Could not find psfgen executable, or no execute permissions are given. "
                "Run `conda install psfgen -c acellera`."
            )
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if _clean:
        _cleanOutDir(outdir)
    if topo is None:
        topo = defaultTopo()
    if param is None:
        param = defaultParam()
    if stream is None:
        stream = defaultStream()
    if caps is None:
        caps = _defaultCaps(mol)
    # patches that must _not_ be regenerated
    if noregen is None:
        noregen = ["FHEM", "PHEM", "PLOH", "PLO2", "PLIG", "PSUL"]

    alltopo = topo.copy()
    allparam = param.copy()

    # Splitting the stream files and adding them to the list of parameter and topology files
    charmmdir = htmdCharmmHome()
    for s in stream:
        if s[0] != "." and path.isfile(path.join(charmmdir, s)):
            s = path.join(charmmdir, s)
        outrtf, outprm = _prepareStream(s)
        alltopo.append(outrtf)
        allparam.append(outprm)

    # _missingChain(mol)
    # _checkProteinGaps(mol)
    if patches is None:
        patches = []
    if isinstance(patches, str):
        patches = [patches]
    allpatches = []
    allpatches += patches
    # Find protonated residues and add patches for them
    allpatches += _protonationPatches(mol)

    f = open(path.join(outdir, "build.vmd"), "w")
    f.write("# psfgen file generated by charmm.build\n")
    f.write("package require psfgen;\n")
    f.write("psfcontext reset;\n\n")

    # Copying and printing out the topologies
    if not path.exists(path.join(outdir, "topologies")):
        os.makedirs(path.join(outdir, "topologies"))
    for i in range(len(alltopo)):
        if alltopo[i][0] != "." and path.isfile(path.join(charmmdir, alltopo[i])):
            alltopo[i] = path.join(charmmdir, alltopo[i])
        localname = "{}.".format(i) + path.basename(alltopo[i])
        shutil.copy(alltopo[i], path.join(outdir, "topologies", localname))
        f.write("topology " + path.join("topologies", localname) + "\n")
    f.write("\n")

    _printAliases(f)
    if aliasresidues is not None:  # User defined aliases
        for key, val in aliasresidues.items():
            mol.resname[mol.resname == key] = val
            f.write("        pdbalias residue {} {}\n".format(val, key))

    # Printing out segments
    if not path.exists(path.join(outdir, "segments")):
        os.makedirs(path.join(outdir, "segments"))
    logger.info("Writing out segments.")
    segments = _getSegments(mol)
    wateratoms = mol.atomselect("water")
    for seg in segments:
        pdbname = "segment" + seg + ".pdb"
        segatoms = mol.segid == seg
        mol.write(path.join(outdir, "segments", pdbname), sel=segatoms)

        segwater = wateratoms & segatoms
        f.write("segment " + seg + " {\n")
        if np.all(
            segatoms == segwater
        ):  # If segment only contains waters, set: auto none
            f.write("\tauto none\n")
        f.write("\tpdb " + path.join("segments", pdbname) + "\n")
        if caps is not None and seg in caps:
            for c in caps[seg]:
                f.write("\t" + c + "\n")
        f.write("}\n")
        f.write("coordpdb " + path.join("segments", pdbname) + " " + seg + "\n\n")

    # Printing out patches for the disulfide bridges
    # TODO: Remove this once we deprecate the class
    from htmd.builder.builder import DisulfideBridge
    from moleculekit.molecule import UniqueResidueID

    if (
        disulfide is not None
        and len(disulfide) != 0
        and isinstance(disulfide[0], DisulfideBridge)
    ):
        newdisu = []
        for d in disulfide:
            r1 = UniqueResidueID.fromMolecule(
                mol, "resid {} and segname {}".format(d.resid1, d.segid1)
            )
            r2 = UniqueResidueID.fromMolecule(
                mol, "resid {} and segname {}".format(d.resid2, d.segid2)
            )
            newdisu.append([r1, r2])
        disulfide = newdisu
    # TODO: Remove up to here ----------------------

    if (
        disulfide is not None
        and len(disulfide) != 0
        and isinstance(disulfide[0][0], str)
    ):
        disulfide = convertDisulfide(mol, disulfide)

    if disulfide is None:
        disulfide = detectDisulfideBonds(mol)

    if len(disulfide) != 0:
        for d in disulfide:
            str0 = "{}:{}{}".format(d[0].segid, d[0].resid, d[0].insertion)
            str1 = "{}:{}{}".format(d[1].segid, d[1].resid, d[1].insertion)
            f.write("patch DISU {} {}\n".format(str0, str1))
        f.write("\n")

    noregenpatches = [p for p in allpatches if p.split()[1] in noregen]
    regenpatches = [p for p in allpatches if p.split()[1] not in noregen]

    # Printing regenerable patches
    if len(regenpatches) != 0:
        for p in regenpatches:
            f.write(p + "\n")
        f.write("\n")

    # Regenerate angles and dihedrals
    if regenerate is not None:
        f.write("regenerate {}\n".format(" ".join(regenerate)))
        f.write("\n")

    # Printing non-regenerable patches
    if len(noregenpatches) != 0:
        for p in noregenpatches:
            f.write(p + "\n")
        f.write("\n")

    f.write("guesscoord\n")
    f.write("writepsf " + prefix + ".psf\n")
    f.write("writepdb " + prefix + ".pdb\n")
    # f.write('quit\n')
    f.close()

    if allparam is not None:
        combine(allparam, path.join(outdir, "parameters"))

    molbuilt = None
    if execute:
        logpath = os.path.abspath("{}/log.txt".format(outdir))
        logger.info("Starting the build.")
        currdir = os.getcwd()
        os.chdir(outdir)
        f = open(logpath, "w")
        # call([vmd, '-dispdev', 'text', '-e', './build.vmd'], stdout=f)
        my_env = os.environ.copy()
        my_env["LC_ALL"] = "C"
        call([psfgen, "./build.vmd"], stdout=f, stderr=f, env=my_env)
        f.close()
        errors = _logParser(logpath)
        os.chdir(currdir)
        if errors:
            raise BuildError(
                errors
                + [
                    "Check {} for further information on errors in building.".format(
                        logpath
                    )
                ]
            )
        logger.info("Finished building.")

        if path.isfile(path.join(outdir, "structure.pdb")) and path.isfile(
            path.join(outdir, "structure.psf")
        ):
            molbuilt = Molecule(path.join(outdir, "structure.pdb"))
            molbuilt.read(path.join(outdir, "structure.psf"))
        else:
            raise BuildError(
                "No structure pdb/psf file was generated. Check {} for errors in building.".format(
                    logpath
                )
            )

        if ionize:
            os.makedirs(path.join(outdir, "pre-ionize"))
            data = glob(path.join(outdir, "*"))
            for f in data:
                shutil.move(f, path.join(outdir, "pre-ionize"))
            totalcharge = np.sum(molbuilt.charge)
            nwater = np.sum(molbuilt.atomselect("water and noh"))
            anion, cation, anionatom, cationatom, nanion, ncation = ionizef(
                totalcharge,
                nwater,
                saltconc=saltconc,
                anion=saltanion,
                cation=saltcation,
            )
            newmol = ionizePlace(
                mol, anion, cation, anionatom, cationatom, nanion, ncation
            )
            # Redo the whole build but now with ions included
            return build(
                newmol,
                topo=alltopo,
                param=allparam,
                stream=[],
                prefix=prefix,
                outdir=outdir,
                ionize=False,
                caps=caps,
                execute=execute,
                saltconc=saltconc,
                disulfide=disulfide,
                regenerate=regenerate,
                patches=patches,
                noregen=noregen,
                aliasresidues=aliasresidues,
                psfgen=psfgen,
                _clean=False,
            )
    _checkFailedAtoms(molbuilt)
    _recoverProtonations(molbuilt)
    detectCisPeptideBonds(molbuilt)  # Warn in case of cis bonds
    return molbuilt


def _cleanOutDir(outdir):
    from glob import glob

    files = glob(os.path.join(outdir, "structure.*"))
    files += glob(os.path.join(outdir, "log.*"))
    files += glob(os.path.join(outdir, "*.log"))
    files += glob(os.path.join(outdir, "*.vmd"))
    files += glob(os.path.join(outdir, "parameters"))
    for f in files:
        os.remove(f)
    folders = glob(os.path.join(outdir, "segments"))
    folders += glob(os.path.join(outdir, "topologies"))
    folders += glob(os.path.join(outdir, "pre-ionize"))
    for f in folders:
        shutil.rmtree(f)


def _getSegments(mol):
    # Calculate unique segments but keep sorting
    indexes = np.unique(mol.segid, return_index=True)[1]
    uqseg = [mol.segid[index] for index in sorted(indexes)]
    return uqseg


def _printAliases(f):
    lines = """
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
        pdbalias residue FE2 Fe2p
        pdbalias atom Fe2p FE Fe2p 

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
        pdbalias residue NI Ni2p
        pdbalias residue HG Hg2p

        # Aliases for Maestro residues
        pdbalias residue AR0 ARG
        pdbalias residue GLH GLU
        pdbalias residue ASH ASP
        pdbalias residue LYN LYS
        pdbalias residue HIE HSE
        pdbalias residue HID HSD
        pdbalias residue HIP HSP
        pdbalias residue CYX CYS
        pdbalias residue CYM CYS
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

        # Aliases for carbohydrates
        pdbalias residue NAG BGLCNA
        pdbalias residue BMA BMAN
        pdbalias residue GLA AGAL
        pdbalias residue GAL BGAL
        pdbalias residue MAN AMAN
        pdbalias residue FUC AFUC
        pdbalias residue FUL BFUC

    """
    f.write(textwrap.dedent(lines))
    f.write("\n\n")


def _defaultCaps(mol):
    # neutral for protein, nothing for any other segment
    # of course this might not be ideal for protein which require charged terminals
    prot = mol.atomselect("protein")
    segsProt = np.unique(mol.segid[prot])
    segsNonProt = np.unique(mol.segid[~prot])
    caps = dict()
    for s in segsProt:
        if len(np.unique(mol.resid[mol.segid == s])) < 10:
            logger.warning(
                "Segment {} consists of a peptide with less than 10 residues. It will not be capped by "
                "default. If you want to cap it use the caps argument of charmm.build to manually define "
                "caps for all segments".format(s)
            )
            continue
        nter, cter = _removeCappedResidues(mol, s)
        caps[s] = ["first {}".format(nter), "last {}".format(cter)]
    for s in segsNonProt:
        caps[s] = ["first none", "last none"]
    return caps


def _removeCappedResidues(mol, seg):
    # Default caps for charmm
    nter = "ACE"
    cter = "CT3"

    # Mapping from various residue caps to charmm patches
    """
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
    """

    ntercaps = {"ACE": "ACE"}
    ctercaps = {"NMA": "CT3"}

    # If caps already exist, remove them and convert them to charmm caps
    for n in ntercaps:
        asel = (mol.segid == seg) & (mol.resname == n)
        if np.sum(asel) != 0:
            mol.remove(asel, _logger=False)
            nter = ntercaps[n]
            break  # I expect only one capped n-terminal residue in one segment!
    for c in ctercaps:
        asel = (mol.segid == seg) & (mol.resname == c)
        if np.sum(asel) != 0:
            mol.remove(asel, _logger=False)
            cter = ctercaps[c]
            break  # I expect only one capped c-terminal residue in one segment!
    return nter, cter


# Mapping Maestro protonated residue names to CHARMM patches
def _protonationPatches(mol):
    protonations = {
        "GLH": "GLUP",
        "ASH": "ASPP",
        "LYN": "LSN",
        "AR0": "RN1",
        "CYM": "CYSD",
    }
    aliases = {}  # Some protonations don't exist in CHARMM
    # TODO: Remember to alias all of them before applying patches
    patches = []

    for pro in sorted(protonations.keys()):
        sel = (mol.resname == pro) & (mol.name == "CA")
        pseg = mol.segid[sel]
        pres = mol.resid[sel]
        if len(pseg) == 0:
            continue
        for r in range(len(pseg)):
            # from IPython.core.debugger import Tracer
            # Tracer()()
            patch = "patch {} {}:{}".format(protonations[pro], pseg[r], pres[r])
            patches.append(patch)

    """for pro in aliases:
        sel = mol.atomselect('resname {}'.format(pro))
        if np.sum(sel) != 0:
            logger.warning('Found resname {}. This protonation state does not exist in CHARMM '
                           'and will be reverted to {}.'.format(pro, aliases[pro]))
            mol.set('resname', aliases[pro], sel=sel)"""
    for pro in aliases:
        sel = mol.resname == pro
        if np.any(sel):
            raise RuntimeError(
                "Found resname {}. This protonation state does not exist in CHARMM. Cannot build."
            )

    return patches


# Recover protonation states of residues after building (CHARMM renames protonated residues to standard names)
def _recoverProtonations(mol):
    mol.set("resname", "ASH", sel="same residue as (resname ASP and name HD2)")
    mol.set("resname", "GLH", sel="same residue as (resname GLU and name HE2)")
    mol.set(
        "resname",
        "LYN",
        sel="resname LYS and not (same residue as (resname LYS and name HZ3))",
    )  # The LYN patch removes the HZ3 proton
    mol.set(
        "resname",
        "AR0",
        sel="resname ARG and not (same residue as (resname ARG and name HE))",
    )  # The AR0 patch removes the HE proton
    # Histidine protonations keep their names in CHARMM. No need to rename them


def combine(prmlist, outfile):
    """ Combines CHARMM parameter files
    Take a list of parameters files and combine them into a single file (useful for acemd)

    Parameters
    ----------
    prmlist: list
        List of parameter files to combine
    outfile: str
        Output filename of combined parameter files

    """
    # Process parameter files
    prm_list = [
        "!COMMENTS\n",
        "!ATOMS\n",
        "BONDS\n",
        "ANGLES\n",
        "DIHEDRALS\n",
        "IMPROPER\n",
        "CMAP\n",
        "NONBONDED\n",
        "NBFIX\n",
        "HBOND\n",
    ]

    charmmdir = htmdCharmmHome()
    for myfile in prmlist:
        if myfile[0] != "." and path.isfile(path.join(charmmdir, myfile)):
            myfile = path.join(charmmdir, myfile)
        if not path.isfile(myfile):
            raise FileNotFoundError(
                myfile + " file does not exist. Cannot create combined parameter file."
            )
        fn = os.path.basename(myfile)
        with open(myfile, "r") as fh:
            context = 0
            for line in fh:
                if re.search(r"^ATOMS", line):
                    context = 1
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^BOND", line):
                    context = 2
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^ANGL", line):
                    context = 3
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^DIHE", line) or re.search(r"^THET", line):
                    context = 4
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^IMPR", line) or re.search(r"^IMPH", line):
                    context = 5
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^CMAP", line) or re.search(r"^NBON", line):
                    context = 6
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^NONB", line):
                    context = 7
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^NBFI", line):
                    context = 8
                    prm_list[context] += _sec_name(fn)
                elif re.search(r"^HBON", line):
                    context = 9
                    prm_list[context] += _sec_name(fn)
                else:
                    if context == 0:  # COMMENTS
                        if re.search(r"^\s*\!+", line, re.I) or re.search(
                            r"^\s*\*+", line, re.I
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 1:  # ATOMS
                        if re.search(r"^\s*\!+", line, re.I):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 2:  # BONDS
                        if (
                            re.search(r"^\s*\!+", line)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 3:  # ANGLES
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 4:  # DIHEDRALS
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 5:  # IMPROPER
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 6:  # CMAP
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+",
                                line,
                            )
                            or re.search(
                                r"^\s*[0-9\.-]+\s+[0-9\.-]+\s+[0-9\.-]+\s+[0-9\.-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 7:  # NONBONDED
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[_A-Z0-9a-z\.]+\s+[0-9\.\-]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 8:  # NBFIX
                        if (
                            re.search(r"^\s*\!+", line, re.I)
                            or re.search(r"^\s*$", line)
                            or re.search(
                                r"^\s*[A-Z0-9a-z_]+\s+[A-Z0-9a-z_]+\s+[0-9\.\-]+\s+",
                                line,
                            )
                        ):
                            prm_list[context] += line
                        else:
                            prm_list[context] += "!" + line
                    elif context == 9:  # HBOND
                        if not re.search(r"^END", line, re.I):
                            prm_list[context] += line
                    else:
                        continue

    prm = "".join(map(str, prm_list)) + "END"
    prmfh = open(outfile, "w")
    prmfh.write(prm)
    prmfh.close()


# Create string to indicate source of section
def _sec_name(filename):
    return "!Following lines added from %s\n" % (filename)


def split(filename, outdir):
    """ Splits a stream file into an rtf and prm file.

    Parameters
    ----------
    filename : str
        Stream file name
    """
    regex = re.compile("^(toppar_)?(.*)\.str$")
    base = os.path.basename(os.path.normpath(filename))
    base = regex.findall(base)[0][1]
    outrtf = os.path.join(outdir, "top_{}.rtf".format(base))
    outprm = os.path.join(outdir, "par_{}.prm".format(base))

    startrtf = re.compile("^read rtf card", flags=re.IGNORECASE)
    startprm = re.compile("^read para\w* card", flags=re.IGNORECASE)
    endsection = re.compile("^end", flags=re.IGNORECASE)

    rtfsection = 0
    prmsection = 0
    section = "junk"

    rtfstr = ""
    prmstr = ""

    f = open(filename, "r")
    for line in f:
        if startrtf.match(line):
            rtfsection += 1
            if rtfsection > 1:
                rtfstr += "! WARNING -- ANOTHER rtf SECTION FOUND\n"
            section = "rtf"
        elif startprm.match(line):
            prmsection += 1
            if prmsection > 1:
                prmstr += "! WARNING -- ANOTHER para SECTION FOUND\n"
            section = "prm"
        elif endsection.match(line):
            section = "junk"
        elif section == "rtf":
            rtfstr += line
        elif section == "prm":
            prmstr += line
    f.close()

    if rtfsection > 1:
        raise BuildError(
            "Multiple ({}) rtf topology sections found in {} stream file.".format(
                rtfsection, filename
            )
        )
    if prmsection > 1:
        raise BuildError(
            "Multiple ({}) prm parameter sections found in {} stream file.".format(
                prmsection, filename
            )
        )

    f = open(outrtf, "w")
    f.write(rtfstr + "END\n")
    f.close()
    f = open(outprm, "w")
    f.write(prmstr + "END\n")
    f.close()
    return outrtf, outprm


def _prepareStream(filename):
    from htmd.util import tempname

    tmpout = tempname()
    os.makedirs(tmpout)
    return split(filename, tmpout)


def _logParser(fname):
    import re

    unknownres_regex = re.compile("unknown residue type (\w+)")
    failedcoor = re.compile("Warning: failed to set coordinate for atom")
    failedangle = re.compile("Warning: failed to guess coordinate due to bad angle")
    poorlycoor = re.compile("Warning: poorly guessed coordinate(s?)")
    otherwarn = re.compile("Warning")

    failedcoorcount = 0
    failedanglecount = 0
    poorlycoorcount = -1  # Discount the summary report message in the log
    otherwarncount = 0
    unknownres = []
    with open(fname, "r") as f:
        for line in f:
            if failedcoor.search(line):
                failedcoorcount += 1
            elif failedangle.search(line):
                failedanglecount += 1
            elif poorlycoor.search(line):
                poorlycoorcount += 1
            elif otherwarn.search(line):
                otherwarncount += 1
            elif unknownres_regex.search(line):
                unknownres.append(unknownres_regex.findall(line)[0])

    warnings = False
    errors = []
    if failedcoorcount > 0:
        warnings = True
        logger.warning(
            "Failed to set coordinates for {} atoms.".format(failedcoorcount)
        )
    if failedanglecount > 0:
        warnings = True
        logger.warning(
            "Failed to guess coordinates for {} atoms due to bad angles.".format(
                failedanglecount
            )
        )
    if poorlycoorcount > 0:
        warnings = True
        logger.warning(
            "Poorly guessed coordinates for {} atoms.".format(poorlycoorcount)
        )
    if otherwarncount > 0:
        warnings = True
        logger.warning(
            "{} undefined warnings were produced during building.".format(
                otherwarncount
            )
        )
    if len(unknownres):
        errors.append(
            UnknownResidueError(
                "Unknown residue(s) {} found in the input structure. "
                "You are either missing a topology definition for the residue or you need to "
                "rename it to the correct residue name".format(np.unique(unknownres))
            )
        )
    if warnings:
        logger.warning("Please check {} for further information.".format(fname))

    return errors


def _checkFailedAtoms(mol):
    if mol is None:
        return
    idx = np.where(np.sum(mol.coords == 0, axis=1) == 3)[0]
    if len(idx) != 0:
        logger.critical(
            "Atoms with indexes {} are positioned at [0,0,0]. This can cause simulations to crash. "
            "Check log file for more details.".format(idx)
        )


class _TestCharmmBuild(TestCase):
    def test_build(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdbids = ["3PTB", "1A25", "1GZM", "1U5U"]
        for pdb in pdbids:
            with self.subTest(pdb=pdb):
                print("Building {}".format(pdb))
                inFile = os.path.join(
                    preparedInputDir, pdb, "{}-prepared.pdb".format(pdb)
                )
                mol = Molecule(inFile)
                mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

                np.random.seed(1)  # Needed for ions
                smol = solvate(mol)
                topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]
                params = ["par/par_all36_prot.prm", "par/par_water_ions.prm"]
                tmpdir = tempname()
                _ = build(smol, topo=topos, param=params, outdir=tmpdir)

                compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
                assertSameAsReferenceDir(compareDir, tmpdir)

                # shutil.rmtree(tmpdir)

    def test_customDisulfideBonds(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdb = "1GZM"
        inFile = os.path.join(preparedInputDir, pdb, "{}-prepared.pdb".format(pdb))
        mol = Molecule(inFile)
        mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

        np.random.seed(1)  # Needed for ions
        smol = solvate(mol)
        topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]
        params = ["par/par_all36_prot.prm", "par/par_water_ions.prm"]
        disu = [
            ["segid 1 and resid 110", "segid 1 and resid 187"],
            ["segid 0 and resid 110", "segid 0 and resid 187"],
        ]
        tmpdir = tempname()
        _ = build(smol, topo=topos, param=params, outdir=tmpdir, disulfide=disu)

        compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
        assertSameAsReferenceDir(compareDir, tmpdir)

        # TODO: Remove this after deprecation
        from htmd.builder.builder import DisulfideBridge

        np.random.seed(1)  # Needed for ions
        disu = [
            DisulfideBridge("1", 110, "1", 187),
            DisulfideBridge("0", 110, "0", 187),
        ]
        tmpdir = tempname()
        _ = build(smol, topo=topos, param=params, outdir=tmpdir, disulfide=disu)
        compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
        assertSameAsReferenceDir(compareDir, tmpdir)
        # TODO: Remove up to here -----------

    def test_disulfideWithInsertion(self):
        from moleculekit.molecule import Molecule
        from htmd.builder.solvate import solvate
        from htmd.home import home
        from htmd.util import tempname, assertSameAsReferenceDir
        import os
        import numpy as np

        # Use pre-prepared files so we can tell whether the error is in prepare or in build
        # Inputs are reference outputs of proteinprepare.
        preparedInputDir = home(dataDir="test-proteinprepare")

        pdb = "3PTB"

        print("Building {}".format(pdb))
        inFile = os.path.join(preparedInputDir, pdb, "{}-prepared.pdb".format(pdb))
        mol = Molecule(inFile)
        mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing

        np.random.seed(1)  # Needed for ions
        smol = solvate(mol)
        topos = ["top/top_all36_prot.rtf", "top/top_water_ions.rtf"]
        params = ["par/par_all36_prot.prm", "par/par_water_ions.prm"]

        smol.insertion[
            smol.resid == 42
        ] = "A"  # Adding an insertion to test that disulfide bonds with insertions work
        tmpdir = tempname()
        _ = build(smol, topo=topos, param=params, outdir=tmpdir)
        compareDir = home(dataDir=os.path.join("test-charmm-build", "3PTB_insertion"))
        assertSameAsReferenceDir(compareDir, tmpdir)


if __name__ == "__main__":
    import unittest

    unittest.main()

    import doctest

    doctest.testmod()

