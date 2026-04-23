# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function

import numpy as np
import os.path as path
import os
import shutil
import textwrap

from subprocess import call
from htmd.home import home
from moleculekit.molecule import Molecule
from moleculekit.util import _missingSegID
from htmd.builder.builder import (
    detectDisulfideBonds,
    convertDisulfide,
    detectCisPeptideBonds,
    _checkMixedSegment,
    _checkLongResnames,
    MissingResidueError,
    BuildError,
)
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from glob import glob
import unittest


import logging

logger = logging.getLogger(__name__)

_psfgen_exists = shutil.which("psfgen", mode=os.X_OK) is not None


def htmdCharmmHome():
    """Returns the location of the CHARMM files distributed with HTMD"""
    return os.path.abspath(
        os.path.join(home(shareDir=True), "builder", "charmmfiles", "")
    )


def listFiles():
    """Lists all available Charmm topology, parameter and stream files.

    charmm.build only consumes topology/stream files, but parameter files are
    listed too since the user typically needs to copy them (or equivalent
    stream files) into their simulation directory to run MD afterwards.

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
    """Searches for CHARMM files containing a given definition.

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
    """Returns the default topology/stream files used by charmm.build."""
    return [
        "top/top_all36_prot.rtf",
        "top/top_all36_lipid.rtf",
        "top/top_water_ions.rtf",
        "top/top_all36_cgenff.rtf",
        "str/prot/toppar_all36_prot_arg0.str",
    ]


def build(
    mol,
    topo=None,
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
    """Builds a system for CHARMM

    Uses VMD and psfgen to build a system for CHARMM. Additionally it allows for ionization and adding of disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object containing the system
    topo : list of str
        A list of topology files (``.rtf`` or CHARMM ``.str`` stream files).
        Stream files are passed straight to psfgen's ``topology`` command —
        psfgen parses the RTF portion and ignores the ``PARAMETER`` sections.
        Use :func:`charmm.listFiles <htmd.builder.charmm.listFiles>` to get a list of available topology/stream files.
        Default: :func:`defaultTopo <htmd.builder.charmm.defaultTopo>`
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
    >>> topos = ['top/top_all36_prot.rtf', './BEN.rtf', 'top/top_water_ions.rtf']
    >>> disu = [['segid P and resid 157', 'segid P and resid 13'], ['segid K and resid 1', 'segid K and resid 25']]
    >>> ar = {'SAPI24': 'SP24'}  # Alias large resnames to a short-hand version
    >>> molbuilt = charmm.build(mol, topo=topos, outdir='/tmp/build', saltconc=0.15, disulfide=disu, aliasresidues=ar)  # doctest: +SKIP
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
    if caps is None:
        caps = _defaultCaps(mol)
    if noregen is None:
        noregen = ["FHEM", "PHEM", "PLOH", "PLO2", "PLIG", "PSUL"]
    if patches is None:
        patches = []
    if isinstance(patches, str):
        patches = [patches]

    # Stage topology/stream files for psfgen. Stream files are handed to
    # psfgen directly via its `topology` command — psfgen parses the RTF
    # section and ignores the parameter section. Parameter files are *not*
    # staged here; the user is responsible for copying the parameter/stream
    # files they need into their simulation directory afterwards.
    topology_local_paths = _stage_files(list(topo), outdir, "topologies")

    # Apply user-requested residue aliases up-front so every subsequent step
    # (segment writing, psfgen script, protonation-patch detection) sees the
    # aliased names.
    if aliasresidues is not None:
        for key, val in aliasresidues.items():
            mol.resname[mol.resname == key] = val

    allpatches = list(patches) + _protonationPatches(mol)

    # Resolve disulfide bonds once against the full mol (residues that take
    # part in disulfides are always in the solute).
    if (
        disulfide is not None
        and len(disulfide) != 0
        and isinstance(disulfide[0][0], str)
    ):
        disulfide = convertDisulfide(mol, disulfide)
    if disulfide is None:
        disulfide = detectDisulfideBonds(mol)

    water_sel = mol.atomselect("water")
    has_water = bool(np.any(water_sel))
    do_ionize = ionize and execute and has_water

    # Fast solute-only psfgen pass to compute the total charge, mirroring the
    # amber.build flow. This avoids running the expensive full build twice
    # (once for charge, once for the ionized system).
    if do_ionize:
        logger.info("Running solute-only psfgen pass to compute total charge.")
        solute_mol = mol.copy(sel=~water_sel)
        solvent_mol = mol.copy(sel=water_sel)

        _write_segments(solute_mol, outdir)
        _write_psfgen_script(
            outdir,
            "build_solute.vmd",
            topology_local_paths,
            solute_mol,
            caps=caps,
            disulfide=disulfide,
            allpatches=allpatches,
            noregen=noregen,
            regenerate=regenerate,
            aliasresidues=aliasresidues,
            prefix="solute_charge",
        )
        molbuilt_solute = _run_psfgen(
            psfgen, outdir, "build_solute.vmd", "solute_charge"
        )

        totalcharge = float(np.sum(molbuilt_solute.charge))
        nwater = int(np.sum(solvent_mol.atomselect("water and noh")))

        # Remove solute-only build artifacts so the final build writes clean.
        for ext in ("psf", "pdb"):
            fpath = path.join(outdir, f"solute_charge.{ext}")
            if path.exists(fpath):
                os.remove(fpath)

        anion, cation, anionatom, cationatom, nanion, ncation = ionizef(
            molbuilt_solute,
            totalcharge,
            nwater,
            saltconc=saltconc,
            anion=saltanion,
            cation=saltcation,
        )
        solvent_mol = ionizePlace(
            solvent_mol,
            solute_mol,
            anion,
            cation,
            anionatom,
            cationatom,
            nanion,
            ncation,
        )
        mol = solute_mol.copy()
        mol.append(solvent_mol)

    # Write segments and the full build script for the (possibly ionized) mol.
    _write_segments(mol, outdir)
    _write_psfgen_script(
        outdir,
        "build.vmd",
        topology_local_paths,
        mol,
        caps=caps,
        disulfide=disulfide,
        allpatches=allpatches,
        noregen=noregen,
        regenerate=regenerate,
        aliasresidues=aliasresidues,
        prefix=prefix,
    )

    if not execute:
        return None

    molbuilt = _run_psfgen(psfgen, outdir, "build.vmd", prefix)
    _checkFailedAtoms(molbuilt)
    _recoverProtonations(molbuilt)
    detectCisPeptideBonds(molbuilt, respect_bonds=True)
    return molbuilt


def _stage_files(files, outdir, subdir):
    """Copy a list of CHARMM topology/parameter/stream files into ``outdir/subdir``.

    Each file is resolved against ``htmdCharmmHome()`` when it isn't an absolute
    path / doesn't start with ".", and copied with a numeric-prefixed basename
    so that order is preserved and collisions avoided.

    Returns a list of paths relative to ``outdir`` (e.g. ``topologies/0.top_all36_prot.rtf``).
    """
    charmmdir = htmdCharmmHome()
    target = path.join(outdir, subdir)
    if not path.exists(target):
        os.makedirs(target)

    local_paths = []
    for i, f in enumerate(files):
        if f[0] != "." and path.isfile(path.join(charmmdir, f)):
            f = path.join(charmmdir, f)
        if not path.isfile(f):
            raise FileNotFoundError(
                f"File {f} does not exist. Cannot stage to {target}."
            )
        localname = f"{i}.{path.basename(f)}"
        shutil.copy(f, path.join(target, localname))
        local_paths.append(path.join(subdir, localname))
    return local_paths


def _write_segments(mol, outdir):
    """Write one segment PDB per unique segid into ``outdir/segments/``."""
    segdir = path.join(outdir, "segments")
    if not path.exists(segdir):
        os.makedirs(segdir)
    logger.info("Writing out segments.")
    for seg in _getSegments(mol):
        pdbname = f"segment{seg}.pdb"
        mol.write(path.join(segdir, pdbname), sel=mol.segid == seg)


def _write_psfgen_script(
    outdir,
    script_name,
    topology_local_paths,
    mol,
    caps,
    disulfide,
    allpatches,
    noregen,
    regenerate,
    aliasresidues,
    prefix,
):
    """Write a psfgen input script that rebuilds the structure contained in ``mol``."""
    wateratoms = mol.atomselect("water")
    with open(path.join(outdir, script_name), "w") as f:
        f.write("# psfgen file generated by charmm.build\n")
        f.write("package require psfgen;\n")
        f.write("psfcontext reset;\n\n")

        for p in topology_local_paths:
            f.write(f"topology {p}\n")
        f.write("\n")

        _printAliases(f)
        if aliasresidues is not None:
            for key, val in aliasresidues.items():
                f.write(f"        pdbalias residue {val} {key}\n")

        for seg in _getSegments(mol):
            pdbname = f"segment{seg}.pdb"
            segatoms = mol.segid == seg
            segwater = wateratoms & segatoms
            f.write(f"segment {seg} {{\n")
            if np.all(segatoms == segwater):
                # Pure-water segments disable psfgen's angle/dihedral autogen
                # (TIP3 residues carry their own topology).
                f.write("\tauto none\n")
            f.write(f"\tpdb {path.join('segments', pdbname)}\n")
            if caps is not None and seg in caps:
                for c in caps[seg]:
                    f.write(f"\t{c}\n")
            f.write("}\n")
            f.write(f"coordpdb {path.join('segments', pdbname)} {seg}\n\n")

        if disulfide is not None and len(disulfide):
            for d in sorted(disulfide, key=lambda x: x[0].segid):
                str0 = f"{d[0].segid}:{d[0].resid}{d[0].insertion}"
                str1 = f"{d[1].segid}:{d[1].resid}{d[1].insertion}"
                f.write(f"patch DISU {str0} {str1}\n")
            f.write("\n")

        noregenpatches = [p for p in allpatches if p.split()[1] in noregen]
        regenpatches = [p for p in allpatches if p.split()[1] not in noregen]

        if regenpatches:
            for p in regenpatches:
                f.write(p + "\n")
            f.write("\n")

        if regenerate is not None:
            f.write("regenerate {}\n\n".format(" ".join(regenerate)))

        if noregenpatches:
            for p in noregenpatches:
                f.write(p + "\n")
            f.write("\n")

        f.write("guesscoord\n")
        f.write(f"writepsf {prefix}.psf\n")
        f.write(f"writepdb {prefix}.pdb\n")


def _run_psfgen(psfgen, outdir, script_name, prefix):
    """Invoke psfgen on ``outdir/script_name`` and return the resulting Molecule."""
    logpath = os.path.abspath(path.join(outdir, "log.txt"))
    logger.info("Starting the build.")
    currdir = os.getcwd()
    os.chdir(outdir)
    try:
        my_env = os.environ.copy()
        my_env["LC_ALL"] = "C"
        with open(logpath, "w") as f:
            call([psfgen, f"./{script_name}"], stdout=f, stderr=f, env=my_env)
        errors = _logParser(logpath)
    finally:
        os.chdir(currdir)

    if errors:
        raise BuildError(
            errors
            + [f"Check {logpath} for further information on errors in building."],
            errors,
        )
    logger.info("Finished building.")

    pdb_path = path.join(outdir, f"{prefix}.pdb")
    psf_path = path.join(outdir, f"{prefix}.psf")
    if not (path.isfile(pdb_path) and path.isfile(psf_path)):
        raise BuildError(
            f"No {prefix} pdb/psf file was generated. Check {logpath} for errors in building."
        )

    molbuilt = Molecule(pdb_path)
    molbuilt.read(psf_path)
    return molbuilt


def _cleanOutDir(outdir):
    from glob import glob

    files = glob(os.path.join(outdir, "structure.*"))
    files += glob(os.path.join(outdir, "log.*"))
    files += glob(os.path.join(outdir, "*.log"))
    files += glob(os.path.join(outdir, "*.vmd"))
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


def _logParser(fname):
    import re

    unknownres_regex = re.compile(r"unknown residue type (\w+)")
    failedcoor = re.compile(r"Warning: failed to set coordinate for atom")
    failedangle = re.compile(r"Warning: failed to guess coordinate due to bad angle")
    poorlycoor = re.compile(r"Warning: poorly guessed coordinate(s?)")
    otherwarn = re.compile(r"Warning")

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
        logger.warning(f"Failed to set coordinates for {failedcoorcount} atoms.")
    if failedanglecount > 0:
        warnings = True
        logger.warning(
            f"Failed to guess coordinates for {failedanglecount} atoms due to bad angles."
        )
    if poorlycoorcount > 0:
        warnings = True
        logger.warning(f"Poorly guessed coordinates for {poorlycoorcount} atoms.")
    if otherwarncount > 0:
        warnings = True
        logger.warning(
            f"{otherwarncount} undefined warnings were produced during building."
        )
    if len(unknownres):
        errors.append(
            MissingResidueError(
                f"Unknown residue(s) {np.unique(unknownres)} found in the input structure. "
                "You are either missing a topology definition for the residue or you need to "
                "rename it to the correct residue name",
                unknownres,
            )
        )
    if warnings:
        logger.warning(f"Please check {fname} for further information.")

    return errors


def _checkFailedAtoms(mol):
    if mol is None:
        return
    idx = np.where(np.sum(mol.coords == 0, axis=1) == 3)[0]
    if len(idx) != 0:
        logger.critical(
            f"Atoms with indexes {idx} are positioned at [0,0,0]. This can cause simulations to crash. "
            "Check log file for more details."
        )


class _TestCharmmBuild(unittest.TestCase):
    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
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
                # 1U5U contains deprotonated arginines (AR0 → patch RN1),
                # which are only defined in the arg0 stream file.
                if pdb == "1U5U":
                    topos += [
                        "str/prot/toppar_all36_prot_arg0.str",
                        "str/misc/toppar_ions_won.str",
                    ]
                tmpdir = tempname()
                _ = build(smol, topo=topos, outdir=tmpdir)

                compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
                assertSameAsReferenceDir(compareDir, tmpdir)

                # shutil.rmtree(tmpdir)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
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
        disu = [
            ["segid 1 and resid 110", "segid 1 and resid 187"],
            ["segid 0 and resid 110", "segid 0 and resid 187"],
        ]
        tmpdir = tempname()
        _ = build(smol, topo=topos, outdir=tmpdir, disulfide=disu)

        compareDir = home(dataDir=os.path.join("test-charmm-build", pdb))
        assertSameAsReferenceDir(compareDir, tmpdir)

    @unittest.skipUnless(_psfgen_exists, "Requires psfgen")
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

        smol.insertion[smol.resid == 42] = (
            "A"  # Adding an insertion to test that disulfide bonds with insertions work
        )
        tmpdir = tempname()
        _ = build(smol, topo=topos, outdir=tmpdir)
        compareDir = home(dataDir=os.path.join("test-charmm-build", "3PTB_insertion"))
        assertSameAsReferenceDir(compareDir, tmpdir)


if __name__ == "__main__":
    import unittest

    unittest.main()
