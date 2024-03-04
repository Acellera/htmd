# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function


from htmd.home import home
import numpy as np
import os
from os.path import join
from glob import glob
from moleculekit.util import _missingSegID, sequenceID
import shutil
from htmd.builder.builder import (
    detectDisulfideBonds,
    detectCisPeptideBonds,
    convertDisulfide,
    _checkMixedSegment,
    BuildError,
    MissingResidueError,
    MissingAngleError,
    MissingBondError,
    MissingParameterError,
    MissingTorsionError,
    MissingAtomTypeError,
)
from subprocess import call
from moleculekit.molecule import Molecule
from moleculekit.tools.sequencestructuralalignment import sequenceStructureAlignment
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.util import ensurelist
import unittest
import sys
import logging

logger = logging.getLogger(__name__)

# All atom types reserved by GAFF2. We will use all other low caps for parameterize
gaff2_at = [
    "h1",
    "h2",
    "h3",
    "h4",
    "h5",
    "ha",
    "hc",
    "hn",
    "ho",
    "hp",
    "hs",
    "hw",
    "hx",
    "o",
    "o2",
    "oh",
    "op",
    "oq",
    "os",
    "ow",
    "c",
    "c1",
    "c2",
    "c3",
    "ca",
    "cc",
    "cd",
    "ce",
    "cf",
    "cg",
    "ch",
    "cp",
    "cs",
    "cq",
    "cu",
    "cv",
    "cx",
    "cy",
    "cz",
    "n",
    "n1",
    "n2",
    "n3",
    "n4",
    "n5",
    "n6",
    "n7",
    "n8",
    "n9",
    "na",
    "nb",
    "nc",
    "nd",
    "ne",
    "nf",
    "nh",
    "ni",
    "nj",
    "nk",
    "nl",
    "nm",
    "nn",
    "no",
    "np",
    "nq",
    "ns",
    "nt",
    "nu",
    "nv",
    "nx",
    "ny",
    "nz",
    "n+",
    "s",
    "s2",
    "s3",
    "s4",
    "s6",
    "sh",
    "sp",
    "sq",
    "ss",
    "sx",
    "sy",
    "p2",
    "p3",
    "p4",
    "p5",
    "pb",
    "pc",
    "pd",
    "pe",
    "pf",
    "px",
    "py",
    "f",
    "cl",
    "br",
    "i",
]


def _get_parameterize_at():
    """The GAFF and GAFF2 atom types are in lower case.
    This is the mechanism by which the GAFF force field is kept independent from the macromolecular AMBER force fields.
    All of the traditional AMBER force fields use uppercase atom types.
    In this way the GAFF and traditional force fields can be mixed in the same calculation.
    Taken from: http://ambermd.org/tutorials/basic/tutorial4b/
    """
    import string

    for x in string.ascii_lowercase[::-1]:
        for y in string.ascii_lowercase + string.digits:
            if x + y not in gaff2_at:
                yield x + y

    for x in string.ascii_lowercase[::-1]:
        if x not in gaff2_at:
            yield x


parameterize_types = list(_get_parameterize_at())


def _findTeLeap():
    teleap = shutil.which("teLeap", mode=os.X_OK)
    if not teleap:
        raise FileNotFoundError(
            "teLeap not found. You should have AmberTools installed "
            "(to install AmberTools do: conda install ambertools -c conda-forge)"
        )
    if os.path.islink(teleap):
        if os.path.isabs(os.readlink(teleap)):
            teleap = os.readlink(teleap)
        else:
            teleap = os.path.join(os.path.dirname(teleap), os.readlink(teleap))
    return teleap


def defaultAmberHome(teleap=None):
    """Returns the default AMBERHOME as defined by the location of teLeap binary

    Parameters:
    -----------
    teleap : str
        Path to teLeap executable used to build the system for AMBER
    """
    if teleap is None:
        teleap = _findTeLeap()
    else:
        if shutil.which(teleap) is None:
            raise NameError(f"Could not find executable: `{teleap}` in the PATH.")

    return os.path.normpath(os.path.join(os.path.dirname(teleap), "../"))


_defaultAmberSearchPaths = {
    "ff": join("dat", "leap", "cmd"),
    "topo": join("dat", "leap", "prep"),
    "param": join("dat", "leap", "parm"),
    "lib": join("dat", "leap", "lib"),
}


def _convert_in_file_to_prepis(in_file, outdir):
    "This was used to convert the .in files in cofactors/ff-ncaa/ff-ptm used in _cofactors_noncanonical_ptm_params to prepi files"
    with open(in_file, "r") as f:
        mollines = []
        name = None
        for line in f.readlines():
            mollines.append(line)
            if line.strip().startswith("CORR"):
                name = mollines[-2].strip().split()[0]
                if len(mollines) != 6:
                    # Add empty lines to start if they don't exist
                    mollines = ["\n"] * (6 - len(mollines)) + mollines
            if line.strip() == "DONE":
                with open(os.path.join(outdir, f"{name.upper()}.prepi"), "w") as fout:
                    for ll in mollines:
                        fout.write(ll)
                    fout.write("STOP")

                name = None
                mollines = []


def _cofactors_ncaa_ptm_params():
    res_dict = {}
    htmdamberdir = htmdAmberHome()

    dirs = [os.path.join(htmdamberdir, dd) for dd in ("cofactors", "ff-ncaa", "ff-ptm")]
    for dir in dirs:
        for cc in glob(os.path.join(dir, "*.frcmod")):
            basename = os.path.splitext(os.path.basename(cc))[0]
            if basename in res_dict:
                raise AssertionError(
                    f"Duplicate residue parameters detected for {basename} in {res_dict[basename][0]} and {cc}"
                )
            prepi = os.path.join(dir, f"{basename}.prepi")
            if not os.path.exists(prepi):
                raise AssertionError(
                    f"prepi file should exist for frcmod {cc} in folder {os.path.dirname(cc)}"
                )
            res_dict[basename] = (cc, prepi, os.path.basename(dir))
    return res_dict


def _detect_cofactors_ncaa_ptm(mol, param, topo):
    import itertools
    import string

    names = {
        "cofactors": "Cofactor(s)",
        "ff-ncaa": "Non-canonical amino-acid(s)",
        "ff-ptm": "Post-translational modification(s)",
    }

    segid_gen = itertools.product(*[["c"], string.digits + string.ascii_lowercase])

    # Detect known cofactors, non-canonical AA and post-translational modifications and add the required files to the build
    cofactors_etc = _cofactors_ncaa_ptm_params()
    param_basenames = [os.path.splitext(os.path.basename(ff))[0] for ff in param]
    uqresid = sequenceID((mol.resid, mol.insertion, mol.segid))
    detected = {k: [] for k in names}
    for ncres in sorted(cofactors_etc):
        cof = cofactors_etc[ncres]
        if any(mol.resname == ncres) and ncres not in param_basenames:
            detected[cof[2]].append(ncres)
            param.append(cof[0])
            topo.append(cof[1])
            if cof[2] in ("cofactors", "ff-ptm"):
                # Remove hydrogens to avoid naming issues
                h_sel = (mol.resname == ncres) & (mol.element == "H")

            if cof[2] == "cofactors":
                uqresid = uqresid[~h_sel]
                mol.remove(h_sel, _logger=False)
                # Fix segids of cofactors
                uq_rid = np.unique(uqresid[mol.resname == ncres])
                for rid in uq_rid:
                    new_segid = "".join(next(segid_gen))
                    while new_segid in mol.segid:
                        new_segid = "".join(next(segid_gen))
                    mol.segid[(mol.resname == ncres) & (uqresid == rid)] = new_segid

    for key in names:
        if len(detected[key]):
            logger.info(
                f"{names[key]} {', '.join(detected[key])} detected in system. Automatically adding parameters and topology for AMBER."
            )


def htmdAmberHome():
    """Returns the location of the AMBER files distributed with HTMD"""

    return os.path.abspath(os.path.join(home(shareDir=True), "builder", "amberfiles"))


def listFiles():
    """Lists all AMBER forcefield files available in HTMD

    Example
    -------
    >>> from htmd.builder import amber
    >>> amber.listFiles()             # doctest: +ELLIPSIS
    ---- Forcefield files list: ... ----
    leaprc.amberdyes
    leaprc.conste
    leaprc.constph
    leaprc.DNA.bsc1
    ...
    """

    amberhome = defaultAmberHome()

    # Original AMBER FFs
    ffdir = join(amberhome, _defaultAmberSearchPaths["ff"])
    ffs = glob(join(ffdir, "*"))
    print("---- Forcefield files list: " + join(ffdir, "") + " ----")
    for f in sorted(ffs, key=str.lower):
        if os.path.isdir(f):
            continue
        print(f.replace(join(ffdir, ""), ""))

    oldffdir = join(amberhome, _defaultAmberSearchPaths["ff"], "oldff")
    ffs = glob(join(oldffdir, "*"))
    print("---- OLD Forcefield files list: " + join(ffdir, "") + " ----")
    for f in sorted(ffs, key=str.lower):
        print(f.replace(join(ffdir, ""), ""))

    topodir = os.path.join(amberhome, _defaultAmberSearchPaths["topo"])
    topos = glob(join(topodir, "*"))
    print("---- Topology files list: " + join(topodir, "") + " ----")
    for f in sorted(topos, key=str.lower):
        if os.path.isdir(f):
            continue
        print(f.replace(join(topodir, ""), ""))

    # FRCMOD files
    frcmoddir = os.path.join(amberhome, _defaultAmberSearchPaths["param"])
    ffs = glob(join(frcmoddir, "frcmod.*"))
    print("---- Parameter files list: " + join(frcmoddir, "") + " ----")
    for f in sorted(ffs, key=str.lower):
        print(f.replace(join(frcmoddir, ""), ""))

    htmdamberdir = htmdAmberHome()

    # Extra AMBER FFs on HTMD
    extraffs = glob(join(htmdamberdir, "*", "leaprc.*"))
    print("---- Extra forcefield files list: " + join(htmdamberdir, "") + " ----")
    for f in sorted(extraffs, key=str.lower):
        print(f.replace(join(htmdamberdir, ""), ""))

    # Extra AMBER FFs on HTMD (*.frcmod, *.in) @cuzzo87
    extratopos = glob(join(htmdamberdir, "*", "*.in"))
    print("---- Extra topology files list: " + join(htmdamberdir, "") + " ----")
    for f in sorted(extratopos, key=str.lower):
        print(f.replace(join(htmdamberdir, ""), ""))

    extraparams = glob(join(htmdamberdir, "*", "*.frcmod"))
    print("---- Extra parameter files list: " + join(htmdamberdir, "") + " ----")
    for f in sorted(extraparams, key=str.lower):
        print(f.replace(join(htmdamberdir, ""), ""))


def _locateFile(fname, ftype, teleap):
    amberhome = defaultAmberHome(teleap=teleap)
    htmdamberdir = htmdAmberHome()
    searchdir = os.path.join(amberhome, _defaultAmberSearchPaths[ftype])
    foundfile = glob(os.path.join(searchdir, fname))
    if len(foundfile) != 0:
        return foundfile[0]
    foundfile = glob(os.path.join(htmdamberdir, fname))
    if len(foundfile) != 0:
        return foundfile[0]
    logger.warning(f"Was not able to find {ftype} file {fname}")


def _getTeLeapImportFlags(teleap=None, ff=None, teleapimports=()):
    if not teleapimports:
        teleapimports = []
        # Source default Amber (i.e. the same paths tleap imports)
        amberhome = defaultAmberHome(teleap=teleap)
        teleapimports += [
            os.path.join(amberhome, s) for s in _defaultAmberSearchPaths.values()
        ]
        if len(teleapimports) == 0:
            raise RuntimeWarning(
                f"No default Amber force-field found. Check teLeap location: {teleap}"
            )
        # Source HTMD Amber paths that contain ffs
        htmdamberdir = htmdAmberHome()
        teleapimports += [
            os.path.join(htmdamberdir, os.path.dirname(f))
            for f in ff
            if os.path.isfile(os.path.join(htmdamberdir, f))
        ]
        if len(teleapimports) == 0:
            raise RuntimeError(
                "No default Amber force-field imports found. Check "
                "`htmd.builder.amber.defaultAmberHome()` and `htmd.builder.amber.htmdAmberHome()`"
            )

    # Set import flags for teLeap
    teleapimportflags = []
    for p in teleapimports:
        teleapimportflags.append("-I")
        teleapimportflags.append(str(p))
    return teleapimportflags


def _getResiduesInAllFFs():
    amberhome = defaultAmberHome()
    # Original AMBER FFs
    ffdir = join(amberhome, _defaultAmberSearchPaths["ff"])
    ffs = [
        os.path.relpath(ff, ffdir)
        for ff in glob(join(ffdir, "*"))
        if os.path.isfile(ff)
    ]

    htmdamberdir = htmdAmberHome()
    # Extra AMBER FFs on HTMD
    ffs += [
        os.path.relpath(ff, htmdamberdir)
        for ff in glob(join(htmdamberdir, "*", "leaprc.*"))
        if os.path.isfile(ff)
    ]

    residues = {}
    for ff in ffs:
        residues[ff] = _getResiduesInFF(ff)

    return residues


def _getResiduesInFF(ff, teleap=None):
    import tempfile

    if teleap is None:
        teleap = _findTeLeap()

    with tempfile.TemporaryDirectory() as tmpdir:
        ff = _locateFile(ff, "ff", teleap)
        shutil.copy(ff, tmpdir)

        importflgs = _getTeLeapImportFlags(teleap, [ff])

        with open(os.path.join(tmpdir, "tleap.in"), "w") as f:
            f.write(f"logFile leap.log\nsource {os.path.basename(ff)}\nlist\nquit")

        logpath = os.path.join(tmpdir, "log.txt")
        outlog = os.path.join(tmpdir, "leap.log")
        with open(logpath, "w") as f:
            try:
                cmd = [teleap, "-f", "./tleap.in"]
                cmd[1:1] = importflgs
                logger.debug(cmd)
                call(cmd, stdout=f, cwd=tmpdir)
            except Exception:
                raise NameError("teLeap failed at execution")

        residues = []
        starting = False

        with open(outlog, "r") as f:
            for line in f:
                if line.strip() == "> list":
                    starting = True
                    continue
                if line.strip() == "quit":
                    break
                if starting:
                    residues += line.strip().split()
    return residues


def defaultFf():
    """Returns the default leaprc forcefield files used by amber.build"""
    return [
        "leaprc.protein.ff14SB",
        "leaprc.lipid21",
        "leaprc.gaff2",
        "leaprc.RNA.Shaw",
        "leaprc.DNA.bsc1",
        "leaprc.water.tip3p",
    ]


def defaultTopo():
    """Returns the default topology `prepi` files used by amber.build"""
    return []


def defaultParam():
    """Returns the default parameter `frcmod` files used by amber.build"""
    return []


def _prepareMolecule(mol: Molecule, caps, disulfide):
    # Remove pdb bonds as they can be regenerated by tleap
    mol.deleteBonds("all")

    # Check for missing segids or mixed protein / non-protein segments
    _missingSegID(mol)
    _checkMixedSegment(mol)

    # Convert lipids to AMBER naming
    mol = _charmmLipid2Amber(mol)

    cyclic = _detect_cyclic_segments(mol)

    # Add caps to termini
    if caps is None:
        caps = _defaultProteinCaps(mol)

    for cc, _, _ in cyclic:
        if cc in caps:
            logger.info(f"Found cyclic segment {cc}. Disabling capping on it.")
            del caps[cc]

    _add_caps(mol, caps)

    # Before modifying the resids copy the molecule to map back the disulfide bonds
    mol_orig_resid = mol.copy()

    # We need to renumber residues to unique numbers for teLeap. It goes up to 9999 and restarts from 0
    mol.resid[:] = (sequenceID((mol.resid, mol.insertion, mol.segid)) + 1) % 10000

    # Map the old disulfide UniqueResidueIDs to the new ones
    if disulfide is not None and len(disulfide):
        if not isinstance(disulfide[0][0], str):
            raise RuntimeError("Can only accept disulfide bond strings")
        disulfide = convertDisulfide(mol_orig_resid, disulfide)
        disulfide = _mapDisulfide(disulfide, mol, mol_orig_resid)
    else:
        logger.info("Detecting disulfide bonds.")
        disulfide = detectDisulfideBonds(mol)

    # Fix structure to match the disulfide patching
    if len(disulfide) != 0:
        torem = np.zeros(mol.numAtoms, dtype=bool)
        for d1, d2 in disulfide:
            # Rename the residues to CYX if there is a disulfide bond
            atoms1 = d1.selectAtoms(mol, indexes=False)
            atoms2 = d2.selectAtoms(mol, indexes=False)
            d1.resname = "CYX"
            d2.resname = "CYX"
            mol.resname[atoms1] = "CYX"
            mol.resname[atoms2] = "CYX"
            # Remove (eventual) HG hydrogens on these CYS (from systemPrepare)
            torem |= (atoms1 & (mol.name == "HG")) | (atoms2 & (mol.name == "HG"))
        mol.remove(torem, _logger=False)

    return disulfide


def _write_residue_mapping(molbuilt, mol_orig, outdir):
    # Align with original molecule to rename residues back to original resids
    try:
        _, mapping = sequenceStructureAlignment(molbuilt, mol_orig, maxalignments=1)
        if len(mapping):
            mol_map, ref_map = mapping[0]  # Top alignment
            idx_mol = np.where(mol_map)[0]
            idx_ref = np.where(ref_map)[0]
            residmap = []
            for im, ir in zip(idx_mol, idx_ref):
                residmap.append(
                    [
                        str(molbuilt.resid[im]),
                        str(mol_orig.resid[ir]),
                        mol_orig.insertion[ir],
                        mol_orig.chain[ir],
                    ]
                )
            with open(os.path.join(outdir, "residue_mapping.csv"), "w") as fcsv:
                fcsv.write("new_resid,old_resid,old_insertion,old_chain\n")
                for mm in residmap:
                    fcsv.write(",".join(mm) + "\n")
    except Exception:
        pass


def build(
    mol,
    ff=None,
    topo=None,
    param=None,
    prefix="structure",
    outdir="./build",
    caps=None,
    ionize=True,
    saltconc=0,
    saltanion=None,
    saltcation=None,
    disulfide=None,
    teleap=None,
    teleapimports=None,
    execute=True,
    atomtypes=None,
    offlibraries=None,
    gbsa=False,
    igb=2,
):
    """Builds a system for AMBER

    Uses tleap to build a system for AMBER. Additionally it allows the user to ionize and add disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object containing the system
    ff : list of str
        A list of leaprc forcefield files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available forcefield files.
        Default: :func:`amber.defaultFf <htmd.builder.amber.defaultFf>`
    topo : list of str
        A list of topology `prepi/prep/in/mol2/cif` files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available topology files.
        Default: :func:`amber.defaultTopo <htmd.builder.amber.defaultTopo>`
        When passing residues parameterized with the `parameterize` tool, please pass the .cif file.
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
        A dictionary specifying the caps. It accepts two different formats.
        One with keys segids and values lists of strings describing the caps for a particular protein segment.
        e.g. caps['P'] = ['ACE', 'NME'] or caps['P'] = ['none', 'none']. Default: will apply ACE and NME caps to every
        protein segment.
        The second format is more manual and uses as keys an atomselection of a residue and as values the cap
        to apply to that residue. e.g. caps = {"chain A and resid 5": "ACE", "chain A and resid 10": "NME"}
    ionize : bool
        Enable or disable ionization
    saltconc : float
        Salt concentration to add to the system after neutralization.
    saltanion : {'Cl-'}
        The anion type. Please use only AMBER ion atom names.
    saltcation : {'Na+', 'K+', 'Cs+'}
        The cation type. Please use only AMBER ion atom names.
    disulfide : list of pairs of atomselection strings
        If None it will guess disulfide bonds. Otherwise provide a list pairs of atomselection strings for each pair of
        residues forming the disulfide bridge.
    teleap : str
        Path to teLeap executable used to build the system for AMBER
    teleapimports : list
        A list of paths to pass to teLeap '-I' flag, i.e. directories to be searched
        Default: determined from :func:`amber.defaultAmberHome <htmd.builder.amber.defaultAmberHome>` and
        :func:`amber.htmdAmberHome <htmd.builder.amber.htmdAmberHome>`
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
    molbuilt : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The built system in a Molecule object

    Example
    -------
    >>> from htmd.ui import *  # doctest: +SKIP
    >>> mol = Molecule("3PTB")
    >>> molbuilt = amber.build(mol, outdir='/tmp/build')  # doctest: +SKIP
    >>> # More complex example
    >>> disu = [['segid P and resid 157', 'segid P and resid 13'], ['segid K and resid 1', 'segid K and resid 25']]
    >>> molbuilt = amber.build(mol, outdir='/tmp/build', saltconc=0.15, disulfide=disu)  # doctest: +SKIP
    """
    mol_orig = mol
    mol = mol.copy()

    if ff is None:
        ff = defaultFf()
    if topo is None:
        topo = defaultTopo()
    if param is None:
        param = defaultParam()

    # Warning for bad FF combination
    ff19tip3p = [False, False]
    for x in ff:
        if "ff19SB" in x:
            ff19tip3p[0] = True
        if "tip3p" in x:
            ff19tip3p[1] = True
    if all(ff19tip3p):
        logger.warning(
            "CAUTION: AMBER Forcefield ff19SB is NOT compatible with TIP3P water model."
            " Consider using ff14SB instead or using the OPC water model."
        )

    disulfide = _prepareMolecule(mol, caps, disulfide)
    _detect_cofactors_ncaa_ptm(mol, param, topo)

    if ionize:
        molbuilt = _build(
            mol,
            ff=ff,
            topo=topo,
            param=param,
            prefix=prefix,
            outdir=outdir,
            disulfide=disulfide,
            teleap=teleap,
            teleapimports=teleapimports,
            execute=execute,
            atomtypes=atomtypes,
            offlibraries=offlibraries,
            gbsa=gbsa,
            igb=igb,
        )
        shutil.move(
            os.path.join(outdir, "structure.crd"),
            os.path.join(outdir, "structure.noions.crd"),
        )
        shutil.move(
            os.path.join(outdir, "structure.prmtop"),
            os.path.join(outdir, "structure.noions.prmtop"),
        )
        totalcharge = np.sum(molbuilt.charge)
        nwater = np.sum(molbuilt.atomselect("water and noh"))
        anion, cation, anionatom, cationatom, nanion, ncation = ionizef(
            molbuilt,
            totalcharge,
            nwater,
            saltconc=saltconc,
            anion=saltanion,
            cation=saltcation,
        )
        mol = ionizePlace(mol, anion, cation, anionatom, cationatom, nanion, ncation)

    molbuilt = _build(
        mol,
        ff=ff,
        topo=topo,
        param=param,
        prefix=prefix,
        outdir=outdir,
        disulfide=disulfide,
        teleap=teleap,
        teleapimports=teleapimports,
        execute=execute,
        atomtypes=atomtypes,
        offlibraries=offlibraries,
        gbsa=gbsa,
        igb=igb,
    )
    _write_residue_mapping(molbuilt, mol_orig, outdir)
    return molbuilt


def _create_atomtype_section(mol):
    atypes_str = "addAtomTypes {\n"
    for at in sorted(np.unique(mol.atomtype)):
        el = mol.element[mol.atomtype == at][0]
        atypes_str += f'	{{ "{at}"  "{el}" "sp3" }}\n'
    atypes_str += "}\n"
    return atypes_str


def _build(
    mol,
    ff=None,
    topo=None,
    param=None,
    prefix="structure",
    outdir="./build",
    disulfide=None,
    teleap=None,
    teleapimports=None,
    execute=True,
    atomtypes=None,
    offlibraries=None,
    gbsa=False,
    igb=2,
):
    if teleap is None:
        teleap = _findTeLeap()
    else:
        if shutil.which(teleap) is None:
            raise NameError(
                f"Could not find executable: `{teleap}` in the PATH. Cannot build for AMBER. Please install it with `conda install ambertools -c conda-forge`"
            )

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    _cleanOutDir(outdir)

    with open(os.path.join(outdir, "tleap.in"), "w") as f:
        f.write("# tleap file generated by amber.build\n")

        # Printing out the forcefields
        for i, force in enumerate(ensurelist(ff)):
            if not os.path.isfile(force):
                force = _locateFile(force, "ff", teleap)
                if force is None:
                    continue
            newname = f"ff{i}_{os.path.basename(force)}"
            shutil.copy(force, os.path.join(outdir, newname))
            f.write(f"source {newname}\n")
        f.write("\n")

        if gbsa:
            gbmodels = {
                1: "mbondi",
                2: "mbondi2",
                5: "mbondi2",
                7: "bondi",
                8: "mbondi3",
            }
            f.write(f"set default PBradii {gbmodels[igb]}\n\n")

        # Adding custom atom types
        if atomtypes is not None:
            atomtypes = ensurelist(tocheck=atomtypes[0], tomod=atomtypes)
            f.write("addAtomTypes {\n")
            for at in atomtypes:
                if len(at) != 3:
                    raise RuntimeError(
                        "Atom type definitions have to be triplets. Check the AMBER documentation."
                    )
                f.write(f'    {{ "{at[0]}" "{at[1]}" "{at[2]}" }}\n')
            f.write("}\n\n")

        # Loading OFF libraries
        if offlibraries is not None:
            offlibraries = ensurelist(offlibraries)
            for i, off in enumerate(offlibraries):
                if not os.path.isfile(off):
                    raise RuntimeError(f"Could not find off-library in location {off}")
                newname = f"offlib{i}_{os.path.basename(off)}"
                shutil.copy(off, os.path.join(outdir, newname))
                f.write(f"loadoff {newname}\n")

        # Copy param and topo files to output folder and rename to unique names
        newparam = []
        for i, fname in enumerate(sorted(param, key=lambda x: os.path.basename(x))):
            if not os.path.isfile(fname):
                fname = _locateFile(fname, "param", teleap)
                if fname is None:
                    continue
            newname = os.path.join(outdir, f"param{i}_{os.path.basename(fname)}")
            shutil.copy(fname, newname)
            newparam.append(newname)

        newtopo = []
        for i, fname in enumerate(sorted(topo, key=lambda x: os.path.basename(x))):
            if not os.path.isfile(fname):
                fname = _locateFile(fname, "topo", teleap)
                if fname is None:
                    continue

            bn = os.path.basename(fname)
            newname = os.path.join(outdir, f"topo{i}_{bn}")
            shutil.copy(fname, newname)
            newtopo.append(newname)

        _fix_parameterize_atomtype_collisions(mol, newparam, newtopo)

        # Loading prepi topologies
        f.write("# Loading topologies\n")
        for fname in newtopo:
            if fname.lower().endswith(".mol2"):
                logger.warning(
                    "Please refrain from using MOL2 files as topologies as they are missing crucial information. Use `parameterize` .cif files or .prepi files instead."
                )
                resname = Molecule(fname).resname[0]
                f.write(f"{resname} = loadmol2 {os.path.basename(fname)}\n")
            elif fname.lower().endswith(".cif"):
                _resmol = Molecule(fname)
                _mol2f = os.path.join(fname.replace(".cif", ".mol2"))
                _resmol.write(_mol2f)
                f.write(_create_atomtype_section(_resmol))
                f.write(f"{_resmol.resname[0]} = loadmol2 {os.path.basename(_mol2f)}\n")
            else:
                f.write(f"loadamberprep {os.path.basename(fname)}\n")
        f.write("\n")

        # Loading frcmod parameters
        f.write("# Loading parameter files\n")
        for fname in newparam:
            f.write(f"loadamberparams {os.path.basename(fname)}\n")
        f.write("\n")

        cyclic = _detect_cyclic_segments(mol)
        cyclic_segs = [c[0] for c in cyclic]

        f.write("# Loading the system\n")
        nonc_mol = mol.copy()
        if len(cyclic):
            nonc_mol.remove(np.isin(mol.segid, cyclic_segs), _logger=False)
        if nonc_mol.numAtoms != 0:
            nonc_mol.write(os.path.join(outdir, "input.pdb"))
            f.write("mol = loadpdb input.pdb\n\n")

        if len(cyclic):
            f.write(
                "# Writing cyclic peptide segments. clearPdbResMap stops terminal patching\n"
            )
            f.write("clearPdbResMap\n")
            cyc_mol = mol.copy()
            cyc_mol.filter(np.isin(mol.segid, cyclic_segs), _logger=False)
            cyc_mol.write(os.path.join(outdir, "cyclic.pdb"))
            f.write("cyc = loadpdb cyclic.pdb\n")
            for seg, res_start, res_end in cyclic:
                f.write(f"bond cyc.{res_start}.N cyc.{res_end}.C\n")
            f.write("mol = combine {mol cyc}\n\n")

        # Write disulfide patches
        if len(disulfide) != 0:
            f.write("# Adding disulfide bonds\n")
            for d in disulfide:
                atoms1 = d[0].selectAtoms(mol, indexes=False)
                atoms2 = d[1].selectAtoms(mol, indexes=False)
                uqres1 = int(np.unique(mol.resid[atoms1]))
                uqres2 = int(np.unique(mol.resid[atoms2]))
                f.write(f"bond mol.{uqres1}.SG mol.{uqres2}.SG\n")
            f.write("\n")

        # Calculate the bounding box and store it in the CRD file
        f.write('setBox mol "vdw"\n\n')

        f.write("# Writing out the results\n")
        f.write(f"saveamberparm mol {prefix}.prmtop {prefix}.crd\n")
        f.write("quit")

    if execute:
        teleapimportflags = _getTeLeapImportFlags(teleap, ff, teleapimports)
        logpath = os.path.abspath(os.path.join(outdir, "log.txt"))
        logger.info("Starting the build.")
        currdir = os.getcwd()
        os.chdir(outdir)
        f = open(logpath, "w")
        try:
            cmd = [teleap, "-f", "./tleap.in"]
            cmd[1:1] = teleapimportflags
            logger.debug(cmd)
            call(cmd, stdout=f)
        except Exception:
            raise NameError("teLeap failed at execution")
        f.close()
        errors = _logParser(logpath)
        os.chdir(currdir)
        if errors:
            raise BuildError(
                errors
                + [f"Check {logpath} for further information on errors in building."],
                errors,
            )
        logger.info("Finished building.")

        if (
            os.path.exists(os.path.join(outdir, "structure.crd"))
            and os.path.getsize(os.path.join(outdir, "structure.crd")) != 0
            and os.path.getsize(os.path.join(outdir, "structure.prmtop")) != 0
        ):
            try:
                molbuilt = Molecule(
                    os.path.join(outdir, "structure.prmtop"), validateElements=False
                )
                molbuilt.read(os.path.join(outdir, "structure.crd"))
            except Exception as e:
                raise RuntimeError(
                    f"Failed at reading structure.prmtop/structure.crd due to error: {e}"
                )
        else:
            raise BuildError(
                f"No structure pdb/prmtop file was generated. Check {logpath} for errors in building."
            )

        molbuilt.write(os.path.join(outdir, "structure.pdb"), writebonds=False)
        detectCisPeptideBonds(molbuilt, respect_bonds=True)  # Warn in case of cis bonds
        return molbuilt


def _detect_cyclic_segments(mol: Molecule):
    cyclic = []
    prot = mol.atomselect("protein")
    for seg in np.unique(mol.segid):
        segatoms = mol.segid == seg
        if not np.all(prot[segatoms]):
            # Only for full-protein segments
            continue
        first_resid = mol.resid[segatoms][0]
        first_inser = mol.insertion[segatoms][0]
        last_resid = mol.resid[segatoms][-1]
        last_inser = mol.insertion[segatoms][-1]
        first_mask = (
            (mol.segid == seg)
            & (mol.resid == first_resid)
            & (mol.insertion == first_inser)
            & (mol.name == "N")
        )
        last_mask = (
            (mol.segid == seg)
            & (mol.resid == last_resid)
            & (mol.insertion == last_inser)
            & (mol.name == "C")
        )
        if first_mask.sum() > 1:
            raise RuntimeError(
                f"More than one atom named N in first residue of segid {seg}"
            )
        if last_mask.sum() > 1:
            raise RuntimeError(
                f"More than one atom named C in last residue of segid {seg}"
            )
        if first_mask.sum() == 0 or last_mask.sum() == 0:
            logger.warning(
                "Could not find N atom in N-terminal or C atom in C-terminal. Skipping cyclic peptide detection."
            )
        dist = np.linalg.norm(
            mol.coords[first_mask, :, 0] - mol.coords[last_mask, :, 0]
        )
        # Amide bond distance ranges between 1.325 - 1.346
        if dist < 1.35:
            cyclic.append((seg, first_resid, last_resid))
    return cyclic


def _add_caps(mol: Molecule, caps: dict):
    from moleculekit.molecule import UniqueResidueID

    capdir = os.path.join(home(shareDir=True), "builder", "caps")

    remove_atoms = {
        "ACE": ["H1", "H2", "H3", "HT1", "HT2", "HT3", "H"],
        "NME": ["OXT", "OT1", "O", "HXT"],
        "NHE": ["OXT", "OT1", "O", "HXT"],
    }

    uqres_caps = []
    if (
        caps is not None
        and len(caps)
        and isinstance(caps[list(caps.keys())[0]], (list, tuple))
    ):
        # Converting from segment caps to residue caps
        for segid in caps:
            # Keep first and last indexes of the segment
            segidx = np.where(mol.segid == segid)[0][[0, -1]]
            for i in range(2):
                cap = caps[segid][i]
                if cap is not None and cap != "none":
                    uqres_caps.append(
                        (UniqueResidueID.fromMolecule(mol, idx=segidx[i]), cap)
                    )
    else:
        # Converting from cap residue selections to unique residue IDs
        for sel, cap in caps.items():
            uqres_caps.append((UniqueResidueID.fromMolecule(mol, sel), cap))

    # For each caps definition
    for uqres, cap in uqres_caps:
        # Get info on segment and its terminals
        mask = uqres.selectAtoms(mol, indexes=False)

        if uqres.resname == cap:  # Cap already exists
            continue

        capmol = Molecule(os.path.join(capdir, f"{cap}.pdb"))
        align_names = capmol.name[capmol.beta == 1]
        mol_idx = np.where(np.isin(mol.name, align_names) & mask)[0]
        # Ensuring the atom order matches for the alignment
        cap_idx = [
            np.where((capmol.name == nn) & (capmol.resname == "XXX"))[0][0]
            for nn in mol.name[mol_idx]
        ]

        if len(cap_idx) != len(mol_idx):
            raise RuntimeError("Could not find all required backbone atoms!")

        capmol.align(cap_idx, refmol=mol, refsel=mol_idx)
        capmol.resname[capmol.resname == "XXX"] = uqres.resname
        capmol.segid[:] = uqres.segid
        capmol.chain[:] = uqres.chain
        capmol.resid[:] += uqres.resid
        capmol.remove(cap_idx, _logger=False)  # Remove atoms which were aligned

        # Remove the atoms defined in remove_atoms to clean up the termini
        mol.remove(np.isin(mol.name, remove_atoms[cap]) & mask, _logger=False)

        # Insert the cap into the mol
        res_idx = uqres.selectAtoms(mol, indexes=True)
        mol.insert(capmol, res_idx[0] if cap == "ACE" else res_idx[-1] + 1)

    # Remove terminal hydrogens regardless of caps or no caps. tleap confuses itself when it makes residues into terminals.
    # For example HID has an atom H which in NHID should become H[123] and crashes tleap.
    for seg in np.unique(mol.get("segid", "protein")):
        seg_mask = mol.segid == seg
        resids = mol.resid[seg_mask]
        mask = (
            ((mol.resid == resids[0]) | (mol.resid == resids[-1]))
            & seg_mask
            & (mol.element == "H")
        )
        mol.remove(mask, _logger=False)


def _defaultProteinCaps(mol):
    # Defines ACE and NME (neutral terminals) as default for protein segments
    # Of course, this might not be ideal for proteins that require charged terminals
    segsProt = np.unique(mol.get("segid", sel="protein"))
    caps = dict()
    for s in segsProt:
        if len(np.unique(mol.resid[mol.segid == s])) < 10:
            logger.warning(
                f"Segment {s} consists of a peptide with less than 10 residues. It will not be capped by "
                "default. If you want to cap it use the caps argument of amber.build to manually define caps"
                "for all segments."
            )
            continue
        caps[s] = ["ACE", "NME"]
    return caps


def _cleanOutDir(outdir):
    from glob import glob

    files = glob(os.path.join(outdir, "*"))
    for f in files:
        os.remove(f)


def _charmmLipid2Amber(mol):
    """Convert a CHARMM lipid membrane to AMBER format

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The Molecule object containing the membrane

    Returns
    -------
    newmol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        A new Molecule object with the membrane converted to AMBER
    """

    resdict = _readcsvdict(
        os.path.join(home(shareDir=True), "builder", "charmmlipid2amber.csv")
    )

    natoms = mol.numAtoms
    neworder = np.array(
        list(range(natoms))
    )  # After renaming the atoms and residues I have to reorder them

    begs = np.zeros(natoms, dtype=bool)
    fins = np.zeros(natoms, dtype=bool)
    begters = np.zeros(natoms, dtype=bool)
    finters = np.zeros(natoms, dtype=bool)

    # Iterate over the translation dictionary
    incrresids = sequenceID((mol.resid, mol.insertion, mol.segid))

    for res in resdict.keys():
        molresidx = mol.resname == res
        if not np.any(molresidx):
            continue
        names = (
            mol.name.copy()
        )  # Need to make a copy or I accidentally double-modify atoms

        atommap = resdict[res]
        for atom in atommap.keys():
            rule = atommap[atom]

            molatomidx = np.zeros(len(names), dtype=bool)
            molatomidx[molresidx] = names[molresidx] == atom

            mol.set("resname", rule.replaceresname, sel=molatomidx)
            mol.set("name", rule.replaceatom, sel=molatomidx)
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

    # if len(begters) > 999:
    #    raise NameError('More than 999 lipids. Cannot define separate segments for all of them.')

    for i in range(len(begters)):
        mapping = np.zeros(len(mol.resid), dtype=bool)
        mapping[begters[i] : finters[i] + 1] = True
        mol.set("resid", sequenceID(mol.get("resname", sel=mapping)), sel=mapping)
        mol.set("segid", f"L{i % 2}", sel=mapping)

    return mol


def _reorderMol(mol, order):
    for k in mol._atom_and_coord_fields:
        if mol.__dict__[k] is not None and np.size(mol.__dict__[k]) != 0:
            if k == "coords":
                mol.__dict__[k] = mol.__dict__[k][order, :, :]
            else:
                mol.__dict__[k] = mol.__dict__[k][order]


def _readcsvdict(filename):
    import csv
    from collections import namedtuple

    if os.path.isfile(filename):
        csvfile = open(filename, "r")
    else:
        raise NameError("File " + filename + " does not exist")

    resdict = dict()

    Rule = namedtuple(
        "Rule", ["replaceresname", "replaceatom", "order", "natoms", "ter"]
    )

    # Skip header line of csv file. Line 2 contains dictionary keys:
    csvfile.readline()
    csvreader = csv.DictReader(csvfile)
    for line in csvreader:
        searchres = line["search"].split()[1]
        searchatm = line["search"].split()[0]
        if searchres not in resdict:
            resdict[searchres] = dict()
        resdict[searchres][searchatm] = Rule(
            line["replace"].split()[1],
            line["replace"].split()[0],
            int(line["order"]),
            int(line["num_atom"]),
            line["TER"] == "True",
        )
    csvfile.close()

    return resdict


def _mapDisulfide(disulfide, mol, mol_orig):
    from moleculekit.molecule import UniqueResidueID

    disulfide_new = []
    for dd1, dd2 in disulfide:
        idx1 = dd1.selectAtoms(mol_orig)
        idx2 = dd2.selectAtoms(mol_orig)
        uqnew1 = UniqueResidueID.fromMolecule(mol, idx=idx1)
        uqnew2 = UniqueResidueID.fromMolecule(mol, idx=idx2)
        disulfide_new.append([uqnew1, uqnew2])
    return disulfide_new


def _logParser(fname):
    import re

    unknownres_regex = re.compile(r"Unknown residue:\s+(\w+)")
    missingparam_regex = re.compile(
        r"For atom: (.*) Could not find vdW \(or other\) parameters for type:\s+(\w+)"
    )
    missingtorsion_regex = re.compile(r"No torsion terms for\s+(.*)$")
    missingbond_regex = re.compile(r"Could not find bond parameter for:\s+(.*)$")
    missingangle_regex = re.compile(r"Could not find angle parameter:\s+(.*)$")
    missingatomtype_regex = re.compile(r"FATAL:\s+Atom (.*) does not have a type")

    unknownres = []
    missingparam = []
    missingtorsion = []
    missingbond = []
    missingangle = []
    missingatomtype = []
    with open(fname, "r") as f:
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
        errors.append(
            MissingResidueError(
                f"Unknown residue(s) {np.unique(unknownres)} found in the input structure. "
                "You are either missing a topology definition for the residue or you need to "
                "rename it to the correct residue name",
                np.unique(unknownres),
            )
        )
    if len(missingparam):
        errors.append(
            MissingParameterError(
                f"Missing parameters for atom {missingparam} and type {missingparam}",
                missingparam,
            )
        )
    if len(missingtorsion):
        errors.append(
            MissingTorsionError(
                f"Missing torsion terms for {missingtorsion}", missingtorsion
            )
        )
    if len(missingbond):
        errors.append(
            MissingBondError(f"Missing bond parameters for {missingbond}", missingbond)
        )
    if len(missingangle):
        errors.append(
            MissingAngleError(
                f"Missing angle parameters for {missingangle}", missingangle
            )
        )
    if len(missingatomtype):
        errors.append(
            MissingAtomTypeError(
                f"Missing atom type for {missingatomtype}", missingatomtype
            )
        )

    return errors


def _resname_from_fname(ff):
    return (
        os.path.splitext(os.path.basename(ff))[0]
        .split("_", maxsplit=1)[1]
        .replace("-orig", "")
    )


def _fix_parameterize_atomtype_collisions(mol, params, topos):
    def _gen_atomtypes():
        for x in parameterize_types:
            yield x

    gen = _gen_atomtypes()

    def _fix_frcmod(fname, repl):
        with open(fname, "r") as f:
            lines = f.readlines()

        const = ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB")
        with open(fname, "w") as f:
            f.write("Created by HTMD\n")
            for line in lines[1:]:
                if line.strip() == "" or any([line.startswith(x) for x in const]):
                    f.write(line)
                    continue
                # Split on two white spaces. Some atom types have a white space H - N -zs for example
                types, rest = line.split(sep="  ", maxsplit=1)
                types = types.split("-")
                for i in range(len(types)):
                    if types[i] in repl:
                        types[i] = repl[types[i]]
                f.write("-".join(types) + "  " + rest)

    def _fix_prepi(fname, repl):
        with open(fname, "r") as f:
            lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].strip().startswith("CORR"):
                break

        for i in range(i + 2, len(lines)):
            if lines[i].strip() == "":
                break
            old_at = lines[i][12:14]
            if old_at in repl:
                lines[i] = lines[i][:12] + repl[old_at] + lines[i][14:]

        with open(fname, "w") as f:
            for line in lines:
                f.write(line)

    def _fix_mol2(fname, repl):
        mm = Molecule(fname)
        for orig, new in repl.items():
            mm.atomtype[mm.atomtype == orig] = new
        mm.write(fname)

    def _fix_mol(bn, mol, repl):
        idx = np.where(mol.resname == bn)[0]
        for i in idx:
            if mol.atomtype[i] in repl:
                mol.atomtype[i] = repl[mol.atomtype[i]]

    # Match topos to frcmods
    topos_bn = {_resname_from_fname(x): x for x in topos}
    frcmd_bn = {_resname_from_fname(x): x for x in params}
    replacements = {}
    atomtypes = []
    for bn in frcmd_bn:
        replacements[bn] = {}

        # Find all parameterize atom types in the frcmod file
        with open(frcmd_bn[bn], "r") as f:
            file_at = []
            for line in f.readlines()[2:]:
                if line.strip() == "":
                    break
                at = line.split()[0]
                if at not in parameterize_types:
                    continue
                file_at.append(at)

        # Check for collisions with previous frcmod files and invent new type
        for at in file_at:
            if at in atomtypes:
                # Invent new atom type and rename in frcmod
                new_at = next(gen)
                while new_at in atomtypes + file_at:
                    new_at = next(gen)

                atomtypes.append(new_at)
                replacements[bn][at] = new_at
                logger.info(
                    f"Converting {bn} atom type {at} to {new_at} to avoid collisions."
                )
            else:
                atomtypes.append(at)

    # Correct all atom types to the new ones
    for bn in replacements:
        if len(replacements[bn]) == 0:
            continue
        _fix_frcmod(frcmd_bn[bn], replacements[bn])
        if bn in topos_bn:
            ext = os.path.splitext(topos_bn[bn])[1].lower()
            if ext in (".mol2", ".mol", ".cif"):
                _fix_mol2(topos_bn[bn], replacements[bn])
            elif ext in (".prepi", ".prep"):
                _fix_prepi(topos_bn[bn], replacements[bn])
            else:
                RuntimeError(
                    f"Don't know how to repair topologies with extension {ext}"
                )

        if bn in mol.resname:
            _fix_mol(bn, mol, replacements[bn])
        else:
            raise RuntimeError(
                f"Could not find residue {bn} in the input structure. Please ensure that the frcmod / prepi files are named as RES.frcmod / RES.prepi where RES the residue name to which they correspond."
            )


try:
    _findTeLeap()
    tleap_installed = True
except Exception:
    tleap_installed = False


@unittest.skipIf(
    not tleap_installed, "teLeap is not installed. Cannot test amber.build"
)
class _TestAmberBuild(unittest.TestCase):
    currentResult = None  # holds last result object passed to run method

    def setUp(self):
        from htmd.util import tempname

        self.testDir = os.environ.get("TESTDIR", tempname())
        print(f"Running tests in {self.testDir}")

    def run(self, result=None):
        self.currentResult = result  # remember result for use in tearDown
        unittest.TestCase.run(self, result)  # call superclass run method

    def tearDown(self):
        pass
        # if self.currentResult.wasSuccessful():
        #     shutil.rmtree(self.testDir)

    @staticmethod
    def _compareResultFolders(
        compare, tmpdir, pid, ignore_ftypes=(".log", ".txt", ".frcmod")
    ):
        import filecmp

        def _cutfirstline(infile, outfile):
            # Cut out the first line of prmtop which has a build date in it
            with open(infile, "r") as fin:
                data = fin.read().splitlines(True)
            with open(outfile, "w") as fout:
                fout.writelines(data[1:])

        try:
            from ffevaluation.ffevaluate import FFEvaluate, loadParameters

            prm2 = loadParameters(os.path.join(compare, "structure.prmtop"))
            mol2 = Molecule(os.path.join(compare, "structure.prmtop"))
            mol2.read(os.path.join(compare, "structure.pdb"))
            ffev2 = FFEvaluate(mol2, prm2)
            energies2, _, _ = ffev2.calculate(mol2.coords)

            prm = loadParameters(os.path.join(tmpdir, "structure.prmtop"))
            mol = Molecule(os.path.join(tmpdir, "structure.prmtop"))
            mol.read(os.path.join(tmpdir, "structure.pdb"))
            ffev = FFEvaluate(mol, prm)
            energies, _, _ = ffev.calculate(mol.coords)
            ene_diff = np.abs(energies - energies2)
            if any(ene_diff > 1e-3):
                print(f"ENERGY DIFF:\n{ene_diff}")
            else:
                print("SIMILAR ENERGIES")
        except Exception as e:
            print(f"Could not compare energies... {e}")

        files = []
        deletefiles = []
        for f in glob(join(compare, "*")):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith("prmtop"):
                _cutfirstline(f, join(compare, fname + ".mod"))
                _cutfirstline(join(tmpdir, fname), os.path.join(tmpdir, fname + ".mod"))
                files.append(os.path.basename(f) + ".mod")
                deletefiles.append(join(compare, fname + ".mod"))
            else:
                files.append(os.path.basename(f))

        match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError(
                f"Different results produced by amber.build for test {pid} between {compare} and {tmpdir} in files {mismatch}."
            )

        for f in deletefiles:
            os.remove(f)

    def test_with_protein_prepare(self):
        from moleculekit.tools.preparation import systemPrepare
        from htmd.builder.solvate import solvate
        from moleculekit.tools.autosegment import autoSegment

        # Test with systemPrepare
        pdbids = ["3PTB"]
        # pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)
        for pid in pdbids:
            np.random.seed(1)
            mol = Molecule(pid)
            mol.filter("protein")
            mol = systemPrepare(mol, pH=7.0)
            mol.filter("protein")  # Fix for bad systemPrepare hydrogen placing
            mol = autoSegment(mol)
            smol = solvate(mol)
            ffs = defaultFf()
            tmpdir = os.path.join(self.testDir, "withProtPrep", pid)
            _ = build(smol, ff=ffs, outdir=tmpdir)

            refdir = home(dataDir=join("test-amber-build", "pp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

    def test_without_protein_prepare(self):
        from htmd.builder.solvate import solvate

        # Test without systemPrepare
        pdbids = ["3PTB"]
        # pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
        for pid in pdbids:
            np.random.seed(1)
            mol = Molecule(pid)
            mol.filter("protein")
            smol = solvate(mol)
            ffs = defaultFf()
            tmpdir = os.path.join(self.testDir, "withoutProtPrep", pid)
            _ = build(smol, ff=ffs, outdir=tmpdir)

            refdir = home(dataDir=join("test-amber-build", "nopp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

    def test_protein_ligand(self):
        # Test protein ligand building with parametrized ligand
        refdir = home(dataDir=join("test-amber-build", "protLig"))
        tmpdir = os.path.join(self.testDir, "protLig")
        os.makedirs(tmpdir)

        mol = Molecule(join(refdir, "3ptb_mod.pdb"))
        lig = Molecule(join(refdir, "BEN.mol2"))
        lig.segid[:] = "L"

        newmol = Molecule()
        newmol.append(mol)
        newmol.append(lig)

        params = defaultParam() + [join(refdir, "BEN.frcmod")]
        topos = defaultTopo() + [join(refdir, "BEN.mol2")]

        _ = build(newmol, outdir=tmpdir, param=params, topo=topos, ionize=False)

        resdir = home(dataDir=join("test-amber-build", "protLig", "results"))
        _TestAmberBuild._compareResultFolders(resdir, tmpdir, "3PTB with mol2")

    def test_custom_disulfide_bonds(self):
        from htmd.builder.solvate import solvate

        pdbids = ["1GZM"]
        for pid in pdbids:
            np.random.seed(1)
            mol = Molecule(pid)
            mol.filter("protein")
            smol = solvate(mol)
            ffs = defaultFf()
            disu = [
                ["segid 0 and resid 110", "segid 0 and resid 187"],
                ["segid 1 and resid 110", "segid 1 and resid 187"],
            ]
            tmpdir = os.path.join(self.testDir, "withoutProtPrep", pid)
            _ = build(smol, ff=ffs, outdir=tmpdir, disulfide=disu)

            refdir = home(dataDir=join("test-amber-build", "nopp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

            np.random.seed(1)
            tmpdir = os.path.join(self.testDir, "withoutProtPrep", pid)
            _ = build(smol, ff=ffs, outdir=tmpdir)

            refdir = home(dataDir=join("test-amber-build", "nopp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

    def test_custom_teleap_imports(self):
        from htmd.builder.solvate import solvate

        pdbids = ["3PTB"]
        # pdbids = ['3PTB', '1A25', '1GZM', '1U5U']
        for pid in pdbids:
            np.random.seed(1)
            mol = Molecule(pid)
            mol.filter("protein")
            smol = solvate(mol)
            ffs = defaultFf()
            tmpdir = os.path.join(self.testDir, "withoutProtPrep", pid)

            amberhome = defaultAmberHome()
            teleapimports = [
                os.path.join(amberhome, _defaultAmberSearchPaths["ff"]),
                os.path.join(amberhome, _defaultAmberSearchPaths["lib"]),
                os.path.join(amberhome, _defaultAmberSearchPaths["param"]),
            ]

            _ = build(smol, ff=ffs, outdir=tmpdir, teleapimports=teleapimports)

            refdir = home(dataDir=join("test-amber-build", "nopp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

    def test_rna(self):
        from htmd.builder.solvate import solvate

        np.random.seed(1)

        mol = Molecule("6VA1")
        smol = solvate(mol)

        tmpdir = os.path.join(self.testDir, "rna", "6VA1")
        _ = build(smol, outdir=tmpdir)

        refdir = home(dataDir=join("test-amber-build", "rna", "6VA1"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "6VA1")

    def test_dna(self):
        from htmd.builder.solvate import solvate

        np.random.seed(1)

        mol = Molecule("1BNA")
        smol = solvate(mol)

        tmpdir = os.path.join(self.testDir, "dna", "1BNA")
        _ = build(smol, outdir=tmpdir)

        refdir = home(dataDir=join("test-amber-build", "dna", "1BNA"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "1BNA")

    def test_protein_rna(self):
        from htmd.builder.solvate import solvate
        from moleculekit.tools.preparation import systemPrepare
        from moleculekit.tools.autosegment import autoSegment

        np.random.seed(1)

        mol = Molecule("3WBM")
        mol.filter("not water")
        mol = autoSegment(mol, field="both")
        pmol = systemPrepare(mol, pH=7.0)
        smol = solvate(pmol)

        tmpdir = os.path.join(self.testDir, "protein-rna", "3WBM")
        _ = build(smol, outdir=tmpdir)

        refdir = home(dataDir=join("test-amber-build", "protein-rna", "3WBM"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "3WBM")

    def test_caps(self):
        from htmd.builder.solvate import solvate
        from moleculekit.tools.preparation import systemPrepare

        np.random.seed(1)

        mol = Molecule("6A5J")
        pmol = systemPrepare(mol, pH=7.0)
        smol = solvate(pmol)

        tmpdir = os.path.join(self.testDir, "peptide-cap", "6A5J")
        _ = build(smol, outdir=tmpdir, ionize=False)

        refdir = home(dataDir=join("test-amber-build", "peptide-cap", "6A5J"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "6A5J")

        np.random.seed(1)

        mol = Molecule("6A5J")
        pmol = systemPrepare(mol, pH=7.0)
        pmol.remove("(resid 1 13 and not backbone) or (resid 13 and name OXT)")
        smol = solvate(pmol)

        tmpdir = os.path.join(self.testDir, "peptide-cap-only-backbone", "6A5J")
        _ = build(smol, outdir=tmpdir, ionize=False)

        refdir = home(
            dataDir=join("test-amber-build", "peptide-cap-only-backbone", "6A5J")
        )

        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "6A5J")

        tmpdir = os.path.join(
            self.testDir, "peptide-cap-only-backbone-manual-caps", "6A5J"
        )
        _ = build(
            smol,
            outdir=tmpdir,
            ionize=False,
            caps={"resid 1": "ACE", "resid 13": "NME"},
        )

        refdir = home(
            dataDir=join("test-amber-build", "peptide-cap-only-backbone", "6A5J")
        )
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "6A5J")

    def test_non_standard_residue_building(self):
        homedir = home(dataDir=join("test-amber-build", "non-standard"))
        pdbids = ["5VBL", "1AWF"]

        for pdbid in pdbids:
            protdir = os.path.join(homedir, pdbid)
            with self.subTest(pdbid=pdbid):
                tmpdir = os.path.join(self.testDir, f"ns-res-{pdbid}")
                pmol = Molecule(os.path.join(protdir, f"{pdbid}_nolig.pdb"))
                _ = build(
                    pmol,
                    topo=glob(os.path.join(protdir, "*.prepi")),
                    param=glob(os.path.join(protdir, "*.frcmod")),
                    ionize=False,
                    outdir=tmpdir,
                )
                refdir = home(
                    dataDir=join("test-amber-build", "non-standard", pdbid, "build")
                )
                _TestAmberBuild._compareResultFolders(refdir, tmpdir, pdbid)

    @unittest.skipIf(sys.platform.startswith("darwin"), "Fails on OSX")
    def test_cofactor_building(self):
        homedir = home(dataDir=join("test-amber-build", "cofactors"))
        tmpdir = os.path.join(self.testDir, "cofactor")

        mol = Molecule(os.path.join(homedir, "cofactors.pdb"))
        _ = build(mol, ionize=False, outdir=tmpdir)
        refdir = home(dataDir=join("test-amber-build", "cofactors", "build"))

        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "Cofactors")

    def test_post_translational_modifications_building(self):
        homedir = home(dataDir=join("test-amber-build", "post-translational"))
        tmpdir = os.path.join(self.testDir, "posttransmod")

        mol = Molecule(os.path.join(homedir, "4EFP_nolig.pdb"))
        _ = build(mol, ionize=False, outdir=tmpdir)
        refdir = home(dataDir=join("test-amber-build", "post-translational", "build"))

        _TestAmberBuild._compareResultFolders(
            refdir, tmpdir, "post-translational modifications"
        )

    def test_cyclic_peptide_building(self):
        np.random.seed(1)

        mol = Molecule("5VAV")

        tmpdir = os.path.join(self.testDir, "cyclic-peptide", "5VAV")
        _ = build(mol, outdir=tmpdir, ionize=False, caps={"0": ("none", "none")})

        refdir = home(dataDir=join("test-amber-build", "cyclic-peptide", "5VAV"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "5VAV")

    def test_membrane(self):
        np.random.seed(1)

        homedir = home(dataDir=join("test-amber-build", "membrane"))

        mol = Molecule(os.path.join(homedir, "structure.psf"))
        mol.read(os.path.join(homedir, "structure.pdb"))

        tmpdir = os.path.join(self.testDir, "membrane")
        _ = build(mol, outdir=tmpdir)

        refdir = home(dataDir=join("test-amber-build", "membrane", "build"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "Membrane")

    def test_cif_building(self):
        # Tests that mol2 and cif building produce same results
        np.random.seed(1)

        homedir = home(dataDir=join("test-amber-build", "parameterize-cif"))
        ciff = os.path.join(homedir, "MOL-orig.cif")
        fmod = os.path.join(homedir, "MOL.frcmod")

        mol = Molecule("3ptb")
        mol.remove("resname BEN")
        mol2 = Molecule(ciff)
        mol2.segid[:] = "L"
        mol.append(mol2)

        tmpdir = os.path.join(self.testDir, "parameterize-cif", "3PTB")
        refdir = home(dataDir=join("test-amber-build", "parameterize-cif", "build"))

        build(mol, ionize=False, outdir=tmpdir, topo=[ciff], param=[fmod])
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "3PTB")

    def test_atomtype_decollisioning(self):
        # Tests that mol2 and cif building produce same results
        np.random.seed(1)

        homedir = home(dataDir=join("test-amber-build", "atomtype-decollisioning"))
        mol1f = os.path.join(homedir, "m42.cif")
        mol2f = os.path.join(homedir, "m48.cif")
        fmod1 = os.path.join(homedir, "m42.frcmod")
        fmod2 = os.path.join(homedir, "m48.frcmod")

        mol1 = Molecule(mol1f)
        mol1.resname[:] = "m42"
        mol1.segid[:] = "L1"
        mol2 = Molecule(mol2f)
        mol2.resname[:] = "m48"
        mol2.segid[:] = "L2"
        mol1.append(mol2)

        tmpdir = os.path.join(self.testDir, "atomtype-decollisioning", "decollisioning")
        refdir = home(
            dataDir=join("test-amber-build", "atomtype-decollisioning", "build")
        )

        build(
            mol1, ionize=False, outdir=tmpdir, topo=[mol1f, mol2f], param=[fmod1, fmod2]
        )
        _TestAmberBuild._compareResultFolders(
            refdir,
            tmpdir,
            "decollisioning",
            ignore_ftypes=(".log", ".txt", ".frcmod", ".in", ".mol2", ".cif"),
        )


if __name__ == "__main__":
    import doctest
    import unittest

    failure_count, _ = doctest.testmod()
    unittest.main(verbosity=2)
    # if failure_count != 0:
    #     raise Exception('Doctests failed')
