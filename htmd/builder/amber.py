# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
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
import subprocess
from subprocess import call
from moleculekit.molecule import Molecule
from moleculekit.tools.sequencestructuralalignment import sequenceStructureAlignment
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.util import ensurelist
import logging

logger = logging.getLogger(__name__)


# All atom types reserved by GAFF2. We will use all other low caps for parameterize
gaff2_at = [
    "c",
    "cs",
    "c1",
    "c2",
    "c3",
    "ca",
    "cp",
    "cq",
    "cc",
    "cd",
    "ce",
    "cf",
    "cg",
    "ch",
    "cx",
    "cy",
    "c5",
    "c6",
    "cu",
    "cv",
    "cz",
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
    "f",
    "cl",
    "br",
    "i",
    "n",
    "n1",
    "n2",
    "n3",
    "n4",
    "na",
    "nb",
    "nc",
    "nd",
    "ne",
    "nf",
    "nh",
    "no",
    "ns",
    "nt",
    "nx",
    "ny",
    "nz",
    "n+",
    "nu",
    "nv",
    "n7",
    "n8",
    "n9",
    "ni",
    "nj",
    "nk",
    "nl",
    "nm",
    "nn",
    "np",
    "nq",
    "n5",
    "n6",
    "o",
    "oh",
    "op",
    "oq",
    "os",
    "ow",
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
    "s",
    "s2",
    "s4",
    "s6",
    "sh",
    "ss",
    "sx",
    "sy",
    "sp",
    "sq",
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
        return None
    if os.path.islink(teleap):
        if os.path.isabs(os.readlink(teleap)):
            teleap = os.readlink(teleap)
        else:
            teleap = os.path.join(os.path.dirname(teleap), os.readlink(teleap))
    return teleap


def _pyodideAmberHome():
    """Returns the AMBERHOME-equivalent path from the tleap_pyodide package.

    The tleap_pyodide package bundles dat/leap/{cmd,prep,lib,parm} under its
    package directory, mirroring the $AMBERHOME layout.
    """
    try:
        import tleap_pyodide

        return tleap_pyodide._PACKAGE_DIR
    except ImportError:
        return None


def _resolve_backend(teleap=None):
    """Determine whether to use native teLeap or tleap_pyodide.

    Returns
    -------
    tuple of (str, str)
        ("native", teleap_path) or ("pyodide", amberhome_path).

    Raises
    ------
    FileNotFoundError
        If neither native teLeap nor tleap_pyodide is available.
    """
    if teleap is not None:
        resolved = shutil.which(teleap)
        if resolved is None:
            raise FileNotFoundError(
                f"Could not find executable: `{teleap}` in the PATH."
            )
        return ("native", resolved)

    native = _findTeLeap()
    if native is not None:
        return ("native", native)

    pyodide_home = _pyodideAmberHome()
    if pyodide_home is not None:
        return ("pyodide", pyodide_home)

    raise FileNotFoundError(
        "Neither native teLeap nor tleap_pyodide found. "
        "Install AmberTools (conda install ambertools -c conda-forge) "
        "or install tleap_pyodide for browser/pyodide environments."
    )


def defaultAmberHome(teleap: str | None = None) -> str:
    """Return the default AMBERHOME directory.

    Determined from the location of the teLeap binary, or the tleap_pyodide
    package directory when running in a pyodide environment.

    Parameters
    ----------
    teleap : str, optional
        Path to the teLeap executable. If None, the executable is located
        automatically from PATH.

    Returns
    -------
    amberhome : str
        The AMBERHOME directory path.
    """
    backend, value = _resolve_backend(teleap)
    if backend == "pyodide":
        return value

    return os.path.normpath(os.path.join(os.path.dirname(value), "../"))


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


# Modified amino acids that AMBER parameterizes only in mod_amino.lib (the
# ff*SB_modAA forcefield), NOT the base ff14SB/ff19SB libraries. Mirrors the
# unit list + PDB aliases in leaprc.protein.ff14SB_modAA.
_MODAA_RESNAMES = frozenset(
    {"ALY", "AZF", "CYF", "CNX", "MSE", "4II", "4CF", "MTN"}
)


def _detect_modaa_residues(mol, ff):
    """Auto-load AMBER's modified-amino-acid forcefield when such a residue is
    present. MSE, ALY, AZF, CYF, CNX live in ``mod_amino.lib`` (loaded by
    ``leaprc.protein.ff*SB_modAA``), not the base ff14SB/ff19SB libraries, so
    without this tleap reports "Unknown residue". The modAA leaprc is purely
    additive (atom types + frcmod + lib + PDB name map), designed to load
    alongside the base protein ff; the variant is matched to the loaded protein
    ff (ff19SB vs ff14SB). Appends the leaprc to ``ff`` in place and returns the
    detected resnames.
    """
    present = sorted(
        {str(r) for r in np.unique(mol.resname)} & _MODAA_RESNAMES
    )
    if not present:
        return []
    variant = "ff19SB" if any("ff19SB" in str(f) for f in ff) else "ff14SB"
    leaprc = f"leaprc.protein.{variant}_modAA"
    if leaprc not in ff:
        ff.append(leaprc)
        logger.info(
            f"Modified amino acid(s) {', '.join(present)} detected in system. "
            f"Automatically loading {leaprc} for AMBER."
        )
    return present


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
                uqresid = uqresid[~h_sel]
                mol.remove(h_sel, _logger=False)

            if cof[2] == "cofactors":
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


def _locateFile(fname, ftype, amberhome):
    htmdamberdir = htmdAmberHome()
    searchdir = os.path.join(amberhome, _defaultAmberSearchPaths[ftype])
    foundfile = glob(os.path.join(searchdir, fname))
    if len(foundfile) != 0:
        return foundfile[0]
    foundfile = glob(os.path.join(htmdamberdir, fname))
    if len(foundfile) != 0:
        return foundfile[0]
    logger.warning(f"Was not able to find {ftype} file {fname}")


def _getTeLeapImportFlags(amberhome=None, ff=None, teleapimports=()):
    if not teleapimports:
        if amberhome is None:
            amberhome = defaultAmberHome()
        teleapimports = []
        teleapimports += [
            os.path.join(amberhome, s) for s in _defaultAmberSearchPaths.values()
        ]
        if len(teleapimports) == 0:
            raise RuntimeError(
                "No default Amber force-field found. Check AMBERHOME location."
            )
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

    teleapimportflags = []
    for p in teleapimports:
        teleapimportflags.append("-I")
        teleapimportflags.append(str(p))
    return teleapimportflags


def _getHtmdExtraImports(ff=None):
    """Return extra import directories from htmdAmberHome for tleap_pyodide."""
    htmdamberdir = htmdAmberHome()
    extra = []
    if ff is not None:
        for f in ff:
            dirpath = os.path.join(htmdamberdir, os.path.dirname(f))
            if os.path.isdir(dirpath) and dirpath not in extra:
                extra.append(dirpath)
    return extra


def _getResiduesInAllFFs():
    amberhome = defaultAmberHome()
    ffdir = join(amberhome, _defaultAmberSearchPaths["ff"])
    ffs = [
        os.path.relpath(ff, ffdir)
        for ff in glob(join(ffdir, "*"))
        if os.path.isfile(ff)
    ]

    htmdamberdir = htmdAmberHome()
    ffs += [
        os.path.relpath(ff, htmdamberdir)
        for ff in glob(join(htmdamberdir, "*", "leaprc.*"))
        if os.path.isfile(ff)
    ]

    residues = {}
    for ff in ffs:
        residues[ff] = _getResiduesInFF(ff)

    return residues


def _parseResidueList(log_text):
    """Parse the output of tleap 'list' command to extract residue names."""
    residues = []
    starting = False
    for line in log_text.splitlines():
        stripped = line.strip()
        if stripped == "> list":
            starting = True
            continue
        if stripped == "quit":
            break
        if starting:
            residues += stripped.split()
    return residues


def _getResiduesInFF(ff, teleap=None):
    import tempfile

    backend, value = _resolve_backend(teleap)
    if backend == "native":
        amberhome = os.path.normpath(
            os.path.join(os.path.dirname(value), "../")
        )
    else:
        amberhome = value

    with tempfile.TemporaryDirectory() as tmpdir:
        ff = _locateFile(ff, "ff", amberhome)
        if ff is None:
            raise FileNotFoundError(
                "Could not find forcefield file. Check that AmberTools "
                "or tleap_pyodide is installed correctly."
            )
        shutil.copy(ff, tmpdir)

        script = f"logFile leap.log\nsource {os.path.basename(ff)}\nlist\nquit"

        if backend == "pyodide":
            from tleap_pyodide import run_tleap as _pyodide_run

            with open(os.path.join(tmpdir, "tleap.in"), "w") as f:
                f.write(script)
            result = _pyodide_run(work_dir=tmpdir, script_file="tleap.in")
            outlog = os.path.join(tmpdir, "leap.log")
            if not os.path.exists(outlog):
                if "leap.log" in result:
                    with open(outlog, "wb") as f:
                        f.write(result["leap.log"])
                else:
                    raise BuildError("tleap_pyodide did not produce leap.log")
            with open(outlog, "r") as f:
                log_text = f.read()
        else:
            importflgs = _getTeLeapImportFlags(amberhome, [ff])
            with open(os.path.join(tmpdir, "tleap.in"), "w") as f:
                f.write(script)
            logpath = os.path.join(tmpdir, "log.txt")
            with open(logpath, "w") as f:
                try:
                    cmd = [value, "-f", "./tleap.in"]
                    cmd[1:1] = importflgs
                    logger.debug(cmd)
                    call(cmd, stdout=f, cwd=tmpdir)
                except Exception:
                    raise BuildError("teLeap failed at execution")
            outlog = os.path.join(tmpdir, "leap.log")
            with open(outlog, "r") as f:
                log_text = f.read()

        return _parseResidueList(log_text)


def defaultFf():
    """Returns the default leaprc forcefield files used by amber.build"""
    return [
        "leaprc.protein.ff14SB",
        "leaprc.lipid21",
        "leaprc.gaff2",
        "leaprc.RNA.OL3",
        "leaprc.DNA.bsc1",
        "leaprc.water.tip3p",
    ]


def defaultTopo():
    """Returns the default topology `prepi` files used by amber.build"""
    return []


def defaultParam():
    """Returns the default parameter `frcmod` files used by amber.build"""
    return []


def _prepareMolecule(mol: Molecule, caps, disulfide, custombonds, remove):
    from moleculekit.molecule import UniqueAtomID, UniqueResidueID

    # Check for missing segids or mixed protein / non-protein segments
    _missingSegID(mol)
    _checkMixedSegment(mol)

    # Detect cyclic segments while the input bonds are still present, so an
    # explicit head-to-tail closure bond is honored: a closure modeled longer
    # than the 1.35 A distance heuristic (7BTI's phalloidin at ~1.47 A,
    # microcystin's at ~1.37 A) is still recognised. Must run before deleteBonds
    # and before _add_caps (which inserts atoms and would invalidate the
    # index-based bond list).
    cyclic_segids = [c[0] for c in _detect_cyclic_segments(mol)]

    # Remove pdb bonds as they can be regenerated by tleap
    mol.deleteBonds("all")

    # Convert lipids to AMBER naming
    mol = _charmmLipid2Amber(mol)

    # PDB stores the C-terminal amide cap as "NH2" but AMBER ff14SB calls it
    # "NHE"; the residue is otherwise the same. Rename so tLeap recognises it.
    if np.any(mol.resname == "NH2"):
        mol.resname[mol.resname == "NH2"] = "NHE"

    # Add caps to termini
    if caps is None:
        caps = _defaultProteinCaps(mol)

    for cc in cyclic_segids:
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

    # Map the old remove and bond selections to new ones
    if custombonds is not None:
        newcustombonds = []
        for bb in custombonds:
            if not isinstance(bb[0], str) or not isinstance(bb[1], str):
                raise RuntimeError("All custombonds selections should be strings")
            a1 = UniqueAtomID.fromMolecule(mol_orig_resid, bb[0])
            a2 = UniqueAtomID.fromMolecule(mol_orig_resid, bb[1])
            idx1 = a1.selectAtom(mol_orig_resid)
            idx2 = a2.selectAtom(mol_orig_resid)
            a1 = UniqueAtomID.fromMolecule(mol, idx=idx1)
            a2 = UniqueAtomID.fromMolecule(mol, idx=idx2)
            newcustombonds.append([a1, a2])
        custombonds = newcustombonds

    if remove is not None:
        newremove = []
        for rr in remove:
            if not isinstance(rr, str):
                raise RuntimeError("All remove selections should be strings")
            natoms = mol_orig_resid.atomselect(rr).sum()
            if natoms == 0:
                raise RuntimeError(f"No atoms selected for removal in selection {rr}")
            if natoms == 1:
                a = UniqueAtomID.fromMolecule(mol_orig_resid, rr)
                idx = a.selectAtom(mol_orig_resid)
                a = UniqueAtomID.fromMolecule(mol, idx=idx)
                newremove.append(a)
            else:
                r = UniqueResidueID.fromMolecule(mol_orig_resid, rr)
                idx = r.selectAtoms(mol_orig_resid)
                r = UniqueResidueID.fromMolecule(mol, idx=idx)
                newremove.append(r)
        remove = newremove

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

    # Apply anchor renames + H-drop for atoms that participate in custombonds
    # (defense-in-depth for scaffolded peptides; idempotent with
    # systemPrepare's force_protonation path). For each endpoint of a
    # custombond, look up its (resname, atom_name) in ANCHOR_VARIANTS; if a
    # variant residue is defined, rename the residue and drop the displaced
    # hydrogens. No-op for endpoints that aren't canonical-residue anchors
    # (e.g. the scaffold side of a scaffolded-peptide bond).
    if custombonds is not None and len(custombonds):
        from moleculekit.tools._anchor_variants import lookup_anchor

        torem = np.zeros(mol.numAtoms, dtype=bool)
        for a1, a2 in custombonds:
            for atom_id in (a1, a2):
                idx = atom_id.selectAtom(mol)
                resname = str(mol.resname[idx])
                atom_name = str(mol.name[idx])
                entry = lookup_anchor(resname, atom_name)
                if entry is None:
                    continue
                # Mask of all atoms in this residue.
                res_mask = (
                    (mol.resid == mol.resid[idx])
                    & (mol.segid == mol.segid[idx])
                    & (mol.insertion == mol.insertion[idx])
                    & (mol.chain == mol.chain[idx])
                )
                if entry["ff_variant"] is not None and resname != entry["ff_variant"]:
                    mol.resname[res_mask] = entry["ff_variant"]
                    atom_id.resname = entry["ff_variant"]
                for h_name in entry["drop_h"]:
                    torem |= res_mask & (mol.name == h_name)
        if np.any(torem):
            mol.remove(torem, _logger=False)

    return disulfide, custombonds, remove, cyclic_segids


def _write_residue_mapping(molbuilt, mol_orig, outdir):
    # Align with the original molecule to map built (renumbered) resids back to
    # the originals. sequenceStructureAlignment aligns one residue type at a
    # time, so handle protein and nucleic separately and combine - a system with
    # both (e.g. a protein plus a nucleic-like ligand such as ADP, 7BTI) would
    # otherwise raise "both protein and nucleic residues".
    residmap = []
    for seltype in ("protein", "nucleic"):
        if not np.any(molbuilt.atomselect(seltype)) or not np.any(
            mol_orig.atomselect(seltype)
        ):
            continue
        try:
            _, mapping = sequenceStructureAlignment(
                molbuilt,
                mol_orig,
                molsel=seltype,
                refsel=seltype,
                maxalignments=1,
            )
        except Exception as e:
            logger.error(f"Error aligning {seltype} for residue mapping: {e}")
            continue
        if len(mapping):
            mol_map, ref_map = mapping[0]  # Top alignment
            idx_mol = np.where(mol_map)[0]
            idx_ref = np.where(ref_map)[0]
            for im, ir in zip(idx_mol, idx_ref):
                residmap.append(
                    [
                        str(molbuilt.resid[im]),
                        str(mol_orig.resid[ir]),
                        mol_orig.insertion[ir],
                        mol_orig.chain[ir],
                    ]
                )
    if residmap:
        with open(os.path.join(outdir, "residue_mapping.csv"), "w") as fcsv:
            fcsv.write("new_resid,old_resid,old_insertion,old_chain\n")
            for mm in residmap:
                fcsv.write(",".join(mm) + "\n")


def build(
    mol: Molecule,
    ff: list | None = None,
    topo: list | None = None,
    param: list | None = None,
    prefix: str = "structure",
    outdir: str = "./build",
    caps: dict | None = None,
    ionize: bool = True,
    saltconc: float = 0,
    saltanion: str | None = None,
    saltcation: str | None = None,
    disulfide: list | None = None,
    teleap: str | None = None,
    teleapimports: list | None = None,
    execute: bool = True,
    atomtypes: list | None = None,
    offlibraries: list | str | None = None,
    gbsa: bool = False,
    igb: int = 2,
    custombonds: list | None = None,
    remove: list | None = None,
) -> Molecule:
    """Build a system for AMBER.

    Uses tleap to build a system for AMBER. Additionally allows for ionization
    and adding disulfide bridges.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The Molecule object containing the system.
    ff : list, optional
        A list of leaprc forcefield files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available forcefield files.
        If None, uses :func:`amber.defaultFf <htmd.builder.amber.defaultFf>`.
    topo : list, optional
        A list of topology ``prepi/prep/in/cif`` files.
        CIF and MOL2 files are automatically converted to mol2 format.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available topology files.
        If None, uses :func:`amber.defaultTopo <htmd.builder.amber.defaultTopo>`.
        When passing residues parameterized with the ``parameterize`` tool, pass the .cif file.
    param : list, optional
        A list of parameter ``frcmod`` files.
        Use :func:`amber.listFiles <htmd.builder.amber.listFiles>` to get a list of available parameter files.
        If None, uses :func:`amber.defaultParam <htmd.builder.amber.defaultParam>`.
    prefix : str
        The prefix for the generated pdb and prmtop files.
    outdir : str
        The path to the output directory.
    caps : dict, optional
        A dictionary specifying the caps. Accepts two formats.
        First format: keys are segids, values are lists of cap strings for that segment.
        e.g. ``caps['P'] = ['ACE', 'NME']`` or ``caps['P'] = ['none', 'none']``.
        If None, ACE and NME caps are applied to every protein segment.
        Second format: keys are atom selection strings for a residue, values are the cap name.
        e.g. ``caps = {"chain A and resid 5": "ACE", "chain A and resid 10": "NME"}``.
    ionize : bool
        Enable or disable ionization.
    saltconc : float
        Salt concentration (in Molar) to add to the system after neutralization.
    saltanion : str, optional
        The anion type. Available: ``'Cl-'``. Also accepts ``'CL'``, ``'chloride'``, ``'CLA'``.
    saltcation : str, optional
        The cation type. Available: ``'Na+'``, ``'K+'``, ``'Cs+'``, ``'Mg2+'``, ``'Ca2+'``, ``'Zn2+'``.
        Also accepts formats such as ``'NA'``, ``'sodium'``, ``'SOD'``.
    disulfide : list, optional
        If None, disulfide bonds are guessed automatically. Otherwise provide a list of pairs
        of atom selection strings for each pair of residues forming a disulfide bridge.
    teleap : str, optional
        Path to the teLeap executable. If None, located automatically from PATH.
    teleapimports : list, optional
        A list of directory paths to pass to teLeap via the ``-I`` flag.
        If None, determined from :func:`amber.defaultAmberHome <htmd.builder.amber.defaultAmberHome>`
        and :func:`amber.htmdAmberHome <htmd.builder.amber.htmdAmberHome>`.
    execute : bool
        If True, run the full build. If False, only write the tleap input script without building.
        Ionization is skipped when False.
    atomtypes : list, optional
        Custom atom types as a list of ``('type', 'element', 'hybrid')`` triplets,
        e.g. ``[('C1', 'C', 'sp2'), ('CI', 'C', 'sp3')]``. See ``addAtomTypes`` in AmberTools docs.
    offlibraries : list or str, optional
        A path or list of paths to OFF library files. See ``loadOFF`` in AmberTools docs.
    gbsa : bool
        Modify radii for the GBSA implicit water model.
    igb : int
        GB model: 1 for mbondi, 2 and 5 for mbondi2, 7 for bondi, 8 for mbondi3.
        See section 4 of the AMBER manual on the Generalized Born/Surface Area model.
    custombonds : list, optional
        A list of pairs of atom selection strings specifying custom bonds to add.
    remove : list, optional
        A list of atom selection strings. Matching atoms or residues are removed from the
        system. Useful for removing atoms added during building.

    Returns
    -------
    molbuilt : :class:`Molecule <moleculekit.molecule.Molecule>`
        The built system as a Molecule object.

    Examples
    --------
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

    disulfide, custombonds, remove, cyclic_segids = _prepareMolecule(
        mol, caps, disulfide, custombonds, remove
    )
    _detect_cofactors_ncaa_ptm(mol, param, topo)
    _detect_modaa_residues(mol, ff)

    backend, backend_value = _prepare_build(
        mol,
        ff=ff,
        topo=topo,
        param=param,
        prefix=prefix,
        outdir=outdir,
        disulfide=disulfide,
        teleap=teleap,
        atomtypes=atomtypes,
        offlibraries=offlibraries,
        gbsa=gbsa,
        igb=igb,
        custombonds=custombonds,
        remove=remove,
        cyclic_segids=cyclic_segids,
    )

    if ionize and execute and np.any(mol.atomselect("water")):
        molbuilt_solute = _run_tleap(
            outdir,
            "solute_charge",
            backend=backend,
            backend_value=backend_value,
            ff=ff,
            teleapimports=teleapimports,
            script="tleap_solute.in",
        )
        totalcharge = np.sum(molbuilt_solute.charge)
        nwater = np.sum(mol.atomselect("water and noh"))

        for ext in (".crd", ".pdb", ".prmtop"):
            fpath = os.path.join(outdir, f"solute_charge{ext}")
            if os.path.exists(fpath):
                os.remove(fpath)

        anion, cation, anionatom, cationatom, nanion, ncation = ionizef(
            molbuilt_solute,
            totalcharge,
            nwater,
            saltconc=saltconc,
            anion=saltanion,
            cation=saltcation,
        )

        water_sel = mol.atomselect("water")
        solvent_mol = mol.copy(sel=water_sel)
        solute_mol = mol.copy(sel=~water_sel)
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
        solvent_mol.write(os.path.join(outdir, "solvent.pdb"))

    if execute:
        molbuilt = _run_tleap(
            outdir,
            prefix,
            backend=backend,
            backend_value=backend_value,
            ff=ff,
            teleapimports=teleapimports,
        )
    else:
        molbuilt = None

    if molbuilt is not None:
        _write_residue_mapping(molbuilt, mol_orig, outdir)
    return molbuilt


def _create_atomtype_section(mol):
    atypes_str = "addAtomTypes {\n"
    for at in sorted(np.unique(mol.atomtype)):
        el = mol.element[mol.atomtype == at][0]
        atypes_str += f'	{{ "{at}"  "{el}" "sp3" }}\n'
    atypes_str += "}\n"
    return atypes_str


def _build_load_mol2_commands(mol, mol2_path, unit_name="_lig"):
    """Generate tleap commands to load a ligand from mol2 correctly.

    Produces:
      1. `addAtomTypes` block mapping every unique atom type to its element
         (needed so parameter lookups resolve the correct element).
      2. `<unit> = loadmol2 <mol2_path>`.
      3. `set <unit>.1.<name> element "<Element>"` for every atom.

    The explicit `set element` commands are required because tleap's
    `loadmol2` reader derives the atom element from only the first
    character of the atom type string (see `tripos.c` line 209-211 in
    AmberTools).  For standard GAFF2 types (c3, ca, h1, oh, …) this
    happens to work, but for custom types from parameterize (za, zb …)
    or OpenFF outputs it would misclassify hydrogens and produce an
    empty BONDS_INC_HYDROGEN list, breaking SHAKE.  Writing
    `set .element` after `loadmol2` unconditionally corrects this.
    """
    # addAtomTypes block
    lines = ["addAtomTypes {"]
    for at in sorted(np.unique(mol.atomtype)):
        el = mol.element[mol.atomtype == at][0]
        # pick a sensible hybridization default - tleap only uses this
        # for coordinate-building, which we never invoke for a ligand
        # whose geometry is already given.
        lines.append(f'    {{ "{at}"  "{el}" "sp3" }}')
    lines.append("}")

    # load mol2
    lines.append(f"{unit_name} = loadmol2 {os.path.basename(mol2_path)}")

    # per-atom element override
    for i in range(mol.numAtoms):
        name = mol.name[i]
        el = mol.element[i]
        if not el:
            raise ValueError(
                f"Atom {i} ({name}) has no element; cannot build mol2 load script."
            )
        lines.append(f'set {unit_name}.1.{name} element "{el}"')

    lines.append("")
    return "\n".join(lines)


def _tleap_residue_positions(mol, cyc_info, include_water):
    """Return a per-atom 1-based residue position in the final combined
    tLeap unit, matching the load-and-combine order:

        mol = loadpdb input.pdb       # solute residues 1..N_solute
        wat = loadpdb solvent.pdb     # waters appended via combine -> N_solute+1..
        cyc_X = loadpdb cyclic_X.pdb  # each cyclic seg appended in turn
        mol = combine {mol wat cyc_X ...}

    The combined-unit position - NOT the PDB ``resid`` field - is what
    tLeap's ``unit.<N>.<atom>`` selector resolves against. Using
    ``mol.resid`` directly only works when no residues are split out of
    input.pdb between mol's load and the final combine; the moment
    waters or cyclic segments are written to separate files, the
    htmd-internal resid (which counts those split-out residues) drifts
    past the combined-unit's sequential index.
    """
    n = mol.numAtoms
    pos = np.full(n, -1, dtype=int)

    cyclic_segs = [ci[0].replace("cyc_", "", 1) for ci in cyc_info] if cyc_info else []
    cyclic_mask = (
        np.isin(mol.segid, list(cyclic_segs)) if cyclic_segs else np.zeros(n, bool)
    )
    water_mask = mol.atomselect("water") & ~cyclic_mask
    solute_mask = ~water_mask & ~cyclic_mask

    offset = 0
    if solute_mask.any():
        sub = sequenceID(
            (mol.resid[solute_mask], mol.insertion[solute_mask], mol.segid[solute_mask])
        )
        pos[solute_mask] = sub + 1
        offset = int(sub.max()) + 1

    if water_mask.any() and include_water:
        sub = sequenceID(
            (mol.resid[water_mask], mol.insertion[water_mask], mol.segid[water_mask])
        )
        pos[water_mask] = sub + 1 + offset
        offset += int(sub.max()) + 1

    for seg in cyclic_segs:
        seg_mask = mol.segid == seg
        if not seg_mask.any():
            continue
        sub = sequenceID(
            (mol.resid[seg_mask], mol.insertion[seg_mask], mol.segid[seg_mask])
        )
        pos[seg_mask] = sub + 1 + offset
        offset += int(sub.max()) + 1

    return pos


def _write_tleap_script(
    filepath,
    ff_sources,
    gbsa,
    igb,
    atomtypes,
    offlib_names,
    topo_cmds,
    param_names,
    has_solute,
    has_water,
    include_solvent,
    cyc_info,
    mol,
    remove,
    disulfide,
    custombonds,
    prefix,
):
    """Write a tleap input script.

    When include_solvent is False, the solvent.pdb loading is omitted,
    producing a solute-only build (useful for fast charge calculation).
    """
    res_pos = _tleap_residue_positions(
        mol, cyc_info, include_water=(has_water and include_solvent)
    )
    with open(filepath, "w") as f:
        f.write("# tleap file generated by amber.build\n")

        for name in ff_sources:
            f.write(f"source {name}\n")
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

        for name in offlib_names:
            f.write(f"loadoff {name}\n")

        f.write("# Loading topologies\n")
        for cmd in topo_cmds:
            f.write(cmd)
        f.write("\n")

        f.write("# Loading parameter files\n")
        for name in param_names:
            f.write(f"loadamberparams {name}\n")
        f.write("\n")

        f.write("# Loading the system\n")
        if has_solute:
            f.write("mol = loadpdb input.pdb\n")

        if include_solvent and has_water:
            f.write("wat = loadpdb solvent.pdb\n")
            if has_solute:
                f.write("mol = combine {mol wat}\n")
            else:
                f.write("mol = wat\n")

        f.write("\n")

        if len(cyc_info):
            f.write(
                "# Writing cyclic peptide segments. clearPdbResMap stops terminal patching\n"
            )
            f.write("clearPdbResMap\n")
            cyc_vars = []
            for cyc_var, fname, res_start, res_end in cyc_info:
                cyc_vars.append(cyc_var)
                f.write(f"{cyc_var} = loadpdb {fname}\n")
                f.write(f"bond {cyc_var}.{res_start}.N {cyc_var}.{res_end}.C\n")
            cyc_list = " ".join(cyc_vars)
            has_non_cyclic = has_solute or (include_solvent and has_water)
            if has_non_cyclic:
                f.write(f"mol = combine {{mol {cyc_list}}}\n\n")
            elif len(cyc_vars) == 1:
                f.write(f"mol = {cyc_vars[0]}\n\n")
            else:
                f.write(f"mol = combine {{{cyc_list}}}\n\n")

        if remove is not None and len(remove) != 0:
            from moleculekit.molecule import UniqueAtomID, UniqueResidueID

            f.write("# Removing custom atoms/residues\n")
            for d in remove:
                if isinstance(d, UniqueAtomID):
                    atom = d.selectAtom(mol, indexes=False)
                    uqres = int(res_pos[atom])
                    name = mol.name[atom][0]
                    f.write(f"remove mol.{uqres} mol.{uqres}.{name}\n")
                elif isinstance(d, UniqueResidueID):
                    atoms = d.selectAtoms(mol, indexes=False)
                    uqres = int(res_pos[atoms][0])
                    f.write(f"remove mol mol.{uqres}\n")
            f.write("\n")

        if disulfide is not None and len(disulfide) != 0:
            f.write("# Adding disulfide bonds\n")
            for d in disulfide:
                atoms1 = d[0].selectAtoms(mol, indexes=False)
                atoms2 = d[1].selectAtoms(mol, indexes=False)
                uqres1 = int(np.unique(res_pos[atoms1])[0])
                uqres2 = int(np.unique(res_pos[atoms2])[0])
                f.write(f"bond mol.{uqres1}.SG mol.{uqres2}.SG\n")
            f.write("\n")

        if custombonds is not None and len(custombonds) != 0:
            f.write("# Adding custom bonds\n")
            for d in custombonds:
                atom1 = d[0].selectAtom(mol, indexes=False)
                atom2 = d[1].selectAtom(mol, indexes=False)
                uqres1 = int(res_pos[atom1][0])
                uqres2 = int(res_pos[atom2][0])
                name1 = mol.name[atom1][0]
                name2 = mol.name[atom2][0]
                f.write(f"bond mol.{uqres1}.{name1} mol.{uqres2}.{name2}\n")
            f.write("\n")

        f.write('setBox mol "vdw"\n\n')

        f.write("# Writing out the results\n")
        f.write(f"saveamberparm mol {prefix}.prmtop {prefix}.crd\n")
        f.write("quit")


def _prepare_build(
    mol,
    ff=None,
    topo=None,
    param=None,
    prefix="structure",
    outdir="./build",
    disulfide=None,
    teleap=None,
    atomtypes=None,
    offlibraries=None,
    gbsa=False,
    igb=2,
    custombonds=None,
    remove=None,
    cyclic_segids=None,
):
    """Copy forcefield/param/topo files, write PDB files, and generate tleap scripts.

    Splits the molecule into up to three PDB files:
    - input.pdb: non-water, non-cyclic solute atoms
    - solvent.pdb: water molecules
    - cyclic_<segid>.pdb: cyclic peptide segments (if any)

    Generates two tleap scripts:
    - tleap.in: full build (solute + solvent + cyclic)
    - tleap_solute.in: solute-only build for fast charge calculation (only when water is present)

    Returns
    -------
    tuple of (str, str)
        The (backend, value) tuple from _resolve_backend.
    """
    backend, value = _resolve_backend(teleap)
    if backend == "native":
        amberhome = os.path.normpath(os.path.join(os.path.dirname(value), "../"))
    else:
        amberhome = value

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # --- Copy force field files ---
    ff_sources = []
    for i, force in enumerate(ensurelist(ff)):
        if not os.path.isfile(force):
            force = _locateFile(force, "ff", amberhome)
            if force is None:
                continue
        newname = f"ff{i}_{os.path.basename(force)}"
        shutil.copy(force, os.path.join(outdir, newname))
        ff_sources.append(newname)

    if gbsa:
        if igb not in (1, 2, 5, 7, 8):
            raise ValueError(
                f"Invalid igb value {igb}. Supported values are: [1, 2, 5, 7, 8]"
            )

    # --- Copy OFF libraries ---
    offlib_names = []
    if offlibraries is not None:
        offlibraries = ensurelist(offlibraries)
        for i, off in enumerate(offlibraries):
            if not os.path.isfile(off):
                raise RuntimeError(f"Could not find off-library in location {off}")
            newname = f"offlib{i}_{os.path.basename(off)}"
            shutil.copy(off, os.path.join(outdir, newname))
            offlib_names.append(newname)

    # --- Copy param files ---
    newparam = []
    for i, fname in enumerate(sorted(param, key=lambda x: os.path.basename(x))):
        if not os.path.isfile(fname):
            fname = _locateFile(fname, "param", amberhome)
            if fname is None:
                continue
        newname = os.path.join(outdir, f"param{i}_{os.path.basename(fname)}")
        shutil.copy(fname, newname)
        newparam.append(newname)

    # --- Copy topo files ---
    newtopo = []
    for i, fname in enumerate(sorted(topo, key=lambda x: os.path.basename(x))):
        if not os.path.isfile(fname):
            fname = _locateFile(fname, "topo", amberhome)
            if fname is None:
                continue
        bn = os.path.basename(fname)
        newname = os.path.join(outdir, f"topo{i}_{bn}")
        shutil.copy(fname, newname)
        newtopo.append(newname)

    _fix_parameterize_atomtype_collisions(mol, newparam, newtopo)

    # --- Generate topology loading commands ---
    topo_cmds = []
    for fname in newtopo:
        ext = os.path.splitext(fname)[1].lower()
        if ext in (".mol2", ".cif"):
            # Load the ligand via mol2 (AMBER's officially recommended format
            # for small molecules).  CIFs are converted to mol2 first since
            # tleap only reads mol2/mol3.  The loader emits explicit
            # `set .element` commands for every atom so that custom atom
            # types (OpenFF, parameterize collision-avoidance types, etc.)
            # don't break the BONDS_INC_HYDROGEN split.
            _resmol = Molecule(fname)
            mol2_path = f"{os.path.splitext(fname)[0]}.mol2"
            if ext == ".cif":
                _resmol.write(mol2_path)
            # Use the residue name as the tleap variable name - tleap's
            # `loadpdb` looks up residue templates in the global variable
            # dictionary by residue name, so the unit must be registered
            # under that exact name.
            unit_name = _resmol.resname[0]
            topo_cmds.append(
                _build_load_mol2_commands(_resmol, mol2_path, unit_name=unit_name)
            )
        else:
            topo_cmds.append(f"loadamberprep {os.path.basename(fname)}\n")

    param_names = [os.path.basename(p) for p in newparam]

    # --- Cyclic segments (decided earlier, while bonds were present) and write
    #     PDB files. Break auto-sequencing at custom-bond (e.g. isopeptide)
    #     junctions on the per-PDB copies only - the chain reassignment must not
    #     touch the molecule used below for residue-position / custombond
    #     resolution (those UniqueAtomIDs match on chain). ---
    cyclic = _cyclic_segment_endpoints(mol, cyclic_segids or [])
    cyclic_segs = [c[0] for c in cyclic]
    break_points = _custombond_break_points(mol, custombonds)

    nonc_mol = mol.copy()
    if len(cyclic):
        nonc_mol.remove(np.isin(mol.segid, cyclic_segs), _logger=False)

    water_sel = nonc_mol.atomselect("water")
    has_water = np.any(water_sel)
    has_solute = np.any(~water_sel)

    if has_solute:
        solute_mol = nonc_mol.copy(sel=~water_sel)
        solute_mol = _apply_chain_breaks(solute_mol, break_points)
        solute_mol.write(os.path.join(outdir, "input.pdb"))

    if has_water:
        water_mol = nonc_mol.copy(sel=water_sel)
        water_mol.write(os.path.join(outdir, "solvent.pdb"))

    cyc_info = []
    if len(cyclic):
        for seg, res_start, res_end in cyclic:
            seg_mol = mol.copy(sel=mol.segid == seg)
            seg_mol = _apply_chain_breaks(seg_mol, break_points)
            fname = f"cyclic_{seg}.pdb"
            seg_mol.write(os.path.join(outdir, fname))
            cyc_var = f"cyc_{seg}"
            cyc_info.append((cyc_var, fname, res_start, res_end))

    # --- Write tleap scripts ---
    script_kwargs = dict(
        ff_sources=ff_sources,
        gbsa=gbsa,
        igb=igb,
        atomtypes=atomtypes,
        offlib_names=offlib_names,
        topo_cmds=topo_cmds,
        param_names=param_names,
        has_solute=has_solute,
        has_water=has_water,
        cyc_info=cyc_info,
        mol=mol,
        remove=remove,
        disulfide=disulfide,
        custombonds=custombonds,
    )

    _write_tleap_script(
        os.path.join(outdir, "tleap.in"),
        include_solvent=True,
        prefix=prefix,
        **script_kwargs,
    )

    if has_water:
        _write_tleap_script(
            os.path.join(outdir, "tleap_solute.in"),
            include_solvent=False,
            prefix="solute_charge",
            **script_kwargs,
        )

    return backend, value


def _read_tleap_output(outdir, prefix, logpath):
    """Shared post-processing: parse log, read prmtop/crd, return Molecule."""
    if not os.path.exists(logpath) or os.path.getsize(logpath) == 0:
        raise BuildError(
            "teLeap produced no output log. It likely crashed before initialization."
        )

    errors = _logParser(logpath)
    if errors:
        raise BuildError(
            errors
            + [f"Check {logpath} for further information on errors in building."],
            errors,
        )
    logger.info("Finished building.")

    crd_path = os.path.join(outdir, f"{prefix}.crd")
    prmtop_path = os.path.join(outdir, f"{prefix}.prmtop")
    if (
        os.path.exists(crd_path)
        and os.path.getsize(crd_path) != 0
        and os.path.exists(prmtop_path)
        and os.path.getsize(prmtop_path) != 0
    ):
        try:
            molbuilt = Molecule(prmtop_path, validateElements=False)
            # Force the AMBER inpcrd reader - tLeap's `saveamberparm` produces
            # an inpcrd-format file with a .crd extension, but moleculekit
            # defaults to the CHARMM .crd reader on that extension and loses
            # the periodic box stored on the file's last line.
            molbuilt.read(crd_path, type="inpcrd")
        except Exception as e:
            raise RuntimeError(
                f"Failed at reading {prefix}.prmtop/{prefix}.crd due to error: {e}"
            )
    else:
        raise BuildError(
            f"No {prefix} pdb/prmtop file was generated. Check {logpath} for errors in building."
        )

    molbuilt.write(os.path.join(outdir, f"{prefix}.pdb"), writebonds=False)
    detectCisPeptideBonds(molbuilt, respect_bonds=True)
    return molbuilt


def _run_tleap_native(
    outdir, prefix, teleap, amberhome, ff=None, teleapimports=None, script="tleap.in"
):
    """Execute native teLeap binary as a subprocess."""
    teleapimportflags = _getTeLeapImportFlags(amberhome, ff, teleapimports)
    logpath = os.path.abspath(os.path.join(outdir, "log.txt"))
    logger.info("Starting the build.")

    with open(logpath, "w") as f:
        cmd = [teleap, "-f", f"./{script}"]
        cmd[1:1] = teleapimportflags
        logger.debug(cmd)
        result = subprocess.run(cmd, stdout=f, stderr=f, cwd=outdir)

    if result.returncode != 0:
        errors = _logParser(logpath)
        if errors:
            raise BuildError(
                errors
                + [f"Check {logpath} for further information on errors in building."],
                errors,
            )
        raise BuildError(
            f"teLeap exited with code {result.returncode}.\n"
            f"Check {logpath} for details."
        )

    return _read_tleap_output(outdir, prefix, logpath)


def _run_tleap_pyodide(outdir, prefix, ff=None, script="tleap.in"):
    """Execute tleap via tleap_pyodide.run_tleap using the work_dir parameter."""
    from tleap_pyodide import run_tleap as _pyodide_run

    logger.info("Starting the build (pyodide).")

    extra_imports = _getHtmdExtraImports(ff)

    result = _pyodide_run(
        work_dir=outdir,
        script_file=script,
        extra_imports=extra_imports if extra_imports else None,
    )

    logpath = os.path.abspath(os.path.join(outdir, "log.txt"))
    with open(logpath, "wb") as f:
        f.write(result.get("stdout", b""))

    return _read_tleap_output(outdir, prefix, logpath)


def _run_tleap(
    outdir,
    prefix,
    backend="native",
    backend_value=None,
    ff=None,
    teleapimports=None,
    script="tleap.in",
):
    """Dispatch tleap execution to native or pyodide backend."""
    if backend == "pyodide":
        return _run_tleap_pyodide(outdir, prefix, ff=ff, script=script)
    else:
        teleap = backend_value
        if teleap is None:
            teleap = _findTeLeap()
        amberhome = os.path.normpath(
            os.path.join(os.path.dirname(teleap), "../")
        )
        return _run_tleap_native(
            outdir,
            prefix,
            teleap,
            amberhome,
            ff=ff,
            teleapimports=teleapimports,
            script=script,
        )


# Longest first-N to last-C distance still accepted as a real head-to-tail
# amide closure when the input declares it as an explicit bond. A real amide
# C-N is ~1.33 A; even poorly-modeled ones stay well under this. Beyond it an
# explicit closure bond is treated as a misassigned LINK / CONECT record.
MAX_CYCLIC_CLOSURE_DIST = 2.0


def _detect_cyclic_segments(mol: Molecule):
    cyclic = []
    prot = mol.atomselect("protein")
    for seg in np.unique(mol.segid):
        segatoms = mol.segid == seg
        if not np.any(prot[segatoms]):
            # Skip segments with no peptide character at all (water, ions, free
            # ligands). A cyclic peptide is still admitted when some of its
            # residues are non-canonical - e.g. microcystin's beta-amino-acid
            # Adda fails the "protein" selection, which previously (np.all)
            # silently dropped the whole macrocycle. The first-N / last-C
            # closure test below is the actual cyclic criterion.
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
            continue
        dist = np.linalg.norm(
            mol.coords[first_mask, :, 0] - mol.coords[last_mask, :, 0]
        )
        # An explicit head-to-tail amide bond between the first residue's N
        # and the last residue's C is authoritative: the peptide is cyclic no
        # matter how long that bond was modeled, as long as the distance stays
        # physically plausible. Crystallographic amides routinely exceed the
        # 1.35 A heuristic below (e.g. 7BTI's phalloidin closure is 1.468 A),
        # so honor the input bond up to MAX_CYCLIC_CLOSURE_DIST. Beyond that
        # the "bond" is almost certainly a misassigned LINK / CONECT record;
        # ignore it rather than ask tLeap to close a ring across a nonsensical
        # gap. Fall back to the distance heuristic when no such bond exists.
        has_closure_bond = False
        if len(mol.bonds):
            pair = {int(np.where(first_mask)[0][0]), int(np.where(last_mask)[0][0])}
            has_closure_bond = any(
                {int(b[0]), int(b[1])} == pair for b in mol.bonds
            )
        if has_closure_bond and dist > MAX_CYCLIC_CLOSURE_DIST:
            logger.warning(
                f"Segment {seg} carries an explicit first-N to last-C bond, but "
                f"it spans {dist:.2f} A - too long for a real amide. Ignoring it "
                "for cyclic-peptide detection (likely a misassigned bond record)."
            )
            has_closure_bond = False
        # Amide bond distance ranges between 1.325 - 1.346
        if has_closure_bond or dist < 1.35:
            cyclic.append((seg, first_resid, last_resid))
    return cyclic


def _cyclic_segment_endpoints(mol: Molecule, cyclic_segids):
    """Return ``(segid, first_resid, last_resid)`` for segments already decided
    cyclic (by :func:`_detect_cyclic_segments` run earlier, while the input
    bonds were still present). The endpoints are read positionally from the
    current molecule - cyclic segments are never capped, so cap-insertion does
    not move their first/last residue."""
    out = []
    for seg in cyclic_segids:
        segatoms = mol.segid == seg
        if not np.any(segatoms):
            continue
        resids = mol.resid[segatoms]
        out.append((seg, int(resids[0]), int(resids[-1])))
    return out


_CHAIN_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _custombond_break_points(mol: Molecule, custombonds):
    """Find consecutive-residue junctions joined by a custom (e.g. isopeptide)
    bond rather than a standard backbone amide. Returns a set of
    ``(segid, lower_resid)`` "break-after" points. Resolved on the full molecule
    where the custombond ``UniqueAtomID``s match."""
    points = set()
    if not custombonds:
        return points
    for a1, a2 in custombonds:
        i1 = int(np.atleast_1d(a1.selectAtom(mol))[0])
        i2 = int(np.atleast_1d(a2.selectAtom(mol))[0])
        if str(mol.segid[i1]) != str(mol.segid[i2]):
            continue
        r1, r2 = int(mol.resid[i1]), int(mol.resid[i2])
        if abs(r1 - r2) == 1:
            points.add((str(mol.segid[i1]), min(r1, r2)))
    return points


def _apply_chain_breaks(mol: Molecule, break_points):
    """Reassign chains so tLeap writes a TER (and does not head-to-tail
    auto-sequence) at each break-after point. A custom-bonded junction's real
    bond is emitted explicitly; the spurious backbone auto-bond must not be
    created (e.g. it would bond microcystin's ACB.CA to the next residue's N and
    corrupt CA's chirality). Only the chain field changes - resids and atom order
    are preserved, so residue numbering and bond references are unaffected. Apply
    only to the per-PDB copies written for tLeap, never to the molecule used for
    residue-position / custombond resolution (those ``UniqueAtomID``s match on
    chain)."""
    from collections import defaultdict

    byseg = defaultdict(list)
    for seg, resid in break_points:
        if np.any(mol.segid == seg):
            byseg[seg].append(resid)
    if not byseg:
        return mol
    mol = mol.copy()
    used = set(np.unique(mol.chain).tolist())
    pool = [c for c in _CHAIN_POOL if c not in used]
    pi = 0
    for seg, brks in byseg.items():
        segmask = mol.segid == seg
        prev = None
        for b in sorted(brks):
            if prev is not None:
                if pi >= len(pool):
                    raise RuntimeError(
                        "Ran out of free chain IDs while breaking auto-sequencing "
                        "at custom-bond junctions."
                    )
                mol.chain[segmask & (mol.resid > prev) & (mol.resid <= b)] = pool[pi]
                pi += 1
            prev = b
        if pi >= len(pool):
            raise RuntimeError(
                "Ran out of free chain IDs while breaking auto-sequencing at "
                "custom-bond junctions."
            )
        mol.chain[segmask & (mol.resid > prev)] = pool[pi]
        pi += 1
    return mol


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

    # Remove backbone terminal hydrogens regardless of caps or no caps.
    # tleap confuses itself when it makes residues into terminals.
    # For example HID has an atom H which in NHID should become H[123] and crashes tleap.
    sidechain_h = ("HG", "HB2", "HB3")
    for seg in np.unique(mol.get("segid", "protein")):
        seg_mask = mol.segid == seg
        resids = mol.resid[seg_mask]
        mask = (
            ((mol.resid == resids[0]) | (mol.resid == resids[-1]))
            & seg_mask
            & (mol.element == "H")
            & (~np.isin(mol.name, sidechain_h))  # Don't remove sidechain hydrogens
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
                "default. If you want to cap it use the caps argument of amber.build to manually define caps "
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

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File {filename} does not exist")

    resdict = dict()

    Rule = namedtuple(
        "Rule", ["replaceresname", "replaceatom", "order", "natoms", "ter"]
    )

    with open(filename, "r") as csvfile:
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
                if line.strip() == "" or any(line.startswith(x) for x in const):
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
                raise RuntimeError(
                    f"Don't know how to repair topologies with extension {ext}"
                )

        if bn in mol.resname:
            _fix_mol(bn, mol, replacements[bn])
        else:
            raise RuntimeError(
                f"Could not find residue {bn} in the input structure. Please ensure that the frcmod / prepi files are named as RES.frcmod / RES.prepi where RES the residue name to which they correspond."
            )
