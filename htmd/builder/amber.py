# (c) 2015-2018 Acellera Ltd http://www.acellera.com
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
    UnknownResidueError,
    MissingAngleError,
    MissingBondError,
    MissingParameterError,
    MissingTorsionError,
    MissingAtomTypeError,
)
from subprocess import call, check_output, DEVNULL
from moleculekit.molecule import Molecule
from htmd.builder.ionize import ionize as ionizef, ionizePlace
from htmd.util import ensurelist
import unittest

import logging

logger = logging.getLogger(__name__)


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


def htmdAmberHome():
    """ Returns the location of the AMBER files distributed with HTMD"""

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


def defaultFf():
    """ Returns the default leaprc forcefield files used by amber.build """
    return [
        "leaprc.protein.ff19SB",
        "leaprc.lipid17",
        "leaprc.gaff2",
        "leaprc.RNA.Shaw",
        "leaprc.DNA.bsc1",
        "leaprc.water.tip3p",
    ]


def defaultTopo():
    """ Returns the default topology `prepi` files used by amber.build """
    return []


def defaultParam():
    """ Returns the default parameter `frcmod` files used by amber.build """
    return []


def _prepareMolecule(mol, caps, disulfide):
    # Remove pdb protein bonds as they can be regenerated by tleap. Keep non-protein bonds i.e. for ligands
    _removeProteinBonds(mol)

    # Check for missing segids or mixed protein / non-protein segments
    _missingSegID(mol)
    _checkMixedSegment(mol)

    # Convert lipids to AMBER naming
    mol = _charmmLipid2Amber(mol)

    # Add caps to termini
    if caps is None:
        caps = _defaultProteinCaps(mol)
    _applyProteinCaps(mol, caps)

    # Before modifying the resids copy the molecule to map back the disulfide bonds
    mol_orig_resid = mol.copy()

    # We need to renumber residues to unique numbers for teLeap. It goes up to 9999 and restarts from 0
    mol.resid[:] = (sequenceID((mol.resid, mol.insertion, mol.segid)) + 1) % 10000

    # curr_resid = 0
    # prev_combo = (None, None, None)
    # for i in range(mol.numAtoms):
    #     curr_combo = (mol.resid[i], mol.insertion[i], mol.segid[i])
    #     if prev_combo != curr_combo:
    #         curr_resid += 1
    #         if prev_combo[0] is not None and curr_combo[0] - prev_combo[0] > 1:
    #             curr_resid += 1  # Insert a resid gap if a gap existed before
    #         prev_combo = curr_combo

    #     mol.resid[i] = curr_resid % 10000

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
            # Remove (eventual) HG hydrogens on these CYS (from proteinPrepare)
            torem |= (atoms1 & (mol.name == "HG")) | (atoms2 & (mol.name == "HG"))
        mol.remove(torem, _logger=False)

    return disulfide


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
    mol = mol.copy()

    disulfide = _prepareMolecule(mol, caps, disulfide)

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
    return molbuilt


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
    if ff is None:
        ff = defaultFf()
    if topo is None:
        topo = defaultTopo()
    if param is None:
        param = defaultParam()

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

        # Loading frcmod parameters
        f.write("# Loading parameter files\n")
        for i, p in enumerate(param):
            if not os.path.isfile(p):
                p = _locateFile(p, "param", teleap)
                if p is None:
                    continue
            newname = f"param{i}_{os.path.basename(p)}"
            shutil.copy(p, os.path.join(outdir, newname))
            f.write(f"loadamberparams {newname}\n")
        f.write("\n")

        # Loading prepi topologies
        f.write("# Loading prepi topologies\n")
        for i, t in enumerate(topo):
            if not os.path.isfile(t):
                t = _locateFile(t, "topo", teleap)
                if t is None:
                    continue
            newname = f"topo{i}_{os.path.basename(t)}"
            shutil.copy(t, os.path.join(outdir, newname))
            f.write(f"loadamberprep {newname}\n")
        f.write("\n")

        f.write("# Loading the system\n")
        f.write("mol = loadpdb input.pdb\n\n")

        if np.sum(mol.atomtype != "") != 0:
            f.write("# Loading the ligands\n")
            segs = np.unique(mol.segid[mol.atomtype != ""])

            # teLeap crashes if you try to combine too many molecules in a single command so we will do them by 10s
            for k in range(0, len(segs), 10):
                segments_string = ""
                for seg in segs[k : min(k + 10, len(segs))]:
                    name = f"segment{seg}"
                    segments_string += f" {name}"

                    mol2name = os.path.join(outdir, f"{name}.mol2")
                    mol.write(mol2name, (mol.atomtype != "") & (mol.segid == seg))
                    if not os.path.isfile(mol2name):
                        raise NameError("Failed writing ligand mol2 file.")

                    f.write(f"{name} = loadmol2 {name}.mol2\n")
                f.write(f"mol = combine {{mol{segments_string}}}\n\n")

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

    # Printing and loading the PDB file. AMBER can work with a single PDB file if the segments are separate by TER
    logger.debug("Writing PDB file for input to tleap.")
    pdbname = os.path.join(outdir, "input.pdb")

    # mol2 files have atomtype, here we only write parts not coming from mol2
    # We need to write the input.pdb at the end since we modify the resname for disulfide bridges in mol
    mol.write(pdbname, mol.atomtype == "")
    if not os.path.isfile(pdbname):
        raise NameError("Could not write a PDB file out of the given Molecule.")

    if execute:
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
        except:
            raise NameError("teLeap failed at execution")
        f.close()
        errors = _logParser(logpath)
        os.chdir(currdir)
        if errors:
            raise BuildError(
                errors
                + [f"Check {logpath} for further information on errors in building."]
            )
        logger.info("Finished building.")

        if (
            os.path.exists(os.path.join(outdir, "structure.crd"))
            and os.path.getsize(os.path.join(outdir, "structure.crd")) != 0
            and os.path.getsize(os.path.join(outdir, "structure.prmtop")) != 0
        ):
            try:
                molbuilt = Molecule(os.path.join(outdir, "structure.prmtop"))
                molbuilt.read(os.path.join(outdir, "structure.crd"))
            except Exception as e:
                raise RuntimeError(
                    f"Failed at reading structure.prmtop/structure.crd due to error: {e}"
                )
        else:
            raise BuildError(
                f"No structure pdb/prmtop file was generated. Check {logpath} for errors in building."
            )

        tmpbonds = molbuilt.bonds
        molbuilt.bonds = []  # Removing the bonds to speed up writing
        molbuilt.write(os.path.join(outdir, "structure.pdb"))
        molbuilt.bonds = tmpbonds  # Restoring the bonds
        detectCisPeptideBonds(molbuilt, respect_bonds=True)  # Warn in case of cis bonds
        return molbuilt


def _applyProteinCaps(mol, caps):

    # AMBER capping
    # =============
    # This is the (horrible) way of adding caps in tleap:
    # For now, this is hardwired for ACE and NME
    # 1. Change one of the hydrogens of the N terminal (H[T]?[123]) to the ACE C atom, giving it a new resid
    # 1a. If no hydrogen present, create the ACE C atom.
    # 1b. Add the ACE CH3 atom to define the direction of the cap
    # 2. Change one of the oxygens of the C terminal ({O,OT1,OXT}) to the NME N atom, giving it a new resid
    # 2a. If no oxygen present, create the NME N atom.
    # 3. Reorder to put the new atoms first and last
    # 4. Remove the lingering hydrogens of the N terminal and oxygens of the C terminal.

    # Define the atoms to be replaced (0 and 1 corresponds to N- and C-terminal caps)
    terminalatoms = {
        "ACE": ["H1", "H2", "H3", "HT1", "HT2", "HT3"],
        "NME": ["OXT", "OT1", "O"],
    }  # XPLOR names for H[123] and OXT are HT[123]
    # and OT1, respectively.
    capresname = ["ACE", "NME"]
    capatomtype = ["C", "N"]

    # For each caps definition
    for seg in caps:
        # Can't move this out since we remove atoms in this loop
        prot = mol.atomselect("protein")
        # Get the segment
        segment = np.where(mol.segid == seg)[0]
        # Test segment
        if len(segment) == 0:
            raise RuntimeError(f"There is no segment {seg} in the molecule.")
        if not np.any(prot & (mol.segid == seg)):
            raise RuntimeError(
                f"Segment {seg} is not protein. Capping for non-protein segments is not supported."
            )
        # For each cap
        passed = False
        for i, cap in enumerate(caps[seg]):
            if cap is None or (isinstance(cap, str) and cap == "none"):
                continue
            # Get info on segment and its terminals
            segidm = mol.segid == seg  # Mask for segid
            segididx = np.where(segidm)[0]
            terminalids = [segididx[0], segididx[-1]]
            terminalresids = [mol.resid[segididx][0], mol.resid[segididx][-1]]
            terminalchains = [mol.chain[segididx][0], mol.chain[segididx][-1]]
            residm = mol.resid == terminalresids[i]  # Mask for resid

            if not passed:
                orig_terminalresids = terminalresids
                passed = True

            if cap is None or cap == "":  # In case there is no cap defined
                logger.warning(
                    f"No cap provided for resid {terminalresids[i]} on segment {seg}. Did not apply it."
                )
                continue
            elif cap not in capresname:  # If it is defined, test if supported
                raise RuntimeError(
                    f"In segment {seg}, the {cap} cap is not supported. Try using {capresname} instead."
                )

            # Test if cap is already applied
            testcap = np.where(segidm & residm & (mol.resname == cap))[0]
            if len(testcap) != 0:
                logger.warning(
                    f"Cap {cap} already exists on segment {seg}. Did not re-apply it."
                )
                continue

            # Test if the atom to change exists
            termatomsids = np.isin(mol.name, terminalatoms[cap])
            termatomsids = np.where(termatomsids & segidm & residm)[0]

            # Create new atom if none of the replace atoms exist
            if len(termatomsids) == 0:
                termcaid = np.where(segidm & residm & (mol.name == "CA"))[0]
                termcenterid = np.where(
                    segidm & residm & (mol.name == capatomtype[1 - i])
                )[0]
                atom = Molecule()
                atom.empty(1)
                atom.coords = mol.coords[termcenterid] + 0.33 * np.subtract(
                    mol.coords[termcenterid], mol.coords[termcaid]
                )
                new_idx = terminalids[i] if i == 0 else terminalids[i] + 1
                mol.insert(atom, new_idx)
                termatomsids = [new_idx]

            # Select atom to change, do changes to cap, and change resid
            replace_atom = np.max(termatomsids)
            mol.record[replace_atom] = "ATOM"
            mol.name[replace_atom] = capatomtype[i]
            mol.resname[replace_atom] = cap
            # if i=0 => resid-1; i=1 => resid+1
            # TODO: Increase all following resids???
            mol.resid[replace_atom] = terminalresids[i] - 1 + 2 * i
            mol.element[replace_atom] = capatomtype[i]
            mol.segid[replace_atom] = seg
            mol.chain[replace_atom] = terminalchains[i]

            # Swap positions of cap atom with the first/last atom of the segment
            neworder = np.arange(mol.numAtoms)
            neworder[replace_atom] = terminalids[i]
            neworder[terminalids[i]] = replace_atom
            _reorderMol(mol, neworder)

        # For each cap
        for i, cap in enumerate(caps[seg]):
            if cap is None or (isinstance(cap, str) and cap == "none"):
                continue
            # Remove lingering hydrogens or oxygens in terminals
            mol.remove(
                f'segid {seg} and resid "{orig_terminalresids[i]}" and name {" ".join(terminalatoms[cap])}',
                _logger=False,
            )

    # Remove terminal hydrogens regardless of caps or no caps. tleap confuses itself when it makes residues into terminals.
    # For example HID has an atom H which in NHID should become H[123] and crashes tleap.
    for seg in np.unique(mol.get("segid", "protein")):
        segidm = mol.segid == seg  # Mask for segid
        segididx = np.where(segidm)[0]
        resids = mol.resid[segididx]
        mol.remove(
            f'(resid "{resids[0]}" "{resids[-1]}") and segid {seg} and hydrogen',
            _logger=False,
        )


def _removeProteinBonds(mol):
    # Keeping just bonds with atomtypes to catch ligand cases. Bad heuristic
    segs = np.unique(mol.segid[mol.atomtype != ""])
    mol.deleteBonds(~np.in1d(mol.segid, segs))


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

    files = glob(os.path.join(outdir, "structure.*"))
    files += glob(os.path.join(outdir, "log.*"))
    files += glob(os.path.join(outdir, "*.log"))
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

    unknownres_regex = re.compile("Unknown residue:\s+(\w+)")
    missingparam_regex = re.compile(
        "For atom: (.*) Could not find vdW \(or other\) parameters for type:\s+(\w+)"
    )
    missingtorsion_regex = re.compile("No torsion terms for\s+(.*)$")
    missingbond_regex = re.compile("Could not find bond parameter for:\s+(.*)$")
    missingangle_regex = re.compile("Could not find angle parameter:\s+(.*)$")
    missingatomtype_regex = re.compile("FATAL:\s+Atom (.*) does not have a type")

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
            UnknownResidueError(
                f"Unknown residue(s) {np.unique(unknownres)} found in the input structure. "
                "You are either missing a topology definition for the residue or you need to "
                "rename it to the correct residue name"
            )
        )
    if len(missingparam):
        errors.append(
            MissingParameterError(
                f"Missing parameters for atom {missingparam} and type {missingparam}"
            )
        )
    if len(missingtorsion):
        errors.append(
            MissingTorsionError(f"Missing torsion terms for {missingtorsion}")
        )
    if len(missingbond):
        errors.append(MissingBondError(f"Missing bond parameters for {missingbond}"))
    if len(missingangle):
        errors.append(MissingAngleError(f"Missing angle parameters for {missingangle}"))
    if len(missingatomtype):
        errors.append(MissingAtomTypeError(f"Missing atom type for {missingatomtype}"))

    return errors


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
        from moleculekit.tools.preparation import proteinPrepare
        from htmd.builder.solvate import solvate
        from moleculekit.tools.autosegment import autoSegment

        # Test with proteinPrepare
        pdbids = ["3PTB"]
        # pdbids = ['3PTB', '1A25', '1GZM']  # '1U5U' out because it has AR0 (no parameters)
        for pid in pdbids:
            np.random.seed(1)
            mol = Molecule(pid)
            mol.filter("protein")
            mol = proteinPrepare(mol)
            mol.filter("protein")  # Fix for bad proteinPrepare hydrogen placing
            mol = autoSegment(mol)
            smol = solvate(mol)
            ffs = defaultFf()
            tmpdir = os.path.join(self.testDir, "withProtPrep", pid)
            _ = build(smol, ff=ffs, outdir=tmpdir)

            refdir = home(dataDir=join("test-amber-build", "pp", pid))
            _TestAmberBuild._compareResultFolders(refdir, tmpdir, pid)

    def test_without_protein_prepare(self):
        from htmd.builder.solvate import solvate

        # Test without proteinPrepare
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
        from htmd.builder.solvate import solvate

        # Test protein ligand building with parametrized ligand
        refdir = home(dataDir=join("test-amber-build", "protLig"))
        tmpdir = os.path.join(self.testDir, "protLig")
        os.makedirs(tmpdir)

        mol = Molecule(join(refdir, "3ptb_mod.pdb"))
        lig = Molecule(join(refdir, "benzamidine.mol2"))
        lig.segid[:] = "L"

        # params =
        newmol = Molecule()
        newmol.append(mol)
        newmol.append(lig)
        smol = solvate(newmol)

        params = defaultParam() + [
            join(refdir, "benzamidine.frcmod"),
        ]

        _ = build(smol, outdir=tmpdir, param=params, ionize=False)

        refdir = home(dataDir=join("test-amber-build", "protLig", "results"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "3PTB")

    def test_custom_disulfide_bonds(self):
        from htmd.builder.solvate import solvate

        # Test without proteinPrepare
        pdbids = [
            "1GZM",
        ]
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

        # Test without proteinPrepare
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
        from moleculekit.tools.preparation import proteinPrepare
        from moleculekit.tools.autosegment import autoSegment

        np.random.seed(1)

        mol = Molecule("3WBM")
        mol.filter("not water")
        mol = autoSegment(mol, field="both")
        pmol = proteinPrepare(mol)
        smol = solvate(pmol)

        tmpdir = os.path.join(self.testDir, "protein-rna", "3WBM")
        _ = build(smol, outdir=tmpdir)

        refdir = home(dataDir=join("test-amber-build", "protein-rna", "3WBM"))
        _TestAmberBuild._compareResultFolders(refdir, tmpdir, "3WBM")

    def test_caps(self):
        # TODO: Add test where I remove all side chain from the termini and force build to reconstruct ACE NME from scratch
        pass


if __name__ == "__main__":
    import doctest
    import unittest

    failure_count, _ = doctest.testmod()
    unittest.main(verbosity=2)
    # if failure_count != 0:
    #     raise Exception('Doctests failed')
