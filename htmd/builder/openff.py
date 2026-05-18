# (c) 2015-2025 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
"""OpenMM/OpenFF-based system builder for molecular dynamics.

Provides an alternative to the tleap-based AMBER builder, using OpenMM
force fields for parameterization and OpenFF Interchange (or ParmEd)
for AMBER prmtop/inpcrd export.

Dependencies (imported lazily at call time):
    openmm, openmmforcefields, openff-toolkit, openff-interchange
"""

import contextlib
import os
import numpy as np
import logging

from moleculekit.molecule import Molecule
from moleculekit.util import _missingSegID, sequenceID
from htmd.builder.builder import (
    detectDisulfideBonds,
    detectCisPeptideBonds,
    convertDisulfide,
    _checkMixedSegment,
    BuildError,
)

logger = logging.getLogger(__name__)


# ====================================================================
# Public API
# ====================================================================


def defaultFf():
    """Returns the default OpenMM XML force field files.

    The ``amber14/`` prefix ships with OpenMM itself.  If the
    ``openmmforcefields`` package is installed, force fields under the
    ``amber/`` or ``amber19/`` prefixes are also available.

    Returns
    -------
    list of str
        Force field XML file names loadable by ``openmm.app.ForceField``.
    """
    return [
        "amber14/protein.ff14SB.xml",
        "amber14/lipid17.xml",
        "amber14/DNA.bsc1.xml",
        "amber14/RNA.OL3.xml",
        "amber14/tip3p.xml",
    ]


def build(
    mol,
    ff=None,
    extra_xml=None,
    small_molecule_ff=None,
    molecules=None,
    prefix="structure",
    outdir="./build",
    caps=None,
    ionize=True,
    saltconc=0,
    saltanion=None,
    saltcation=None,
    disulfide=None,
    custombonds=None,
    solvate=True,
    padding=10.0,
    water_model="tip3p",
    boxsize=None,
    gbsa=False,
):
    """Build a system using OpenMM force fields and export to AMBER format.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The input molecular system.
    ff : list of str, optional
        OpenMM XML force field file names (see ``openmm.app.ForceField``).
        Default: :func:`defaultFf`.
    extra_xml : str or list of str, optional
        Paths to additional OpenMM XML files for non-standard residues.
    small_molecule_ff : str, optional
        Small-molecule force field for the template generator, e.g.
        ``"gaff-2.2.20"`` or ``"openff-2.3.0"``.  Requires *molecules*.
    molecules : list, optional
        ``openff.toolkit.Molecule`` objects **or** paths to SDF files
        describing the small molecules present in the system.
    prefix : str
        Prefix for output files.
    outdir : str
        Output directory path.
    caps : dict, optional
        Terminal capping specification.  Accepts two formats:

        * ``{segid: [nterm_cap, cterm_cap]}`` — e.g. ``{"P": ["ACE","NME"]}``
        * ``{atomsel: cap}`` — e.g. ``{"chain A and resid 5": "ACE"}``

        Default: ACE/NME on every protein segment with >= 10 residues.
    ionize : bool
        Neutralise the system (and add salt if *saltconc* > 0).
    saltconc : float
        Salt concentration in molar, added on top of neutralisation.
    saltanion : str, optional
        Anion type.  Accepts ``"Cl-"``, ``"CL"``, ``"chloride"`` etc.
    saltcation : str, optional
        Cation type.  Accepts ``"Na+"``, ``"K+"``, ``"Cs+"`` and
        divalent ``"Mg2+"``, ``"Ca2+"``, ``"Zn2+"``.
    disulfide : list of pairs of str, optional
        Manual disulfide bonds as pairs of atom-selection strings.
        ``None`` triggers automatic detection.
    custombonds : list of pairs of str, optional
        Extra bonds as pairs of atom-selection strings.
    solvate : bool
        Add explicit water via ``Modeller.addSolvent()``.
    padding : float
        Box padding in **Angstroms** (used when *solvate* is True).
    water_model : str
        Water model name for ``Modeller.addSolvent()``.
    boxsize : float or list of float, optional
        Explicit box dimensions ``[x, y, z]`` in Angstroms.  Overrides
        *padding*.
    gbsa : bool
        Use GBSA implicit solvent (OBC2 model).

    Returns
    -------
    molbuilt : :class:`Molecule <moleculekit.molecule.Molecule>`
        The fully built system.
    system : ``openmm.System``
        The parameterised OpenMM System object.
    """
    import openmm
    import openmm.app as app
    import openmm.unit as unit

    mol = mol.copy()
    if ff is None:
        ff = defaultFf()

    disulfide_ids, custombond_ids, cyclic = _prepare_molecule(
        mol, caps, disulfide, custombonds
    )

    # Extract bond specs as topology-independent identifiers (chain, resid,
    # insertion, atom_name).  This must happen *before* the OpenMM conversion
    # because Modeller.addHydrogens() may change atom indices.
    bond_specs = _extract_bond_specs(mol, cyclic, disulfide_ids, custombond_ids)

    os.makedirs(outdir, exist_ok=True)

    topology, positions = _mol_to_openmm(mol, outdir, extra_xml=extra_xml)

    forcefield = _setup_forcefield(ff, extra_xml, small_molecule_ff, molecules)

    # Cyclic, disulfide, and custom bonds must be registered in the
    # topology *before* addHydrogens so that OpenMM's template matcher
    # sees the correct external-bond connectivity (e.g. cyclic peptide
    # terminals are not free, CYX pairs are already bonded).
    _apply_bond_specs(topology, bond_specs)

    # Add missing hydrogens.  This can fail for unusual topologies
    # (cyclic peptides, stripped sidechains, non-standard residues)
    # where the template matcher cannot find a match.  In those cases
    # the input structure must already contain all required H atoms.
    try:
        modeller = app.Modeller(topology, positions)
        modeller.addHydrogens(forcefield)
        topology = modeller.topology
        positions = modeller.positions
        # Re-apply after addHydrogens (it rebuilds the topology,
        # possibly dropping our custom bonds).
        _apply_bond_specs(topology, bond_specs)
    except ValueError as exc:
        logger.warning(
            f"Modeller.addHydrogens() failed ({exc}).  "
            "Proceeding with the input hydrogen atoms as-is."
        )

    has_water = np.any(mol.atomselect("water"))
    did_solvate = False

    if solvate and not has_water:
        logger.info("Solvating system with Modeller.addSolvent()...")
        modeller = app.Modeller(topology, positions)
        solvent_kw = _build_solvent_kwargs(
            water_model,
            padding,
            boxsize,
            ionize,
            saltconc,
            saltanion,
            saltcation,
            unit,
        )
        modeller.addSolvent(forcefield, **solvent_kw)
        topology = modeller.topology
        positions = modeller.positions
        did_solvate = True
    elif solvate and has_water:
        logger.warning(
            "Molecule already contains water — skipping solvation.  "
            "Pass solvate=False to silence this warning."
        )

    if has_water and ionize and not did_solvate:
        logger.info(
            "Pre-solvated system detected.  Building temporary system "
            "to calculate charge for ionisation..."
        )
        _ensure_box_vectors(topology, positions, unit, mol)
        tmp_system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=app.HBonds,
            rigidWater=True,
        )
        mol = _ionize_presolvated(mol, tmp_system, saltconc, saltanion, saltcation)
        topology, positions = _mol_to_openmm(mol, outdir, extra_xml=extra_xml)
        modeller = app.Modeller(topology, positions)
        modeller.addHydrogens(forcefield)
        topology = modeller.topology
        positions = modeller.positions
        _apply_bond_specs(topology, bond_specs)

    has_water_final = has_water or did_solvate
    # Always assign periodic box vectors to the topology unless GBSA
    # (which is intrinsically non-periodic). ``amber.build`` writes
    # ``setBox mol "vdw"`` in its tLeap script and so its prmtops
    # always carry a box, even without solvation. Mirroring that
    # behaviour here keeps the two output prmtops comparable AND
    # satisfies ``Interchange.from_openmm``, which requires the
    # NonbondedForce to use PME (and therefore the topology to be
    # periodic) at export time.
    if not gbsa:
        _ensure_box_vectors(topology, positions, unit, mol)

    logger.info("Creating OpenMM System...")

    # For GBSA we need the implicit solvent XML loaded into the FF used
    # to create the simulation System.  But ParmEd cannot export GBSA
    # systems, so we keep a clean copy of the FF for the AMBER export.
    export_ff = forcefield
    if gbsa:
        import copy

        export_ff = copy.deepcopy(forcefield)
        forcefield.loadFile("implicit/obc2.xml")

    sys_kw = {
        "nonbondedCutoff": 1.0 * unit.nanometers,
        "constraints": app.HBonds,
    }
    if gbsa:
        sys_kw["nonbondedMethod"] = app.NoCutoff
    elif has_water_final:
        sys_kw["nonbondedMethod"] = app.PME
        sys_kw["rigidWater"] = True
    else:
        sys_kw["nonbondedMethod"] = app.NoCutoff

    system = forcefield.createSystem(topology, **sys_kw)

    # Serialize the OpenMM System as XML alongside the prmtop. The
    # XML preserves the System exactly as the ForceField built it
    # (per-instance parameters, full force-type info), so it can be
    # round-tripped through ``openmm.XmlSerializer`` without going
    # through ParmEd's atom-type bucketing. Useful for downstream
    # energy comparisons and reproducibility.
    system_xml_path = os.path.join(outdir, f"{prefix}.system.xml")
    with open(system_xml_path, "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))
    logger.info(f"Wrote {system_xml_path}")

    _export_amber(topology, system, positions, outdir, prefix, export_ff)

    pdb_path = os.path.join(outdir, f"{prefix}.pdb")
    with open(pdb_path, "w") as fh:
        app.PDBFile.writeFile(topology, positions, fh)
    logger.info(f"Wrote {pdb_path}")

    molbuilt = _read_built_molecule(outdir, prefix)
    detectCisPeptideBonds(molbuilt, respect_bonds=True)

    logger.info("Finished building.")
    return molbuilt, system


# ====================================================================
# Molecule preparation
# ====================================================================


def _prepare_molecule(mol, caps, disulfide, custombonds):
    """Prepare *mol* in-place: caps, disulfide, lipid naming, etc.

    Returns ``(disulfide_ids, custombond_ids, cyclic)`` where the first
    two are lists of ``UniqueResidueID`` / ``UniqueAtomID`` pairs and
    *cyclic* is the list returned by ``_detect_cyclic_segments``.
    """
    from htmd.builder.amber import _detect_cyclic_segments

    _missingSegID(mol)
    _checkMixedSegment(mol)

    _fix_water_naming(mol)
    _fix_nucleic_naming(mol)

    cyclic = _detect_cyclic_segments(mol)

    if caps is None:
        caps = _defaultProteinCaps(mol)
    for seg, _, _ in cyclic:
        if seg in caps:
            logger.info(f"Found cyclic segment {seg}.  Disabling capping on it.")
            del caps[seg]
    _add_caps(mol, caps)

    # --- Disulfide bonds ---
    if disulfide is not None and len(disulfide):
        if not isinstance(disulfide[0][0], str):
            raise RuntimeError("Can only accept disulfide bond strings")
        disulfide = convertDisulfide(mol, disulfide)
    else:
        logger.info("Detecting disulfide bonds.")
        disulfide = detectDisulfideBonds(mol)

    # --- Custom bonds ---
    if custombonds is not None:
        from moleculekit.molecule import UniqueAtomID

        new_cb = []
        for bb in custombonds:
            if not isinstance(bb[0], str) or not isinstance(bb[1], str):
                raise RuntimeError("All custombonds selections should be strings")
            new_cb.append(
                [
                    UniqueAtomID.fromMolecule(mol, bb[0]),
                    UniqueAtomID.fromMolecule(mol, bb[1]),
                ]
            )
        custombonds = new_cb

    # --- Apply CYS→CYX rename and HG removal for disulfide bonds ---
    if len(disulfide):
        torem = np.zeros(mol.numAtoms, dtype=bool)
        for d1, d2 in disulfide:
            a1 = d1.selectAtoms(mol, indexes=False)
            a2 = d2.selectAtoms(mol, indexes=False)
            d1.resname = "CYX"
            d2.resname = "CYX"
            mol.resname[a1] = "CYX"
            mol.resname[a2] = "CYX"
            torem |= (a1 & (mol.name == "HG")) | (a2 & (mol.name == "HG"))
        mol.remove(torem, _logger=False)

    return disulfide, custombonds, cyclic


def _add_caps(mol, caps):
    """Insert ACE / NME cap residues into *mol* in-place.

    Re-implemented from ``amber._add_caps`` **without** the tleap-specific
    terminal-hydrogen stripping that follows cap insertion in the original.
    OpenMM infers correct terminal chemistry from its residue templates.
    """
    from moleculekit.molecule import UniqueResidueID
    from htmd.home import home

    if caps is None or len(caps) == 0:
        return

    capdir = os.path.join(home(shareDir=True), "builder", "caps")
    _remove_atoms = {
        "ACE": ["H1", "H2", "H3", "HT1", "HT2", "HT3", "H"],
        "NME": ["OXT", "OT1", "O", "HXT"],
        "NHE": ["OXT", "OT1", "O", "HXT"],
    }

    uqres_caps = []
    first_val = caps[next(iter(caps))]
    if isinstance(first_val, (list, tuple)):
        for segid, pair in caps.items():
            segidx = np.where(mol.segid == segid)[0][[0, -1]]
            for i, cap in enumerate(pair):
                if cap is not None and cap != "none":
                    uqres_caps.append(
                        (UniqueResidueID.fromMolecule(mol, idx=segidx[i]), cap)
                    )
    else:
        for sel, cap in caps.items():
            uqres_caps.append((UniqueResidueID.fromMolecule(mol, sel), cap))

    for uqres, cap in uqres_caps:
        mask = uqres.selectAtoms(mol, indexes=False)
        if uqres.resname == cap:
            continue

        capmol = Molecule(os.path.join(capdir, f"{cap}.pdb"))
        align_names = capmol.name[capmol.beta == 1]
        mol_idx = np.where(np.isin(mol.name, align_names) & mask)[0]
        cap_idx = [
            np.where((capmol.name == n) & (capmol.resname == "XXX"))[0][0]
            for n in mol.name[mol_idx]
        ]
        if len(cap_idx) != len(mol_idx):
            raise RuntimeError("Could not find all required backbone atoms!")

        capmol.align(cap_idx, refmol=mol, refsel=mol_idx)
        capmol.resname[capmol.resname == "XXX"] = uqres.resname
        capmol.segid[:] = uqres.segid
        capmol.chain[:] = uqres.chain
        capmol.resid[:] += uqres.resid
        capmol.remove(cap_idx, _logger=False)

        mol.remove(np.isin(mol.name, _remove_atoms[cap]) & mask, _logger=False)

        res_idx = uqres.selectAtoms(mol, indexes=True)
        insert_pos = res_idx[0] if cap == "ACE" else res_idx[-1] + 1
        mol.insert(capmol, insert_pos)


def _defaultProteinCaps(mol):
    segs = np.unique(mol.get("segid", sel="protein"))
    caps = {}
    for s in segs:
        nres = len(np.unique(mol.resid[mol.segid == s]))
        if nres < 10:
            logger.warning(
                f"Segment {s} has fewer than 10 residues — not capped by "
                "default.  Use the caps argument to override."
            )
            continue
        caps[s] = ["ACE", "NME"]
    return caps


def _fix_water_naming(mol):
    """Rename water residues / atoms to match OpenMM AMBER templates."""
    water = mol.atomselect("water")
    if not np.any(water):
        return
    for old in ("TIP3", "WAT", "SOL", "T3P", "TP3"):
        mol.resname[water & (mol.resname == old)] = "HOH"
    mol.name[water & (mol.name == "OH2")] = "O"
    mol.name[water & (mol.name == "OW")] = "O"


_AMBER_TO_OPENMM_RNA = {
    "RG5": "G5",
    "RA5": "A5",
    "RC5": "C5",
    "RU5": "U5",
    "RG3": "G3",
    "RA3": "A3",
    "RC3": "C3",
    "RU3": "U3",
    "RG": "G",
    "RA": "A",
    "RC": "C",
    "RU": "U",
}

_AMBER_TO_OPENMM_DNA = {
    "DG5": "DG5",
    "DA5": "DA5",
    "DC5": "DC5",
    "DT5": "DT5",
    "DG3": "DG3",
    "DA3": "DA3",
    "DC3": "DC3",
    "DT3": "DT3",
    "DG": "DG",
    "DA": "DA",
    "DC": "DC",
    "DT": "DT",
}

_AMBER_TO_OPENMM_NUCLEIC_ATOMS = {
    "H2'1": "H2'",
    "H5'1": "H5'",
    "H5'2": "H5''",
    "HO'2": "HO2'",
    "HO'3": "HO3'",
    "H5T": "HO5'",
    "H3T": "HO3'",
}


def _fix_nucleic_naming(mol):
    """Convert AMBER-style nucleic acid residue/atom names to OpenMM naming."""
    resname_map = {**_AMBER_TO_OPENMM_RNA, **_AMBER_TO_OPENMM_DNA}
    nucleic_mask = np.zeros(mol.numAtoms, dtype=bool)
    for amber_name, omm_name in resname_map.items():
        mask = mol.resname == amber_name
        if np.any(mask):
            nucleic_mask |= mask
            mol.resname[mask] = omm_name

    if not np.any(nucleic_mask):
        return

    for amber_aname, omm_aname in _AMBER_TO_OPENMM_NUCLEIC_ATOMS.items():
        mask = nucleic_mask & (mol.name == amber_aname)
        if np.any(mask):
            mol.name[mask] = omm_aname


# ====================================================================
# OpenMM conversion & force-field setup
# ====================================================================


def _mol_to_openmm(mol, outdir, extra_xml=None):
    """Write *mol* to PDB and read back with OpenMM.

    Before reading the PDB we temporarily register the residue bond
    definitions from *extra_xml* (the OpenMM force-field XML written by
    :func:`htmd.builder.nonstandard.parameterizeFromSpecs`) into
    :class:`openmm.app.Topology`'s ``_standardBonds`` dict. PDBFile's
    ``createStandardBonds`` then auto-inserts the NCAA / ligand intra-
    residue bonds AND the peptide bonds to their canonical neighbours,
    so the topology arrives at the ForceField matcher with the right
    bond graph. The registration is scoped to the PDBFile read - we
    unregister immediately after so subsequent builds in the same
    Python process don't see leftover bond defs (the dict is
    class-level / process-wide).

    Returns ``(topology, positions)`` ready for ``ForceField.createSystem``.
    """
    import openmm.app as app

    _register_amber_variant_bond_defs()

    pdb_path = os.path.join(outdir, "input.pdb")
    mol.write(pdb_path, writebonds=False)
    with _temporary_residue_bond_defs(extra_xml):
        pdb = app.PDBFile(pdb_path)

    topology = pdb.topology
    _add_missing_bonds(mol, topology)

    return topology, pdb.positions


@contextlib.contextmanager
def _temporary_residue_bond_defs(extra_xml):
    """Register the ``<Residue>``/``<Bond>`` entries from each OpenMM
    force-field XML in *extra_xml* into ``Topology._standardBonds`` for
    the duration of the ``with`` block.

    The force-field XML format used by :func:`parameterizeFromSpecs`
    spells bonds as ``<Bond atomName1="..." atomName2="..."/>`` inside
    each ``<Residue>``. ``Topology._standardBonds`` expects
    ``(from, to)`` tuples (the same shape ``Topology.loadBondDefinitions``
    populates). We translate one to the other here.

    On exit we restore the snapshot of every entry we touched so other
    builds in the same Python process see a pristine dict.
    """
    from openmm.app import Topology
    import xml.etree.ElementTree as ET

    if extra_xml is None:
        yield
        return
    paths = [extra_xml] if isinstance(extra_xml, str) else list(extra_xml)

    snapshot = {}
    try:
        for path in paths:
            if not path or not os.path.isfile(path):
                continue
            root = ET.parse(path).getroot()
            for res in root.iter("Residue"):
                name = res.attrib.get("name")
                if not name:
                    continue
                if name not in snapshot:
                    snapshot[name] = Topology._standardBonds.get(name)
                bonds = []
                for bond in res.findall("Bond"):
                    a = bond.attrib.get("atomName1") or bond.attrib.get("from")
                    b = bond.attrib.get("atomName2") or bond.attrib.get("to")
                    if a and b:
                        bonds.append((a, b))
                # The force-field XML's ``<ExternalBond atomName="N"/>``
                # tag declares that the backbone N has a bond outside
                # the residue - i.e. the peptide bond to the previous
                # residue's C. Encode it as ``(-C, N)`` so PDBFile's
                # createStandardBonds inserts the peptide N-C bond at
                # the residue boundary. We only emit the pull-from-
                # previous direction (matching OpenMM's own
                # ``residues.xml`` convention); the push-to-next side
                # is created by the next residue's own ``(-C, N)``
                # entry. Emitting both would duplicate the same bond.
                # ``ExternalBond C`` and other ExternalBond atoms
                # (sidechain crosslink anchors, ligand junctions) are
                # handled via ``custombonds`` since their partner is
                # not known a-priori.
                atom_names = {a.attrib.get("name") for a in res.findall("Atom")}
                for ext in res.findall("ExternalBond"):
                    aname = ext.attrib.get("atomName")
                    if aname == "N" and "N" in atom_names:
                        bonds.append(("-C", "N"))
                Topology._standardBonds[name] = bonds
        yield
    finally:
        for name, prev in snapshot.items():
            if prev is None:
                Topology._standardBonds.pop(name, None)
            else:
                Topology._standardBonds[name] = prev


_AMBER_VARIANT_BONDS_REGISTERED = False


def _register_amber_variant_bond_defs():
    """Teach :class:`openmm.app.Topology` about the AMBER ff14SB
    protonation-state variants of canonical residues so PDBFile can
    auto-infer their intra-residue and peptide N-C bonds.

    OpenMM's built-in ``residues.xml`` (consumed by
    ``Topology.createStandardBonds``) only knows the standard PDB
    amino-acid names (CYS, HIS, LYS, ASP, GLU, TYR, ARG). The AMBER
    variant resnames moleculekit emits (CYM, CYX, HID, HIE, HIP, LYN,
    ASH, GLH, TYM, AR0) are absent. Without these definitions PDBFile
    can't insert the peptide bonds at the boundary of a variant
    residue, so OpenMM's ``ForceField`` template matcher rejects the
    NEIGHBOURING canonical residue with ``missing 1 C atom``.

    Each variant is the base residue's bond set minus the H atoms that
    don't exist for that protonation state (e.g. ``HG`` removed for
    CYM/CYX, ``HZ3`` for LYN, ``HE2`` for HID, ``HD1`` for HIE). The
    backbone and side-chain heavy-atom bonds carry over unchanged.

    Idempotent - the first call registers; subsequent calls are no-ops.
    """
    global _AMBER_VARIANT_BONDS_REGISTERED
    if _AMBER_VARIANT_BONDS_REGISTERED:
        return
    from openmm.app import Topology

    # Force the one-shot load of OpenMM's built-in ``residues.xml`` so
    # ``Topology._standardBonds`` is populated with the canonical
    # residues we base our variants on. ``createStandardBonds`` is the
    # only public path that flips ``_hasLoadedStandardBonds`` after
    # loading the bundled file; calling it on an empty Topology is a
    # cheap no-op apart from the load.
    Topology().createStandardBonds()

    # ``_standardBonds[name]`` is a list of (from_atom, to_atom) tuples.
    # ``-X`` / ``+X`` denote the previous / next residue's atom (i.e.
    # the peptide N-C bonds at the residue boundary).
    base_to_variants = {
        "CYS": {
            "CYM": {"HG"},        # deprotonated thiolate
            "CYX": {"HG"},        # disulfide-bonded
        },
        "HIS": {
            "HID": {"HE2"},       # delta-protonated
            "HIE": {"HD1"},       # epsilon-protonated
            "HIP": set(),         # doubly protonated
        },
        "LYS": {
            "LYN": {"HZ3"},       # neutral amine
        },
        "ASP": {
            "ASH": set(),         # protonated carboxyl (HD2 already in base?)
        },
        "GLU": {
            "GLH": set(),         # protonated carboxyl
        },
        "TYR": {
            "TYM": {"HH"},        # deprotonated phenolate
        },
        "ARG": {
            "AR0": {"HE"},        # neutral arginine (one variant)
        },
    }

    for base, variants in base_to_variants.items():
        base_bonds = Topology._standardBonds.get(base)
        if base_bonds is None:
            continue
        for variant, dropped_atoms in variants.items():
            if variant in Topology._standardBonds:
                continue
            kept = [
                (a, b) for a, b in base_bonds
                if a.lstrip("-+") not in dropped_atoms
                and b.lstrip("-+") not in dropped_atoms
            ]
            Topology._standardBonds[variant] = kept

    _AMBER_VARIANT_BONDS_REGISTERED = True


def _add_missing_bonds(mol, topology):
    """Add the cross-residue covalent bonds that PDBFile / Topology's
    standard-bond inference cannot reconstruct on its own.

    PDBFile already covers everything else:

    - Intra-residue bonds for any residue whose name is in
      ``Topology._standardBonds`` (standard amino acids and nucleotides
      out of the box; AMBER protonation variants CYM / CYX / HID / HIE
      / HIP / LYN / ASH / GLH / TYM / AR0 are registered here via
      :func:`_register_amber_variant_bond_defs` so PDBFile knows them
      too).
    - Peptide N-C bonds between adjacent residues whose names are in
      that same set.

    What remains is cross-residue covalent connectivity: disulfide
    SG-SG, isopeptide CD-NZ, scaffold thioether SG-C{ligand}, etc.
    These are exactly what ``amber.build`` materializes via
    ``custombonds``. We mirror that minimal set here by walking
    ``mol.bonds`` and keeping only SG-SG pairs (the disulfide case),
    matching ``amber.build``'s auto-detection. Spec-driven custom
    bonds for other crosslinks are added later via the OpenMM
    ``Modeller.addBond`` calls in ``build()`` from ``custombonds``.

    Bonds to single-atom ion residues (Zn, Mg, Ca, ...) are skipped:
    AMBER/ff14SB models metal coordination as non-bonded point
    charges and ``amber.build`` doesn't add explicit Metal-S /
    Metal-N bonds either.
    """
    if not hasattr(mol, "bonds") or len(mol.bonds) == 0:
        return

    atoms_list = list(topology.atoms())
    n_atoms = len(atoms_list)
    existing = {
        frozenset((b[0].index, b[1].index)) for b in topology.bonds()
    }
    for row in mol.bonds:
        i, j = int(row[0]), int(row[1])
        if i >= n_atoms or j >= n_atoms:
            continue
        ai, aj = atoms_list[i], atoms_list[j]
        # Only auto-add SG-SG disulfides (the one bond category
        # PDBFile won't infer between two standard CYS-like residues
        # but that's an unconditional part of the system topology).
        # All other cross-residue crosslinks come in through
        # ``custombonds`` in ``build()``.
        if str(ai.name) != "SG" or str(aj.name) != "SG":
            continue
        if ai.residue.index == aj.residue.index:
            continue
        key = frozenset((i, j))
        if key in existing:
            continue
        topology.addBond(ai, aj)
        existing.add(key)


def _setup_forcefield(ff, extra_xml, small_molecule_ff, molecules):
    """Create an ``openmm.app.ForceField`` with optional template generators.

    Parameters
    ----------
    ff : list of str
        Base force-field XML file names.
    extra_xml : str or list of str or None
        Additional XML file paths for custom residues.
    small_molecule_ff : str or None
        Name passed to GAFFTemplateGenerator / SMIRNOFFTemplateGenerator.
    molecules : list or None
        OpenFF ``Molecule`` objects or SDF file paths.
    """
    import openmm.app as app

    ff_args = list(ff)
    if extra_xml:
        if isinstance(extra_xml, str):
            extra_xml = [extra_xml]
        ff_args.extend(extra_xml)

    forcefield = app.ForceField(*ff_args)

    if small_molecule_ff and molecules:
        if not isinstance(molecules, (list, tuple)):
            molecules = [molecules]

        off_mols = []
        for entry in molecules:
            if isinstance(entry, str):
                from openff.toolkit import Molecule as OFFMolecule

                off_mols.append(OFFMolecule.from_file(entry))
            else:
                off_mols.append(entry)

        if "gaff" in small_molecule_ff.lower():
            from openmmforcefields.generators import GAFFTemplateGenerator

            gen = GAFFTemplateGenerator(
                molecules=off_mols, forcefield=small_molecule_ff
            )
        else:
            from openmmforcefields.generators import SMIRNOFFTemplateGenerator

            gen = SMIRNOFFTemplateGenerator(
                molecules=off_mols, forcefield=small_molecule_ff
            )
        forcefield.registerTemplateGenerator(gen.generator)

    return forcefield


# ====================================================================
# Topology bond patching (topology-independent)
# ====================================================================


def _extract_bond_specs(mol, cyclic, disulfide, custombonds):
    """Encode special bonds as ``(chain, resid_str, insertion, atom_name)`` tuples.

    This representation is independent of atom indices so it survives
    ``Modeller.addHydrogens()`` which may insert new atoms.
    """
    specs = {"cyclic": [], "disulfide": [], "custom": []}

    def _s(val):
        """Ensure plain Python str (numpy strings break == with OpenMM)."""
        return str(val) if val is not None else ""

    for seg, first_resid, last_resid in cyclic:
        seg_mask = mol.segid == seg
        chain = _s(mol.chain[seg_mask][0])
        first_ins = _s(mol.insertion[seg_mask & (mol.resid == first_resid)][0])
        last_ins = _s(mol.insertion[seg_mask & (mol.resid == last_resid)][0])
        specs["cyclic"].append(
            (
                chain,
                str(int(first_resid)),
                first_ins,
                "N",
                chain,
                str(int(last_resid)),
                last_ins,
                "C",
            )
        )

    for d1, d2 in disulfide:
        specs["disulfide"].append(
            (
                _s(d1.chain),
                str(int(d1.resid)),
                _s(d1.insertion),
                "SG",
                _s(d2.chain),
                str(int(d2.resid)),
                _s(d2.insertion),
                "SG",
            )
        )

    if custombonds:
        for a1, a2 in custombonds:
            i1 = a1.selectAtom(mol)
            i2 = a2.selectAtom(mol)
            specs["custom"].append(
                (
                    _s(mol.chain[i1]),
                    str(int(mol.resid[i1])),
                    _s(mol.insertion[i1]),
                    _s(mol.name[i1]),
                    _s(mol.chain[i2]),
                    str(int(mol.resid[i2])),
                    _s(mol.insertion[i2]),
                    _s(mol.name[i2]),
                )
            )

    return specs


def _find_atom_in_topology(topology, chain_id, resid_str, insertion, atom_name):
    """Locate a single atom in *topology* by chain / residue / name.

    OpenMM uses ``' '`` for empty insertion codes while moleculekit uses
    ``''``, so both sides are stripped before comparison.
    """
    insertion = insertion.strip()
    for chain in topology.chains():
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain.residues():
            if residue.id == resid_str and residue.insertionCode.strip() == insertion:
                for atom in residue.atoms():
                    if atom.name == atom_name:
                        return atom
    return None


def _apply_bond_specs(topology, specs):
    """Add bonds encoded by :func:`_extract_bond_specs` to *topology*.

    Skips bonds that already exist (e.g. disulfide SG-SG bonds that
    OpenMM creates automatically from CYX templates).
    """
    existing = {frozenset((b[0].index, b[1].index)) for b in topology.bonds()}

    for kind in ("cyclic", "disulfide", "custom"):
        for spec in specs[kind]:
            c1, r1, i1, n1, c2, r2, i2, n2 = spec
            a1 = _find_atom_in_topology(topology, c1, r1, i1, n1)
            a2 = _find_atom_in_topology(topology, c2, r2, i2, n2)
            if a1 is None or a2 is None:
                logger.warning(f"Could not locate atoms for {kind} bond: {spec}")
                continue
            key = frozenset((a1.index, a2.index))
            if key in existing:
                logger.debug(
                    f"{kind} bond already present: "
                    f"{a1.residue.name}:{a1.name} — "
                    f"{a2.residue.name}:{a2.name}"
                )
                continue
            topology.addBond(a1, a2)
            existing.add(key)
            if kind == "cyclic":
                logger.info(
                    f"Added cyclic N-C bond: "
                    f"{a1.residue.name}:{a1.name} — "
                    f"{a2.residue.name}:{a2.name}"
                )


# ====================================================================
# Solvation helpers
# ====================================================================


_OPENMM_ION_NAMES = {
    "NA": "Na+",
    "SODIUM": "Na+",
    "NA+": "Na+",
    "SOD": "Na+",
    "K": "K+",
    "POTASSIUM": "K+",
    "K+": "K+",
    "POT": "K+",
    "CS": "Cs+",
    "CESIUM": "Cs+",
    "CS+": "Cs+",
    "CES": "Cs+",
    "CL": "Cl-",
    "CHLORIDE": "Cl-",
    "CL-": "Cl-",
    "CLA": "Cl-",
    "MG": "Mg2+",
    "MAGNESIUM": "Mg2+",
    "MG2+": "Mg2+",
    "CA": "Ca2+",
    "CALCIUM": "Ca2+",
    "CA2+": "Ca2+",
    "CAL": "Ca2+",
    "ZN": "Zn2+",
    "ZINC": "Zn2+",
    "ZN2+": "Zn2+",
}

_ION_CHARGES = {
    "Na+": 1,
    "K+": 1,
    "Cs+": 1,
    "Cl-": -1,
    "Mg2+": 2,
    "Ca2+": 2,
    "Zn2+": 2,
}

_ION_ELEMENTS = {
    "Na+": "Na",
    "K+": "K",
    "Cs+": "Cs",
    "Cl-": "Cl",
    "Mg2+": "Mg",
    "Ca2+": "Ca",
    "Zn2+": "Zn",
}


def _resolve_ion_name(name):
    """Map any common ion alias to the OpenMM residue-template name."""
    if name is None:
        return None
    key = name.upper().replace(" ", "")
    result = _OPENMM_ION_NAMES.get(key)
    if result is None:
        raise ValueError(
            f"Unknown ion type '{name}'.  Accepted values: "
            + ", ".join(sorted(set(_OPENMM_ION_NAMES.values())))
        )
    return result


def _build_solvent_kwargs(
    water_model,
    padding,
    boxsize,
    ionize,
    saltconc,
    saltanion,
    saltcation,
    unit,
):
    """Construct keyword dict for ``Modeller.addSolvent``."""
    kw = {"model": water_model}

    if boxsize is not None:
        if isinstance(boxsize, (int, float)):
            boxsize = [boxsize] * 3
        from openmm import Vec3

        kw["boxSize"] = (
            Vec3(boxsize[0] / 10.0, boxsize[1] / 10.0, boxsize[2] / 10.0)
            * unit.nanometers
        )
    else:
        kw["padding"] = (padding / 10.0) * unit.nanometers

    if ionize:
        kw["neutralize"] = True
        if saltconc > 0:
            kw["ionicStrength"] = saltconc * unit.molar
        if saltcation:
            kw["positiveIon"] = _resolve_ion_name(saltcation)
        if saltanion:
            kw["negativeIon"] = _resolve_ion_name(saltanion)
    else:
        kw["neutralize"] = False

    return kw


def _ensure_box_vectors(topology, positions, unit, mol=None):
    """Set periodic box vectors on *topology* if missing.

    If *mol* has non-zero ``box`` (Angstroms) those dimensions are used.
    Otherwise the bounding box of *positions* (nanometers) is used.
    """
    if topology.getPeriodicBoxVectors() is not None:
        return
    from openmm import Vec3

    bx = by = bz = 0.0
    if mol is not None and mol.box is not None and mol.box.shape[0] >= 3:
        bx = float(mol.box[0, 0]) / 10.0
        by = float(mol.box[1, 0]) / 10.0
        bz = float(mol.box[2, 0]) / 10.0

    if bx <= 0 or by <= 0 or bz <= 0:
        coords = np.array([[p.x, p.y, p.z] for p in positions])
        bx, by, bz = (coords.max(axis=0) - coords.min(axis=0)).tolist()

    if bx <= 0 or by <= 0 or bz <= 0:
        raise ValueError("Cannot determine periodic box dimensions.")

    vectors = [
        Vec3(bx, 0, 0),
        Vec3(0, by, 0),
        Vec3(0, 0, bz),
    ] * unit.nanometers
    topology.setPeriodicBoxVectors(vectors)


# ====================================================================
# Pre-solvated ionisation
# ====================================================================


def _get_total_charge(system):
    """Sum particle charges from the NonbondedForce of *system*."""
    import openmm
    import openmm.unit as unit

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            total = 0.0
            for i in range(force.getNumParticles()):
                q, _, _ = force.getParticleParameters(i)
                total += q.value_in_unit(unit.elementary_charge)
            return total
    raise RuntimeError("No NonbondedForce found in system")


def _ionize_presolvated(mol, system, saltconc, saltanion, saltcation):
    """Replace water molecules with ions for a pre-solvated system.

    Uses farthest-point sampling to spread ions throughout the solvent.
    Returns the modified *mol* (in-place).
    """
    from scipy.spatial import cKDTree

    total_q = round(_get_total_charge(system))
    anion = _resolve_ion_name(saltanion) if saltanion else "Cl-"
    cation = _resolve_ion_name(saltcation) if saltcation else "Na+"
    cat_q = _ION_CHARGES[cation]
    an_q = _ION_CHARGES[anion]

    # --- Neutralisation ---
    ncation = nanion = 0
    if total_q < 0:
        ncation = int(abs(total_q) // cat_q)
        leftover = abs(total_q) - ncation * cat_q
        if leftover:
            nanion = int(leftover // abs(an_q))
    elif total_q > 0:
        nanion = int(total_q // abs(an_q))
        leftover = total_q - nanion * abs(an_q)
        if leftover:
            ncation = int(leftover // cat_q)

    # --- Additional salt ---
    if saltconc and saltconc > 0:
        nwater = int(np.sum(mol.atomselect("water and noh")))
        nsalt = max(0, int(round(0.0187 * saltconc * nwater)))
        ncation += nsalt
        nanion += nsalt

    ntotal = ncation + nanion
    if ntotal == 0:
        logger.info("System already neutral and no salt requested.")
        return mol

    logger.info(f"Adding {ncation} {cation} and {nanion} {anion}")

    # --- Pick water oxygens far from solute (>5 A) ---
    water_o = mol.atomselect("water and noh")
    solute = ~mol.atomselect("water")
    water_o_idx = np.where(water_o)[0]
    solute_coords = mol.coords[solute, :, 0]
    water_coords = mol.coords[water_o_idx, :, 0]

    tree = cKDTree(solute_coords)
    dists, _ = tree.query(water_coords)

    far_mask = dists > 5.0
    far_idx = water_o_idx[far_mask]
    far_coords = water_coords[far_mask]

    if len(far_idx) < ntotal:
        raise BuildError(
            f"Only {len(far_idx)} water molecules are >5 A from solute, "
            f"but {ntotal} ions are needed."
        )

    # Farthest-point sampling
    picked = [int(np.argmax(dists[far_mask]))]
    for _ in range(ntotal - 1):
        d_to_picked = cKDTree(far_coords[picked]).query(far_coords)[0]
        d_to_picked[picked] = -np.inf
        picked.append(int(np.argmax(d_to_picked)))

    picked_o = far_idx[np.array(picked)]

    # Assign cation/anion labels (shuffled for spatial mixing)
    labels = [cation] * ncation + [anion] * nanion
    rng = np.random.default_rng(42)
    rng.shuffle(labels)

    # Replace waters with ions
    uqresid = sequenceID((mol.resid, mol.insertion, mol.segid))
    torem = np.zeros(mol.numAtoms, dtype=bool)
    for o_idx, ion in zip(picked_o, labels):
        res_mask = uqresid == uqresid[o_idx]
        h_of_water = res_mask & (mol.element == "H")
        torem |= h_of_water
        mol.resname[o_idx] = _ION_ELEMENTS[ion]
        mol.name[o_idx] = ion
        mol.element[o_idx] = _ION_ELEMENTS[ion]
        mol.segid[o_idx] = "I"
        # Remove any extra atoms in the residue that are not the ion
        extra = res_mask & (np.arange(mol.numAtoms) != o_idx) & ~h_of_water
        torem |= extra

    mol.remove(torem, _logger=False)
    return mol


# ====================================================================
# Export & I/O
# ====================================================================


def _export_amber(topology, system, positions, outdir, prefix, forcefield=None):
    """Write AMBER prmtop / inpcrd via ParmEd.

    A second System is built **without** constraints so that every bond
    carries an explicit force constant (AMBER's prmtop format requires
    that). ParmEd's ``openmm.load_topology`` then turns the OpenMM
    System into an ``AmberParm`` and writes prmtop / inpcrd.

    Note: ``parmed.openmm.load_topology`` synthesises AMBER atom types
    by bucketing on ``(element, sigma, epsilon)``, which collapses
    ff14SB's ``CT`` / ``CX`` into a single prmtop type. The resulting
    prmtop is still well-formed for OpenMM (per-instance bond / angle
    / dihedral params), so loading it back via ``parmed.amber.AmberParm``
    and computing energy through OpenMM is fine. Tools that consolidate
    parameters per atom-type pair (e.g. ``ffevaluation.loadParameters``)
    cannot consume this output. ``openff.interchange.Interchange``
    preserves the original atom-type names but its ``to_prmtop`` is
    prohibitively slow (pint unit-string re-parsing per parameter) on
    protein-sized systems, so it's not used here.
    """
    import openmm.app as app
    import openmm.unit as unit

    prmtop = os.path.join(outdir, f"{prefix}.prmtop")
    inpcrd = os.path.join(outdir, f"{prefix}.inpcrd")

    try:
        import parmed

        export_system = system
        if forcefield is not None:
            nb_method = app.NoCutoff
            if topology.getPeriodicBoxVectors() is not None:
                nb_method = app.PME
            export_system = forcefield.createSystem(
                topology,
                nonbondedMethod=nb_method,
                nonbondedCutoff=1.0 * unit.nanometers,
                constraints=None,
                rigidWater=False,
            )

        # ``condense_atom_types=False`` keeps a unique atom type per
        # atom (rather than bucketing on (element, sigma, epsilon)).
        # The default bucketing collapses ff14SB's ``CT`` / ``CX``
        # into a single prmtop type, which produces multiple distinct
        # bond / angle / dihedral parameter values for the same
        # collapsed type pair - breaking ``ffevaluation.loadParameters``
        # which assumes one parameter set per type pair. With
        # condense disabled each bond / angle / dihedral instance has
        # a unique type tuple, so the prmtop satisfies that invariant.
        struct = parmed.openmm.load_topology(
            topology,
            export_system,
            xyz=positions,
            condense_atom_types=False,
        )
        struct.save(prmtop, overwrite=True)
        struct.save(inpcrd, overwrite=True)
        logger.info(f"Wrote {prmtop} and {inpcrd} via ParmEd")
    except Exception as exc2:
        logger.error(f"ParmEd AMBER export failed: {exc2}")
        raise


def _read_built_molecule(outdir, prefix):
    """Read back the built system as a ``Molecule``.

    Reads topology from prmtop and coordinates from the PDB file
    (ParmEd writes restart-format inpcrd that moleculekit cannot
    always parse, so PDB is the reliable coordinate source).
    """
    prmtop = os.path.join(outdir, f"{prefix}.prmtop")
    pdb = os.path.join(outdir, f"{prefix}.pdb")

    if os.path.exists(prmtop) and os.path.getsize(prmtop) > 0:
        try:
            molbuilt = Molecule(prmtop, validateElements=False)
            if os.path.exists(pdb):
                pdb_mol = Molecule(pdb)
                molbuilt.coords = pdb_mol.coords.copy()
                molbuilt.box = pdb_mol.box.copy()
            return molbuilt
        except Exception:
            logger.warning("Could not read prmtop — falling back to PDB.")

    return Molecule(pdb)
