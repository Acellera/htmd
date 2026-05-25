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
# OpenFF Interchange parameterisation (free ligands)
# ====================================================================


def parameterizeLigandsOpenFF(
    mol,
    ligand_ff="openff_unconstrained-2.3.0.offxml",
    charge_method="gasteiger",
    resnames=None,
):
    """Parameterise free ligand residues via OpenFF Interchange.

    Slices each ligand residue out of *mol*, assigns partial charges
    (default RDKit Gasteiger, which honours the net formal charge unlike
    antechamber's ``-c gas``), applies the chosen SMIRNOFF force field
    via :func:`openff.interchange.Interchange.from_smirnoff`, and returns
    one :class:`Interchange` per resname. The returned objects can be
    combined with other Interchanges via ``ic.combine(other)`` and
    exported to OpenMM via ``ic.to_openmm()``.

    Mirrors the OpenFF protein-ligand tutorial pattern::

        ligand_ic = Interchange.from_smirnoff(
            force_field=ForceField("openff_unconstrained-2.3.0.offxml"),
            topology=[ligand_offmol],
        )

    Parameters
    ----------
    mol : :class:`moleculekit.molecule.Molecule`
        Molecule containing one or more free ligand residues. Each ligand
        must already have explicit hydrogens and explicit integer bond
        orders, e.g. via
        :meth:`moleculekit.molecule.Molecule.templateResidueFromSmiles`.
    ligand_ff : str, optional
        Name of an OpenFF force field offxml. ``"openff_unconstrained-*"``
        variants are recommended when the final system will be
        energy-minimised; constrained variants are only correct for
        rigid-water MD.
    charge_method : str or None, optional
        ``"gasteiger"`` (default) - compute Gasteiger PEOE charges with
        RDKit before SMIRNOFF application and pass them through
        ``charge_from_molecules`` so SMIRNOFF does not overwrite them.
        Honours the net formal charge.
        ``None`` - let SMIRNOFF assign its own charges (e.g. AM1-BCC for
        Sage). Requires sqm from AmberTools at runtime.
    resnames : list[str] or None, optional
        Subset of resnames to parameterise. ``None`` parameterises every
        unique resname in *mol*.

    Returns
    -------
    dict[str, openff.interchange.Interchange]
        One Interchange per resname.
    """
    from openff.toolkit import ForceField
    from openff.interchange import Interchange
    from htmd.builder._ambertools import _assign_rdkit_gasteiger_charges

    if resnames is None:
        resnames = [str(r) for r in np.unique(mol.resname) if str(r)]
    if not resnames:
        raise ValueError("No ligand resnames found in mol")

    ff = ForceField(ligand_ff)
    out = {}
    for resname in resnames:
        sub = mol.copy()
        sub.filter(f'resname "{resname}"', _logger=False)
        if sub.numAtoms == 0:
            raise ValueError(f"resname {resname!r} matched no atoms in mol")

        if charge_method == "gasteiger":
            _assign_rdkit_gasteiger_charges(sub)
        elif charge_method is None:
            pass
        else:
            raise ValueError(
                f"Unsupported charge_method {charge_method!r}. "
                f"Supported: 'gasteiger', None."
            )

        off_mol = sub.toOpenFFMolecule(sanitize=True, assignStereo=True)

        kwargs = {"force_field": ff, "topology": [off_mol]}
        if charge_method is not None:
            kwargs["charge_from_molecules"] = [off_mol]
        ic = Interchange.from_smirnoff(**kwargs)
        out[resname] = ic
        logger.info(
            f"Parameterised ligand {resname!r} with {ligand_ff} ({charge_method=}): "
            f"{off_mol.n_atoms} atoms, {off_mol.n_bonds} bonds"
        )
    return out


# OpenMM XML expects kJ/mol, nm, and radians. openff.units handles the
# conversion natively via Quantity.m_as(target_unit); call sites just
# ask for the right target. No manual conversion factors needed.


def _emit_openmm_xml_from_interchange(
    ic,
    resname,
    external_bond_atom_names,
    xml_path,
    *,
    class_prefix=None,
):
    """Free-residue (single-template) wrapper around the multi-residue
    XML emitter. Builds a one-residue ``residue_assignment`` and forwards.

    See :func:`_emit_openmm_xml_from_cluster_interchange` for the general
    multi-residue API and full parameter notes.
    """
    n_atoms = ic.topology.molecule(0).n_atoms
    residue_assignment = [(resname, list(range(n_atoms)))]
    cap_atoms = set()
    external_bonds = {resname: list(external_bond_atom_names)}
    return _emit_openmm_xml_from_cluster_interchange(
        ic=ic,
        residue_assignment=residue_assignment,
        cap_atoms=cap_atoms,
        external_bonds=external_bonds,
        xml_path=xml_path,
        class_prefix=class_prefix,
    )


def _emit_openmm_xml_from_cluster_interchange(
    ic,
    residue_assignment,
    cap_atoms,
    external_bonds,
    xml_path,
    *,
    class_prefix=None,
    atom_names=None,
    class_overrides=None,
    type_overrides=None,
):
    """Multi-residue OpenMM ``<ForceField>`` XML emitter.

    Splits one cluster :class:`Interchange` (covering several residues
    plus optional ACE/NME cap atoms) into a single XML document
    containing one ``<Residue>`` per cluster residue and global force
    sections that reference synthesised per-atom classes.

    Cap atoms (``cap_atoms``) are excluded from the emitted residue
    templates AND from all force-section entries. Cap atoms typically
    parameterise the cluster boundary in context (ACE/NME on chain
    residues, neighbour-residue stubs for cross-linked clusters); after
    parameterisation they are no longer relevant in the final peptide,
    so their parameters get dropped at emission time. The peptide-bond
    parameters at the canonical-NCAA junction in the final peptide come
    from ff14SB's protein XML (loaded alongside our output).

    Parameters
    ----------
    ic : :class:`openff.interchange.Interchange`
        Interchange resolved via :func:`Interchange.from_smirnoff` on a
        cluster topology (full residues + caps).
    residue_assignment : list[tuple[str, list[int]]]
        Per-cluster-residue ``(resname, atom_indices)`` pairs. Atom
        indices index into the Interchange topology. The order
        determines the order of ``<Residue>`` elements in the output.
    cap_atoms : set[int]
        Topology-atom indices for ACE/NME-style cap atoms - dropped at
        emit time.
    external_bonds : dict[str, list[str]]
        Per-resname list of atom NAMES that need ``<ExternalBond>``
        markers. For chain residues this is typically ``["N", "C"]``
        (or one of them on a terminal residue).
    xml_path : str
        File path to write the XML to.
    class_prefix : str or None
        Class-name prefix for synthesised atom classes. Defaults to
        ``OFF_<first_resname>`` for backward compatibility with the
        single-residue path. Pass an explicit prefix to disambiguate
        across multiple clusters that share a resname.

    Returns
    -------
    str
        ``xml_path`` (the file written).
    """
    import xml.etree.ElementTree as ET
    from openff.units import unit

    if class_prefix is None:
        first_resname = residue_assignment[0][0]
        class_prefix = f"OFF_{first_resname}"

    off_mol = ic.topology.molecule(0)
    n_atoms = off_mol.n_atoms

    if class_overrides is None:
        class_overrides = {}
    if type_overrides is None:
        # Default: type == class (true for our synthesised OFF_ classes
        # and for backward compatibility callers that pass only
        # class_overrides). For ff14SB overrides callers should pass
        # type_overrides separately: ff14SB's "protein-CT" is a TYPE
        # whose CLASS is "CT" - force entries match by class, residue
        # templates carry the type, so the two strings differ.
        type_overrides = dict(class_overrides)

    # Per-atom class names (used in force-section entries). Default to
    # synthesised unique ``OFF_<prefix>_<i>``; atoms in ``class_overrides``
    # get the ff14SB class (e.g. "CT" for backbone CA).
    classes = [
        class_overrides.get(i, f"{class_prefix}_{i}") for i in range(n_atoms)
    ]
    # Per-atom type names (used in <Residue><Atom type=...> attributes).
    # Default to the same as the class; backbone atoms get the ff14SB
    # type (e.g. "protein-CT" for CA).
    types_per_atom = [
        type_overrides.get(i, f"{class_prefix}_{i}") for i in range(n_atoms)
    ]
    is_overridden = [i in class_overrides for i in range(n_atoms)]

    if atom_names is None:
        atom_names = [a.name for a in off_mol.atoms]
    else:
        assert len(atom_names) == n_atoms, (
            f"atom_names length {len(atom_names)} != topology n_atoms "
            f"{n_atoms}"
        )
    elements = [_element_symbol_from_atomic_number(a.atomic_number) for a in off_mol.atoms]
    masses = [float(a.mass.m_as(unit.dalton)) for a in off_mol.atoms]
    is_kept = [i not in cap_atoms for i in range(n_atoms)]

    # Per-atom partial charges from the Electrostatics collection.
    charges = [0.0] * n_atoms
    for key, q in ic["Electrostatics"].charges.items():
        idx = key.atom_indices[0]
        if idx < n_atoms:
            charges[idx] = float(q.m_as(unit.elementary_charge))

    # Build the root <ForceField> element.
    root = ET.Element("ForceField")
    resnames = ",".join(r for r, _ in residue_assignment)
    root.append(ET.Comment(
        f" Auto-generated from OpenFF Interchange. Residues: {resnames}. "
        f"Cap atoms dropped: {len(cap_atoms)}. Do not edit by hand. "
    ))

    # <AtomTypes>: one entry per kept atom. Atoms with a class override
    # are skipped because their class is declared by a different XML
    # (typically ff14SB's protein.ff14SB.xml).
    atom_types_el = ET.SubElement(root, "AtomTypes")
    for i, (cls, elem, mass) in enumerate(zip(classes, elements, masses)):
        if not is_kept[i] or is_overridden[i]:
            continue
        ET.SubElement(
            atom_types_el,
            "Type",
            {"name": cls, "class": cls, "element": elem, "mass": f"{mass:.6f}"},
        )

    # <Residues>: one <Residue> per cluster residue. Bonds within the
    # residue go in the residue template; bonds crossing the residue
    # boundary become implicit (recognized via <ExternalBond> markers).
    residues_el = ET.SubElement(root, "Residues")
    atom_residue = {}  # atom_idx -> position in residue_assignment
    for ri, (_resname, atom_indices) in enumerate(residue_assignment):
        for ai in atom_indices:
            atom_residue[ai] = ri

    for ri, (resname, atom_indices) in enumerate(residue_assignment):
        residue_el = ET.SubElement(residues_el, "Residue", {"name": resname})
        for ai in atom_indices:
            ET.SubElement(
                residue_el,
                "Atom",
                {
                    "name": atom_names[ai],
                    "type": types_per_atom[ai],
                    "charge": f"{charges[ai]:.6f}",
                },
            )
        # Intra-residue bonds.
        for bond in off_mol.bonds:
            a, b = bond.atom1_index, bond.atom2_index
            if a in cap_atoms or b in cap_atoms:
                continue
            if atom_residue.get(a) == ri and atom_residue.get(b) == ri:
                ET.SubElement(
                    residue_el,
                    "Bond",
                    {
                        "atomName1": atom_names[a],
                        "atomName2": atom_names[b],
                    },
                )
        # External bonds: from the per-resname mapping the caller supplied
        # (chain N/C, scaffolded-peptide cross bonds, etc).
        for ext_name in external_bonds.get(resname, []):
            ET.SubElement(residue_el, "ExternalBond", {"atomName": ext_name})

    # Global force sections - dropping any entry that touches a cap atom,
    # and any entry where ALL atoms have a class override (ff14SB owns
    # those parameters, our XML should not redeclare them).
    overridden_atoms = {i for i in range(n_atoms) if is_overridden[i]}
    _write_bond_force(
        root, ic, classes,
        drop_atoms=cap_atoms, all_overridden_drop=overridden_atoms,
    )
    _write_angle_force(
        root, ic, classes,
        drop_atoms=cap_atoms, all_overridden_drop=overridden_atoms,
    )
    _write_torsion_force(
        root, ic, classes,
        drop_atoms=cap_atoms, all_overridden_drop=overridden_atoms,
    )
    _write_nonbonded_force(
        root, ic, classes,
        drop_atoms=cap_atoms | overridden_atoms,
    )

    ET.indent(root, space="  ")
    tree = ET.ElementTree(root)
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)
    n_kept = sum(is_kept)
    logger.info(
        f"Wrote OpenMM XML fragment to {xml_path} - "
        f"{len(residue_assignment)} residue(s) [{resnames}], "
        f"{n_kept}/{n_atoms} atoms kept ({len(cap_atoms)} cap atoms dropped)"
    )
    return xml_path


def _element_symbol_from_atomic_number(atomic_number):
    """Atomic number -> element symbol via RDKit's PeriodicTable.

    RDKit is already a hard moleculekit dep so this stays self-contained.
    """
    from rdkit.Chem import GetPeriodicTable

    return GetPeriodicTable().GetElementSymbol(int(atomic_number))


def _write_bond_force(root, ic, classes, drop_atoms=frozenset(),
                      all_overridden_drop=frozenset()):
    import xml.etree.ElementTree as ET
    from openff.units import unit

    bonds_el = ET.SubElement(root, "HarmonicBondForce")
    for top_key, pot_key in ic["Bonds"].key_map.items():
        i, j = top_key.atom_indices
        if i in drop_atoms or j in drop_atoms:
            continue
        # Drop if ALL atoms are class-overridden (ff14SB owns the param).
        if i in all_overridden_drop and j in all_overridden_drop:
            continue
        pot = ic["Bonds"].potentials[pot_key]
        length_nm = float(pot.parameters["length"].m_as(unit.nanometer))
        k_kj = float(pot.parameters["k"].m_as(
            unit.kilojoule_per_mole / unit.nanometer**2,
        ))
        ET.SubElement(
            bonds_el,
            "Bond",
            {
                "class1": classes[i],
                "class2": classes[j],
                "length": f"{length_nm:.8f}",
                "k": f"{k_kj:.6f}",
            },
        )


def _write_angle_force(root, ic, classes, drop_atoms=frozenset(),
                       all_overridden_drop=frozenset()):
    import xml.etree.ElementTree as ET
    from openff.units import unit

    angles_el = ET.SubElement(root, "HarmonicAngleForce")
    for top_key, pot_key in ic["Angles"].key_map.items():
        i, j, k = top_key.atom_indices
        if i in drop_atoms or j in drop_atoms or k in drop_atoms:
            continue
        if (
            i in all_overridden_drop
            and j in all_overridden_drop
            and k in all_overridden_drop
        ):
            continue
        pot = ic["Angles"].potentials[pot_key]
        angle_rad = float(pot.parameters["angle"].m_as(unit.radian))
        k_kj = float(pot.parameters["k"].m_as(
            unit.kilojoule_per_mole / unit.radian**2,
        ))
        ET.SubElement(
            angles_el,
            "Angle",
            {
                "class1": classes[i],
                "class2": classes[j],
                "class3": classes[k],
                "angle": f"{angle_rad:.8f}",
                "k": f"{k_kj:.6f}",
            },
        )


def _write_torsion_force(root, ic, classes, drop_atoms=frozenset(),
                         all_overridden_drop=frozenset()):
    """Write proper + improper torsions. SMIRNOFF spreads multi-term
    dihedrals across `mult`-indexed entries with the same atom tuple; we
    merge them back into one ``<Proper>``/``<Improper>`` element with
    periodicity1/k1/phase1, periodicity2/..., etc. so the OpenMM XML
    schema is compact.
    """
    import xml.etree.ElementTree as ET
    from collections import defaultdict
    from openff.units import unit

    # OpenMM expands SMIRNOFF impropers in trefoil symmetry (3 terms per
    # <Improper> entry, one per permutation of the non-central atoms);
    # the AMBER ordering applies each <Improper> as a single term. Since
    # the Interchange's stored impropers come from SMIRNOFF's SMARTS-based
    # assignment, we want SMIRNOFF ordering to match Interchange.to_openmm.
    pt_el = ET.SubElement(root, "PeriodicTorsionForce", {"ordering": "smirnoff"})

    def _emit(collection_name, xml_tag):
        # Group multi-term entries by atom tuple.
        grouped = defaultdict(list)
        for top_key, pot_key in ic[collection_name].key_map.items():
            atoms = tuple(top_key.atom_indices)
            if any(a in drop_atoms for a in atoms):
                continue
            if all(a in all_overridden_drop for a in atoms):
                continue
            pot = ic[collection_name].potentials[pot_key]
            periodicity = int(pot.parameters["periodicity"].m_as(unit.dimensionless))
            phase_rad = float(pot.parameters["phase"].m_as(unit.radian))
            # SMIRNOFF stores improper k undivided: the per-term
            # contribution is k/idivf (idivf=3 by SMIRNOFF convention so
            # the three trefoil-permuted entries sum to k). OpenMM applies
            # each <Improper> as a single term, so we bake the divide in
            # here. Propers carry idivf=1 (no division).
            idivf = 1.0
            if "idivf" in pot.parameters:
                idivf = float(pot.parameters["idivf"].m_as(unit.dimensionless))
            k_kj = float(
                pot.parameters["k"].m_as(unit.kilojoule_per_mole) / idivf
            )
            grouped[atoms].append((periodicity, phase_rad, k_kj))

        for atoms, terms in grouped.items():
            attrs = {
                "class1": classes[atoms[0]],
                "class2": classes[atoms[1]],
                "class3": classes[atoms[2]],
                "class4": classes[atoms[3]],
            }
            for idx, (per, phase, k) in enumerate(terms, start=1):
                attrs[f"periodicity{idx}"] = str(per)
                attrs[f"phase{idx}"] = f"{phase:.8f}"
                attrs[f"k{idx}"] = f"{k:.6f}"
            ET.SubElement(pt_el, xml_tag, attrs)

    _emit("ProperTorsions", "Proper")
    _emit("ImproperTorsions", "Improper")


def _write_nonbonded_force(root, ic, classes, drop_atoms=frozenset()):
    """Emit ``<NonbondedForce>`` with per-class sigma/epsilon. Charges
    live on the ``<Residue><Atom>`` entries so we don't repeat them here.
    """
    import xml.etree.ElementTree as ET
    from openff.units import unit

    # SMIRNOFF 1-4 scaling. OpenFF Sage uses 1.0 / (1/0.833333) for
    # coulomb (which is 0.833...) and 0.5 for LJ. Read them from the
    # vdW collection's metadata when available, fall back to AMBER
    # defaults otherwise.
    coulomb14_scale = 0.833333
    lj14_scale = 0.5
    if hasattr(ic.collections.get("vdW", None), "scale_14"):
        lj14_scale = float(ic["vdW"].scale_14)
    if hasattr(ic.collections.get("Electrostatics", None), "scale_14"):
        coulomb14_scale = float(ic["Electrostatics"].scale_14)

    nb_el = ET.SubElement(
        root,
        "NonbondedForce",
        {
            "coulomb14scale": f"{coulomb14_scale:.6f}",
            "lj14scale": f"{lj14_scale:.6f}",
        },
    )
    # UseAttributeFromResidue: tell OpenMM that per-atom charges come
    # from the <Residue> template, not from the NonbondedForce entries.
    ET.SubElement(nb_el, "UseAttributeFromResidue", {"name": "charge"})

    n_atoms = ic.topology.n_atoms
    sigmas = [0.0] * n_atoms
    epsilons = [0.0] * n_atoms
    for top_key, pot_key in ic["vdW"].key_map.items():
        (idx,) = top_key.atom_indices
        pot = ic["vdW"].potentials[pot_key]
        sigmas[idx] = float(pot.parameters["sigma"].m_as(unit.nanometer))
        epsilons[idx] = float(
            pot.parameters["epsilon"].m_as(unit.kilojoule_per_mole)
        )
    for i, (cls, sigma, epsilon) in enumerate(zip(classes, sigmas, epsilons)):
        if i in drop_atoms:
            continue
        ET.SubElement(
            nb_el,
            "Atom",
            {
                "type": cls,
                "sigma": f"{sigma:.8f}",
                "epsilon": f"{epsilon:.6f}",
            },
        )


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

    # CHARMM PSFs encode the rigid-water SHAKE constraint as an explicit
    # H1-H2 bond. OpenMM's HOH template only has O-H bonds, so emitting
    # the H-H bond via CONECT trips "1 H-H bond too many" at template
    # matching. Drop them here.
    if len(mol.bonds) == 0:
        return
    is_h = mol.element == "H"
    i, j = mol.bonds[:, 0].astype(int), mol.bonds[:, 1].astype(int)
    hh_water = water[i] & water[j] & is_h[i] & is_h[j]
    if np.any(hh_water):
        mol.bonds = mol.bonds[~hh_water]
        if mol.bondtype is not None and len(mol.bondtype) == len(hh_water):
            mol.bondtype = mol.bondtype[~hh_water]


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
    mol.write(pdb_path, writebonds=True)
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
            "CYM": {"HG"},  # deprotonated thiolate
            "CYX": {"HG"},  # disulfide-bonded
        },
        "HIS": {
            "HID": {"HE2"},  # delta-protonated
            "HIE": {"HD1"},  # epsilon-protonated
            "HIP": set(),  # doubly protonated
        },
        "LYS": {
            # Neutral lysine. ff14SB / systemPrepare emit LYN with
            # HZ2 + HZ3 on NZ (dropping HZ1), not HZ1 + HZ2 (dropping
            # HZ3). Keep both Hs that ff14SB defines.
            "LYN": {"HZ1"},
        },
        "ASP": {
            "ASH": set(),  # protonated carboxyl (HD2 already in base?)
        },
        "GLU": {
            "GLH": set(),  # protonated carboxyl
        },
        "TYR": {
            "TYM": {"HH"},  # deprotonated phenolate
        },
        "ARG": {
            "AR0": {"HE"},  # neutral arginine (one variant)
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
                (a, b)
                for a, b in base_bonds
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
    existing = {frozenset((b[0].index, b[1].index)) for b in topology.bonds()}
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


def _register_extra_xml_template_matcher(forcefield, extra_xml):
    """Register a template matcher that prefers the template whose NAME
    matches the topology residue's name, whenever such a same-named
    template exists and actually matches the residue's atoms.

    Why this is safe and necessary:

    The default OpenMM matcher selects templates by atom-element +
    bond-graph signature, which can return multiple matches when two
    templates are topologically indistinguishable. Two failure modes
    we hit in this codebase:

    1. *Same atom set, different bonded partners on an ExternalBond.*
       Our XX1 (Cys-renamed, S-C scaffold thio-ether) vs ff14SB's NCYX
       (Cys-disulfide-to-another-CYX): identical atom set, identical
       ExternalBond count on SG, but chemically different. The default
       matcher reports both and errors.
    2. *Stereoisomer NCAAs.* DAL (D-Ala) and ff14SB ALA (L-Ala) have
       identical atom topology because chirality isn't in the bond
       graph. A topology residue named "ALA" gets matched by both.

    Resolving by NAME is the right call because:
    - The user / pipeline emitted the topology residue with a specific
      NAME for a reason. If they wanted ff14SB's ALA they named the
      residue ALA; if they wanted our DAL they named it DAL.
    - We only return the same-named template when it ACTUALLY matches
      the residue's atoms - otherwise we fall through to default logic.

    The matcher does NOT replace e.g. NCYX residues with our XX1
    parameters: it looks up ``_templates[residue.name]``, so a real
    NCYX residue still gets ff14SB's NCYX (assuming no other XML
    redefined it).
    """
    import xml.etree.ElementTree as ET
    from openmm.app.internal import compiled

    if isinstance(extra_xml, str):
        extra_xml = [extra_xml]

    # We only need to ARM the matcher when extra_xml exists - if there's
    # no extra XML the default OpenMM logic already produces the right
    # answer everywhere. Collecting the resname set is just a sanity
    # check that something was loaded.
    extra_resnames = set()
    for path in extra_xml:
        try:
            root = ET.parse(path).getroot()
        except (ET.ParseError, FileNotFoundError):
            continue
        for res in root.findall("./Residues/Residue"):
            name = res.attrib.get("name")
            if name:
                extra_resnames.add(name)

    if not extra_resnames:
        return

    def _matcher(forcefield, residue, bonded_to_atom, ignore_external, ignore_extra):
        template = forcefield._templates.get(residue.name)
        if template is None:
            return None
        # Verify the same-named template actually matches before
        # returning it - registerTemplateMatcher raises if a returned
        # template doesn't match. This guards against unusual cases
        # where a residue is renamed (e.g. HIS->HID at apply time) but
        # the topology still carries the input name.
        match = compiled.matchResidueToTemplate(
            residue, template, bonded_to_atom, ignore_external, ignore_extra,
        )
        if match is None:
            return None
        return template

    forcefield.registerTemplateMatcher(_matcher)


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

    # Templates from ``extra_xml`` take priority over the base ff14SB / GAFF2
    # templates when the topology residue NAME matches the extra_xml
    # template name. This disambiguates structurally-equivalent residues
    # whose chemistry differs only in *what's on the other end of an
    # ExternalBond* - the classic case being CYS-renamed XX1 (S-C scaffold
    # bond) vs ff14SB's NCYX (S-S disulfide bond): identical atom set, same
    # ExternalBond count on SG, but different bond partners. OpenMM's
    # default signature-based matcher can't tell them apart and raises
    # "Multiple non-identical matching templates"; the matcher below short-
    # circuits the choice by name.
    if extra_xml:
        _register_extra_xml_template_matcher(forcefield, extra_xml)

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

        struct = parmed.openmm.load_topology(topology, export_system, xyz=positions)
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
