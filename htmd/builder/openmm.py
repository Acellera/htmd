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
import sys
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


def defaultFf() -> list:
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


def _phosaa_ffxml_namespaced(outdir):
    """Write a copy of ``openmmforcefields``' ``phosaa14SB.xml`` made compatible
    with OpenMM's bundled ``amber14/protein.ff14SB.xml`` base, and return its path.

    Two adjustments. (1) amber14 names its protein atom TYPES ``protein-*`` while
    phosaa references them bare (``N``, ``CX``, ...), so the residue templates'
    standard-type references are prefixed with ``protein-`` (phosaa's own
    phosphate types OP/OQ/.../CG and the phosphorus ``P`` are left bare).
    (2) the protein-only base defines no phosphorus, so the ``P`` atom type and
    its nonbonded term are copied in from the monolithic ``amber/ff14SB.xml``.
    Classes stay bare in both, so all bonded/nonbonded terms match unchanged.
    Generated from the installed package each build, so it tracks its version."""
    import re
    import glob
    import openmmforcefields

    ffdir = os.path.join(os.path.dirname(openmmforcefields.__file__), "ffxml")

    def _find(fn):
        hits = [
            m
            for m in glob.glob(os.path.join(ffdir, "**", fn), recursive=True)
            if os.path.basename(m) == fn
        ]
        if not hits:
            raise RuntimeError(
                f"openmmforcefields is missing {fn}; cannot build phosphorylated "
                "residues (PTR/SEP/TPO) with the OpenMM builder."
            )
        return hits[0]

    phos = open(_find("phosaa14SB.xml")).read()
    mono = open(_find("ff14SB.xml")).read()

    keep_bare = set(re.findall(r'<Type\b[^>]*\bname="([^"]+)"', phos)) | {"P"}

    def _fix_residue(match):
        def _fix_atom(atom):
            whole, typ = atom.group(0), atom.group(1)
            if typ in keep_bare:
                return whole
            return whole.replace(f'type="{typ}"', f'type="protein-{typ}"')

        return re.sub(r'<Atom\b[^>]*\btype="([^"]+)"[^>]*/>', _fix_atom, match.group(0))

    phos = re.sub(r"<Residue\b.*?</Residue>", _fix_residue, phos, flags=re.S)

    p_type = next(
        line for line in re.findall(r"<Type\b[^>]*>", mono) if re.search(r'\bname="P"', line)
    )
    p_nonbonded = next(
        line
        for line in re.findall(
            r"<Atom\b[^>]*>",
            re.search(r"<NonbondedForce\b.*?</NonbondedForce>", mono, re.S).group(0),
        )
        if re.search(r'class="P"', line)
    )
    phos = phos.replace("</AtomTypes>", f"  {p_type}\n </AtomTypes>")
    phos = re.sub(r"(<NonbondedForce\b[^>]*>)", rf"\1\n  {p_nonbonded}", phos, count=1)

    out_path = os.path.join(outdir, "phosaa14SB_ff14SBnamespaced.xml")
    with open(out_path, "w") as fh:
        fh.write(phos)
    return out_path


def _maybe_add_phosaa(mol, outdir, extra_xml):
    """If ``mol`` carries a phosphorylated amino acid, generate the amber14-
    compatible phosaa XML into ``outdir`` and append it to ``extra_xml`` (so it
    reaches both the topology bond-definitions and the ForceField). Returns the
    (possibly-augmented) ``extra_xml``."""
    from htmd.builder.amber import _PHOSAA_RESNAMES

    if not np.any(np.isin(mol.resname, list(_PHOSAA_RESNAMES))):
        return extra_xml
    phos_xml = _phosaa_ffxml_namespaced(outdir)
    logger.info(
        "Phosphorylated residue(s) detected; auto-loading namespaced "
        "phosaa14SB for the OpenMM builder."
    )
    if not extra_xml:
        return [phos_xml]
    if isinstance(extra_xml, str):
        return [extra_xml, phos_xml]
    return list(extra_xml) + [phos_xml]


# AMBER modified-residue libraries that openmmforcefields never converted to
# OpenMM XML. Each is converted at build time (ParmEd) to a minimal
# amber14-compatible ffxml. ``params`` / ``lib`` are relative to AMBER's
# dat/leap/{parm,lib}; ``prefix`` is the base ff's type namespace (e.g. amber14
# names protein types ``protein-*``).
def _amber_modres_libs():
    from moleculekit.residues import MODIFIED_NUCLEIC_RESIDUE_NAMES

    return (
        {
            "params": ("parm10.dat", "frcmod.ff14SB", "frcmod.ff14SBmodAA"),
            "lib": "mod_amino.lib",
            "residues": frozenset({"ALY", "AZF", "CYF", "CNX", "MSE"}),
            "base_xml": "amber14/protein.ff14SB.xml",
            "prefix": "protein",
        },
        {
            "params": ("parm10.dat", "all_modrna08.frcmod"),
            "lib": "all_modrna08.lib",
            "residues": frozenset(MODIFIED_NUCLEIC_RESIDUE_NAMES),
            "base_xml": "amber14/RNA.OL3.xml",
            "prefix": "RNA",
        },
    )


def _amber14_base_classes(base_xml):
    import re
    import openmm.app as app

    txt = open(os.path.join(os.path.dirname(app.__file__), "data", base_xml)).read()
    return set(re.findall(r'<Type[^>]*\bclass="([^"]+)"', txt))


def _amber_modres_ffxml(present, outdir):
    """Convert every AMBER modified-residue library that provides a resname in
    ``present`` to a MINIMAL amber14-compatible OpenMM ffxml, and return the
    written paths.

    ParmEd converts the library (params + OFF residue templates); from that we
    keep only the target residue templates, the atom types they introduce that
    the base ff lacks (a "new" class), and the force terms touching a new class.
    Standard-type references in the residue templates are namespaced to the base
    convention (``protein-CT`` ...) so the modified residue's peptide junction to
    ordinary residues resolves; charge comes from the residue template
    (``UseAttributeFromResidue``), matching amber14."""
    import xml.etree.ElementTree as ET
    from parmed.amber import AmberParameterSet, AmberOFFLibrary
    from parmed.openmm import OpenMMParameterSet
    from htmd.builder.amber import defaultAmberHome, _defaultAmberSearchPaths

    present = set(present)
    home = defaultAmberHome()
    pdir = os.path.join(home, _defaultAmberSearchPaths["param"])
    ldir = os.path.join(home, _defaultAmberSearchPaths["lib"])

    paths = []
    for lib in _amber_modres_libs():
        want = sorted(present & lib["residues"])
        if not want:
            continue
        aps = AmberParameterSet(*[os.path.join(pdir, f) for f in lib["params"]])
        reslib = AmberOFFLibrary.parse(os.path.join(ldir, lib["lib"]))
        omm = OpenMMParameterSet.from_parameterset(aps)
        for r in want:
            omm.residues[r] = reslib[r]
        full = os.path.join(outdir, f"_full_{lib['lib']}.xml")
        omm.write(full, provenance=None)
        root = ET.parse(full).getroot()
        os.remove(full)

        base_cls = _amber14_base_classes(lib["base_xml"])
        prefix = lib["prefix"]
        type_def = {t.get("name"): t for t in root.findall(".//AtomTypes/Type")}
        res_elems = [
            res for res in root.findall(".//Residues/Residue") if res.get("name") in want
        ]
        used = {a.get("type") for res in res_elems for a in res.findall("Atom")}
        new_types = {t for t in used if type_def[t].get("class") not in base_cls}
        new_classes = {type_def[t].get("class") for t in new_types}

        ff = ET.Element("ForceField")
        atblock = ET.SubElement(ff, "AtomTypes")
        for t in sorted(new_types):
            atblock.append(type_def[t])
        resblock = ET.SubElement(ff, "Residues")
        for res in res_elems:
            for a in res.findall("Atom"):
                if a.get("type") not in new_types:
                    a.set("type", f"{prefix}-{a.get('type')}")
            resblock.append(res)

        def _touches_new(el):
            vals = set()
            for k, v in el.attrib.items():
                if k.startswith("class"):
                    vals.add(v)
                elif k == "type":
                    vals.add(type_def[v].get("class") if v in type_def else v)
            return bool(vals & new_classes)

        for tag in (
            "HarmonicBondForce", "HarmonicAngleForce",
            "PeriodicTorsionForce", "NonbondedForce",
        ):
            src = root.find(f".//{tag}")
            if src is None:
                continue
            keep = [c for c in list(src) if _touches_new(c)]
            if tag == "NonbondedForce" and not keep:
                continue
            fe = ET.SubElement(ff, tag, src.attrib)
            if tag == "NonbondedForce":
                # amber14 supplies charge via the residue template.
                ET.SubElement(fe, "UseAttributeFromResidue", {"name": "charge"})
            for c in keep:
                if tag == "NonbondedForce":
                    c.attrib.pop("charge", None)
                fe.append(c)

        out_path = os.path.join(outdir, f"{lib['lib'].split('.')[0]}_ff14SBnamespaced.xml")
        ET.ElementTree(ff).write(out_path)
        paths.append(out_path)
    return paths


def _maybe_add_amber_modres(mol, outdir, extra_xml):
    """Auto-convert + load any AMBER-library modified residue present in ``mol``
    (MSE / ALY / ... from ff14SB_modAA). Returns the augmented ``extra_xml``."""
    paths = _amber_modres_ffxml({str(r) for r in np.unique(mol.resname)}, outdir)
    if not paths:
        return extra_xml
    logger.info(
        "AMBER-library modified residue(s) detected; auto-loading converted "
        f"OpenMM XML: {[os.path.basename(p) for p in paths]}."
    )
    if not extra_xml:
        return list(paths)
    if isinstance(extra_xml, str):
        return [extra_xml, *paths]
    return list(extra_xml) + list(paths)


def _ffptm_prepi_residues():
    """Resnames htmd ships as ff-PTM prepi units (post-translational
    modifications AMBER parameterizes through ``ff-ptm/<RES>.prepi`` +
    ``<RES>.frcmod``, auto-loaded by :func:`amber._detect_cofactors_ncaa_ptm`).

    These are amino-acid modifications built on the standard ff14SB atom types
    (the prepi carries custom charges + a handful of extra bonded terms in the
    frcmod; no new atom types), so their OpenMM counterpart is a residue
    template over ``protein-*`` types. Phospho residues were removed from
    ff-ptm in favour of phosaa, so there is no overlap with
    :data:`amber._PHOSAA_RESNAMES`."""
    from htmd.builder.amber import htmdAmberHome

    ffptm = os.path.join(htmdAmberHome(), "ff-ptm")
    if not os.path.isdir(ffptm):
        return frozenset()
    return frozenset(
        os.path.splitext(f)[0]
        for f in os.listdir(ffptm)
        if f.endswith(".prepi")
    )


def _parse_prepi(path):
    """Parse an AMBER amino-acid ``prepi`` unit tleap-free.

    Returns ``(atoms, bonds)`` where ``atoms`` is a list of
    ``(name, amber_type, charge)`` for the real (non-dummy) atoms and ``bonds``
    is a list of ``(name1, name2)`` intra-residue covalent bonds. Bonds are
    read from the internal-coordinate tree (each atom's ``NA`` parent column)
    plus the explicit ``LOOP`` ring-closure section; bonds to the three leading
    ``DUMM`` atoms (the backbone-N chain entry) are dropped as those become
    inter-residue external bonds.

    Parameters
    ----------
    path : str
        Path to the ``.prepi`` file.

    Returns
    -------
    atoms : list of tuple
        ``(name, amber_type, charge)`` per real atom, in prepi order.
    bonds : list of tuple
        ``(name1, name2)`` covalent bonds.
    """
    idx_name = {}
    idx_isdummy = {}
    atoms = []
    bonds = []
    in_loop = False
    for ln in open(path).read().splitlines():
        f = ln.split()
        if not f:
            continue
        if f[0] == "LOOP":
            in_loop = True
            continue
        if f[0] in ("IMPROPER", "DONE", "STOP", "CHARGE"):
            in_loop = False
            continue
        if in_loop:
            if len(f) == 2:
                bonds.append((f[0], f[1]))
            continue
        # atom row: "I NAME TYPE TOPO NA NB NC R THETA PHI CHARGE"
        if len(f) >= 11 and f[0].isdigit() and f[3] in (
            "M", "S", "B", "E", "3", "4", "5", "6",
        ):
            i = int(f[0])
            name, atype = f[1], f[2]
            idx_name[i] = name
            idx_isdummy[i] = atype == "DU"
            if atype == "DU":
                continue
            atoms.append((name, atype, float(f[10])))
            na = int(f[4])
            if na in idx_name and not idx_isdummy[na]:
                bonds.append((name, idx_name[na]))
    return atoms, bonds


def _ffptm_prepi_ffxml(present, outdir):
    """Convert each ff-PTM prepi residue in ``present`` to a minimal
    amber14-compatible OpenMM ffxml plus a hydrogen-definition XML.

    The residue template is read straight from the prepi (atoms with their
    ff-PTM charges over ``protein-*`` types, the tree/LOOP bonds, and backbone
    N/C external bonds). The frcmod's extra bonded terms (a handful of
    bond/angle/dihedral/improper terms over standard classes that ff14SB lacks)
    are converted to OpenMM units by ParmEd - the frcmod carries no MASS
    section, so it defines no new atom types and ParmEd emits exactly those
    terms. A companion hydrogen-definition XML (H name -> parent, from the
    prepi) lets ``Modeller.addHydrogens`` rebuild the residue's hydrogens in
    the template's naming, sidestepping the prepi-vs-PDBv3 hydrogen-name
    mismatch. Returns a list of ``(resname, ffxml_path, hdef_path)``."""
    import xml.etree.ElementTree as ET
    from parmed.amber import AmberParameterSet
    from parmed.openmm import OpenMMParameterSet
    from htmd.builder.amber import htmdAmberHome

    ffptm = os.path.join(htmdAmberHome(), "ff-ptm")
    want = sorted(set(present) & _ffptm_prepi_residues())
    results = []
    for res in want:
        atoms, bonds = _parse_prepi(os.path.join(ffptm, f"{res}.prepi"))
        names = {a[0] for a in atoms}

        ff = ET.Element("ForceField")
        resblock = ET.SubElement(ff, "Residues")
        rel = ET.SubElement(resblock, "Residue", {"name": res})
        for name, atype, charge in atoms:
            ET.SubElement(rel, "Atom", {
                "name": name, "type": f"protein-{atype}", "charge": f"{charge:.6f}",
            })
        for n1, n2 in bonds:
            ET.SubElement(rel, "Bond", {"atomName1": n1, "atomName2": n2})
        for bb in ("N", "C"):
            if bb in names:
                ET.SubElement(rel, "ExternalBond", {"atomName": bb})

        # Extra bonded terms from the frcmod, converted to OpenMM units.
        frcmod = os.path.join(ffptm, f"{res}.frcmod")
        aps = AmberParameterSet(frcmod)
        if any((aps.bond_types, aps.angle_types, aps.dihedral_types,
                getattr(aps, "improper_periodic_types", {}))):
            tmp = os.path.join(outdir, f"_frcmod_{res}.xml")
            OpenMMParameterSet.from_parameterset(aps).write(tmp, provenance=None)
            froot = ET.parse(tmp).getroot()
            os.remove(tmp)
            for tag in ("HarmonicBondForce", "HarmonicAngleForce",
                        "PeriodicTorsionForce"):
                src = froot.find(f".//{tag}")
                if src is None or len(src) == 0:
                    continue
                dst = ET.SubElement(ff, tag, src.attrib)
                for c in src:
                    # AMBER wildcard "X" -> OpenMM empty class.
                    for k, v in list(c.attrib.items()):
                        if k.startswith("class") and v == "X":
                            c.set(k, "")
                    dst.append(c)

        ffxml_path = os.path.join(outdir, f"{res}_ff14SBnamespaced.xml")
        ET.ElementTree(ff).write(ffxml_path)

        # Hydrogen definition: parent heavy atom of each H, from the prepi bonds.
        parent = {}
        for a, b in bonds:
            if a.startswith("H") ^ b.startswith("H"):
                h, p = (a, b) if a.startswith("H") else (b, a)
                parent[h] = p
        hroot = ET.Element("Residues")
        hres = ET.SubElement(hroot, "Residue", {"name": res})
        for name, _, _ in atoms:
            if name.startswith("H") and name in parent:
                ET.SubElement(hres, "H", {"name": name, "parent": parent[name]})
        hdef_path = os.path.join(outdir, f"{res}_hydrogens.xml")
        ET.ElementTree(hroot).write(hdef_path)

        results.append((res, ffxml_path, hdef_path))
    return results


def _maybe_add_ffptm_prepi(mol, outdir, extra_xml):
    """Auto-convert + load any ff-PTM prepi residue present in ``mol`` (CGU,
    CSO, ...). Registers each residue's hydrogen definition and strips its
    input hydrogens **in place** so ``Modeller.addHydrogens`` rebuilds them in
    the template's naming (mirroring how amber.build strips + rebuilds these
    from the prepi via tleap). Returns the augmented ``extra_xml``."""
    import openmm.app as app

    results = _ffptm_prepi_ffxml({str(r) for r in np.unique(mol.resname)}, outdir)
    if not results:
        return extra_xml
    present = [r for r, _, _ in results]
    logger.info(
        f"ff-PTM modified residue(s) {', '.join(present)} detected; auto-loading "
        "converted OpenMM XML and rebuilding their hydrogens from the prepi."
    )
    paths = []
    for res, ffxml_path, hdef_path in results:
        app.Modeller.loadHydrogenDefinitions(hdef_path)
        mol.remove(
            (mol.resname == res) & (mol.element == "H"), _logger=False
        )
        paths.append(ffxml_path)
    if not extra_xml:
        return paths
    if isinstance(extra_xml, str):
        return [extra_xml, *paths]
    return list(extra_xml) + paths


def build(
    mol: Molecule,
    ff: list | None = None,
    extra_xml: str | list | None = None,
    small_molecule_ff: str | None = None,
    molecules: list | None = None,
    prefix: str = "structure",
    outdir: str = "./build",
    caps: dict | None = None,
    ionize: bool = True,
    saltconc: float = 0,
    saltanion: str | None = None,
    saltcation: str | None = None,
    disulfide: list | None = None,
    custombonds: list | None = None,
    solvate: bool = True,
    padding: float = 10.0,
    water_model: str = "tip3p",
    boxsize: float | list | None = None,
    gbsa: bool = False,
):
    """Build a system using OpenMM force fields and export to AMBER format.

    Parameters
    ----------
    mol : Molecule
        The input molecular system.
    ff : list of str, optional
        OpenMM XML force field file names (see ``openmm.app.ForceField``).
        If None, uses :func:`defaultFf`.
    extra_xml : str or list of str, optional
        Paths to additional OpenMM XML files for non-standard residues.
    small_molecule_ff : str, optional
        Small-molecule force field for the template generator, e.g.
        ``"gaff-2.2.20"`` or ``"openff-2.3.0"``.  Requires *molecules*.
    molecules : list, optional
        ``openff.toolkit.Molecule`` objects or paths to SDF files
        describing the small molecules present in the system.
    prefix : str
        Prefix for output files.
    outdir : str
        Output directory path.
    caps : dict, optional
        Terminal capping specification.  Accepts two formats:

        * ``{segid: [nterm_cap, cterm_cap]}`` - e.g. ``{"P": ["ACE","NME"]}``
        * ``{atomsel: cap}`` - e.g. ``{"chain A and resid 5": "ACE"}``

        If None, ACE/NME caps are added to every protein segment with >= 10 residues.
    ionize : bool
        Neutralise the system (and add salt if *saltconc* > 0).
    saltconc : float
        Salt concentration in molar, added on top of neutralisation.
    saltanion : str, optional
        Anion type.  Accepts ``"Cl-"``, ``"CL"``, ``"chloride"`` etc.
    saltcation : str, optional
        Cation type.  Accepts ``"Na+"``, ``"K+"``, ``"Cs+"`` and
        divalent ``"Mg2+"``, ``"Ca2+"``, ``"Zn2+"``.
    disulfide : list, optional
        Manual disulfide bonds as pairs of atom-selection strings.
        If None, automatic detection is performed.
    custombonds : list, optional
        Extra bonds as pairs of atom-selection strings.
    solvate : bool
        Add explicit water via ``Modeller.addSolvent()``.
    padding : float
        Box padding in Angstroms (used when *solvate* is True).
    water_model : str
        Water model name for ``Modeller.addSolvent()``.
    boxsize : float or list of float, optional
        Explicit box dimensions ``[x, y, z]`` in Angstroms.  Overrides
        *padding*.
    gbsa : bool
        Use GBSA implicit solvent (OBC2 model).

    Returns
    -------
    molbuilt : Molecule
        The fully built system.
    system : openmm.System
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

    # Resnames whose backbone N is a custombond (isopeptide / side-chain
    # crosslink) endpoint. createStandardBonds must not ALSO add a (-C, N)
    # backbone peptide bond for them (see _isopeptide_acceptor_resnames).
    skip_peptide_n = _isopeptide_acceptor_resnames(mol, custombond_ids)

    os.makedirs(outdir, exist_ok=True)

    # Phosphorylated residues (PTR/SEP/TPO) need AMBER's phosaa parameters; the
    # amber14 base has none. Generate an amber14-compatible phosaa XML and load
    # it alongside (both for topology bond defs and the ForceField).
    extra_xml = _maybe_add_phosaa(mol, outdir, extra_xml)
    # Same for AMBER-library modified residues openmmforcefields never converted
    # (MSE / ALY / ... from ff14SB_modAA): convert them on the fly.
    extra_xml = _maybe_add_amber_modres(mol, outdir, extra_xml)
    # ff-PTM prepi residues (CGU, CSO, ...): convert the prepi + frcmod straight
    # to an OpenMM template, register a hydrogen definition, and strip their
    # input Hs so addHydrogens rebuilds them in the template's naming.
    extra_xml = _maybe_add_ffptm_prepi(mol, outdir, extra_xml)

    topology, positions = _mol_to_openmm(
        mol, outdir, extra_xml=extra_xml, skip_peptide_n=skip_peptide_n
    )

    forcefield, smallmol_ffxml = _setup_forcefield(
        ff, extra_xml, small_molecule_ff, molecules
    )

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
            "Molecule already contains water - skipping solvation.  "
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
        topology, positions = _mol_to_openmm(
            mol, outdir, extra_xml=extra_xml, skip_peptide_n=skip_peptide_n
        )
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

    _export_amber(topology, system, positions, outdir, prefix, export_ff)

    # Topology + bonds from the prmtop; coordinates straight from the OpenMM
    # positions (full precision, no PDB round-trip); box from the topology.
    molbuilt = _read_built_molecule(outdir, prefix, topology, positions)
    detectCisPeptideBonds(molbuilt, respect_bonds=True)

    # ForceField-XML handoff for ACEMD: moleculekit writes the mmCIF (always)
    # and the PDB (when it fits the serial limit), plus the parameter set and
    # system.yaml. Replaces the old system.xml.
    _write_ff_handoff(molbuilt, outdir, prefix, ff, extra_xml, smallmol_ffxml)

    logger.info("Finished building.")
    return molbuilt, system


PDB_SERIAL_LIMIT = 99999  # PDB 5-digit atom-serial / CONECT cap


def _write_ff_handoff(molbuilt, outdir, prefix, ff, extra_xml, smallmol_ffxml):
    """Write the self-contained ForceField-XML handoff for ACEMD.

    Emits an *ordered* parameter set (stock FF names first, then copied
    ``extra_xml`` fragments, then the frozen small-molecule ffxml), an mmCIF
    topology (always; preserves all connectivity), a PDB (only when the system
    fits PDB's serial limit - above that CONECT records overflow so we ship
    only the mmCIF), and a partial ``system.yaml`` config. Returns the ordered
    ``parameters`` list.
    """
    import shutil
    import yaml

    parameters = [ff] if isinstance(ff, str) else list(ff)  # stock names, ordered first

    if extra_xml:
        extra_xml = [extra_xml] if isinstance(extra_xml, str) else list(extra_xml)
        for path in extra_xml:
            dst = os.path.join(outdir, os.path.basename(path))
            # An extra_xml generated straight into outdir (e.g. the auto-added
            # phosaa XML) is already in place - copying it onto itself errors.
            if os.path.abspath(path) != os.path.abspath(dst):
                shutil.copy(path, dst)
            parameters.append(os.path.basename(path))

    for i, xml in enumerate(smallmol_ffxml):
        fname = f"{prefix}.smallmol_{i}.xml"
        with open(os.path.join(outdir, fname), "w") as fh:
            fh.write(xml)
        parameters.append(fname)

    cif = f"{prefix}.cif"
    molbuilt.write(os.path.join(outdir, cif))

    # The PDB is a convenience carrier for other tools; moleculekit writes it
    # (its CONECT writer caps at the PDB serial limit). Above the limit the
    # bond records overflow, so we ship only the mmCIF (no such cap).
    if molbuilt.numAtoms <= PDB_SERIAL_LIMIT:
        molbuilt.write(os.path.join(outdir, f"{prefix}.pdb"))
    else:
        logger.info(
            f"System has {molbuilt.numAtoms} atoms (> {PDB_SERIAL_LIMIT}); "
            "PDB CONECT serials overflow - shipping only the mmCIF."
        )

    system_yaml = {
        "structure": cif,
        "coordinates": cif,
        "parameters": parameters,
    }
    # Periodic builds carry a box; vacuum / implicit-solvent builds do not.
    # Omit boxsize when there is no box, so downstream (e.g. setup_equilibration)
    # derives one from the coordinates instead of choking on an empty box.
    if molbuilt.box is not None and molbuilt.box.shape[1] > 0:
        box = [float(b) for b in molbuilt.box[:, molbuilt.frame]]
        if any(box):
            system_yaml["boxsize"] = box
    with open(os.path.join(outdir, "system.yaml"), "w") as fh:
        yaml.safe_dump(system_yaml, fh, sort_keys=False)
    logger.info(f"Wrote ForceField-XML handoff to {outdir}")
    return parameters


# ====================================================================
# OpenFF Interchange parameterisation (free ligands)
# ====================================================================


_MIN_INTERCHANGE_PYTHON = (3, 12)


def _require_openff_python():
    """Raise a clear error if the running Python is too old for
    openff-interchange. From version 0.5 onward openff-interchange uses
    PEP 695 ``type`` statement syntax (e.g. ``type Array = Any``) which
    is a parse error on Python < 3.12. Trying to import it on 3.10 or
    3.11 produces ``SyntaxError`` from inside the package - this guard
    raises a friendlier message before the import is attempted."""
    if sys.version_info < _MIN_INTERCHANGE_PYTHON:
        raise RuntimeError(
            f"openff-interchange requires Python "
            f">= {_MIN_INTERCHANGE_PYTHON[0]}.{_MIN_INTERCHANGE_PYTHON[1]} "
            "(uses PEP 695 type-statement syntax). You have Python "
            f"{sys.version_info.major}.{sys.version_info.minor}."
        )


# NAGL / RESP / Gasteiger charge helpers live in ``_charge_helpers`` so
# the GAFF backend (``_ambertools.py``) and the SMIRNOFF backend
# (this file) can both reach them without a circular import.
from htmd.builder._charge_helpers import (  # noqa: E402
    _assign_rdkit_gasteiger_charges,
    _assign_nagl_charges,
    _assign_resp_charges,
)


def parameterizeLigandsOpenFF(
    mol: Molecule,
    ligand_ff: str = "openff_unconstrained-2.3.0.offxml",
    charge_method: str | None = "nagl",
    resnames: list | None = None,
):
    """Parameterize free ligand residues via OpenFF Interchange.

    Slices each ligand residue out of *mol*, assigns partial charges,
    applies the chosen SMIRNOFF force field via
    :func:`openff.interchange.Interchange.from_smirnoff`, and returns
    one Interchange per resname. The returned objects can be
    combined with other Interchanges via ``ic.combine(other)`` and
    exported to OpenMM via ``ic.to_openmm()``.

    The default charge model is NAGL, which is the OpenFF-recommended
    fast surrogate for AM1-BCC (Sage 2.0-2.2) / AshGC (Sage 2.3.0).
    Using a different charge model logs a warning: Sage's vdW + torsion
    parameters were fit alongside a specific charge model, so other
    charge choices give a thermodynamically inconsistent potential.

    Mirrors the OpenFF protein-ligand tutorial pattern::

        ligand_ic = Interchange.from_smirnoff(
            force_field=ForceField("openff_unconstrained-2.3.0.offxml"),
            topology=[ligand_offmol],
        )

    Parameters
    ----------
    mol : Molecule
        Molecule containing one or more free ligand residues. Each ligand
        must already have explicit hydrogens and explicit integer bond
        orders, e.g. via
        :meth:`moleculekit.molecule.Molecule.templateResidueFromSmiles`.
    ligand_ff : str
        Name of an OpenFF force field offxml. The default
        ``"openff_unconstrained-2.3.0.offxml"`` is the right choice for
        ``openmm.build``: the bond / angle / torsion / vdW / charge
        parameters are byte-identical to the constrained variant
        ``"openff-2.3.0.offxml"`` (per the OpenFF README: "Each mainline
        force field is currently available in two forms - both with and
        without bond constraints to hydrogen"), and ``openmm.build``
        already applies ``constraints=app.HBonds`` at ``createSystem``
        time so the X-H bonds are constrained at the simulator level.
        Pick the constrained variant only if the emitted XML will be
        consumed by a downstream tool that runs ``createSystem`` without
        passing ``constraints``.
    charge_method : str or None, optional
        ``"nagl"`` (recommended) - AM1-BCC-equivalent partial charges
        from the OpenFF NAGL graph neural network. Sage was fit alongside
        this charge model. Needs PyTorch. Other choices log a warning.
        ``"gasteiger"`` - RDKit PEOE. ``"resp"`` / ``"resp-multiconf"``
        - RESP via parameterize + Psi4. ``None`` - let SMIRNOFF assign its
        own charges (e.g. AM1-BCC via ToolkitAM1BCCHandler).
    resnames : list, optional
        Subset of resnames to parameterize. If None, every unique resname
        in *mol* is parameterized.

    Returns
    -------
    dict
        One ``openff.interchange.Interchange`` per resname.
    """
    _require_openff_python()
    from openff.toolkit import ForceField
    from openff.interchange import Interchange
    from htmd.builder._ambertools import _assign_rdkit_gasteiger_charges

    # SMIRNOFF forcefields are fit alongside a specific charge model;
    # warn (don't error) if a different charge model is used since the
    # vdW + torsion parameters are entangled with the charge values.
    if charge_method != "nagl":
        logger.warning(
            f"parameterizeLigandsOpenFF was called with "
            f"charge_method={charge_method!r}; OpenFF Sage was fit "
            "alongside a specific charge model (AM1-BCC for Sage 2.0-"
            "2.2, AshGC for Sage 2.3.0) that NAGL reproduces. Using "
            "any other charge model against Sage's LJ and torsion "
            "parameters gives a thermodynamically inconsistent "
            "potential. Pass charge_method='nagl' for the recommended "
            "behaviour."
        )

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
        elif charge_method == "nagl":
            _assign_nagl_charges(sub)
        elif charge_method == "resp":
            _assign_resp_charges(sub, multi_conf=False)
        elif charge_method == "resp-multiconf":
            _assign_resp_charges(sub, multi_conf=True)
        elif charge_method is None:
            pass
        else:
            raise ValueError(
                f"Unsupported charge_method {charge_method!r}. "
                f"Supported: 'gasteiger', 'nagl', 'resp', "
                f"'resp-multiconf', None."
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

    # PDB stores the C-terminal amide cap as "NH2" but the OpenMM ff14SB
    # template (like AMBER) calls it "NHE"; rename so it is recognised.
    if np.any(mol.resname == "NH2"):
        mol.resname[mol.resname == "NH2"] = "NHE"

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
    from htmd.builder.amber import _PHOSAA_RESNAMES

    # Modified protein residues whose OpenMM template comes from a converted
    # AMBER library (phosaa or ff14SB_modAA) have no OpenMM hydrogen definition,
    # so addHydrogens can't re-add the backbone amide H that capping strips -
    # keep it for them (standard residues are unaffected).
    keep_amide_h = set(_PHOSAA_RESNAMES) | {
        r for lib in _amber_modres_libs() if lib["prefix"] == "protein" for r in lib["residues"]
    }

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

        remove_names = _remove_atoms[cap]
        if cap == "ACE" and uqres.resname in keep_amide_h:
            # Converted-library modified residues have no OpenMM hydrogen
            # definition, so Modeller.addHydrogens cannot re-add the backbone
            # amide H that capping would strip. Keep it (systemPrepare already
            # placed it); the extra N-terminal H2/H3 are still dropped.
            remove_names = [n for n in remove_names if n != "H"]
        mol.remove(np.isin(mol.name, remove_names) & mask, _logger=False)

        res_idx = uqres.selectAtoms(mol, indexes=True)
        insert_pos = res_idx[0] if cap == "ACE" else res_idx[-1] + 1
        mol.insert(capmol, insert_pos)


def _defaultProteinCaps(mol):
    from moleculekit.residues import (
        N_TERMINAL_CAP_RESIDUE_NAMES,
        C_TERMINAL_CAP_RESIDUE_NAMES,
    )

    segs = np.unique(mol.get("segid", sel="protein"))
    caps = {}
    for s in segs:
        segmask = mol.segid == s
        if len(np.unique(mol.resid[segmask])) < 10:
            logger.warning(
                f"Segment {s} has fewer than 10 residues - not capped by "
                "default.  Use the caps argument to override."
            )
            continue
        # A segment already carrying a terminal cap (e.g. a deposited
        # C-terminal amide NH2 / NHE) must not be re-capped: aligning an NME
        # onto an existing cap residue has no backbone to match and would
        # otherwise build a spurious terminus. Mirrors amber._defaultProteinCaps.
        seg_atoms = np.where(segmask)[0]
        first_rn = str(mol.resname[seg_atoms[0]])
        last_rn = str(mol.resname[seg_atoms[-1]])
        nterm = "none" if first_rn in N_TERMINAL_CAP_RESIDUE_NAMES else "ACE"
        cterm = "none" if last_rn in C_TERMINAL_CAP_RESIDUE_NAMES else "NME"
        caps[s] = [nterm, cterm]
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


def _strip_metal_coordination_bonds(mol):
    """Return *mol* (a copy if anything changed) with metal-coordination bonds
    removed.

    moleculekit's structure readers store metal coordination with the dedicated
    bond type ``"mc"``. AMBER / ff14SB models such coordination (Mg2+, Zn2+,
    Ca2+, ...) as non-bonded point charges; ``amber.build`` never carries these
    bonds because it deletes all input bonds and lets tLeap regenerate only
    standard + custom ones. The OpenMM path instead keeps ``mol.bonds`` (written
    as CONECT and read back by ``PDBFile``), so a coordination bond such as ADP
    phosphate-O -> Mg (7BTI) would survive and OpenMM's template matcher would
    then reject the coordinated residue ("externally bonded atoms has 1 O atom
    too many"). Drop the ``"mc"`` bonds to match amber.build's treatment. The
    bonds are needed for the CONECT records of the PDB written next, so a copy is
    returned rather than aliasing the caller's arrays.
    """
    if mol.bonds is None or len(mol.bonds) == 0 or mol.bondtype is None:
        return mol
    bondtype = np.asarray(mol.bondtype, dtype=object)
    if len(bondtype) != len(mol.bonds):
        return mol
    keep = bondtype != "mc"
    if keep.all():
        return mol
    mol = mol.copy()
    mol.bonds = mol.bonds[keep]
    mol.bondtype = mol.bondtype[keep]
    return mol


def _mol_to_openmm(mol, outdir, extra_xml=None, skip_peptide_n=None):
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

    *skip_peptide_n* (from :func:`_isopeptide_acceptor_resnames`) is the set
    of resnames whose backbone N is an isopeptide acceptor and so must NOT
    get a ``(-C, N)`` peptide bond from ``createStandardBonds``.

    Returns ``(topology, positions)`` ready for ``ForceField.createSystem``.
    """
    import openmm.app as app

    _register_amber_variant_bond_defs()
    _register_modified_residue_bond_defs(mol)

    # Metal-ion coordination is non-bonded in AMBER; drop those bonds so they
    # don't reach the topology via CONECT and trip the template matcher.
    mol = _strip_metal_coordination_bonds(mol)

    pdb_path = os.path.join(outdir, "input.pdb")
    mol.write(pdb_path, writebonds=True)
    with _temporary_residue_bond_defs(extra_xml, skip_peptide_n=skip_peptide_n):
        pdb = app.PDBFile(pdb_path)

    topology = pdb.topology
    _add_missing_bonds(mol, topology)

    return topology, pdb.positions


def _isopeptide_acceptor_resnames(mol, custombonds):
    """Resnames whose backbone ``N`` is an endpoint of a custombond - i.e.
    the acceptor of an isopeptide / side-chain crosslink, amide-bonded
    through a side chain rather than by the usual backbone peptide bond.

    ``createStandardBonds`` must NOT add the ``(-C, N)`` peptide bond for
    these residues: the N's external bond is the custombond, and a spurious
    backbone bond would give the *donor* residue (whose backbone C is a free
    acid) one external bond too many and break OpenMM's template match.
    Skipping the acceptor's ``(-C, N)`` is also what keeps the donor's C
    free, so a single skip covers both sides of the junction. This is the
    OpenMM analog of the chain-break the amber builder inserts at isopeptide
    junctions.

    The custombonds (from ``detectNonStandardResidues`` /
    ``parameterizeFromSpecs``) are the authoritative list of non-standard
    inter-residue bonds, so we reuse that classification here rather than
    re-deriving it from atom names. *custombonds* is the resolved list of
    ``[UniqueAtomID, UniqueAtomID]`` pairs returned by
    :func:`_prepare_molecule`.
    """
    resnames = set()
    for a1, a2 in (custombonds or []):
        for aid in (a1, a2):
            idx = aid.selectAtom(mol)
            if str(mol.name[idx]) == "N":
                resnames.add(str(mol.resname[idx]))
    return resnames


@contextlib.contextmanager
def _temporary_residue_bond_defs(extra_xml, skip_peptide_n=None):
    """Register the ``<Residue>``/``<Bond>`` entries from each OpenMM
    force-field XML in *extra_xml* into ``Topology._standardBonds`` for
    the duration of the ``with`` block.

    The force-field XML format used by :func:`parameterizeFromSpecs`
    spells bonds as ``<Bond atomName1="..." atomName2="..."/>`` inside
    each ``<Residue>``. ``Topology._standardBonds`` expects
    ``(from, to)`` tuples (the same shape ``Topology.loadBondDefinitions``
    populates). We translate one to the other here.

    *skip_peptide_n* is the set of resnames whose backbone ``N`` is an
    isopeptide / side-chain acceptor rather than a normal peptide-bond
    acceptor. For those we do NOT emit the ``(-C, N)`` peptide bond (see
    :func:`_isopeptide_acceptor_resnames`).

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
                skip_n = skip_peptide_n is not None and name in skip_peptide_n
                for ext in res.findall("ExternalBond"):
                    aname = ext.attrib.get("atomName")
                    if aname == "N" and "N" in atom_names and not skip_n:
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


def _register_modified_residue_bond_defs(mol):
    """Register intra-residue + peptide bond definitions for the force-field-
    supported MODIFIED residues present in *mol* (HYP, MSE, SEP, ...), so
    ``PDBFile.createStandardBonds`` builds their bonds.

    These residues are parameterized by the stock AMBER FF (``ff14SB.xml`` /
    phosaa) but are absent from OpenMM's ``residues.xml``; and because they are
    canonical, ``systemPrepare``'s bond capture drops their intra-residue bonds
    (tLeap rebuilds them from its libraries - the OpenMM path must do it here).
    Without this the residue reaches the template matcher with no bonds ("the
    set of atoms matches HYP, but the residue has no bonds between its atoms").

    Bond definitions are taken from moleculekit's shipped reference cif (named
    atom pairs). Bonds to terminal atoms (OXT etc.) register harmlessly -
    ``createStandardBonds`` only adds a bond when both atoms are present in the
    residue, so an internal residue ignores them. The ``(-C, N)`` entry encodes
    the peptide bond to the previous residue (OpenMM's own convention).
    """
    from openmm.app import Topology
    from moleculekit.residues import MODIFIED_PROTEIN_RESIDUE_NAMES
    from moleculekit import __share_dir

    present = {str(r) for r in np.unique(mol.resname)} & set(
        MODIFIED_PROTEIN_RESIDUE_NAMES
    )
    if not present:
        return
    # Ensure the bundled residues.xml defs are loaded before we add to the dict.
    Topology().createStandardBonds()
    cif_dir = os.path.join(__share_dir, "residue_cifs")
    for resn in sorted(present):
        if resn in Topology._standardBonds:
            continue
        cif = os.path.join(cif_dir, f"{resn}.cif")
        if not os.path.isfile(cif):
            continue
        ref = Molecule(cif)
        if ref.bonds is None or len(ref.bonds) == 0:
            continue
        defs = [
            (str(ref.name[int(a)]), str(ref.name[int(b)])) for a, b in ref.bonds
        ]
        if "N" in set(ref.name.tolist()):
            defs.append(("-C", "N"))
        Topology._standardBonds[resn] = defs


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

    smallmol_ffxml = []
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

        # Emit a static ffxml per small molecule, with partial charges
        # baked in, so the ACEMD handoff carries the parameters as files
        # instead of a runtime template generator.
        for offmol in off_mols:
            smallmol_ffxml.append(gen.generate_residue_template(offmol))

    return forcefield, smallmol_ffxml


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
                    f"{a1.residue.name}:{a1.name} - "
                    f"{a2.residue.name}:{a2.name}"
                )
                continue
            topology.addBond(a1, a2)
            existing.add(key)
            if kind == "cyclic":
                logger.info(
                    f"Added cyclic N-C bond: "
                    f"{a1.residue.name}:{a1.name} - "
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


def _rigidify_three_point_water(struct):
    """Rewrite 3-point water to AMBER's rigid encoding before the prmtop is saved.

    OpenMM force fields express TIP3P-style water rigidity as constraints, but
    ParmEd's ``load_topology`` only translates Forces, so the export System
    (built flexible so every bond keeps an explicit force constant) yields water
    with two O-H bonds and an H-O-H angle. OpenMM's ``AmberPrmtopFile`` only
    rigidifies water that carries an explicit H-H bond; without it,
    ``constraints=HBonds`` constrains the two O-H bonds but leaves the H-O-H
    angle flexible, which is unstable at timesteps of 2 fs or more. This mirrors
    tleap by adding the H-H bond (length derived from the O-H length and H-O-H
    angle) and removing the H-O-H angle, so the water is fully constrained by any
    downstream consumer. The pass is idempotent: water that already carries an
    H-H bond is left untouched.
    """
    import math

    from parmed.topologyobjects import Bond, BondType

    hh_type = None
    n_rigidified = n_angles_removed = 0
    for res in struct.residues:
        atoms = res.atoms
        if len(atoms) != 3 or sorted(at.atomic_number for at in atoms) != [1, 1, 8]:
            continue
        ox = next(at for at in atoms if at.atomic_number == 8)
        h1, h2 = (at for at in atoms if at.atomic_number == 1)
        if any({b.atom1, b.atom2} == {h1, h2} for b in h1.bonds):
            continue  # already rigid (e.g. a tleap-built water)
        r_oh = next(b.type.req for b in ox.bonds if h1 in (b.atom1, b.atom2))
        theta = next(
            (ang.type.theteq for ang in struct.angles
             if ang.atom2 is ox and {ang.atom1, ang.atom3} == {h1, h2}),
            104.52,
        )
        d_hh = 2.0 * r_oh * math.sin(math.radians(theta) / 2.0)
        if hh_type is None:
            hh_type = BondType(553.0, d_hh)
            struct.bond_types.append(hh_type)
        struct.bonds.append(Bond(h1, h2, type=hh_type))
        n_rigidified += 1
        for ang in [a for a in struct.angles
                    if a.atom2 is ox and {a.atom1, a.atom3} == {h1, h2}]:
            ang.delete()
            struct.angles.remove(ang)
            n_angles_removed += 1
    if n_rigidified:
        logger.info(
            f"Rigidified {n_rigidified} water molecules for the prmtop "
            f"(added H-H bonds, removed {n_angles_removed} H-O-H angles)"
        )


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
        _rigidify_three_point_water(struct)
        struct.save(prmtop, overwrite=True)
        struct.save(inpcrd, overwrite=True)
        logger.info(f"Wrote {prmtop} and {inpcrd} via ParmEd")
    except Exception as exc2:
        logger.error(f"ParmEd AMBER export failed: {exc2}")
        raise


def _read_built_molecule(outdir, prefix, topology=None, positions=None):
    """Read back the built system as a ``Molecule``.

    Topology and bonds come from the prmtop. Coordinates come straight from
    the OpenMM ``positions`` (full precision, no PDB round-trip) and the box
    from the OpenMM ``topology`` when given; otherwise the PDB is used as the
    coordinate source (legacy fallback).
    """
    import numpy as np
    from openmm import unit
    from moleculekit.unitcell import box_vectors_to_lengths_and_angles

    prmtop = os.path.join(outdir, f"{prefix}.prmtop")
    pdb = os.path.join(outdir, f"{prefix}.pdb")

    if os.path.exists(prmtop) and os.path.getsize(prmtop) > 0:
        try:
            molbuilt = Molecule(prmtop, validateElements=False)
            if positions is not None:
                coords = np.array(
                    positions.value_in_unit(unit.angstrom), dtype=np.float32
                )
                molbuilt.coords = coords.reshape(molbuilt.numAtoms, 3, 1)
            elif os.path.exists(pdb):
                molbuilt.coords = Molecule(pdb).coords.copy()

            if topology is not None and topology.getPeriodicBoxVectors() is not None:
                a, b, c = topology.getPeriodicBoxVectors()
                la = box_vectors_to_lengths_and_angles(
                    np.array(a.value_in_unit(unit.angstrom)),
                    np.array(b.value_in_unit(unit.angstrom)),
                    np.array(c.value_in_unit(unit.angstrom)),
                )
                molbuilt.box = np.array(la[:3], dtype=np.float32).reshape(3, 1)
                molbuilt.boxangles = np.array(la[3:], dtype=np.float32).reshape(3, 1)
            elif os.path.exists(pdb):
                molbuilt.box = Molecule(pdb).box.copy()
            return molbuilt
        except Exception:
            logger.warning("Could not read prmtop - falling back to PDB.")

    return Molecule(pdb)
