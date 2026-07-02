"""End-to-end parameterization pipeline for non-canonical residues under
AMBER, driven by the spec list returned from
:func:`moleculekit.tools.nonstandard_residues.detectNonStandardResidues`.

:func:`parameterizeFromSpecs` is the user-facing entry point. It walks
``mol.bonds`` to recover cluster grouping (residues sharing non-peptide
inter-residue bonds), builds a combined antechamber model compound per
cluster (full residues + ACE/NME-style backbone caps), runs antechamber +
parmchk2 once per cluster, and splits the output into per-residue
CIF / frcmod pairs. Free residues (no cluster bonds) are parameterized
standalone via :func:`htmd.builder._ambertools._fftype_antechamber`.

The result :class:`ClusterOutputs` carries the topology paths, frcmod
paths, and ``custombonds`` list in the shape that
:func:`htmd.builder.amber.build` expects.

For canonical residues that the detector renamed (CYS bonded to a
scaffold, ASN glycosylated by a sugar, ...), the per-residue CIF carries
ff14SB atom types pulled from the right AMBER residue template (mid-chain
``CYX`` / N-terminal ``NCYX`` / C-terminal ``CCYX`` and the analogous
forms for LYS/HIS/ASN/...) so that backbone bonds resolve against
ff14SB. Per-atom charges come from the antechamber compute on the
combined model, except the backbone atoms of chain-resident residues,
which are pinned to ff14SB: the whole backbone from the ff14SB libraries
for canonical residues, the charge-class amide charges for NCAAs
(see :func:`_backbone_charge_map`). The frcmod carries cross-FF
junction terms (bond / angle / dihedral entries spanning a
canonical-residue atom and a non-canonical
one) with the canonical-side atom types rewritten from antechamber's
GAFF2 to ff14SB.
"""

import copy
import logging
import math
import os
import tempfile
import warnings
from dataclasses import dataclass, field

import numpy as np
from parmed.amber import AmberParameterSet

from moleculekit.molecule import Molecule, UniqueAtomID, UniqueResidueID
from moleculekit.tools._anchor_variants import (
    ANCHOR_TABLE,
)  # noqa: F401  (re-exported)
from moleculekit.util import calculateAnglesAndDihedrals

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cluster machinery: an internal grouping of residues that share at least one
# non-peptide inter-residue covalent bond. The whole cluster is parameterized
# with a single antechamber compute (combined model compound + ACE/NME caps)
# so charges respect the connectivity, then split per-residue. Cluster
# grouping itself is a parameterizer concern; the user-facing detector in
# moleculekit returns flat per-residue specs and never sees these types.
# ---------------------------------------------------------------------------


@dataclass
class ModelAtom:
    """Per-atom record in a cluster model compound. ``role`` is one of
    ``"residue"`` (atom is part of a cluster residue) or ``"cap"`` (an
    ACE/NME-style backbone cap atom that is dropped at split time).
    ``ff_type`` is unused in the current pipeline and kept for forward
    compatibility."""

    role: str
    ff_type: str | None = None


@dataclass
class ClusterBond:
    """One non-peptide covalent bond between two atoms in a cluster.
    Symmetric (no canonical-side / scaffold-side asymmetry), so it works
    uniformly for NCAA-NCAA crosslinks, canonical-AA-anchored scaffolds,
    and everything in between."""

    atom_a: UniqueAtomID
    atom_b: UniqueAtomID


@dataclass
class ClusterSpec:
    """A connected covalent cluster of residues that share non-peptide
    bonds and need combined parameterization. ``residues`` lists every
    cluster member; the four parallel lists carry per-residue metadata
    (chain residency, canonical/non-canonical, role tag, original
    canonical resname for renamed anchors)."""

    subtype: str
    residues: list  # list[UniqueResidueID]
    is_chain_resident: list  # list[bool]
    is_canonical: list  # list[bool]
    roles: list  # list[str]: "scaffold" | "covalent_ligand" | "stapled_ncaa" | "anchor"
    bonds: list  # list[ClusterBond]
    canonical_resnames: list = field(default_factory=list)
    # Parallel to ``residues``: one of ``""`` (mid-chain or non-canonical),
    # ``"n"`` (N-terminal canonical anchor), or ``"c"`` (C-terminal). Drives
    # the AMBER terminal-template lookup (``CCYX`` / ``NCYX`` / ...) when
    # writing per-residue CIFs for canonical residues.
    canonical_terminus: list = field(default_factory=list)
    # Parallel to ``residues``: per-residue terminus flags for chain-resident
    # residues. Drives whether the per-residue prepi declares HEAD / TAIL
    # connect atoms (omit on the side where the residue is terminal).
    is_n_term: list = field(default_factory=list)
    is_c_term: list = field(default_factory=list)

    def __post_init__(self):
        if not self.canonical_resnames:
            self.canonical_resnames = ["" for _ in self.residues]
        if not self.canonical_terminus:
            self.canonical_terminus = ["" for _ in self.residues]
        if not self.is_n_term:
            self.is_n_term = [False for _ in self.residues]
        if not self.is_c_term:
            self.is_c_term = [False for _ in self.residues]


@dataclass
class ClusterModel:
    """Result of :func:`buildClusterModel`. Carries everything
    :func:`prepareClusterResidues` needs to split antechamber's output
    back into per-residue topology files."""

    spec: ClusterSpec
    cif_path: str
    atom_map: dict  # {final_atom_name: ModelAtom}
    atom_to_residue: dict  # {final_atom_name: cluster_residue_index}
    atom_to_orig_name: dict  # {final_atom_name: original_pdb_atom_name}
    canonical_renames: dict  # {cluster_residue_index: new_resname}


@dataclass
class ClusterOutputs:
    """Aggregated result of :func:`parameterizeFromSpecs` /
    :func:`prepareClusterResidues`. Carries the topology, parameter and
    custombond inputs that the user feeds back into
    :func:`htmd.builder.amber.build`."""

    # Per-residue topology files (one CIF per cluster residue) that the
    # user must register with tLeap via ``topo=``. For NCAA-only
    # clusters: one CIF per NCAA. For canonical-anchored clusters: one
    # CIF per unique bucket resname plus one per non-canonical residue.
    topo_paths: list = field(default_factory=list)
    # Per-residue frcmod files to feed via ``param=``.
    frcmod_paths: list = field(default_factory=list)
    # Custombonds to feed via ``custombonds=``. One pair per cluster bond.
    custombonds: list = field(default_factory=list)

    # OpenMM ``<ForceField>`` XML files for ``htmd.builder.openmm.build``.
    # The GAFF backend (when in use) contributes one combined file named
    # ``gaff_<ff>_combined.xml`` listing every GAFF-typed residue + its
    # parameters; the openff backend contributes one per-cluster fragment
    # named ``cluster_<idx>.xml`` (plus one ``<resname>.xml`` per free
    # ligand). Synthetic atom-type names are globally unique by
    # construction (cluster-index prefix), so the files can be loaded
    # together via ``openmm.app.ForceField(*defaultFf(), *xml_paths)``
    # without name collisions. To know whether the GAFF path ran (and
    # therefore whether ``amber.build`` is also a viable consumer),
    # check ``frcmod_paths``: GAFF populates both ``frcmod_paths`` and
    # an entry in ``xml_paths``; openff populates only ``xml_paths``.
    xml_paths: list = field(default_factory=list)


@dataclass
class _ResidueTemplateData:
    """Per-residue OpenMM ``<Residue>`` template input collected as the
    cluster / free-residue pipeline writes per-residue topology files.

    ``mol`` is the typed-mol slice for this residue (atom name / atom
    type / partial charge / element / intra-residue bonds).
    ``external_bond_atom_names`` lists every atom that should appear as
    an ``<ExternalBond>`` in the emitted template: the N- and C-side
    peptide neighbours for chain-resident NCAAs (skipped on the matching
    side for terminal residues) plus every cluster-bond atom (anchor SG,
    scaffold C, crosslink CG, ...). Deduplication by ``resname`` happens
    at emit time."""

    resname: str
    mol: Molecule
    external_bond_atom_names: list = field(default_factory=list)


def _in_pyodide():
    """Detect Pyodide so AmberTools dispatch can swap subprocess for
    ``antechamber_pyodide.run`` automatically."""
    import sys

    return sys.platform == "emscripten"


def _write_chain_residue_prepi(
    sub,
    prepi_path,
    resname,
    is_n_term=False,
    is_c_term=False,
    use_pyodide=None,
    netcharge=0,
):
    """Convert a single chain-resident residue (already typed and named)
    into an AMBER ``.prepi`` topology so tLeap can splice it into a
    peptide chain. Write mol2, run ``antechamber -fo ac`` to convert it to
    AC format without re-typing, build a mainchain config, then run
    ``prepgen``. ``is_n_term`` / ``is_c_term`` omit the matching
    ``HEAD_NAME`` / ``TAIL_NAME`` lines so tLeap does not try to bond
    that side to a phantom neighbour.

    ``netcharge`` is the residue's integer formal charge, used for the
    antechamber ``-nc`` flag and the prepgen ``CHARGE`` field. It is
    passed in rather than summed from ``sub``: ``sub`` is sliced from an
    antechamber-typed mol2, which cannot round-trip formal charge."""
    import networkx as nx
    from htmd.builder._ambertools import (
        _run_ambertools,
        _fix_prepi_atomname_capitalization,
        _write_pyodide_output,
    )

    if use_pyodide is None:
        use_pyodide = _in_pyodide()

    with tempfile.TemporaryDirectory() as tmpdir:
        mol2_name = f"{resname}.mol2"
        ac_name = f"{resname}_mod.ac"
        prepi_name = f"{resname}.prepi"
        mc_name = f"mainchain.{resname.lower()}"

        mol2_path = os.path.join(tmpdir, mol2_name)
        sub.write(mol2_path)

        input_files = None
        if use_pyodide:
            with open(mol2_path, "rb") as f:
                input_files = {mol2_name: f.read()}

        result = _run_ambertools(
            "antechamber",
            [
                "antechamber",
                "-i",
                mol2_name,
                "-fi",
                "mol2",
                "-o",
                ac_name,
                "-fo",
                "ac",
                "-nc",
                str(int(netcharge)),
                "-pf",
                "y",
                "-dr",
                "n",
                "-j",
                "0",
                "-an",
                "n",
            ],
            cwd=tmpdir,
            use_pyodide=use_pyodide,
            input_files=input_files,
        )

        ac_path = os.path.join(tmpdir, ac_name)
        if use_pyodide:
            _write_pyodide_output(ac_path, result[ac_name])

        n_idx = int(np.where(sub.name == "N")[0][0])
        c_idx = int(np.where(sub.name == "C")[0][0])
        backbone = nx.shortest_path(sub.toGraph(), source=n_idx, target=c_idx)

        mainchain_path = os.path.join(tmpdir, mc_name)
        with open(mainchain_path, "w") as f:
            if not is_n_term:
                f.write("HEAD_NAME N\n")
                f.write("PRE_HEAD_TYPE C\n")
            if not is_c_term:
                f.write("TAIL_NAME C\n")
                f.write("POST_TAIL_TYPE N\n")
            for bb_idx in backbone[1:-1]:
                f.write(f"MAIN_CHAIN {sub.name[bb_idx]}\n")
            f.write(f"CHARGE {float(netcharge):.1f}\n")

        prepgen_inputs = None
        if use_pyodide:
            with open(ac_path, "rb") as f:
                ac_data = f.read()
            with open(mainchain_path, "rb") as f:
                mc_data = f.read()
            prepgen_inputs = {ac_name: ac_data, mc_name: mc_data}

        result = _run_ambertools(
            "prepgen",
            [
                "prepgen",
                "-i",
                ac_name,
                "-o",
                prepi_name,
                "-f",
                "prepi",
                "-m",
                mc_name,
                "-rn",
                resname,
            ],
            cwd=tmpdir,
            use_pyodide=use_pyodide,
            input_files=prepgen_inputs,
        )

        prepi_tmp = os.path.join(tmpdir, prepi_name)
        if use_pyodide:
            _write_pyodide_output(prepi_tmp, result[prepi_name])

        import shutil as _shutil

        _shutil.copy(prepi_tmp, prepi_path)
        _fix_prepi_atomname_capitalization(sub, prepi_path)


def _atom_sel_for_unique(uaid):
    """Build an atomselect string targeting a single :class:`UniqueAtomID`,
    quoted for safety. Mirrors moleculekit's ``_atom_sel_from_id`` so we
    don't have to import a private helper."""
    parts = []
    if uaid.segid:
        parts.append(f'segid "{uaid.segid}"')
    if uaid.chain:
        parts.append(f'chain "{uaid.chain}"')
    parts.append(f"resid {int(uaid.resid)}")
    if uaid.insertion:
        parts.append(f'insertion "{uaid.insertion}"')
    parts.append(f'name "{uaid.name}"')
    return " and ".join(parts)


# Cached ff14SB amino-acid libs. The antechamber backend reads from
# AmberTools' amino*.lib via parmed (canonical source for that path);
# the openff backend reads OpenMM's protein.ff14SB.xml (the canonical
# source the OpenMM ForceField will use at apply time). Each cache is
# populated independently the first time its loader is called.
_FF14SB_AMINO_LIB_AMBER = None
_FF14SB_AMINO_LIB_OPENMM = None


def _load_ff14sb_amino_lib_openmm():
    """Build the ff14SB ``{resname: {atom_name: (type, charge)}}`` dict by
    parsing OpenMM's bundled ``amber14/protein.ff14SB.xml``. This is the
    canonical source for the openff/openmm backend: it matches the data
    the OpenMM ``ForceField`` will actually apply at build time, with no
    AmberTools dependency.

    OpenMM declares types as ``protein-CT``, ``protein-N``, etc.; the
    bare AMBER class names (``CT``, ``N``) are used downstream for type
    comparison, so the ``protein-`` prefix is stripped here.
    """
    global _FF14SB_AMINO_LIB_OPENMM
    if _FF14SB_AMINO_LIB_OPENMM is not None:
        return _FF14SB_AMINO_LIB_OPENMM

    import openmm.app as app
    import xml.etree.ElementTree as ET

    xml_path = os.path.join(
        os.path.dirname(app.__file__),
        "data",
        "amber14",
        "protein.ff14SB.xml",
    )
    if not os.path.isfile(xml_path):
        raise FileNotFoundError(
            f"OpenMM ff14SB XML not found at {xml_path}. "
            "OpenMM may be installed without the amber14 data package."
        )

    prefix = "protein-"
    merged = {}
    root = ET.parse(xml_path).getroot()
    for residue in root.find("Residues").findall("Residue"):
        resname = residue.attrib["name"]
        atoms = {}
        for atom in residue.findall("Atom"):
            atomtype = atom.attrib["type"]
            if atomtype.startswith(prefix):
                atomtype = atomtype[len(prefix):]
            atoms[atom.attrib["name"]] = (
                atomtype,
                float(atom.attrib["charge"]),
            )
        merged[resname] = atoms
    _FF14SB_AMINO_LIB_OPENMM = merged
    return _FF14SB_AMINO_LIB_OPENMM


def _load_ff14sb_amino_lib_amber():
    """Lazy-load the ff14SB amino-acid libraries (mid-chain
    ``amino12.lib``, C-terminal ``aminoct12.lib``, N-terminal
    ``aminont12.lib``) so we can extract canonical atom types per
    residue. Terminal-variant entries (``CCYX``, ``NCYX``, ...) sit
    alongside the mid-chain ones in the returned dict.

    This is the canonical source for the antechamber backend: it matches
    what tleap will pull at build time. Needs AmberTools on PATH (or
    ``tleap_pyodide`` under Pyodide) to locate the data dir.

    Returns a dict ``{resname: {atom_name: (type, charge)}}``."""
    global _FF14SB_AMINO_LIB_AMBER
    if _FF14SB_AMINO_LIB_AMBER is not None:
        return _FF14SB_AMINO_LIB_AMBER
    import shutil

    tleap_path = shutil.which("tleap") or shutil.which("teLeap")
    if tleap_path is not None:
        amber_home = os.path.dirname(os.path.dirname(tleap_path))
    else:
        from htmd.builder.amber import _pyodideAmberHome

        amber_home = _pyodideAmberHome()
        if amber_home is None:
            raise RuntimeError(
                "ff14SB amino lib lookup needs tleap on PATH or the "
                "tleap_pyodide package installed (to locate the AmberTools "
                "data dir)."
            )

    from parmed.amber import AmberOFFLibrary

    libdir = os.path.join(amber_home, "dat", "leap", "lib")
    merged = {}
    main_lib_path = os.path.join(libdir, "amino12.lib")
    if not os.path.isfile(main_lib_path):
        raise FileNotFoundError(
            f"ff14SB amino lib not found at {main_lib_path}. "
            "Set AMBERHOME or check your AmberTools installation."
        )
    for name in ("amino12.lib", "aminoct12.lib", "aminont12.lib"):
        path = os.path.join(libdir, name)
        if not os.path.isfile(path):
            continue
        for resname, residue in AmberOFFLibrary.parse(path).items():
            merged[resname] = {a.name: (a.type, a.charge) for a in residue.atoms}
    _FF14SB_AMINO_LIB_AMBER = merged
    return _FF14SB_AMINO_LIB_AMBER


def _clean_frcmod_params(pset, mol, backbone_types, padding_mask):
    """Drop parameters from ``pset`` (a :class:`parmed.amber.AmberParameterSet`)
    that the kept residues never use, mutating it in place:

    - atom (MASS / NONBON) types absent from ``mol`` once ``padding_mask`` atoms
      (cluster caps, neighbour stubs) are excluded;
    - bond / angle / dihedral types referencing an absent type, built entirely
      from ``backbone_types`` (the protein force field tLeap loads provides
      those), or - for purely non-backbone terms - not realized by ``mol``'s
      connectivity. Terms that touch a ``backbone_types`` entry are kept even
      when not realized in ``mol`` itself: in cluster model compounds the chain
      neighbours are ACE/NME-style caps, so a torsion that spans the
      residue/neighbour boundary in the built system is realized over a cap
      atom here and would otherwise be pruned.
    - reversed-duplicate bond / angle / dihedral entries;
    - improper types referencing only absent or backbone types.

    Mirrors the cleanup the legacy noncanonical builder applied to its frcmod
    files; without it parmchk2's full type-combination dump (caps included) is
    carried into every per-residue frcmod.
    """
    backbone_types = list(backbone_types)
    all_at = np.unique(mol.atomtype[~padding_mask]).tolist() + backbone_types
    angles, dihedrals = calculateAnglesAndDihedrals(mol.bonds)
    realized = {
        "bond_types": mol.atomtype[mol.bonds].tolist(),
        "angle_types": mol.atomtype[angles].tolist(),
        "dihedral_types": mol.atomtype[dihedrals].tolist(),
    }

    for at in [t for t in pset.atom_types if t not in all_at]:
        del pset.atom_types[at]

    for param_t in ("bond_types", "angle_types", "dihedral_types"):
        to_delete, seen = [], []
        for bt in getattr(pset, param_t):
            if (
                not all(np.isin(bt, all_at))
                or all(np.isin(bt, backbone_types))
                or (
                    not any(np.isin(bt, backbone_types))
                    and list(bt) not in realized[param_t]
                    and list(bt)[::-1] not in realized[param_t]
                )
                or list(bt)[::-1] in seen
            ):
                to_delete.append(bt)
            seen.append(list(bt))
        for bt in to_delete:
            del pset.__dict__[param_t][bt]

    for param_t in ("improper_types", "improper_periodic_types"):
        to_delete = [
            bt
            for bt in getattr(pset, param_t)
            if not all(np.isin(bt, all_at)) or all(np.isin(bt, backbone_types))
        ]
        for bt in to_delete:
            del pset.__dict__[param_t][bt]


def _retyped_cluster_mol(typed_mol, model, backbone_at):
    """Return ``(retyped, cap_mask)`` where ``retyped`` is a copy of
    ``typed_mol`` whose ``atomtype`` array matches the per-residue
    topology files written by :func:`prepareClusterResidues`: every
    chain-resident residue's backbone atoms (N, CA, C, O, H, HA) carry
    their ff14SB type, everything else (sidechain atoms, cap atoms,
    free residues) keeps antechamber's GAFF2 type. ``cap_mask`` flags
    the cluster-cap atoms (those not in ``model.atom_to_residue``).

    No per-residue template lookup is needed: systemPrepare guarantees
    canonical backbone atom names, so the mapping is a fixed dict.
    Sidechain atoms intentionally keep GAFF2 types - the frcmod
    duplication pass in :func:`_clean_frcmod_params` extends the
    GAFF2-typed bond / angle / dihedral entries that span the
    sidechain-backbone boundary so they resolve under the
    backbone-rewritten ff14SB types too."""
    retyped = typed_mol.copy()
    cap_mask = np.zeros(retyped.numAtoms, dtype=bool)
    for i, name in enumerate(typed_mol.name):
        name = str(name)
        cidx = model.atom_to_residue.get(name)
        if cidx is None:
            cap_mask[i] = True
            continue
        if not model.spec.is_chain_resident[cidx]:
            continue
        orig_name = model.atom_to_orig_name.get(name, name)
        if orig_name in backbone_at:
            retyped.atomtype[i] = backbone_at[orig_name]
    return retyped, cap_mask


# OpenMM's amber14/protein.ff14SB.xml declares ff14SB AMBER atom classes
# under type names of the form ``protein-<class>`` (e.g. ``protein-N``,
# ``protein-CT``, ``protein-S``). Our cluster pipeline renames NCAA
# backbone atoms and canonical-anchor atoms to those AMBER classes
# (``N``, ``CT``, ``S``, ...), but OpenMM's ``ForceField`` resolves
# ``<Atom type="X">`` via the ``<Type name=...>`` lookup, not by class.
# Rewriting the type strings to ``protein-<class>`` at XML-emit time
# lets OpenMM find the ff14SB ``<Type>`` entries directly so the
# residue templates load cleanly. GAFF2 sidechain types stay as-is.
_FF14SB_OMM_TYPE_PREFIX = "protein-"


def _emit_openmm_xml(residue_templates, parameter_sets, xml_path):
    """Combine cluster + free-residue outputs into one OpenMM ForceField XML.

    Walks every cluster's :class:`AmberParameterSet` (post junction-term
    injection, post backbone-rename duplication, post
    :func:`_clean_frcmod_params`) and unions the bond / angle / dihedral /
    improper / atom-type entries by class tuple. GAFF2 derivations are
    deterministic per type tuple, so collisions across clusters carry the
    same value and last-write-wins is safe. Then builds one parmed
    :class:`ResidueTemplate` per unique resname from the per-residue
    typed Molecule slices and writes the merged set as a single OpenMM
    XML.

    The result is a self-contained ``<ForceField>`` document loadable
    via ``openmm.app.ForceField(*defaultFf(), xml_path)``. Peptide bonds
    resolve through ff14SB (NCAA backbones are renamed to ff14SB classes
    upstream); GAFF2 sidechain parameters and cross-FF junction terms
    come from this file. Per-atom charges live in the residue templates
    so two instances of the same resname could in principle carry
    different antechamber charges, but :func:`parameterizeFromSpecs`
    deduplicates by ``(resname, is_n_term, is_c_term)`` upstream, so the
    first template wins.
    """
    from parmed.amber import AmberParameterSet
    from parmed.modeller import ResidueTemplate
    from parmed.openmm import OpenMMParameterSet
    from parmed import Atom as PmdAtom

    merged = AmberParameterSet()
    for pset in parameter_sets:
        merged.atom_types.update(pset.atom_types)
        merged.bond_types.update(pset.bond_types)
        merged.angle_types.update(pset.angle_types)
        # dihedral_types maps a class tuple to a DihedralTypeList (a list of
        # DihedralType objects, one per Fourier term). dict.update would
        # replace the existing list when two clusters both define the same
        # class tuple with different periodicities, dropping terms. Extend
        # the term list instead, deduping by periodicity so identical
        # entries across clusters don't double-count.
        for key, new_list in pset.dihedral_types.items():
            existing = merged.dihedral_types.get(key)
            if existing is None:
                merged.dihedral_types[key] = new_list
                continue
            seen_per = {dt.per for dt in existing}
            for dt in new_list:
                if dt.per not in seen_per:
                    existing.append(dt)
                    seen_per.add(dt.per)
        merged.improper_types.update(pset.improper_types)
        merged.improper_periodic_types.update(pset.improper_periodic_types)

    # Set of AMBER atom-type names ff14SB owns: every type any residue in
    # the ff14SB amino libraries (mid-chain + N/C-terminal variants) uses.
    # An atomtype string in our cluster output that belongs to this set is
    # an ff14SB type (backbone rename or canonical-anchor template) and
    # has to be rewritten to ``protein-<class>`` so OpenMM finds the
    # corresponding ``<Type>`` declaration via the ff14SB XML. This is
    # the antechamber path, so use the parmed-backed amber library
    # (matches what tleap will apply at build time).
    ff14sb_lib = _load_ff14sb_amino_lib_amber()
    ff14sb_classes = set()
    for resname_atoms in ff14sb_lib.values():
        ff14sb_classes.update(t for t, _charge in resname_atoms.values())

    omm = OpenMMParameterSet.from_parameterset(merged, remediate_residues=False)

    seen_resnames = set()
    for rtd in residue_templates:
        if rtd.resname in seen_resnames:
            continue
        seen_resnames.add(rtd.resname)
        template = ResidueTemplate(name=rtd.resname)
        sub = rtd.mol
        name_to_atom = {}
        for i in range(sub.numAtoms):
            atomtype = str(sub.atomtype[i])
            if atomtype in ff14sb_classes:
                atomtype = f"{_FF14SB_OMM_TYPE_PREFIX}{atomtype}"
            atom = PmdAtom(
                name=str(sub.name[i]), type=atomtype, charge=float(sub.charge[i])
            )
            template.add_atom(atom)
            name_to_atom[atom.name] = atom
        for a, b in sub.bonds:
            template.add_bond(int(a), int(b))
        for ext_name in rtd.external_bond_atom_names:
            atom = name_to_atom.get(ext_name)
            if atom is not None:
                template.connections.append(atom)
        omm.residues[rtd.resname] = template

    # protein-* backbone types live in the ff14SB XML loaded next to this
    # file, not in the local AmberParameterSet, so parmed can't resolve them.
    from parmed.exceptions import ParameterWarning

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Some residue templates are using unavailable AtomTypes",
            category=ParameterWarning,
        )
        omm.write(
            xml_path,
            provenance={"Source": ["htmd parameterizeFromSpecs"]},
        )


# ff14SB backbone amide charges (N, H, C, O), keyed by the residue's
# integer net charge. ff14SB constrains these four equal across all
# standard residues of one charge class - neutral, cationic (+1, e.g.
# LYS / ARG / HIP), anionic (-1, e.g. ASP / GLU) - so they are the safe
# fallback for an NCAA, which is absent from the ff14SB libraries.
# CA / HA are residue specific and cannot be set this way; canonical
# residues instead get their full per-residue backbone from the
# libraries via _backbone_charge_map.
_FF14SB_BACKBONE_CHARGES_BY_CLASS = {
    0: {"N": -0.4157, "H": 0.2719, "C": 0.5973, "O": -0.5679},
    1: {"N": -0.3479, "H": 0.2747, "C": 0.7341, "O": -0.5894},
    -1: {"N": -0.5163, "H": 0.2936, "C": 0.5366, "O": -0.5819},
}

# Atom names treated as backbone when a canonical residue is pinned to
# its ff14SB charges: the peptide backbone plus the terminal-variant
# atoms (N-terminal H1 / H2 / H3, C-terminal OXT).
_BACKBONE_ATOM_NAMES = frozenset(
    {"N", "H", "H1", "H2", "H3", "CA", "HA", "HA2", "HA3", "C", "O", "OXT"}
)


def _backbone_charge_map(
    canonical_resname,
    terminus,
    present_names,
    net_charge,
    n_has_bonded_h=True,
    *,
    ff14sb_lib,
):
    """Return ``{atom_name: ff14SB charge}`` for the backbone atoms to pin
    on one chain-resident cluster residue.

    For a canonical residue (``canonical_resname`` is its standard
    resname, e.g. ``"CYS"``) the ff14SB amino libraries supply the whole
    backbone: residue-specific CA / HA included, and the N-terminal /
    C-terminal variants (looked up as ``N`` / ``C`` + resname), so a
    terminal residue is pinned correctly too. ``terminus`` is ``"n"``,
    ``"c"`` or ``""``.

    For an NCAA (``canonical_resname`` empty, or a canonical residue
    whose name is absent from the libraries, e.g. ``HIS`` - the libraries
    are keyed ``HID`` / ``HIE`` / ``HIP``) the libraries have no entry,
    so only the amide atoms N / H / C / O can be pinned, taken from the
    charge class ``net_charge`` selects (a mid-chain residue's net charge
    equals its sidechain charge class). The map is non-empty only for a
    mid-chain residue carrying all four, and only for a net charge in
    ``{-1, 0, +1}``, the classes ff14SB defines.

    For a **proline-analogue** NCAA (e.g. OIC, pipecolic acid, HYP-style
    residues) the amide N has no bonded H - it is ring-closed to a
    sidechain CD instead - and the universal non-PRO amide-N fallback
    would mis-charge N (which assumes an N-H bond partner). For these
    residues the backbone is pinned to **ff14SB PRO** values for
    ``N / CA / C / O / HA``. PRO is the only canonical residue with this
    topology, so it is the right anchor per the Cornell/Cieplak 1995
    convention (proline is its own equivalence class) and the Betz /
    Ramos practitioner convention (freeze backbone to the
    chemically-matching canonical residue's ff14SB values). CD and any
    extra ring carbons refit because their chemistry varies by ring
    size (5-membered pyrrolidine in PRO vs. fused bicyclic in OIC vs.
    6-membered piperidine in pipecolic, ...).

    ``present_names`` is the set of atom names the residue actually has.
    """
    if canonical_resname:
        prefix = {"n": "N", "c": "C"}.get(terminus, "")
        entry = ff14sb_lib.get(f"{prefix}{canonical_resname}")
        if entry is not None:
            return {
                name: entry[name][1]
                for name in present_names
                if name in _BACKBONE_ATOM_NAMES and name in entry
            }

    if terminus:
        return {}

    # Proline-analogue NCAA: the backbone N has no bonded H (ring-
    # closed to a sidechain CD/CE instead of an amide H). Typical of
    # OIC, pipecolic acid, HYP analogues, ... The caller computes
    # ``n_has_bonded_h`` from the actual bond graph - an N-terminal
    # residue with NH3+ atoms named H1/H2/H3 has bonded H and is
    # *not* proline-like even though "H" is absent from the name
    # set. Pin to ff14SB PRO charges for N / CA / C / O / HA; CD and
    # any extra ring atoms refit per the ff15ipq-m / Betz convention.
    # Only fires for net_charge == 0 (ff14SB has no charged-PRO).
    if "N" in present_names and not n_has_bonded_h and net_charge == 0:
        pro_entry = ff14sb_lib.get("PRO")
        if pro_entry is not None:
            return {
                name: pro_entry[name][1]
                for name in ("N", "CA", "C", "O", "HA")
                if name in present_names and name in pro_entry
            }

    classed = _FF14SB_BACKBONE_CHARGES_BY_CLASS.get(net_charge)
    if classed is None or not classed.keys() <= present_names:
        return {}
    return dict(classed)


def _normalize_residue_charges(sub, net_charge, charge_map=None):
    """Pin atoms named in ``charge_map`` (if provided) to those fixed
    charges, then - if ``net_charge`` is not ``None`` - rebalance
    ``sub.charge`` so it sums to ``net_charge`` by spreading the
    residual equally over the remaining (non-pinned) atoms. Returns a
    boolean mask of which atoms were pinned (useful for downstream
    cluster-level redistributions over the same non-pinned set).

    Two roles in one helper:

    * **Backbone pin** (``charge_map`` non-empty, ``net_charge`` given).
      Atoms named in ``charge_map = {atom_name: charge}`` are first set
      to those exact values - typically the ff14SB backbone partial
      charges from :func:`_backbone_charge_map`. The residual relative
      to ``net_charge`` is then spread equally over the remaining
      atoms. Fixes backbone electrostatics for chain-resident residues
      so the modified residue matches the standard ff14SB residues it
      is peptide-bonded to. This is the Robin Betz / R.E.D. / Carlos
      Ramos tutorial convention.

    * **Per-residue total normalization** (``charge_map`` empty / ``None``,
      ``net_charge`` given). No atoms are pinned; the residual is
      spread equally over all atoms. Brings scaffolds and free
      residues to their integer formal charge after a cluster-wide
      charge computation has left them with a small fractional total
      (the smear at C-S / C-N junctions across the cluster), so each
      emitted ``.prepi`` / ``.cif`` is an integer-charged unit - the
      AMBER convention for tLeap.

    * **Pin-only** (``net_charge=None``). Only the ``charge_map`` atoms
      are written; no rebalancing. Used by the cluster-level
      normalization path, which combines per-residue pinning with a
      single cluster-wide residual distribution after all subs are
      pinned.

    Equal distribution of the residual is the minimum-L2 correction
    and is what tutorials like R.E.D. (Cieplak / Dupradeau) achieve
    implicitly through the constrained RESP fit.

    Mutates ``sub.charge`` in place. ``net_charge`` is the residue's
    integer formal charge; the caller sources it from the cluster CIF,
    since a mol2 cannot round-trip formal charge.
    """
    pinned = np.zeros(sub.numAtoms, dtype=bool)
    if charge_map:
        for name, charge in charge_map.items():
            sel = sub.name == name
            if not sel.any():
                continue
            sub.charge[sel] = charge
            pinned |= sel

    if net_charge is None:
        return pinned

    free = ~pinned
    n_free = int(free.sum())
    if n_free == 0:
        return pinned
    residual = float(net_charge) - float(np.sum(sub.charge))
    sub.charge[free] += residual / n_free
    return pinned


# Valid values for the ``normalize`` argument of :func:`parameterizeFromSpecs`
# and :func:`prepareClusterResidues`. See the parameter docstring.
NORMALIZE_MODES = ("per_residue", "cluster", None)


@dataclass
class _ClusterDispatch:
    """Output of :func:`_detect_clusters`. Bundles the per-residue spec map
    plus the multi-residue cluster groupings the parameterise step iterates
    over.

    ``cluster_membership`` is a list of ``(member_residue_indices, bonds)``
    tuples covering every NCAA cluster (cross-linked or singleton). Free
    residues outside any cluster are referenced by ``standalone_residues``
    (residue indices that need standalone parameterisation), with
    ``in_cluster_set`` recording which residue indices are already covered
    by ``cluster_membership``.
    """
    # All residue groups in mol, indexed 0..N-1 (output of
    # :func:`_residue_groups_with_index`).
    groups: list = field(default_factory=list)
    # Per-atom residue index (output of :func:`_residue_groups_with_index`).
    a2r: object = None
    # Map: residue index in mol -> spec for the residues we'll parameterise.
    spec_by_res_idx: dict = field(default_factory=dict)
    # Per-cluster (member residue indices, inter-residue bond tuples).
    cluster_membership: list = field(default_factory=list)
    # Residue indices already absorbed by ``cluster_membership``.
    in_cluster_set: set = field(default_factory=set)


def _engine_for_forcefield(ff_name):
    """Return the internal engine tag (``"antechamber"`` or ``"openff"``)
    for a given forcefield string. Any name starting with ``gaff``
    (case-insensitive) goes through antechamber + parmchk2; everything
    else is treated as a SMIRNOFF offxml filename and goes through
    OpenFF Interchange."""
    if not isinstance(ff_name, str):
        raise TypeError(
            f"forcefield name must be a string, got {type(ff_name).__name__}"
        )
    return "antechamber" if ff_name.lower().startswith("gaff") else "openff"


def _forcefields_in(forcefield):
    """Yield every forcefield name referenced by the user-supplied
    ``forcefield`` argument, regardless of whether it is a string or a
    per-resname dict. Used by validators that need to inspect the whole
    set up front.
    """
    if isinstance(forcefield, str):
        yield forcefield
        return
    if isinstance(forcefield, dict):
        yield from forcefield.values()


def _warn_if_openff_mismatched_charges(forcefield, charge_method):
    """Log a warning if a SMIRNOFF (openff) forcefield is used with any
    charge method other than ``"nagl"``. The Sage vdW + torsion
    parameters were fit alongside a specific charge model (AM1-BCC for
    Sage 2.0-2.2, AshGC for Sage 2.3.0), and NAGL is the OpenFF-
    recommended fast surrogate that reproduces those charges. Using
    any other charge model against Sage's LJ + torsions gives a
    thermodynamically inconsistent potential, but we leave the choice
    to the caller (warn, don't error).
    """
    if charge_method == "nagl":
        return
    for ff_name in _forcefields_in(forcefield):
        if _engine_for_forcefield(ff_name) != "openff":
            continue
        logger.warning(
            f"forcefield={ff_name!r} (SMIRNOFF) is being used with "
            f"charge_method={charge_method!r}; OpenFF Sage was fit "
            "alongside a specific charge model (AM1-BCC for Sage 2.0-"
            "2.2, AshGC for Sage 2.3.0) that NAGL reproduces. Using "
            "any other charge model against Sage's LJ + torsion "
            "parameters gives a thermodynamically inconsistent "
            "potential. Switch to charge_method='nagl' for the "
            "recommended behaviour."
        )
        return  # one warning per call is enough


def _assert_no_chain_resident_modified_nucleotide(mol, specs):
    """Refuse a non-canonical residue that is covalently spliced into a
    nucleic-acid chain by a phosphodiester bond (its ``P`` to a neighbour's
    ``O3'``, or its ``O3'`` to a neighbour's ``P``).

    Such a residue is a chain-resident MODIFIED NUCLEOTIDE. detect has no nucleic
    chain-residency rule (only peptide bonds make a residue chain-resident), so
    it falls through to the standalone scaffold / ligand path - which GAFF-types
    the whole sugar-phosphate backbone and emits no inter-residue bond, leaving
    the nucleotide DISCONNECTED from the chain in the built topology (a silent
    correctness bug). Parameterizing a modified nucleotide in place needs the
    nucleic analogue of the chain-resident backbone retyping (not yet
    implemented); library-backed modified nucleotides (modrna08, e.g. 5MC/PSU)
    are recognised earlier by detect and never reach here. Fail with an
    actionable message rather than build a broken chain.
    """
    from moleculekit.tools.nonstandard_residues import (
        ScaffoldSpec,
        CovalentLigandSpec,
        LigandSpec,
    )

    if mol.bonds is None or not len(mol.bonds):
        return
    standalone = (ScaffoldSpec, CovalentLigandSpec, LigandSpec)
    spec_resname = {}  # standalone-spec atom index -> resname
    for spec in specs:
        if isinstance(spec, standalone):
            for i in np.where(spec.residue.selectAtoms(mol, indexes=False))[0]:
                spec_resname[int(i)] = spec.resname
    if not spec_resname:
        return

    name, resid, chain, segid = mol.name, mol.resid, mol.chain, mol.segid
    o3 = ("O3'", "O3*")
    for a, b in mol.bonds:
        a, b = int(a), int(b)
        if resid[a] == resid[b] and chain[a] == chain[b] and segid[a] == segid[b]:
            continue  # intra-residue
        for inside, outside in ((a, b), (b, a)):
            if inside not in spec_resname:
                continue
            ni, no = str(name[inside]), str(name[outside])
            if (ni == "P" and no in o3) or (ni in o3 and no == "P"):
                raise RuntimeError(
                    f"Residue {spec_resname[inside]} is a modified nucleotide "
                    f"covalently linked into a nucleic-acid chain (phosphodiester "
                    f"at {ni}). Chain-resident modified nucleotides cannot yet be "
                    f"parameterized de novo - only library-backed ones (modrna08, "
                    f"e.g. 5MC / PSU) are supported. Building it through the "
                    f"standalone path would leave it disconnected from the chain. "
                    f"Supply pre-built parameters (a frcmod + prepi/off named "
                    f"<RESNAME>) via amber.build's `param` / `topo`, swap in a "
                    f"supported modified nucleotide, or remove it before building."
                )


def _warn_if_nonalpha_backbone(mol, specs):
    """Loudly warn for every chain-resident residue whose backbone is non-alpha
    - a shortest N -> C path longer than the standard 3-atom ``N-CA-C`` (beta /
    gamma amino acids such as microcystin's Adda, ``N-CA-C18-C``).

    Such residues build, and their junction atoms ``N`` / ``C`` are still pinned
    to ff14SB so the bonds to neighbouring canonical residues stay consistent.
    But the backbone atoms BEYOND ``N-CA-C`` fall outside the alpha backbone
    atom-name set, so they are treated as sidechain: their partial charges are
    refit by GAFF / OpenFF and their backbone torsions are assigned by analogy.
    That is known to be inaccurate for non-alpha amino-acid backbones - GAFF /
    SMIRNOFF torsions are fit to small molecules, not peptide backbones - so we
    warn rather than silently emit unreliable parameters.
    """
    from collections import deque
    from moleculekit.tools.nonstandard_residues import ChainResidueSpec

    for spec in specs:
        if not isinstance(spec, ChainResidueSpec):
            continue
        rid = spec.residue
        idx = np.where(
            (mol.segid == str(rid.segid))
            & (mol.chain == str(rid.chain))
            & (mol.resid == int(rid.resid))
            & (mol.insertion == str(rid.insertion))
        )[0]
        if not len(idx):
            continue
        idxset = {int(i) for i in idx}
        byname = {}
        for i in idx:
            byname.setdefault(str(mol.name[int(i)]), int(i))
        if "N" not in byname or "C" not in byname:
            continue

        adj = {int(i): [] for i in idx}
        for a, b in mol.bonds:
            a, b = int(a), int(b)
            if a in idxset and b in idxset:
                adj[a].append(b)
                adj[b].append(a)

        # BFS for the shortest N -> C path within the residue.
        src, dst = byname["N"], byname["C"]
        prev = {src: None}
        dq = deque([src])
        while dq:
            x = dq.popleft()
            if x == dst:
                break
            for nb in adj[x]:
                if nb not in prev:
                    prev[nb] = x
                    dq.append(nb)
        if dst not in prev:
            continue
        path = []
        x = dst
        while x is not None:
            path.append(x)
            x = prev[x]
        if len(path) <= 3:
            continue  # standard alpha backbone (N-CA-C)

        path_names = " ".join(str(mol.name[i]) for i in reversed(path))
        logger.warning(
            f"Residue {spec.resname} {rid.chain}:{rid.resid}{rid.insertion} has a "
            f"non-alpha amino-acid backbone (shortest N->C path: {path_names}). It "
            f"will build and its junction atoms N/C are pinned to ff14SB, but the "
            f"backbone atoms beyond the standard N-CA-C are treated as sidechain: "
            f"their partial charges are refit by GAFF/OpenFF and their backbone "
            f"torsions are assigned by analogy. This is KNOWN TO BE INACCURATE for "
            f"non-alpha (beta/gamma) amino-acid backbones - GAFF/SMIRNOFF torsions "
            f"are fit to small molecules, not peptide backbones. For production, "
            f"validate the backbone or supply bespoke QM-derived torsion "
            f"parameters (e.g. OpenFF BespokeFit)."
        )


def _resolve_forcefield(forcefield):
    """Normalise the ``forcefield`` argument to a callable
    ``resname -> (engine, ff_name)``.

    ``engine`` is the internal dispatch tag (``"antechamber"`` for GAFF
    / ``"openff"`` for SMIRNOFF); ``ff_name`` is the user-supplied force-
    field name passed through to the engine's parameter generator
    (antechamber's ``-at`` flag for GAFF, the offxml filename for
    SMIRNOFF). ``forcefield`` may be a single string applied to every
    resname, or a dict mapping ``resname -> name`` with optional
    ``"default"`` key providing the fallback.
    """
    if isinstance(forcefield, str):
        engine = _engine_for_forcefield(forcefield)
        return lambda resname: (engine, forcefield)
    if isinstance(forcefield, dict):
        default = forcefield.get("default", "gaff2")
        # Validate every entry up front so a typo doesn't only surface
        # when the offending resname is dispatched.
        for k, v in forcefield.items():
            _engine_for_forcefield(v)  # raises TypeError if not a string
        def _lookup(resname):
            ff = forcefield.get(resname, default)
            return (_engine_for_forcefield(ff), ff)
        return _lookup
    raise TypeError(
        f"forcefield must be str or dict, got {type(forcefield).__name__}"
    )


def _detect_clusters(mol, specs):
    """Walk ``specs`` and ``mol.bonds`` to group residues into clusters.

    Two responsibilities:
      1. Map every spec to its residue index in ``mol`` and validate that
         the spec'd residues are protonated (delegated to
         :func:`_check_specs_protonated`).
      2. Group spec'd residues into clusters by union-find over non-peptide
         inter-residue bonds. Chain-resident non-canonical residues that
         end up outside any cluster are promoted to singleton clusters so
         the cluster pipeline still runs for them.

    CYS-CYS disulfides (``ChainResidueSpec`` with ``new_resname="CYX"``)
    are filtered out - they have a native ff14SB template and tLeap stitches
    the S-S bond via ``amber.build``'s disulfide detection.

    Returns
    -------
    :class:`_ClusterDispatch`
        Bundles the cluster layout that the parameterise loop iterates
        over.
    """
    from moleculekit.tools.nonstandard_residues import (
        ChainResidueSpec,
        PROTEIN_RESNAMES,
        infer_nonstandard_junction_bonds,
    )

    a2r, groups = _residue_groups_with_index(mol)

    res_key_to_idx = {
        (g["segid"], g["chain"], int(g["resid"]), g["insertion"]): ri
        for ri, g in enumerate(groups)
    }
    spec_by_res_idx = {}
    for s in specs:
        # CYS-CYS disulfides are detected as ChainResidueSpec with
        # ``new_resname="CYX"``. ff14SB has a native CYX template and
        # tLeap stitches the S-S bond via amber.build's disulfide
        # detection, so the cluster pipeline must not see them - they
        # would otherwise be wrongly clustered through their non-peptide
        # SG-SG bond.
        if isinstance(s, ChainResidueSpec) and s.new_resname == "CYX":
            continue
        key = _spec_residue_key(s.residue)
        if key not in res_key_to_idx:
            raise ValueError(
                f"Spec residue {s.residue} does not map to any residue in mol "
                f"(segid={key[0]!r}, chain={key[1]!r}, resid={key[2]}, "
                f"insertion={key[3]!r})"
            )
        spec_by_res_idx[res_key_to_idx[key]] = s

    _check_specs_protonated(mol, spec_by_res_idx, groups)

    # Walk mol.bonds for non-peptide inter-residue bonds among spec'd residues.
    # Also consult inter-residue backbone-continuation bonds the input omits at a
    # non-standard junction (an undeposited side-chain isopeptide), recovered from
    # geometry so the cluster forms without those bonds being stored on mol.bonds.
    spec_indices = set(spec_by_res_idx)
    inferred = infer_nonstandard_junction_bonds(mol)
    inter_bonds = []
    for a1, a2 in list(mol.bonds) + [tuple(p) for p in inferred]:
        a1, a2 = int(a1), int(a2)
        r1, r2 = int(a2r[a1]), int(a2r[a2])
        if r1 == r2 or r1 < 0 or r2 < 0:
            continue
        if r1 not in spec_indices or r2 not in spec_indices:
            continue
        if {str(mol.name[a1]), str(mol.name[a2])} == {"N", "C"}:
            continue
        inter_bonds.append((a1, a2, r1, r2))

    # Union-find clusters over the spec'd residues.
    parent = {r: r for r in spec_indices}

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a, b):
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    for a1, a2, r1, r2 in inter_bonds:
        _union(r1, r2)

    components = {}
    for r in spec_indices:
        components.setdefault(_find(r), set()).add(r)

    in_cluster_set = set()
    cluster_membership = []
    for members in components.values():
        cbonds = [b for b in inter_bonds if b[2] in members and b[3] in members]
        if cbonds:
            cluster_membership.append((sorted(members), cbonds))
            in_cluster_set.update(members)

    # Promote each chain-resident NCAA that didn't end up in a multi-
    # residue cluster into a singleton cluster of its own. The cluster
    # pipeline picks up ACE/NME-style backbone caps from the live mol
    # automatically; that's what gives free NCAAs the same combined-
    # parameterization treatment as crosslinked ones.
    for r_idx, spec in spec_by_res_idx.items():
        if r_idx in in_cluster_set:
            continue
        if isinstance(spec, ChainResidueSpec) and spec.resname not in PROTEIN_RESNAMES:
            cluster_membership.append(([r_idx], []))
            in_cluster_set.add(r_idx)

    return _ClusterDispatch(
        groups=groups,
        a2r=a2r,
        spec_by_res_idx=spec_by_res_idx,
        cluster_membership=cluster_membership,
        in_cluster_set=in_cluster_set,
    )


def _parameterize_cluster_antechamber(
    cluster_model,
    cluster_outdir,
    *,
    forcefield,
    charge_method,
    am1_path_length,
    use_pyodide,
    pin_backbone_charges,
    normalize,
    residue_templates,
    parameter_sets,
):
    """Antechamber backend for one cluster: type the combined model compound
    via :func:`_fftype_antechamber` using the chosen GAFF variant, then
    split per-residue via :func:`prepareClusterResidues`. Mutates
    ``residue_templates`` and ``parameter_sets`` in place with the
    per-residue contributions for the final XML emit.

    Returns the :class:`ClusterOutputs`-shaped artifact bundle for this
    cluster (topo, frcmod, custombonds path lists).
    """
    from htmd.builder._ambertools import _fftype_antechamber

    cluster_mol = Molecule(cluster_model.cif_path)
    netcharge = int(round(float(np.sum(cluster_mol.formalcharge))))
    typed_mol, typed_path, frcmod_path = _fftype_antechamber(
        cluster_mol,
        tmpdir=cluster_outdir,
        forcefield=forcefield,
        netcharge=netcharge,
        charge_method=charge_method,
        am1_path_length=am1_path_length,
        use_pyodide=use_pyodide,
    )
    del typed_mol  # re-loaded inside prepareClusterResidues
    return prepareClusterResidues(
        typed_path,
        frcmod_path,
        cluster_model,
        outdir=cluster_outdir,
        use_pyodide=use_pyodide,
        residue_templates=residue_templates,
        parameter_sets=parameter_sets,
        pin_backbone_charges=pin_backbone_charges,
        normalize=normalize,
    )


def _parameterize_free_residue_antechamber(
    mol,
    group,
    ligand_dir,
    outdir,
    *,
    forcefield,
    charge_method,
    am1_path_length,
    use_pyodide,
    normalize,
    residue_templates,
    parameter_sets,
):
    """Antechamber backend for one free residue: write a per-residue CIF,
    type and charge it via :func:`_fftype_antechamber` using the chosen
    GAFF variant, clean the frcmod, and append to the residue-template /
    parameter-set accumulators.

    Returns ``(out_cif_path, out_frcmod_path)`` for the per-residue
    ``ClusterOutputs.topo_paths`` / ``ClusterOutputs.frcmod_paths`` lists.
    """
    from htmd.builder._ambertools import _fftype_antechamber

    ligand_cif = os.path.join(ligand_dir, f"{group['resname']}.cif")
    sub = _write_free_residue_cif(mol, group, ligand_cif)
    netcharge = int(round(float(np.sum(sub.formalcharge))))
    typed_mol, _, frcmod_path = _fftype_antechamber(
        sub,
        tmpdir=ligand_dir,
        forcefield=forcefield,
        netcharge=netcharge,
        charge_method=charge_method,
        am1_path_length=am1_path_length,
        use_pyodide=use_pyodide,
    )
    out_cif = os.path.join(outdir, f"{group['resname']}.cif")
    out_frcmod = os.path.join(outdir, f"{group['resname']}.frcmod")
    typed_mol.resname[:] = group["resname"]
    # Antechamber's mol2 rounds charges to 4 decimals; with many atoms the
    # sum drifts a few thousandths off the integer net charge. Renormalize
    # the standalone unit unless the user has opted out via
    # ``normalize=None``. (For a single-residue ligand, "cluster" and
    # "per_residue" are identical.)
    if normalize is not None:
        _normalize_residue_charges(typed_mol, netcharge)
    typed_mol.write(out_cif)
    pset = AmberParameterSet(frcmod_path)
    _clean_frcmod_params(
        pset, typed_mol, [], np.zeros(typed_mol.numAtoms, dtype=bool)
    )
    pset.write(out_frcmod, title=f"parameters for {group['resname']}", style="frcmod")
    residue_templates.append(
        _ResidueTemplateData(
            resname=group["resname"],
            mol=typed_mol.copy(),
            external_bond_atom_names=[],
        )
    )
    parameter_sets.append(copy.deepcopy(pset))
    return out_cif, out_frcmod


def _apply_cluster_charge_policy(
    off_mol,
    cluster_model,
    residue_atom_indices,
    cap_atoms,
    original_names,
    pin_backbone_charges,
    normalize,
):
    """Mutate ``off_mol.partial_charges`` in place to apply backbone
    pinning (ff14SB-canonical backbone charges on chain-resident residues)
    and per-residue / cluster normalisation. Run BEFORE
    :func:`Interchange.from_smirnoff` so the charges flow through
    ``charge_from_molecules`` into the resolved Interchange. Mirrors the
    per-residue ``charge_map`` + cluster-residual distribution that
    :func:`prepareClusterResidues` runs on the antechamber-typed mol2.
    """
    from openff.units import unit

    # Openff path: the ff14SB charges we pin come from OpenMM's bundled
    # protein.ff14SB.xml so they match exactly what app.ForceField will
    # apply at build time (no AmberTools / parmed dependency here).
    ff14sb_lib = _load_ff14sb_amino_lib_openmm()

    n_atoms = off_mol.n_atoms
    if off_mol.partial_charges is None:
        raise ValueError(
            "off_mol must have partial_charges assigned before calling "
            "_apply_cluster_charge_policy. Compute Gasteiger or AM1-BCC "
            "first and assign via off_mol.partial_charges = ..."
        )
    charges = np.array(
        [
            float(q.m_as(unit.elementary_charge))
            for q in off_mol.partial_charges
        ],
        dtype=np.float64,
    )

    pinned_mask = np.zeros(n_atoms, dtype=bool)

    # Pass 1: per-residue backbone pin.
    for ri, atom_indices in enumerate(residue_atom_indices):
        if not atom_indices:
            continue
        is_chain = (
            ri < len(cluster_model.spec.is_chain_resident)
            and cluster_model.spec.is_chain_resident[ri]
        )
        if not (pin_backbone_charges and is_chain):
            continue

        # Resolve the residue's atom-name set, terminus, and
        # canonical-resname hint (for canonical anchor residues) so
        # _backbone_charge_map picks the right ff14SB entry.
        cluster_residue = cluster_model.spec.residues[ri]
        present_names = {original_names[ai] for ai in atom_indices}
        # canonical_resnames is only non-empty for canonical anchor residues
        # (CYX/HID-renamed canonicals). For NCAAs it's "" and the charge map
        # falls back to the amide-class table.
        canonical_resname = ""
        if (
            ri < len(cluster_model.spec.canonical_resnames)
            and ri < len(cluster_model.spec.is_canonical)
            and cluster_model.spec.is_canonical[ri]
        ):
            canonical_resname = cluster_model.spec.canonical_resnames[ri]
        # Match the antechamber backend (nonstandard.py:1695-1699): derive
        # terminus from is_n_term / is_c_term, not from canonical_terminus.
        # For a C-terminal NCAA, terminus="c" causes _backbone_charge_map
        # to return {} (no pin) - the universal amide pin would mis-charge
        # the carboxylate C / OXT.
        if ri < len(cluster_model.spec.is_n_term) and cluster_model.spec.is_n_term[ri]:
            terminus = "n"
        elif ri < len(cluster_model.spec.is_c_term) and cluster_model.spec.is_c_term[ri]:
            terminus = "c"
        else:
            terminus = ""
        # Heuristic for proline-like NCAAs (no bonded H on N): look at
        # whether the cluster Interchange has an N-H bond among the
        # residue's atoms.
        n_idx = next(
            (
                ai
                for ai in atom_indices
                if original_names[ai] == "N"
            ),
            None,
        )
        n_has_bonded_h = False
        if n_idx is not None:
            for bond in off_mol.bonds:
                a, b = bond.atom1_index, bond.atom2_index
                if a == n_idx or b == n_idx:
                    other = b if a == n_idx else a
                    elem = off_mol.atoms[other].atomic_number
                    if elem == 1:  # hydrogen
                        n_has_bonded_h = True
                        break

        net_charge = int(round(sum(
            float(off_mol.atoms[ai].formal_charge.m_as(
                unit.elementary_charge
            ))
            for ai in atom_indices
        )))

        charge_map = _backbone_charge_map(
            canonical_resname=canonical_resname,
            terminus=terminus,
            present_names=present_names,
            net_charge=net_charge,
            n_has_bonded_h=n_has_bonded_h,
            ff14sb_lib=ff14sb_lib,
        )
        if not charge_map:
            continue

        for ai in atom_indices:
            name = original_names[ai]
            if name in charge_map:
                charges[ai] = float(charge_map[name])
                pinned_mask[ai] = True

    # Pass 2: residual normalisation.
    if normalize == "per_residue":
        for ri, atom_indices in enumerate(residue_atom_indices):
            if not atom_indices:
                continue
            net_q = int(round(sum(
                float(off_mol.atoms[ai].formal_charge.m_as(
                    unit.elementary_charge
                ))
                for ai in atom_indices
            )))
            non_pinned = [ai for ai in atom_indices if not pinned_mask[ai]]
            if not non_pinned:
                continue
            current = sum(charges[ai] for ai in atom_indices)
            residual = net_q - current
            per_atom = residual / len(non_pinned)
            for ai in non_pinned:
                charges[ai] += per_atom
    elif normalize == "cluster":
        # Distribute the residue residual across RESIDUE non-pinned
        # atoms only (caps untouched). Mirrors the antechamber pipeline
        # at nonstandard.py:1762-1770, where per_residue_data carries
        # only residue atoms - cluster_total / cluster_residual are
        # computed and distributed over residue atoms.
        cluster_atoms = [
            ai for atom_indices in residue_atom_indices for ai in atom_indices
        ]
        if cluster_atoms:
            net_q = int(round(sum(
                float(off_mol.atoms[ai].formal_charge.m_as(
                    unit.elementary_charge
                ))
                for ai in cluster_atoms
            )))
            non_pinned = [ai for ai in cluster_atoms if not pinned_mask[ai]]
            if non_pinned:
                current = sum(charges[ai] for ai in cluster_atoms)
                residual = net_q - current
                per_atom = residual / len(non_pinned)
                for ai in non_pinned:
                    charges[ai] += per_atom
    # normalize is None: leave residual as-is.

    # Interchange.from_smirnoff sanity-checks off_mol's TOTAL charge
    # against the sum of formal charges over ALL atoms (residue + caps).
    # After residue-only normalization, cap atoms still carry their raw
    # Gasteiger values which sum to ~0 but not exactly (RDKit's Gasteiger
    # is iterative and only converges to a small residual). Zero the cap
    # sub-total by spreading the cap-Gasteiger residual equally across
    # cap atoms - residue atoms (which we just normalised) stay untouched.
    cap_list = sorted(cap_atoms)
    if cap_list:
        # Cap formal charge is 0 by construction (ACE/NME methyl).
        cap_residual = sum(charges[ai] for ai in cap_list)
        cap_shift = -cap_residual / len(cap_list)
        for ai in cap_list:
            charges[ai] += cap_shift

    # Write back to off_mol.partial_charges so it flows through
    # Interchange.from_smirnoff via charge_from_molecules=[off_mol].
    off_mol.partial_charges = charges * unit.elementary_charge


def _parameterize_cluster_openff(
    cluster_model,
    cluster_outdir,
    *,
    forcefield,
    charge_method,
    use_pyodide,
    pin_backbone_charges,
    normalize,
    interchange_xml_paths,
    cluster_index=0,
):
    """OpenFF backend for one cluster (one or more NCAA residues plus
    ACE/NME-style caps): parameterise the cluster model compound via
    SMIRNOFF, drop cap atoms, and emit one OpenMM ``<ForceField>`` XML
    fragment covering every cluster residue as a separate ``<Residue>``
    template.

    The cluster Interchange covers the cap atoms during parameterisation
    (they provide chemical context at the cluster boundary), but they are
    NOT emitted into the final XML - the canonical residues on the other
    side of the chain peptide bond come from ff14SB at load time. This
    requires backbone-class renaming on the NCAA residues so ff14SB
    matches the peptide-bond terms across the junction (Phase 2D - not
    yet implemented).
    """
    from htmd.builder.openmm import (
        _emit_openmm_xml_from_cluster_interchange,
        _assign_nagl_charges,
        _assign_resp_charges,
        _require_openff_python,
    )
    from htmd.builder._ambertools import _assign_rdkit_gasteiger_charges

    _require_openff_python()
    from openff.toolkit import ForceField
    from openff.interchange import Interchange

    if use_pyodide:
        raise NotImplementedError(
            "backend='openff' is not yet supported under Pyodide - SMIRNOFF "
            "charge assignment needs RDKit's full chemistry stack which is "
            "not available there. Use backend='antechamber' instead."
        )

    cluster_mol = Molecule(cluster_model.cif_path)

    # Charges. Default to Gasteiger via RDKit if requested, NAGL or RESP
    # if explicitly chosen, otherwise let SMIRNOFF assign (typically AM1-
    # BCC). RESP dispatches to the private Acellera ``parameterize``
    # package, which runs Psi4 and openff-recharge.
    if charge_method == "gasteiger":
        _assign_rdkit_gasteiger_charges(cluster_mol)
    elif charge_method == "nagl":
        _assign_nagl_charges(cluster_mol)
    elif charge_method in ("resp", "resp-multiconf"):
        if charge_method == "resp-multiconf":
            logger.warning(
                "resp-multiconf downgraded to single-conformer RESP for "
                "cluster compounds: RDKit's ETKDG conformer generation is "
                "not appropriate for clusters with ACE/NME caps and "
                "peptide bonds. The existing cluster geometry is used."
            )
        _assign_resp_charges(cluster_mol, multi_conf=False)
    elif charge_method not in (None, "am1-bcc"):
        raise ValueError(
            f"backend='openff' unsupported charge_method {charge_method!r}. "
            f"Use 'gasteiger', 'nagl', 'resp', 'resp-multiconf', "
            f"'am1-bcc', or None."
        )

    off_mol = cluster_mol.toOpenFFMolecule(sanitize=True, assignStereo=True)

    # Map atom index -> cluster residue index (or None for cap atoms).
    atom_names_off = [a.name for a in off_mol.atoms]
    atom_to_residue_idx = [
        cluster_model.atom_to_residue.get(nm) for nm in atom_names_off
    ]
    cap_atoms = {
        i for i, ri in enumerate(atom_to_residue_idx) if ri is None
    }
    original_names = [
        cluster_model.atom_to_orig_name.get(nm, nm) for nm in atom_names_off
    ]

    # Build per-cluster-residue atom-index lists.
    n_cluster_residues_pre = len(cluster_model.spec.residues)
    residue_atom_indices_pre = [[] for _ in range(n_cluster_residues_pre)]
    for i, ri in enumerate(atom_to_residue_idx):
        if ri is None:
            continue
        residue_atom_indices_pre[ri].append(i)

    # Mutate off_mol.partial_charges BEFORE building the Interchange, so
    # backbone pinning / cluster normalisation flow through
    # ``charge_from_molecules`` instead of trying to mutate the resolved
    # Interchange's electrostatics (which exposes ``charges`` as a
    # computed property, not a mutable mapping).
    if (pin_backbone_charges or normalize is not None) and charge_method is not None:
        _apply_cluster_charge_policy(
            off_mol=off_mol,
            cluster_model=cluster_model,
            residue_atom_indices=residue_atom_indices_pre,
            cap_atoms=cap_atoms,
            original_names=original_names,
            pin_backbone_charges=pin_backbone_charges,
            normalize=normalize,
        )

    ff = ForceField(forcefield)
    smirnoff_kwargs = {"force_field": ff, "topology": [off_mol]}
    if charge_method is not None:
        smirnoff_kwargs["charge_from_molecules"] = [off_mol]
    ic = Interchange.from_smirnoff(**smirnoff_kwargs)

    # The per-cluster-residue atom-index lists were already built above
    # (residue_atom_indices_pre) for the charge-policy pass; reuse them.
    residue_atom_indices = residue_atom_indices_pre

    # Build residue_assignment: per-cluster-residue (resname, atom indices).
    residue_assignment = []
    external_bonds = {}
    for ri, atom_indices in enumerate(residue_atom_indices):
        if not atom_indices:
            continue
        cluster_residue = cluster_model.spec.residues[ri]
        # Use the renamed resname if the cluster pipeline applied a
        # canonical rename (e.g. CYS->CYX); otherwise keep the original.
        resname = cluster_model.canonical_renames.get(ri, cluster_residue.resname)
        residue_assignment.append((resname, atom_indices))

        # External-bond markers: N/C for chain residues (one or both,
        # depending on terminal flags). For free / scaffold residues no
        # peptide-bond external bonds; non-peptide cluster crosslinks
        # (disulfides, isopeptides) are emitted separately as custombonds
        # by the main parameterizeFromSpecs loop.
        is_chain = (
            ri < len(cluster_model.spec.is_chain_resident)
            and cluster_model.spec.is_chain_resident[ri]
        )
        if is_chain:
            ext = []
            if not cluster_model.spec.is_n_term[ri]:
                ext.append("N")
            if not cluster_model.spec.is_c_term[ri]:
                ext.append("C")
            external_bonds[resname] = ext
        else:
            external_bonds[resname] = []

    # Backbone class overrides: for each chain-resident NCAA residue,
    # rewrite the backbone atom classes (N/H/CA/HA/C/O and OXT for the
    # C-terminus) to the ff14SB equivalents prefixed with ``protein-``.
    # That lets the peptide bond at the canonical-NCAA junction resolve
    # through ff14SB's tables at apply time - bonds/angles/torsions
    # entirely within the backbone (e.g. N-CA-C-N) are owned by ff14SB
    # and dropped from our XML; entries crossing into the sidechain (e.g.
    # CA-CB) keep ff14SB classes on the backbone side and synthesised
    # OFF classes on the sidechain side, so OpenMM resolves them from
    # our cluster Interchange.
    #
    # Mirrors the antechamber pipeline's backbone_at mapping at
    # ``nonstandard.py:1386``.
    # Per-name (ff14SB TYPE, ff14SB CLASS) pairs. The type goes into the
    # <Residue><Atom type=...> attribute; the class is what force entries
    # use to look up parameters. ff14SB declares e.g. type "protein-CT"
    # WITH class "CT" (no "protein-" prefix on the class), so emitting
    # class1="protein-CT" in HarmonicBondForce would never match. Force
    # entries must use the ff14SB class.
    _BACKBONE_TO_FF14SB = {
        "N":   ("protein-N",   "N"),
        "H":   ("protein-H",   "H"),
        "CA":  ("protein-CT",  "CT"),
        "HA":  ("protein-H1",  "H1"),
        "C":   ("protein-C",   "C"),
        "O":   ("protein-O",   "O"),
        "OXT": ("protein-O2",  "O2"),
    }
    type_overrides = {}
    class_overrides = {}
    for ri, atom_indices in enumerate(residue_atom_indices):
        if not atom_indices:
            continue
        is_chain = (
            ri < len(cluster_model.spec.is_chain_resident)
            and cluster_model.spec.is_chain_resident[ri]
        )
        if not is_chain:
            continue
        for ai in atom_indices:
            orig_name = original_names[ai]
            pair = _BACKBONE_TO_FF14SB.get(orig_name)
            if pair is not None:
                type_overrides[ai] = pair[0]
                class_overrides[ai] = pair[1]

    # Non-peptide crosslinks within the cluster (scaffolded peptides,
    # covalent ligands): each ClusterBond connects two atoms living in
    # different cluster residues. The bond itself is INSIDE the cluster
    # Interchange (parameterised in context), but for OpenMM to apply it
    # at build time we need an ExternalBond marker on each side. The user
    # then has to add the actual bond via topology.addBond - same
    # workflow as amber.build's custombonds=... arg.
    res_key_to_cluster_idx = {
        (r.segid, r.chain, int(r.resid), r.insertion): ri
        for ri, r in enumerate(cluster_model.spec.residues)
    }
    for cluster_bond in cluster_model.spec.bonds:
        for end_atom in (cluster_bond.atom_a, cluster_bond.atom_b):
            key = (
                end_atom.segid,
                end_atom.chain,
                int(end_atom.resid),
                end_atom.insertion,
            )
            ri = res_key_to_cluster_idx.get(key)
            if ri is None:
                # Bond endpoint is outside the cluster - shouldn't happen
                # for ClusterSpec.bonds (those are by construction within
                # the cluster), but we skip defensively.
                continue
            cluster_residue = cluster_model.spec.residues[ri]
            resname = cluster_model.canonical_renames.get(
                ri, cluster_residue.resname
            )
            external_bonds.setdefault(resname, []).append(end_atom.name)

    xml_path = os.path.join(cluster_outdir, f"cluster_{cluster_index:03d}.xml")
    _emit_openmm_xml_from_cluster_interchange(
        ic=ic,
        residue_assignment=residue_assignment,
        cap_atoms=cap_atoms,
        external_bonds=external_bonds,
        xml_path=xml_path,
        class_prefix=f"OFF_C{cluster_index:03d}",
        atom_names=original_names,
        class_overrides=class_overrides,
        type_overrides=type_overrides,
    )
    interchange_xml_paths.append(xml_path)
    return xml_path


def _parameterize_free_residue_openff(
    mol,
    group,
    ligand_dir,
    *,
    forcefield,
    charge_method,
    use_pyodide,
    interchange_xml_paths,
):
    """OpenFF backend for one free residue: build a per-ligand
    :class:`Interchange` via SMIRNOFF and emit it as a self-contained
    OpenMM ``<ForceField>`` XML fragment.

    Unlike the antechamber free-residue path, this does **not** produce
    prepi/cif + frcmod files - the OpenFF backend only supports
    ``htmd.builder.openmm.build``, not ``amber.build``. The emitted XML
    fragment goes into ``ClusterOutputs.xml_paths`` and is consumed via
    ``extra_xml=...`` at build time.
    """
    from htmd.builder.openmm import (
        parameterizeLigandsOpenFF,
        _emit_openmm_xml_from_interchange,
    )

    if use_pyodide:
        raise NotImplementedError(
            "backend='openff' is not yet supported under Pyodide - SMIRNOFF "
            "charge assignment needs RDKit's full chemistry stack which is "
            "not available there. Use backend='antechamber' instead."
        )

    # Slice the residue into its own Molecule (write_free_residue_cif also
    # filters down to the residue's atoms and writes a CIF; we just want
    # the in-memory slice).
    sub_cif = os.path.join(ligand_dir, f"_slice_{group['resname']}.cif")
    sub = _write_free_residue_cif(mol, group, sub_cif)

    ic_map = parameterizeLigandsOpenFF(
        sub,
        ligand_ff=forcefield,
        charge_method=charge_method,
        resnames=[group["resname"]],
    )
    ic = ic_map[group["resname"]]

    xml_path = os.path.join(ligand_dir, f"{group['resname']}.xml")
    _emit_openmm_xml_from_interchange(
        ic,
        resname=group["resname"],
        external_bond_atom_names=[],
        xml_path=xml_path,
    )
    interchange_xml_paths.append(xml_path)
    return xml_path


def prepareClusterResidues(
    typed_path: str,
    frcmod_path: str,
    model: ClusterModel,
    outdir: "str | None" = None,
    use_pyodide: "bool | None" = None,
    residue_templates: "list | None" = None,
    parameter_sets: "list | None" = None,
    pin_backbone_charges: bool = True,
    normalize: "str | None" = "cluster",
) -> ClusterOutputs:
    """Split antechamber output for a cluster model compound into per-
    residue topology files and emit the matching custombonds list.

    For each non-canonical cluster residue the function writes a CIF using
    antechamber's GAFF2 types and the cluster compute's per-atom charges.
    For each canonical anchor the CIF uses the appropriate AMBER residue
    template's ff14SB atom types (mid-chain ``CYX`` / ``NLN`` / ... or
    the matching N- or C-terminal variant ``NCYX`` / ``CCYX`` / ...
    when the residue is at a chain terminus), with per-atom charges from
    the antechamber compute on the combined model. For chain-resident
    residues the backbone charges are pinned to ff14SB by default and
    the residue rebalanced to its integer formal charge (see
    :func:`_backbone_charge_map`); set ``pin_backbone_charges=False``
    to skip the pin and keep the cluster-computed backbone charges. In
    both cases every per-residue file (chain-resident and scaffold) is
    rebalanced to its integer formal charge, so each emitted unit is
    integer-charged. Each canonical residue's bucket resname (assigned
    by detect, e.g. ``CY1``) keeps the residue out of tLeap's built-in
    libraries so our prepi loads instead of the standard template.

    Parameters
    ----------
    typed_path : str
        Antechamber-typed mol2 of the cluster model compound.
    frcmod_path : str
        parmchk2 output for the same model compound.
    model : :class:`ClusterModel`
        Cluster model returned by :func:`buildClusterModel`.
    outdir : str or None, optional
        Output directory; created if missing. If ``None``, a fresh
        temporary directory is used.
    use_pyodide : bool or None, optional
        Force Pyodide dispatch (``True``), native subprocess
        (``False``), or auto-detect (``None``).
    residue_templates : list or None, optional
        If provided, every per-residue typed-mol slice this function
        writes is appended as a :class:`_ResidueTemplateData` entry for
        the downstream OpenMM XML emitter.
    parameter_sets : list or None, optional
        If provided, the cluster's final :class:`AmberParameterSet`
        (post junction-term injection and backbone-rename duplication)
        is appended for the downstream XML emitter.
    pin_backbone_charges : bool, optional
        If ``True``, pin backbone partial charges of chain-resident
        residues to ff14SB values before writing topology files. Set
        ``False`` to keep cluster-computed backbone charges.
    normalize : str or None, optional
        Charge normalization mode applied after backbone pinning.
        ``"cluster"`` distributes any residual across the whole
        cluster; ``"per_residue"`` normalizes each residue
        independently to its integer formal charge; ``None`` leaves
        charges as-is.

    Returns
    -------
    out : :class:`ClusterOutputs`
        Per-residue topology paths, frcmod paths, and custombonds for
        this cluster.
    """
    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="cluster_params_")
    else:
        os.makedirs(outdir, exist_ok=True)

    if use_pyodide is None:
        use_pyodide = _in_pyodide()

    typed_mol = Molecule(typed_path)

    # Per-cluster-residue integer formal charge. The antechamber-typed
    # mol2 carries none (mol2 has no formal-charge field), so read it
    # from the cluster CIF, which does round-trip formal charge, mapping
    # atoms to residues via model.atom_to_residue.
    cif_mol = Molecule(model.cif_path)
    residue_formal_charge = {}
    for i in range(cif_mol.numAtoms):
        cidx = model.atom_to_residue.get(str(cif_mol.name[i]))
        if cidx is None:
            continue
        residue_formal_charge[cidx] = residue_formal_charge.get(cidx, 0.0) + float(
            cif_mol.formalcharge[i]
        )
    residue_formal_charge = {
        cidx: int(round(v)) for cidx, v in residue_formal_charge.items()
    }

    # Group atom indices by cluster_residue_index using model.atom_to_residue.
    # Cap atoms (not in atom_to_residue) are dropped from per-residue files.
    per_residue_atom_idxs = {}
    for i, name in enumerate(typed_mol.name):
        cidx = model.atom_to_residue.get(str(name))
        if cidx is None:
            continue
        per_residue_atom_idxs.setdefault(cidx, []).append(i)

    out = ClusterOutputs()

    # Backbone atom-type rename map applied to chain-resident NCAAs whose
    # backbone atoms picked up GAFF2 types from antechamber. tLeap needs
    # the ff14SB types here to resolve the peptide bond at the chain
    # neighbours.
    backbone_at = {"N": "N", "H": "H", "CA": "CT", "HA": "H1", "C": "C", "O": "O"}
    # Original GAFF2 -> [ff14SB, ...] mapping accumulated across all
    # chain-resident NCAAs whose backbone we retyped. Used to duplicate
    # GAFF2 frcmod entries under the ff14SB-typed backbone so torsions
    # spanning the backbone-sidechain boundary still resolve.
    backbone_type_renames = {}

    # Antechamber path: pin to the ff14SB charges from the AmberTools
    # amino*.lib via parmed - tleap will apply these same values at
    # build time so the pinned charges round-trip exactly.
    ff14sb_lib = _load_ff14sb_amino_lib_amber() if pin_backbone_charges else None

    # Pass 1: build per-residue ``sub`` Molecules and apply the backbone
    # pin (if requested). Defer total-charge normalization to after the
    # loop so the ``normalize=='cluster'`` mode can run a single
    # cluster-wide residual distribution over all the non-pinned atoms.
    per_residue_data = []
    for cidx, residue_id in enumerate(model.spec.residues):
        if cidx not in per_residue_atom_idxs:
            continue
        atom_idxs = np.asarray(per_residue_atom_idxs[cidx], dtype=np.int64)
        sub = typed_mol.copy(sel=atom_idxs)
        new_names = [model.atom_to_orig_name[str(typed_mol.name[i])] for i in atom_idxs]
        sub.name[:] = new_names

        if cidx in model.canonical_renames:
            new_resname = model.canonical_renames[cidx]
        else:
            new_resname = residue_id.resname
        sub.resname[:] = new_resname

        # Chain-resident backbones (canonical and non-canonical alike)
        # must carry ff14SB atom types on N / CA / C / O / H / HA so
        # the peptide bond at the chain neighbour resolves through the
        # ff14SB tables. systemPrepare guarantees the backbone names,
        # so a plain dict suffices - no canonical-template lookup, no
        # sidechain retyping. Sidechain atoms keep antechamber's GAFF2
        # types; the frcmod-duplication pass below extends bond /
        # angle / dihedral entries that span the boundary so they
        # resolve under both the original GAFF2 type and the
        # backbone-rewritten ff14SB type.
        net_charge = residue_formal_charge.get(cidx, 0)
        charge_map = {}
        if model.spec.is_chain_resident[cidx]:
            for name, ff_type in backbone_at.items():
                sel = sub.name == name
                for orig in sub.atomtype[sel]:
                    backbone_type_renames.setdefault(str(orig), set()).add(ff_type)
                sub.atomtype[sel] = ff_type
            # Optionally pin the backbone charges to ff14SB so the
            # modified residue's backbone electrostatics match the
            # standard residues it is peptide-bonded to. Canonical
            # residues get their whole backbone (residue-specific,
            # terminal variants included) from the ff14SB libraries;
            # NCAAs, absent from the libraries, get the universal
            # mid-chain amide charges. Disable via
            # ``pin_backbone_charges=False`` (Forcefield_PTM convention).
            if pin_backbone_charges:
                canonical_resname = (
                    model.spec.canonical_resnames[cidx]
                    if model.spec.is_canonical[cidx]
                    else ""
                )
                terminus = (
                    "n"
                    if model.spec.is_n_term[cidx]
                    else ("c" if model.spec.is_c_term[cidx] else "")
                )
                # Walk sub.bonds for actual H neighbors of N rather
                # than relying on atom-name heuristics (an N-terminal
                # NH3+ has H1/H2/H3, no atom named "H", but it is not
                # proline-like - it has bonded Hs).
                n_idx = np.where(sub.name == "N")[0]
                n_has_bonded_h = False
                if len(n_idx):
                    for nb in sub.getNeighbors(int(n_idx[0])):
                        if str(sub.element[nb]) == "H":
                            n_has_bonded_h = True
                            break
                charge_map = _backbone_charge_map(
                    canonical_resname,
                    terminus,
                    {str(n) for n in sub.name},
                    net_charge,
                    n_has_bonded_h=n_has_bonded_h,
                    ff14sb_lib=ff14sb_lib,
                )

        # Apply the pin but defer the residual rebalance (net_charge=None).
        pinned_mask = _normalize_residue_charges(
            sub, net_charge=None, charge_map=charge_map
        )

        sub.resid[:] = 1
        sub.segid[:] = "A"
        sub.chain[:] = "A"
        sub.insertion[:] = ""

        per_residue_data.append(
            {
                "cidx": cidx,
                "sub": sub,
                "new_resname": new_resname,
                "residue_id": residue_id,
                "pinned_mask": pinned_mask,
                "net_charge": net_charge,
            }
        )

    # Normalization pass. Cluster-wide charge methods (RDKit Gasteiger
    # PEOE, antechamber AM1-BCC) sum the cluster total to the cluster's
    # net formal charge exactly, but two things move per-residue (and
    # potentially cluster) totals off-integer afterwards:
    #   1. The C-S / C-N smear at residue junctions: per-residue slices
    #      get a small fraction of charge that "belongs" to a
    #      neighbouring residue, and vice versa. The cluster sum is
    #      preserved; per-residue sums are not.
    #   2. The backbone pin (when enabled): replacing per-atom Gasteiger
    #      backbone charges with ff14SB values shifts the cluster total
    #      by sum(ff14SB - gasteiger) over all pinned atoms.
    # The "cluster" mode absorbs only (2)'s effect by spreading one
    # residual cluster-wide; "per_residue" mode absorbs both (1) and
    # (2)'s effects into each residue's free atoms; None leaves both.
    #
    #   "cluster" (default): one residual is computed over the whole
    #       cluster and distributed equally across every non-pinned
    #       atom across all residues. The cluster total is integer but
    #       individual residue totals stay at their charge-method
    #       values (modulo a uniform shift). Preserves the natural
    #       per-atom charges better; per-residue totals are fractional.
    #   "per_residue": each residue's free atoms absorb that residue's
    #       residual so the residue sums to its integer formal charge.
    #       AMBER's tLeap convention, used by Betz / R.E.D. / Ramos.
    #       Safer if the same residue might recur in different bonding
    #       contexts.
    #   None: no rebalance at all. Each per-residue file carries the
    #       charge-method's per-atom charges with the backbone-pin
    #       shift left in.
    if normalize == "per_residue":
        for d in per_residue_data:
            sub = d["sub"]
            free = ~d["pinned_mask"]
            n_free = int(free.sum())
            if n_free == 0:
                continue
            residual = float(d["net_charge"]) - float(np.sum(sub.charge))
            sub.charge[free] += residual / n_free
    elif normalize == "cluster":
        cluster_total = sum(float(np.sum(d["sub"].charge)) for d in per_residue_data)
        cluster_target = sum(int(d["net_charge"]) for d in per_residue_data)
        cluster_residual = float(cluster_target) - cluster_total
        n_free_cluster = sum(int((~d["pinned_mask"]).sum()) for d in per_residue_data)
        if n_free_cluster:
            shift = cluster_residual / n_free_cluster
            for d in per_residue_data:
                d["sub"].charge[~d["pinned_mask"]] += shift
    # else normalize is None: leave charges as-is after pin.

    # Pass 2: emit per-residue topology files and (optionally) collect
    # OpenMM XML template data using the now-normalized ``sub``s.
    for d in per_residue_data:
        cidx = d["cidx"]
        sub = d["sub"]
        new_resname = d["new_resname"]
        residue_id = d["residue_id"]

        # Chain-resident residues need a prepi (with HEAD / TAIL / connect
        # atoms) so tLeap can stitch them into their flanking residues
        # via standard peptide bonds. Free / scaffold residues are
        # written as CIF and converted to mol2 internally by amber.build.
        wrote_topo = False
        if model.spec.is_chain_resident[cidx]:
            topo_path = os.path.join(outdir, f"{new_resname}.prepi")
            if not os.path.isfile(topo_path):
                _write_chain_residue_prepi(
                    sub,
                    topo_path,
                    new_resname,
                    is_n_term=model.spec.is_n_term[cidx],
                    is_c_term=model.spec.is_c_term[cidx],
                    use_pyodide=use_pyodide,
                    netcharge=residue_formal_charge.get(cidx, 0),
                )
                out.topo_paths.append(topo_path)
                wrote_topo = True
        else:
            topo_path = os.path.join(outdir, f"{new_resname}.cif")
            if not os.path.isfile(topo_path):
                sub.write(topo_path)
                out.topo_paths.append(topo_path)
                wrote_topo = True

        # Collect a parallel _ResidueTemplateData for the OpenMM XML
        # emitter. ``sub`` carries the antechamber per-atom charges plus
        # the right atom types (ff14SB on the backbone / canonical
        # anchor template; GAFF2 on the sidechain).
        if residue_templates is not None and wrote_topo:
            ext_atoms = []
            if model.spec.is_chain_resident[cidx]:
                if not model.spec.is_n_term[cidx]:
                    ext_atoms.append("N")
                if not model.spec.is_c_term[cidx]:
                    ext_atoms.append("C")
            for cb in model.spec.bonds:
                for atom_id in (cb.atom_a, cb.atom_b):
                    if (
                        str(atom_id.segid) == str(residue_id.segid)
                        and str(atom_id.chain) == str(residue_id.chain)
                        and int(atom_id.resid) == int(residue_id.resid)
                        and str(atom_id.insertion) == str(residue_id.insertion)
                    ):
                        name = str(atom_id.name)
                        if name not in ext_atoms:
                            ext_atoms.append(name)
            residue_templates.append(
                _ResidueTemplateData(
                    resname=new_resname,
                    mol=sub.copy(),
                    external_bond_atom_names=ext_atoms,
                )
            )

    # Load the parmchk2 frcmod, then for canonical-anchored clusters add
    # cross-FF junction terms (bond/angle/dihedral entries spanning a
    # canonical-residue atom and a non-canonical one, with the canonical-
    # side atom types rewritten from antechamber's GAFF2 to ff14SB so
    # tLeap can resolve them when the canonical residue is loaded with
    # ff14SB types) and, for chain-resident NCAAs whose backbone we retyped,
    # duplicate entries under the ff14SB backbone names. Finally prune the
    # parameters none of the kept residues use (parmchk2 emits the full
    # type-combination dump, caps included).
    #
    # Write the result under the basename of EACH per-residue CIF / prepi
    # rather than a single shared name. amber.build's auto-handler for known
    # NCAAs / cofactors / PTMs (``_detect_cofactors_ncaa_ptm``) checks
    # ``param`` basenames against the bundled-resname registry to decide
    # whether to layer in its own prepi - if a cluster residue's resname is
    # in the registry (NLE, MSE, NAG via PTM, ...) but no ``<RESNAME>.frcmod``
    # is in ``param``, the auto-handler kicks in and collides with our
    # cluster prepi. Naming the frcmod ``<RESNAME>.frcmod`` for each cluster
    # residue suppresses the auto-add.
    pset = AmberParameterSet(frcmod_path)
    if backbone_type_renames:
        # Duplicate every frcmod entry that mentions an original GAFF2
        # backbone type under its ff14SB rename so torsions spanning the
        # backbone-sidechain boundary still resolve at build time.
        from htmd.builder._ambertools import _duplicate_parameters

        _duplicate_parameters(
            pset,
            {orig: list(news) for orig, news in backbone_type_renames.items()},
        )
    retyped_mol, cap_mask = _retyped_cluster_mol(typed_mol, model, backbone_at)
    _clean_frcmod_params(pset, retyped_mol, list(backbone_at.values()), cap_mask)

    if parameter_sets is not None:
        parameter_sets.append(copy.deepcopy(pset))

    seen_basenames = set()
    for cif_path in out.topo_paths:
        basename = os.path.splitext(os.path.basename(cif_path))[0]
        if basename in seen_basenames:
            continue
        seen_basenames.add(basename)
        per_res_frcmod = os.path.join(outdir, f"{basename}.frcmod")
        pset.write(
            per_res_frcmod,
            title=f"cluster parameters for {basename}",
            style="frcmod",
        )
        out.frcmod_paths.append(per_res_frcmod)

    # Note: cluster custombonds (one pair per ClusterSpec.bonds entry)
    # are emitted by the parameterizeFromSpecs main loop, not here.
    # That keeps custombond emission per-cluster but independent of
    # whether the cluster's parameterise step ran (the cluster might be
    # deduped if its resname set duplicates an earlier cluster's, but
    # the inter-residue bonds at the junctions still need wiring).

    return out


def _residue_groups_with_index(mol):
    """Walk ``mol`` once and return ``(a2r, groups)`` where ``a2r`` maps
    each atom index to a residue index and ``groups`` is a list of dicts
    with the residue's atoms + identity fields, ordered by residue index."""
    a2r, idx_lists = mol.getResidues(
        fields=("resname", "resid", "insertion", "segid", "chain"),
        return_idx=True,
    )
    groups = []
    for atom_idx in idx_lists:
        first = int(atom_idx[0])
        groups.append(
            {
                "atom_idx": np.asarray(atom_idx, dtype=np.int64),
                "resname": str(mol.resname[first]),
                "resid": int(mol.resid[first]),
                "insertion": str(mol.insertion[first]),
                "segid": str(mol.segid[first]),
                "chain": str(mol.chain[first]),
            }
        )
    return a2r, groups


def _spec_residue_key(uid):
    return (str(uid.segid), str(uid.chain), int(uid.resid), str(uid.insertion))


def _build_internal_cluster_spec(
    member_indices, cluster_bonds, mol, groups, spec_by_res_idx
):
    """Construct a ClusterSpec for the cluster pipeline from the per-residue
    specs that fall inside one connected component of non-peptide bonds.
    A singleton component containing just one chain-resident NCAA is also
    valid - it gets its own 1-residue cluster with peptide-neighbour caps
    drawn from the live mol."""
    from moleculekit.tools.nonstandard_residues import (
        ChainResidueSpec,
        ScaffoldSpec,
        CovalentLigandSpec,
        PROTEIN_RESNAMES,
    )

    residues = []
    is_canonical = []
    is_chain_resident = []
    roles = []
    canonical_resnames = []
    canonical_terminus = []
    is_n_term = []
    is_c_term = []
    has_scaffold = False
    has_canonical = False
    n_chain_resident_nc = 0
    n_nc = 0
    for r_idx in member_indices:
        g = groups[r_idx]
        spec = spec_by_res_idx[r_idx]
        residues.append(
            UniqueResidueID(
                resname=g["resname"],
                chain=g["chain"],
                resid=g["resid"],
                insertion=g["insertion"],
                segid=g["segid"],
            )
        )
        if isinstance(spec, ChainResidueSpec):
            # Role: anchor (canonical at junction) | stapled_ncaa (NCAA at
            # junction) | ncaa (plain NCAA). Canonical-vs-NCAA comes from
            # resname; junction-or-not is on the spec already
            # (anchor_atom is not None).
            is_canonical_origin = spec.resname in PROTEIN_RESNAMES
            if is_canonical_origin:
                is_canonical.append(True)
                is_chain_resident.append(True)
                roles.append("anchor")
                canonical_resnames.append(spec.resname)
                # Trust spec.is_n_term / spec.is_c_term directly: these
                # come from whether the residue has an outgoing peptide bond
                # on each side, which is what matters for the prepi template.
                # An additional atom-presence check (OXT for c-term,
                # H1/H2/H3 for n-term) is too strict because CIF structures
                # often omit OXT from the C-terminal residue, and
                # templateResidueFromSmiles may not have been called on
                # the canonical residue (so H1/H2/H3 might be absent too).
                canonical_terminus.append(
                    "n" if spec.is_n_term else ("c" if spec.is_c_term else "")
                )
                is_n_term.append(spec.is_n_term)
                is_c_term.append(spec.is_c_term)
                has_canonical = True
            else:
                is_canonical.append(False)
                is_chain_resident.append(True)
                roles.append("stapled_ncaa" if spec.anchor_atom is not None else "ncaa")
                canonical_resnames.append("")
                canonical_terminus.append("")
                is_n_term.append(spec.is_n_term)
                is_c_term.append(spec.is_c_term)
                n_chain_resident_nc += 1
                n_nc += 1
        elif isinstance(spec, ScaffoldSpec):
            is_canonical.append(False)
            is_chain_resident.append(False)
            roles.append("scaffold")
            canonical_resnames.append("")
            canonical_terminus.append("")
            is_n_term.append(False)
            is_c_term.append(False)
            has_scaffold = True
            n_nc += 1
        elif isinstance(spec, CovalentLigandSpec):
            is_canonical.append(False)
            is_chain_resident.append(False)
            roles.append("covalent_ligand")
            canonical_resnames.append("")
            canonical_terminus.append("")
            is_n_term.append(False)
            is_c_term.append(False)
            n_nc += 1
        else:
            raise TypeError(
                f"Spec type {type(spec).__name__} is not expected inside a cluster"
            )

    bonds = [
        ClusterBond(
            atom_a=UniqueAtomID.fromMolecule(mol, idx=int(a1)),
            atom_b=UniqueAtomID.fromMolecule(mol, idx=int(a2)),
        )
        for a1, a2, _, _ in cluster_bonds
    ]

    if has_scaffold:
        subtype = "scaffolded_peptide"
    elif has_canonical and n_chain_resident_nc >= 1:
        subtype = "canonical_ncaa_crosslink"
    elif has_canonical:
        subtype = "covalent_ligand"
    elif n_chain_resident_nc == n_nc and n_nc >= 2:
        subtype = "ncaa_crosslink"
    else:
        subtype = "mixed"

    return ClusterSpec(
        subtype=subtype,
        residues=residues,
        is_chain_resident=is_chain_resident,
        is_canonical=is_canonical,
        roles=roles,
        bonds=bonds,
        canonical_resnames=canonical_resnames,
        canonical_terminus=canonical_terminus,
        is_n_term=is_n_term,
        is_c_term=is_c_term,
    )


def _write_free_residue_cif(mol, group, out_path):
    """Slice one residue out of mol and write it as a standalone CIF, with
    resid/segid/chain normalized so antechamber treats it as one molecule."""
    sub = mol.copy(sel=group["atom_idx"])
    sub.resname[:] = group["resname"]
    sub.resid[:] = 1
    sub.segid[:] = "A"
    sub.chain[:] = "A"
    sub.insertion[:] = ""
    sub.write(out_path)
    return sub


# Typical (neutral) valence per element, used for the protonation sanity
# check below. Anything not listed (metals, etc.) contributes no expected
# hydrogens and never trips the check on its own.
_TYPICAL_VALENCE = {
    "H": 1,
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "SI": 4,
    "P": 3,
    "S": 2,
    "CL": 1,
    "AS": 3,
    "SE": 2,
    "BR": 1,
    "I": 1,
}

# Fraction of the all-single-bonds hydrogen count a residue must reach to be
# considered protonated. Even densely fused aromatics keep a real-H /
# sp3-estimate ratio around 0.33 (benzene 0.5, pyrene ~0.38, coronene ~0.33),
# so 0.2 leaves a comfortable margin while still flagging a residue whose
# sidechain hydrogens were never added.
_MIN_HYDROGEN_FRACTION = 0.2

# Don't apply the under-protonation check below this many estimated hydrogens.
# Genuinely H-free small carbon species exist in that range (CO2, CN-, urea,
# hexafluorobenzene, ...) and "stripped of H" is not a realistic failure mode
# for residues that small; real NCAAs / drug-like ligands estimate well above
# it (alanine-sized sidechains already give ~9).
_MIN_EXPECTED_HYDROGENS = 7


def _estimate_sp3_hydrogens(elem_upper, formalcharge, heavy_degree, atom_idx):
    """Upper bound on a residue's hydrogen count, assuming every heavy-heavy
    bond is a single bond: for each heavy atom, ``valence + formal_charge``
    minus the number of bonds it makes to other heavy atoms. Unsaturation
    only lowers the real count, so this stays an over-estimate; the caller
    compares against a small fraction of it."""
    total = 0
    for i in atom_idx:
        el = elem_upper[i]
        if el == "H":
            continue
        val = _TYPICAL_VALENCE.get(el)
        if val is None:
            continue
        eff_val = max(0, val + int(formalcharge[i]))
        total += max(0, eff_val - int(heavy_degree[i]))
    return total


def _check_specs_protonated(mol, spec_by_res_idx, groups):
    """Refuse to parameterize a carbon-bearing residue that is missing its
    hydrogens (entirely, or e.g. a sidechain stripped of H while the backbone
    H survived).

    A residue taken straight from a PDB without explicit hydrogens
    parameterizes into a prepi/frcmod missing those H, which is silently
    wrong. Protonate first - template the residue with
    ``mol.templateResidueFromSmiles(<sel>, <smiles>, addHs=True)`` or
    protonate the whole structure with
    ``moleculekit.tools.preparation.systemPrepare`` - before calling here.

    The test is heuristic: count current H per residue, compare against
    ``_MIN_HYDROGEN_FRACTION`` of the all-single-bonds estimate
    (:func:`_estimate_sp3_hydrogens`). Residues without carbon (monatomic
    ions, oxo-anions such as sulfate/phosphate) and residues whose estimate
    is below ``_MIN_EXPECTED_HYDROGENS`` (small species that are often H-free
    anyway) are left alone.
    """
    if mol.bonds is None or len(mol.bonds) == 0:
        return  # no connectivity to reason about; clustering will fail later

    elem_upper = np.char.upper(np.asarray(mol.element, dtype="U3"))
    is_h = elem_upper == "H"
    bonds = np.asarray(mol.bonds)
    heavy_bonds = bonds[~is_h[bonds[:, 0]] & ~is_h[bonds[:, 1]]]
    heavy_degree = np.zeros(mol.numAtoms, dtype=int)
    np.add.at(heavy_degree, heavy_bonds[:, 0], 1)
    np.add.at(heavy_degree, heavy_bonds[:, 1], 1)
    formalcharge = (
        mol.formalcharge
        if mol.formalcharge is not None and len(mol.formalcharge) == mol.numAtoms
        else np.zeros(mol.numAtoms, dtype=int)
    )

    bad = []
    for r_idx in sorted(spec_by_res_idx):
        g = groups[r_idx]
        atom_idx = np.asarray(g["atom_idx"])
        if "C" not in set(elem_upper[atom_idx]):
            continue
        expected = _estimate_sp3_hydrogens(
            elem_upper, formalcharge, heavy_degree, atom_idx
        )
        if expected < _MIN_EXPECTED_HYDROGENS:
            continue
        h_now = int(is_h[atom_idx].sum())
        if h_now >= math.ceil(_MIN_HYDROGEN_FRACTION * expected):
            continue
        res = (
            f"{g['resname']}:{g['resid']}{g['insertion']}"
            f" (segid {g['segid']!r}, chain {g['chain']!r}): "
            f"{h_now} hydrogen(s) present, ~{expected} expected"
        )
        bad.append(res)
    if bad:
        raise RuntimeError(
            "Residue(s) look under-protonated, refusing to parameterize:\n  "
            + "\n  ".join(bad)
            + "\nAdd hydrogens first - e.g. template the residue with "
            "mol.templateResidueFromSmiles(<sel>, <smiles>, addHs=True), or "
            "protonate the structure with "
            "moleculekit.tools.preparation.systemPrepare - then re-run."
        )


def parameterizeFromSpecs(
    specs: list,
    mol: Molecule,
    outdir: str,
    forcefield: "str | dict" = "gaff2",
    charge_method: str = "am1-bcc",
    am1_path_length: "int | None" = 15,
    pin_backbone_charges: bool = True,
    normalize: "str | None" = "cluster",
    use_pyodide: "bool | None" = None,
) -> ClusterOutputs:
    """Parameterize every non-canonical residue in ``specs`` and return paths
    plus custombonds ready to feed :func:`htmd.builder.amber.build`.

    The function recovers cluster grouping by walking ``mol.bonds`` for
    non-peptide inter-residue bonds and unioning the touching residues.
    Per cluster it builds a combined model compound (full residues +
    ACE/NME-style backbone caps), runs antechamber + parmchk2 once, and
    splits the output into per-residue CIF / frcmod pairs. Free residues
    (no cluster bonds) are parameterized standalone.

    Parameters
    ----------
    specs : list
        Per-residue specs from
        :func:`moleculekit.tools.nonstandard_residues.detectNonStandardResidues`.
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule the specs describe. Must already carry covalent
        bonds (typically the post-``systemPrepare`` molecule).
    outdir : str
        Output directory for all generated CIF / frcmod / XML files.
    forcefield : str or dict, optional
        Force field for the non-canonical atoms. A name starting with
        ``"gaff"`` dispatches through antechamber + parmchk2 and emits
        prepi + frcmod (consumable by :func:`amber.build`) plus a
        combined OpenMM XML; any other string is treated as a SMIRNOFF
        offxml filename (e.g. ``"openff_unconstrained-2.3.0.offxml"``)
        and dispatches through OpenFF Interchange, emitting only
        per-cluster OpenMM XML (consumable by :func:`openmm.build`). A
        dict ``{resname: ff_name, "default": ff_name}`` lets different
        residues use different force fields; mixing within a single
        cluster is not supported (the cluster compound is parameterized
        as one molecule).
    charge_method : str, optional
        Charge model for the non-canonical atoms. Orthogonal to
        ``forcefield`` - every model works with both GAFF and SMIRNOFF
        typing. ``"gasteiger"`` is the recommended choice: fast,
        computed via RDKit, honours the net charge, and is the
        automatic fallback under Pyodide where AM1-BCC's SQM backend is
        unavailable. ``"am1-bcc"`` is more expensive but highly
        accurate. ``"nagl"`` uses the OpenFF NAGL graph neural network
        as an AM1-BCC surrogate (requires PyTorch). ``"resp"`` /
        ``"resp-multiconf"`` fit RESP charges to a Psi4-computed QM
        ESP; requires the private Acellera ``parameterize`` package +
        Psi4. ``"abcg2"`` is AM1-BCC v2, only meaningful with GAFF.
    am1_path_length : int or None, optional
        Maximum path length for AM1-BCC charge equivalence
        determination, passed to antechamber's ``-pl`` flag. Caps
        antechamber's atom-equivalence search so it does not hang on
        cyclic or large molecules. Only used for
        ``charge_method="am1-bcc"`` / ``"abcg2"``; ignored for
        Gasteiger. ``None`` keeps antechamber's own default.
    pin_backbone_charges : bool, optional
        If ``True``, the backbone partial charges of every
        chain-resident residue are pinned to ff14SB (residue-specific
        for canonical residues, charge-class fallback for NCAAs).
        Matches the Robin Betz / R.E.D. / Carlos Ramos tutorial
        convention. Set ``False`` to keep the cluster-computed
        backbone charges (the Forcefield_PTM / Khoury et al. 2014
        convention).
    normalize : str or None, optional
        How to absorb the small per-residue drift left by slicing one
        residue out of a jointly-charged cluster and any shift the
        backbone pin introduces on the cluster total. ``"cluster"``
        (default): only the cluster total is normalized to integer;
        per-residue totals are left at their natural (fractional)
        values. ``"per_residue"``: each emitted unit is
        integer-charged - AMBER's tLeap convention, the safer choice if
        the same residue might recur in different bonding contexts.
        ``None``: no rebalance at all.
    use_pyodide : bool or None, optional
        Force the AmberTools dispatch path. ``True`` dispatches via
        ``antechamber_pyodide.run``; ``False`` uses a native subprocess.
        ``None`` (default) auto-detects Pyodide via ``sys.platform``.

    Returns
    -------
    out : :class:`ClusterOutputs`
        Aggregated topology / parameter files and custombonds for the
        whole system. The ``forcefield`` choice determines what is
        populated:

        - GAFF forcefield: ``out.topo_paths`` (prepi) + ``out.frcmod_paths``
          (frcmod, for :func:`amber.build`) plus one combined OpenMM XML
          (``gaff_combined.xml``) appended to ``out.xml_paths``.
        - SMIRNOFF forcefield: per-cluster XML fragments appended to
          ``out.xml_paths`` only; ``topo_paths`` / ``frcmod_paths`` stay
          empty.
        - Mixed: both above contribute. Synthetic atom-type names are
          globally unique by construction, so the XMLs load together
          via ``openmm.app.ForceField(*defaultFf(), *out.xml_paths)``.

        ``out.custombonds`` is populated in every case and matches the
        ``custombonds=`` argument of :func:`amber.build` /
        :func:`openmm.build`. To tell whether GAFF was involved (and
        therefore whether :func:`amber.build` is also viable), check
        ``frcmod_paths``: it is non-empty iff at least one residue went
        through the GAFF path.

    Examples
    --------
    Build a scaffolded cyclic peptide (3 cysteines thioether-bonded to a
    triazinane scaffold ``LFI``)::

        from moleculekit.molecule import Molecule
        from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
        from moleculekit.tools.preparation import systemPrepare
        from htmd.builder.nonstandard import parameterizeFromSpecs
        from htmd.builder import amber

        mol = Molecule("8QFZ.pdb")
        mol.filter("chain B")
        mol.segid[:] = "P"
        mol.segid[mol.resname == "LFI"] = "L"

        # 1. Inspect the molecule and decide what needs custom params.
        specs = detectNonStandardResidues(mol)

        # 2. Template each non-canonical residue from a SMILES string.
        mol.templateResidueFromSmiles(
            "resname LFI",
            "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr",
            addHs=True,
        )

        # 3. Protonate the canonical part and apply the spec renames /
        #    displaced-H drops in one step.
        pmol, _ = systemPrepare(mol, detect_specs=specs)

        # 4. Run antechamber per cluster and split per-residue using
        #    Gasteiger charges (fast and reliable).
        out = parameterizeFromSpecs(
            specs, pmol, outdir="./params", charge_method="gasteiger"
        )

        # 5. Build.
        built = amber.build(
            pmol,
            outdir="./build",
            custombonds=out.custombonds,
            topo=out.topo_paths,
            param=out.frcmod_paths,
        )
    """
    from moleculekit.tools.nonstandard_residues import (
        ChainResidueSpec,
        ScaffoldSpec,
        CovalentLigandSpec,
        LigandSpec,
        PROTEIN_RESNAMES,
    )

    if normalize not in NORMALIZE_MODES:
        raise ValueError(f"normalize={normalize!r}: expected one of {NORMALIZE_MODES}")
    if pin_backbone_charges and normalize is None:
        raise ValueError(
            "pin_backbone_charges=True with normalize=None is inconsistent: "
            "pinning shifts the cluster total off-integer (by "
            "sum(ff14SB_backbone - charge_method_backbone)), and refusing to "
            "rebalance leaves the system with a non-integer total charge. "
            "Pick one of: pin_backbone_charges=False (raw charges, no shift) "
            "or normalize='cluster' / 'per_residue' (absorb the shift)."
        )

    forcefield_for = _resolve_forcefield(forcefield)

    if use_pyodide is None:
        use_pyodide = _in_pyodide()
    if use_pyodide and str(charge_method).lower() == "am1-bcc":
        logger.info(
            "AM1-BCC needs the SQM backend, unavailable under Pyodide; "
            "falling back to Gasteiger charges."
        )
        charge_method = "gasteiger"

    # SMIRNOFF forcefields are fit alongside a specific charge model
    # (Sage 2.0-2.2 with AM1-BCC, Sage 2.3.0 with AshGC). The vdW + torsion
    # parameters are entangled with the charge values, so using a different
    # charge model gives a thermodynamically inconsistent potential. Warn
    # the user but don't reject - they may have a reason (e.g. comparing
    # backends with the same external charge model for diagnostic tests).
    _warn_if_openff_mismatched_charges(forcefield, charge_method)

    # Non-alpha (beta/gamma) amino-acid backbones build, but their non-alpha
    # backbone atoms get GAFF/OpenFF charges + analogy torsions, which are
    # known-approximate. Warn loudly so the user validates / supplies bespoke
    # torsion params.
    _warn_if_nonalpha_backbone(mol, specs)

    # A modified nucleotide spliced into a nucleic chain has no de-novo
    # parameterization path yet and would silently build disconnected; refuse.
    _assert_no_chain_resident_modified_nucleotide(mol, specs)

    # SMIRNOFF parameters can't be expressed losslessly as AMBER frcmod
    # (e.g. SMIRNOFF supports more torsion periodicities and improper
    # conventions than frcmod can represent), so the openff path emits
    # only OpenMM XML - which means the resulting ClusterOutputs is
    # consumable by ``openmm.build`` but not by ``amber.build``. Log an
    # info reminding the user, since the symptom otherwise ("tleap
    # can't find topology for resname X") is non-obvious.
    if any(
        _engine_for_forcefield(ff) == "openff"
        for ff in _forcefields_in(forcefield)
    ):
        logger.info(
            "SMIRNOFF (openff) forcefield in use; the per-cluster XML "
            "outputs are consumable by htmd.builder.openmm.build only "
            "(via extra_xml=). Use htmd.builder.amber.build only when "
            "every NCAA is parameterised with GAFF/GAFF2."
        )

    os.makedirs(outdir, exist_ok=True)
    dispatch = _detect_clusters(mol, specs)

    out = ClusterOutputs()
    # Accumulators for the OpenMM XML emitter: one entry per per-residue
    # topology file written (resname / typed mol slice / external-bond atom
    # names), one entry per cluster / free-residue AmberParameterSet. The
    # XML emitter dedupes templates by resname and parameter entries by class
    # tuple, so accumulating raw entries here is fine.
    residue_templates = []
    parameter_sets = []

    # Two-pass cluster handling. Pass 1 walks every cluster, extends
    # ``out.custombonds`` with each junction bond, and decides which
    # clusters still need parameterising. Pass 2 runs the actual
    # parameterisation on the deduplicated work list. Splitting the two
    # keeps custombond emission (every junction always counts) cleanly
    # separated from template emission (subject to two dedup gates).
    #
    # Dedup gates filter the Pass 2 work list:
    #
    #   * Singleton chain-NCAA: a chain-resident NCAA cluster with no
    #     cluster bonds and matching (resname, is_n_term, is_c_term) is
    #     the same model compound as a previously-seen one. tLeap loads
    #     one unit per resname anyway, and OpenFF would emit a duplicate
    #     <Residue> entry.
    #
    #   * Multi-residue cluster: a cluster whose resname set is a subset
    #     of any previously-emitted cluster's would clash at
    #     ``app.ForceField`` load with "template <name> with the same
    #     override level already exists". 1R1J's three NAG-Asn
    #     glycosylation sites are the typical case (same {XX*, NAG}
    #     resname set, distinct ND2-C1 junctions each emitted as its
    #     own custombond pair).
    cluster_tasks = []
    seen_singleton_keys = set()
    seen_cluster_resname_sets = []

    for ci, (members, cbonds) in enumerate(dispatch.cluster_membership):
        cluster_spec = _build_internal_cluster_spec(
            members, cbonds, mol, dispatch.groups, dispatch.spec_by_res_idx
        )
        out.custombonds.extend(
            (_atom_sel_for_unique(b.atom_a), _atom_sel_for_unique(b.atom_b))
            for b in cluster_spec.bonds
        )

        if len(members) == 1 and not cbonds:
            spec = dispatch.spec_by_res_idx[members[0]]
            if (
                isinstance(spec, ChainResidueSpec)
                and spec.resname not in PROTEIN_RESNAMES
            ):
                g = dispatch.groups[members[0]]
                key = (
                    g["resname"],
                    bool(spec.is_n_term),
                    bool(spec.is_c_term),
                )
                if key in seen_singleton_keys:
                    continue
                seen_singleton_keys.add(key)

        cluster_resnames = frozenset(
            dispatch.groups[m]["resname"] for m in members
        )
        if any(cluster_resnames <= seen for seen in seen_cluster_resname_sets):
            continue
        seen_cluster_resname_sets.append(cluster_resnames)

        cluster_tasks.append((ci, members, cluster_spec))

    for ci, members, cluster_spec in cluster_tasks:
        cluster_outdir = os.path.join(outdir, f"cluster_{ci:03d}")
        os.makedirs(cluster_outdir, exist_ok=True)

        cluster_model = buildClusterModel(mol, cluster_spec, cluster_outdir)
        # Forcefield is taken from the FIRST member's resname; mixing
        # GAFF and SMIRNOFF within a single multi-residue cluster is not
        # supported (the cluster compound is one molecule, parameterised
        # as a unit).
        cluster_engine, cluster_ff = forcefield_for(
            dispatch.groups[members[0]]["resname"]
        )
        if cluster_engine == "antechamber":
            cluster_out = _parameterize_cluster_antechamber(
                cluster_model,
                cluster_outdir,
                forcefield=cluster_ff,
                charge_method=charge_method,
                am1_path_length=am1_path_length,
                use_pyodide=use_pyodide,
                pin_backbone_charges=pin_backbone_charges,
                normalize=normalize,
                residue_templates=residue_templates,
                parameter_sets=parameter_sets,
            )
            out.topo_paths.extend(cluster_out.topo_paths)
            out.frcmod_paths.extend(cluster_out.frcmod_paths)
        else:  # openff
            _parameterize_cluster_openff(
                cluster_model,
                cluster_outdir,
                forcefield=cluster_ff,
                charge_method=charge_method,
                use_pyodide=use_pyodide,
                pin_backbone_charges=pin_backbone_charges,
                normalize=normalize,
                interchange_xml_paths=out.xml_paths,
                cluster_index=ci,
            )

    # Non-chain residues that ended up outside any cluster: run antechamber
    # standalone and write a per-resname CIF + frcmod pair. ``LigandSpec``
    # is the natural case (no covalent bonds at all); ``ScaffoldSpec`` /
    # ``CovalentLigandSpec`` can also land here when their only inter-
    # residue bonds go to canonical residues that aren't spec'd themselves
    # (e.g. a metal-coordinating ligand whose only neighbours are zinc
    # ions). Chain-resident specs without a cluster shouldn't happen
    # because those are promoted to singleton clusters earlier; if one
    # slips through it's an input bug.
    standalone_types = (LigandSpec, ScaffoldSpec, CovalentLigandSpec)
    seen_ligand_resnames = set()
    for r_idx, spec in dispatch.spec_by_res_idx.items():
        if r_idx in dispatch.in_cluster_set:
            continue
        g = dispatch.groups[r_idx]
        if not isinstance(spec, standalone_types):
            raise RuntimeError(
                f"Spec {type(spec).__name__} for residue {spec.residue} expected "
                f"to sit inside a cluster but no inter-residue bond was found "
                f"in mol.bonds. Check that mol carries the relevant covalent "
                f"bonds."
            )
        # Identical free residues (same resname) yield identical antechamber
        # output, and tLeap only loads one unit per resname. Skip the
        # redundant runs.
        if g["resname"] in seen_ligand_resnames:
            continue
        seen_ligand_resnames.add(g["resname"])
        ligand_dir = os.path.join(outdir, f"ligand_{g['resname']}")
        os.makedirs(ligand_dir, exist_ok=True)
        engine, ff_name = forcefield_for(g["resname"])
        if engine == "antechamber":
            out_cif, out_frcmod = _parameterize_free_residue_antechamber(
                mol,
                g,
                ligand_dir,
                outdir,
                forcefield=ff_name,
                charge_method=charge_method,
                am1_path_length=am1_path_length,
                use_pyodide=use_pyodide,
                normalize=normalize,
                residue_templates=residue_templates,
                parameter_sets=parameter_sets,
            )
            out.topo_paths.append(out_cif)
            out.frcmod_paths.append(out_frcmod)
        else:  # openff
            _parameterize_free_residue_openff(
                mol,
                g,
                ligand_dir,
                forcefield=ff_name,
                charge_method=charge_method,
                use_pyodide=use_pyodide,
                interchange_xml_paths=out.xml_paths,
            )

    if residue_templates:
        # Single combined XML for everything the GAFF/antechamber path
        # produced - same parameters as topo_paths + frcmod_paths, just
        # in OpenMM ForceField format. Named ``gaff_combined.xml`` so the
        # provenance is obvious from the filename; appended to
        # ``xml_paths`` so callers don't have to handle GAFF vs SMIRNOFF
        # XML emission separately (``frcmod_paths`` already tells them
        # whether GAFF was involved).
        gaff_xml_path = os.path.join(outdir, "gaff_combined.xml")
        _emit_openmm_xml(residue_templates, parameter_sets, gaff_xml_path)
        out.xml_paths.append(gaff_xml_path)

    return out


# ---------------------------------------------------------------------------
# Cluster model compound builder.
#
# Combines the cluster's full residues with ACE/NME-equivalent backbone caps
# (trimmed neighbour residues + idealized methyl Hs on the trimmed CA) into
# a single Molecule, ready to feed antechamber + parmchk2. The split back
# into per-residue topology files happens in :func:`prepareClusterResidues`.
# ---------------------------------------------------------------------------


def _idealized_methyl_positions(center, neighbor_pos, bond_length=1.09):
    """Return three idealized H positions around an sp3 ``center`` whose
    only other neighbour is at ``neighbor_pos``. Positions are placed
    tetrahedrally at ``bond_length``, opposite the neighbour direction."""
    direction = neighbor_pos - center
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        direction = np.array([1.0, 0.0, 0.0])
    else:
        direction /= norm
    ref = (
        np.array([1.0, 0.0, 0.0])
        if abs(direction[0]) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )
    x = np.cross(direction, ref)
    x /= np.linalg.norm(x)
    y = np.cross(direction, x)
    cos_t = np.cos(np.deg2rad(180 - 109.471))
    sin_t = np.sin(np.deg2rad(180 - 109.471))
    out = []
    for k in range(3):
        phi = 2 * np.pi * k / 3
        d = cos_t * (-direction) + sin_t * (np.cos(phi) * x + np.sin(phi) * y)
        out.append(center + bond_length * d)
    return out


def _peptide_neighbour_residue_idx(mol, a2r, residue_atoms, side):
    """Return the residue index of the chain neighbour on the given side
    (``"N"`` for the previous residue, ``"C"`` for the next one). Walks
    ``mol.bonds`` first; if no bond is found (sparse CIF inputs frequently
    omit peptide bonds), falls back to a distance-based search for an atom
    of the expected name within 1.6 A of the backbone atom on the relevant
    side. Returns ``None`` if no neighbour can be located."""
    own_atoms = set(int(a) for a in residue_atoms)
    other_name = "C" if side == "N" else "N"
    backbone_atom_idxs = [
        int(a) for a in residue_atoms if str(mol.name[int(a)]) == side
    ]
    for atom_idx in backbone_atom_idxs:
        for nb in mol.getNeighbors(atom_idx):
            nb = int(nb)
            if nb in own_atoms:
                continue
            if str(mol.name[nb]) == other_name:
                return int(a2r[nb])
    for atom_idx in backbone_atom_idxs:
        atom_pos = mol.coords[atom_idx, :, mol.frame]
        for cand_idx in range(mol.numAtoms):
            if cand_idx in own_atoms:
                continue
            if str(mol.name[cand_idx]) != other_name:
                continue
            cand_pos = mol.coords[cand_idx, :, mol.frame]
            if float(np.linalg.norm(cand_pos - atom_pos)) < 1.6:
                return int(a2r[cand_idx])
    return None


def _short_atom_name(element, global_idx, used_names):
    """Build a unique ``<=4`` char antechamber-friendly atom name. Two-char
    elements are title-cased (``Cl``, ``Br``) so antechamber's case-sensitive
    element parser recognises them."""
    el = str(element)[:2]
    el = el[0].upper() + el[1:].lower()
    base = f"{el}{global_idx}"
    if len(base) <= 4 and base not in used_names:
        used_names.add(base)
        return base
    for k in range(10000):
        cand = f"{el[0]}{global_idx}{k:x}"[:4]
        if cand not in used_names:
            used_names.add(cand)
            return cand
    raise RuntimeError(f"Could not generate unique atom name for index {global_idx}")


def _residue_groups_for_model(mol):
    """Return ``(a2r, groups)`` parallel to mol.getResidues, with groups
    carrying ``atom_idx`` plus identity fields. Used inside the cluster
    model builder; mirrors moleculekit's internal helper."""
    a2r, idx_lists = mol.getResidues(
        fields=("resname", "resid", "insertion", "segid", "chain"),
        return_idx=True,
    )
    groups = []
    for atom_idx in idx_lists:
        first = int(atom_idx[0])
        groups.append(
            {
                "atom_idx": np.asarray(atom_idx, dtype=np.int64),
                "resname": str(mol.resname[first]),
                "resid": int(mol.resid[first]),
                "insertion": str(mol.insertion[first]),
                "segid": str(mol.segid[first]),
                "chain": str(mol.chain[first]),
            }
        )
    return a2r, groups


def _build_cluster_model_accurate(mol, spec, cif_path):
    """Build the model compound for combined-antechamber parameterization.

    For each chain-resident cluster residue the immediate peptide
    neighbour on each side is trimmed to ACE/NME-equivalent atoms
    (prev.C, prev.O, prev.CA on the N-side; next.N, next.H, next.CA on
    the C-side), and the trimmed CA is converted into a methyl by adding
    three idealized hydrogens. The cluster bonds are added explicitly if
    the slice didn't preserve them; hydrogens displaced by those bonds
    are assumed to be already removed upstream (canonical anchors via
    :func:`detectNonStandardResidues`, NCAA sidechains via
    :meth:`Molecule.templateResidueFromSmiles`).

    Returns ``(model_mol, atom_map, atom_to_residue, atom_to_orig_name)``.
    """
    a2r, groups = _residue_groups_for_model(mol)

    orig_to_cluster_idx = {}
    for cidx, residue_id in enumerate(spec.residues):
        mask = np.ones(mol.numAtoms, dtype=bool)
        for f in ("resname", "resid", "insertion", "segid", "chain"):
            mask &= getattr(mol, f) == getattr(residue_id, f)
        atom_idxs = np.where(mask)[0]
        if len(atom_idxs) == 0:
            raise ValueError(f"Residue {residue_id} not found in mol")
        for ai in atom_idxs:
            orig_to_cluster_idx[int(ai)] = cidx

    cap_atoms = {}  # {cap_residue_idx_in_mol: {orig_idx: cap_role_str}}

    def _atom_in_residue(residue_atom_idxs, name):
        for ai in residue_atom_idxs:
            if str(mol.name[int(ai)]) == name:
                return int(ai)
        return None

    for cidx, residue_id in enumerate(spec.residues):
        if not spec.is_chain_resident[cidx]:
            continue
        cluster_atom_idxs = np.asarray(
            [ai for ai, c in orig_to_cluster_idx.items() if c == cidx],
            dtype=np.int64,
        )
        prev_r = _peptide_neighbour_residue_idx(mol, a2r, cluster_atom_idxs, "N")
        if prev_r is not None:
            prev_atoms = groups[prev_r]["atom_idx"]
            prev_c = _atom_in_residue(prev_atoms, "C")
            prev_o = _atom_in_residue(prev_atoms, "O")
            prev_ca = _atom_in_residue(prev_atoms, "CA")
            cap = cap_atoms.setdefault(prev_r, {})
            if prev_c is not None:
                cap[prev_c] = "C"
            if prev_o is not None:
                cap[prev_o] = "O"
            if prev_ca is not None:
                cap[prev_ca] = "CH3"
        next_r = _peptide_neighbour_residue_idx(mol, a2r, cluster_atom_idxs, "C")
        if next_r is not None:
            next_atoms = groups[next_r]["atom_idx"]
            next_n = _atom_in_residue(next_atoms, "N")
            next_h = _atom_in_residue(next_atoms, "H")
            next_ca = _atom_in_residue(next_atoms, "CA")
            cap = cap_atoms.setdefault(next_r, {})
            if next_n is not None:
                cap[next_n] = "N"
            if next_h is not None:
                cap[next_h] = "H"
            if next_ca is not None:
                cap[next_ca] = "CH3"

    cap_orig_idxs = {ai for cap in cap_atoms.values() for ai in cap}
    all_orig_idxs = sorted(set(orig_to_cluster_idx) | cap_orig_idxs)
    model = mol.copy(sel=np.asarray(all_orig_idxs, dtype=np.int64))
    orig_to_new = {orig: new for new, orig in enumerate(all_orig_idxs)}

    role_per_pos = []
    cluster_residue_per_pos = []
    orig_name_per_pos = []
    cap_role_per_pos = []
    flat_cap_role = {ai: role for cap in cap_atoms.values() for ai, role in cap.items()}
    for orig in all_orig_idxs:
        if orig in orig_to_cluster_idx:
            role_per_pos.append("residue")
            cluster_residue_per_pos.append(orig_to_cluster_idx[orig])
            orig_name_per_pos.append(str(mol.name[orig]))
            cap_role_per_pos.append(None)
        else:
            role_per_pos.append("cap")
            cluster_residue_per_pos.append(-1)
            orig_name_per_pos.append(str(mol.name[orig]))
            cap_role_per_pos.append(flat_cap_role[orig])

    # Drop intra-cap bonds we don't want; only C-O, C-CH3, N-H, N-CH3 are
    # legitimate edges of an ACE/NME-style cap.
    keep_bonds, keep_bondtypes = [], []
    for bi, (a, b) in enumerate(model.bonds):
        a, b = int(a), int(b)
        if role_per_pos[a] == "cap" and role_per_pos[b] == "cap":
            allowed = {
                ("C", "O"),
                ("O", "C"),
                ("C", "CH3"),
                ("CH3", "C"),
                ("N", "H"),
                ("H", "N"),
                ("N", "CH3"),
                ("CH3", "N"),
            }
            if (cap_role_per_pos[a], cap_role_per_pos[b]) not in allowed:
                continue
        keep_bonds.append([a, b])
        keep_bondtypes.append(model.bondtype[bi])
    model.bonds = (
        np.asarray(keep_bonds, dtype=np.uint32)
        if keep_bonds
        else np.zeros((0, 2), dtype=np.uint32)
    )
    model.bondtype = np.asarray(keep_bondtypes)

    # Add cap topology bonds explicitly from the known role mapping. The
    # slice from mol may not have brought these bonds along - e.g. when
    # the neighbour residue's amide H was added by PDB2PQR after capture
    # so the N-H bond never landed in mol.bonds. The cap topology is
    # deterministic (NME: N-CH3, N-H; ACE: C-CH3, C=O), so adding the
    # bonds by role is not bond guessing.
    cap_bond_pairs = (
        ("N", "H", "1"),
        ("N", "CH3", "1"),
        ("C", "O", "2"),
        ("C", "CH3", "1"),
    )
    existing_pairs = {tuple(sorted((int(a), int(b)))) for a, b in model.bonds}
    new_cap_bonds, new_cap_bondtypes = [], []
    for cap_orig_dict in cap_atoms.values():
        role_to_new = {role: orig_to_new[orig] for orig, role in cap_orig_dict.items()}
        for role_a, role_b, btype in cap_bond_pairs:
            if role_a not in role_to_new or role_b not in role_to_new:
                continue
            a, b = role_to_new[role_a], role_to_new[role_b]
            key = tuple(sorted((int(a), int(b))))
            if key in existing_pairs:
                continue
            new_cap_bonds.append([a, b])
            new_cap_bondtypes.append(btype)
            existing_pairs.add(key)
    if new_cap_bonds:
        model.bonds = np.vstack(
            [model.bonds, np.asarray(new_cap_bonds, dtype=np.uint32)]
        )
        model.bondtype = np.hstack([model.bondtype, np.asarray(new_cap_bondtypes)])

    for bond in spec.bonds:
        a_orig = int(bond.atom_a.selectAtom(mol))
        b_orig = int(bond.atom_b.selectAtom(mol))
        a_new, b_new = orig_to_new[a_orig], orig_to_new[b_orig]
        already = any(
            (int(x) == a_new and int(y) == b_new)
            or (int(x) == b_new and int(y) == a_new)
            for x, y in model.bonds
        )
        if not already:
            model.addBond(a_new, b_new, "1")

    # Add a synthetic amide H on each cap N that lacks one (sparse CIF
    # inputs sometimes come in without backbone Hs).
    extra_amide_h = []
    for i in range(model.numAtoms):
        if cap_role_per_pos[i] != "N":
            continue
        if any(str(model.element[int(nb)]) == "H" for nb in model.getNeighbors(i)):
            continue
        n_pos = model.coords[i, :, model.frame]
        heavy_dirs = []
        for nb in model.getNeighbors(i):
            nb = int(nb)
            if str(model.element[nb]) != "H":
                v = model.coords[nb, :, model.frame] - n_pos
                v /= max(np.linalg.norm(v), 1e-6)
                heavy_dirs.append(v)
        if not heavy_dirs:
            continue
        avg = np.sum(heavy_dirs, axis=0)
        if np.linalg.norm(avg) < 1e-6:
            continue
        h_dir = -avg / np.linalg.norm(avg)
        extra_amide_h.append((i, n_pos + h_dir * 1.01))

    if extra_amide_h:
        synth = Molecule().empty(len(extra_amide_h))
        synth.name[:] = "H"
        synth.element[:] = "H"
        synth.resname[:] = "MOD"
        synth.resid[:] = 1
        synth.segid[:] = "A"
        synth.chain[:] = "A"
        synth.coords = np.array([h_pos for _, h_pos in extra_amide_h])[
            :, :, np.newaxis
        ].astype(np.float32)
        synth.record[:] = "HETATM"
        n_before = model.numAtoms
        model.append(synth, collisions=False)
        for k, (n_idx, _) in enumerate(extra_amide_h):
            model.addBond(n_idx, n_before + k, "1")
            role_per_pos.append("cap")
            cluster_residue_per_pos.append(-1)
            orig_name_per_pos.append("H")
            cap_role_per_pos.append("amideH")

    # Cap each CH3 atom (the trimmed CA) with three methyl Hs.
    methyl_atoms = []
    for i in range(model.numAtoms):
        if cap_role_per_pos[i] != "CH3":
            continue
        ch3_pos = model.coords[i, :, model.frame].copy()
        neighbour_idx = next(
            (
                int(nb)
                for nb in model.getNeighbors(i)
                if str(model.element[int(nb)]) != "H"
            ),
            None,
        )
        if neighbour_idx is None:
            continue
        nb_pos = model.coords[neighbour_idx, :, model.frame]
        for h_pos in _idealized_methyl_positions(ch3_pos, nb_pos):
            methyl_atoms.append((i, h_pos))

    if methyl_atoms:
        synth = Molecule().empty(len(methyl_atoms))
        synth.name[:] = "H"
        synth.element[:] = "H"
        synth.resname[:] = "MOD"
        synth.resid[:] = 1
        synth.segid[:] = "A"
        synth.chain[:] = "A"
        synth.coords = np.array([h_pos for _, h_pos in methyl_atoms])[
            :, :, np.newaxis
        ].astype(np.float32)
        synth.record[:] = "HETATM"
        n_before = model.numAtoms
        model.append(synth, collisions=False)
        for k, (ch3_idx, _) in enumerate(methyl_atoms):
            model.addBond(ch3_idx, n_before + k, "1")
            role_per_pos.append("cap")
            cluster_residue_per_pos.append(-1)
            orig_name_per_pos.append("H")
            cap_role_per_pos.append("methylH")

    model.resname[:] = "MOD"
    model.resid[:] = 1
    model.segid[:] = "A"
    model.chain[:] = "A"
    model.insertion[:] = ""
    model.record[:] = "HETATM"

    used_names = set()
    new_names = []
    atom_map = {}
    atom_to_residue = {}
    atom_to_orig_name = {}
    for i in range(model.numAtoms):
        elem = str(model.element[i]) or str(model.name[i])[:1]
        new_name = _short_atom_name(elem, i + 1, used_names)
        new_names.append(new_name)
        atom_map[new_name] = ModelAtom(role=role_per_pos[i], ff_type=None)
        if role_per_pos[i] == "residue":
            atom_to_residue[new_name] = cluster_residue_per_pos[i]
            atom_to_orig_name[new_name] = orig_name_per_pos[i]
    model.name[:] = new_names

    model.write(cif_path)
    return model, atom_map, atom_to_residue, atom_to_orig_name


def buildClusterModel(mol: Molecule, spec: ClusterSpec, outdir: str) -> ClusterModel:
    """Build the combined model compound for a :class:`ClusterSpec`.

    Assembles full cluster residues plus ACE/NME-style backbone caps
    derived from the live molecule's chain neighbours, then writes the
    result as a CIF ready to feed to antechamber.

    Parameters
    ----------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>`
        The molecule containing the cluster residues and their chain
        neighbours. Must carry covalent bonds.
    spec : :class:`ClusterSpec`
        Specification of the cluster to model (residues, chain
        residency flags, canonical flags, inter-residue bonds).
    outdir : str
        Output directory where the model CIF will be written.

    Returns
    -------
    model : :class:`ClusterModel`
        Cluster model compound with atom-name maps needed by
        :func:`prepareClusterResidues` to split antechamber output back
        into per-residue topology files.
    """
    os.makedirs(outdir, exist_ok=True)
    nc_resnames = "_".join(
        sorted(r.resname for r, c in zip(spec.residues, spec.is_canonical) if not c)
    )
    cif_path = os.path.join(outdir, f"cluster_{spec.subtype}_{nc_resnames}.cif")
    model, atom_map, atom_to_residue, atom_to_orig_name = _build_cluster_model_accurate(
        mol, spec, cif_path
    )
    canonical_renames = {
        cidx: spec.residues[cidx].resname
        for cidx in range(len(spec.residues))
        if spec.is_canonical[cidx]
    }
    return ClusterModel(
        spec=spec,
        cif_path=cif_path,
        atom_map=atom_map,
        atom_to_residue=atom_to_residue,
        atom_to_orig_name=atom_to_orig_name,
        canonical_renames=canonical_renames,
    )
