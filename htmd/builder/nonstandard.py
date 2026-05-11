"""End-to-end parameterization pipeline for non-canonical residues under
AMBER, driven by the spec list returned from
:func:`moleculekit.tools.nonstandard_residues.detectNonStandardResidues`.

:func:`parameterizeFromSpecs` is the user-facing entry point. It walks
``mol.bonds`` to recover cluster grouping (residues sharing non-peptide
inter-residue bonds), builds a combined antechamber model compound per
cluster (full residues + ACE/NME-style backbone caps), runs antechamber +
parmchk2 once per cluster, and splits the output into per-residue
CIF / frcmod pairs. Free residues (no cluster bonds) are parameterized
standalone via :func:`htmd.builder.noncanonical.parameterizeNonCanonicalResidues`
or :func:`htmd.builder.noncanonical._fftype_antechamber`.

The result :class:`ClusterOutputs` carries the topology paths, frcmod
paths, and ``custombonds`` list in the shape that
:func:`htmd.builder.amber.build` expects.

For canonical residues that the detector renamed (CYS bonded to a
scaffold, ASN glycosylated by a sugar, ...), the per-residue CIF carries
ff14SB atom types pulled from the right AMBER residue template (mid-chain
``CYX`` / N-terminal ``NCYX`` / C-terminal ``CCYX`` and the analogous
forms for LYS/HIS/ASN/...) so that backbone bonds resolve against
ff14SB while sidechain charges come from the antechamber compute on the
combined model. The frcmod carries cross-FF junction terms (bond / angle
/ dihedral entries spanning a canonical-residue atom and a non-canonical
one) with the canonical-side atom types rewritten from antechamber's
GAFF2 to ff14SB.
"""

from __future__ import annotations

import copy
import logging
import math
import os
import tempfile
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from parmed.amber import AmberParameterSet

from moleculekit.molecule import Molecule, UniqueAtomID, UniqueResidueID
from moleculekit.tools._anchor_variants import ANCHOR_VARIANTS  # noqa: F401  (re-exported)
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
    ff_type: Optional[str] = None


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
    canonical_original_resnames: list = field(default_factory=list)
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
        if not self.canonical_original_resnames:
            self.canonical_original_resnames = ["" for _ in self.residues]
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


def _live_key(idxs, atom_types, is_stub, ff_types):
    """Build the atom-type tuple a model-compound instance maps to in the
    *live* AMBER system: scaffold atoms keep their GAFF2 types; stub atoms
    are replaced by their canonical-FF type from ANCHOR_VARIANTS. Returns
    ``None`` if any stub atom is missing a canonical-FF mapping."""
    out = []
    for i in idxs:
        if is_stub[i]:
            ff = ff_types[i]
            if not ff:
                return None
            out.append(str(ff))
        else:
            out.append(str(atom_types[i]))
    return tuple(out)


def _splice_into(param_dict, instances, atom_types, is_stub, ff_types):
    """For each instance (an array of atom indices) that spans stub+scaffold,
    look up the parameter value at the GAFF2 key, derive the rewritten
    canonical-FF key, and write the value (copied) under the new key.

    ``param_dict`` is mutated in place. ``instances`` is an ``(N, k)``
    integer array of atom indices (k=2 for bonds, 3 for angles, 4 for
    dihedrals)."""
    for idxs in instances:
        roles = is_stub[idxs]
        if not (roles.any() and (~roles).any()):
            continue  # not a mixed junction instance
        gaff_key = tuple(atom_types[idxs])
        if gaff_key not in param_dict:
            continue
        new_key = _live_key(idxs, atom_types, is_stub, ff_types)
        if new_key is None or new_key in param_dict:
            continue
        # parmed.write deduplicates types by id(); copy so the junction entry
        # is written as a distinct line and not skipped as a "duplicate".
        val_copy = copy.copy(param_dict[gaff_key])
        param_dict[new_key] = val_copy
        param_dict[new_key[::-1]] = val_copy


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


def _in_pyodide():
    """Detect Pyodide so AmberTools dispatch can swap subprocess for
    ``antechamber_pyodide.run`` automatically."""
    import sys

    return sys.platform == "emscripten"


def _write_chain_residue_prepi(
    sub, prepi_path, resname, is_n_term=False, is_c_term=False, use_pyodide=None
):
    """Convert a single chain-resident residue (already typed and named)
    into an AMBER ``.prepi`` topology so tLeap can splice it into a
    peptide chain. Mirrors :func:`htmd.builder.noncanonical._post_process_parameterize`'s
    prepgen step: write mol2, run ``antechamber -fo ac`` to convert it
    to AC format without re-typing, build a mainchain config, then run
    ``prepgen``. ``is_n_term`` / ``is_c_term`` omit the matching
    ``HEAD_NAME`` / ``TAIL_NAME`` lines so tLeap does not try to bond
    that side to a phantom neighbour."""
    import networkx as nx
    from htmd.builder.noncanonical import (
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
                "-i", mol2_name, "-fi", "mol2",
                "-o", ac_name, "-fo", "ac",
                "-nc", str(int(round(float(np.sum(sub.formalcharge))))),
                "-pf", "y", "-dr", "n", "-j", "0", "-an", "n",
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
            f.write(f"CHARGE {float(np.sum(sub.formalcharge)):.1f}\n")

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
                "-i", ac_name,
                "-o", prepi_name, "-f", "prepi",
                "-m", mc_name,
                "-rn", resname,
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



# Cached ff14SB amino-acid lib (atom name -> (type, default charge) per resname).
_FF14SB_AMINO_LIB = None


def _load_ff14sb_amino_lib():
    """Lazy-load the ff14SB amino-acid libraries (mid-chain
    ``amino12.lib``, C-terminal ``aminoct12.lib``, N-terminal
    ``aminont12.lib``) so we can extract canonical atom types per
    residue. Terminal-variant entries (``CCYX``, ``NCYX``, ...) sit
    alongside the mid-chain ones in the returned dict.

    Returns a dict ``{resname: {atom_name: type}}``."""
    global _FF14SB_AMINO_LIB
    if _FF14SB_AMINO_LIB is not None:
        return _FF14SB_AMINO_LIB
    from parmed.amber import AmberOFFLibrary
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
            merged[resname] = {a.name: a.type for a in residue.atoms}
    _FF14SB_AMINO_LIB = merged
    return _FF14SB_AMINO_LIB



def _canonical_variant_template(spec, cidx, ff14sb_amino_lib):
    """Pick the ff14SB atom-name -> type map to use for canonical residue
    ``cidx`` of ``spec``. Looks up the residue's original resname (carried
    on ``spec.canonical_original_resnames``) plus the anchor atom name
    found by walking ``spec.bonds`` against :data:`ANCHOR_VARIANTS` to
    decide which deprotonated variant template applies (CYS+SG -> CYX,
    ASN+ND2 -> NLN, ...). When the residue is at a chain terminus
    (``spec.canonical_terminus[cidx]`` is ``"n"`` or ``"c"``) the
    AMBER terminal variant (``NCYX`` / ``CCYX`` / ...) is preferred so
    OXT and the extra N-terminal hydrogens get the right ff14SB types."""
    from moleculekit.tools._anchor_variants import lookup_anchor_variant

    residue_id = spec.residues[cidx]
    original_resname = (
        spec.canonical_original_resnames[cidx]
        if cidx < len(spec.canonical_original_resnames)
        else ""
    )
    terminus = (
        spec.canonical_terminus[cidx]
        if cidx < len(spec.canonical_terminus)
        else ""
    )
    variant_resname = None
    for bond in spec.bonds:
        for atom_id in (bond.atom_a, bond.atom_b):
            if (
                atom_id.resname == residue_id.resname
                and int(atom_id.resid) == int(residue_id.resid)
                and atom_id.chain == residue_id.chain
                and atom_id.segid == residue_id.segid
            ):
                entry = lookup_anchor_variant(
                    original_resname or atom_id.resname, atom_id.name
                )
                if entry is not None and entry["variant"] is not None:
                    variant_resname = entry["variant"]
                    break
        if variant_resname:
            break

    # Try terminal variants first when the residue is at a chain terminus.
    candidates = []
    if terminus == "n":
        if variant_resname:
            candidates.append("N" + variant_resname)
        if original_resname:
            candidates.append("N" + original_resname)
    elif terminus == "c":
        if variant_resname:
            candidates.append("C" + variant_resname)
        if original_resname:
            candidates.append("C" + original_resname)
    candidates.extend([variant_resname, original_resname, residue_id.resname])
    for cand in candidates:
        if cand and cand in ff14sb_amino_lib:
            return ff14sb_amino_lib[cand]
    return {}


def _cluster_atom_arrays(typed_mol, model, ff14sb_amino_lib):
    """Return ``(gaff2_types, is_canonical_residue, ff14sb_types)`` parallel
    arrays over ``typed_mol`` atoms for cluster-mode junction-frcmod
    rewriting. Cap atoms (not in ``model.atom_to_residue``) get
    ``is_canonical=False`` and empty ``ff14sb_type``. Atoms in canonical
    residues get ``is_canonical=True`` and their ff14SB type from the
    variant template; atoms in non-canonical residues stay GAFF2-typed."""
    spec = model.spec
    n = typed_mol.numAtoms
    gaff2_types = np.asarray(typed_mol.atomtype, dtype=object)
    is_canonical = np.zeros(n, dtype=bool)
    ff14sb_types = np.full(n, "", dtype=object)

    canonical_template = {
        cidx: _canonical_variant_template(spec, cidx, ff14sb_amino_lib)
        for cidx in range(len(spec.residues))
        if spec.is_canonical[cidx]
    }

    for i, name in enumerate(typed_mol.name):
        cidx = model.atom_to_residue.get(str(name))
        if cidx is None:
            continue  # cap atom
        if not spec.is_canonical[cidx]:
            continue  # non-canonical residue
        is_canonical[i] = True
        # Map the model-compound atom name back to the original PDB atom name
        # (via model.atom_to_orig_name) and look up its ff14SB type.
        orig_name = model.atom_to_orig_name.get(str(name))
        if orig_name is None:
            continue
        t = canonical_template[cidx].get(orig_name)
        if t is not None:
            ff14sb_types[i] = t

    return gaff2_types, is_canonical, ff14sb_types


def _add_cluster_junction_terms(pset, typed_mol, model, ff14sb_amino_lib):
    """For every bond/angle/dihedral in ``typed_mol`` that spans both a
    canonical-residue atom and a non-canonical one (or a cap atom),
    rewrite the canonical-side atom types from antechamber GAFF2 to
    ff14SB and add the rewritten entry to the frcmod's parameter dicts.

    This is what makes tLeap able to resolve the canonical-side cluster
    bonds (e.g. the ff14SB ``S`` atom type on CYX SG bonded to a GAFF2
    ``c3`` carbon on the scaffold side)."""
    gaff2_types, is_canonical, ff14sb_types = _cluster_atom_arrays(
        typed_mol, model, ff14sb_amino_lib
    )
    angles, dihedrals = calculateAnglesAndDihedrals(typed_mol.bonds)
    # Reuse the existing splice_into helper - it treats `is_stub` as
    # "atom whose live-system type differs from its antechamber GAFF2
    # type", which is exactly our `is_canonical` semantics here.
    for param_dict, instances in (
        (pset.bond_types, typed_mol.bonds),
        (pset.angle_types, angles),
        (pset.dihedral_types, dihedrals),
    ):
        _splice_into(param_dict, instances, gaff2_types, is_canonical, ff14sb_types)


def prepareClusterResidues(
    typed_path, frcmod_path, model, outdir=None, use_pyodide=None
):
    """Split antechamber output for a cluster model compound into per-
    residue topology files and emit the matching custombonds list.

    For each non-canonical cluster residue the function writes a CIF using
    antechamber's GAFF2 types and the cluster compute's per-atom charges.
    For each canonical anchor the CIF uses the appropriate AMBER residue
    template's ff14SB atom types (mid-chain ``CYX`` / ``NLN`` / ... or
    the matching N- or C-terminal variant ``NCYX`` / ``CCYX`` / ...
    when the residue is at a chain terminus), with charges taken from
    the antechamber compute on the combined model. Each canonical
    residue's bucket resname (assigned by detect, e.g. ``CY1``) keeps
    the residue out of tLeap's built-in libraries so our prepi loads
    instead of the standard template.

    Parameters
    ----------
    typed_path : str
        Antechamber-typed mol2 of the cluster model compound.
    frcmod_path : str
        parmchk2 output for the same model compound.
    model : :class:`ClusterModel`
        Cluster model returned by ``buildClusterModel``.
    outdir : str or None
        Output directory; created if missing. If ``None``, a fresh tempdir
        is used.

    Returns
    -------
    :class:`ClusterOutputs`
    """
    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="cluster_params_")
    else:
        os.makedirs(outdir, exist_ok=True)

    if use_pyodide is None:
        use_pyodide = _in_pyodide()

    typed_mol = Molecule(typed_path)

    # Group atom indices by cluster_residue_index using model.atom_to_residue.
    # Cap atoms (not in atom_to_residue) are dropped from per-residue files.
    per_residue_atom_idxs = {}
    for i, name in enumerate(typed_mol.name):
        cidx = model.atom_to_residue.get(str(name))
        if cidx is None:
            continue
        per_residue_atom_idxs.setdefault(cidx, []).append(i)

    # Lazy-load the ff14SB amino lib only if we have canonical anchors that
    # need ff14SB type rebuilding.
    ff14sb_amino_lib = None
    if model.canonical_renames:
        ff14sb_amino_lib = _load_ff14sb_amino_lib()

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

    for cidx, residue_id in enumerate(model.spec.residues):
        if cidx not in per_residue_atom_idxs:
            continue
        atom_idxs = np.asarray(per_residue_atom_idxs[cidx], dtype=np.int64)
        sub = typed_mol.copy(sel=atom_idxs)
        new_names = [
            model.atom_to_orig_name[str(typed_mol.name[i])] for i in atom_idxs
        ]
        sub.name[:] = new_names

        if cidx in model.canonical_renames:
            # Canonical anchor: pull ff14SB atom types from the matching
            # variant template (CYX / NLN / ... or the N-/C-terminal
            # variant). Charges stay as antechamber computed them on the
            # combined model.
            new_resname = model.canonical_renames[cidx]
            template_atoms = _canonical_variant_template(
                model.spec, cidx, ff14sb_amino_lib
            )
            sub.atomtype[:] = [
                template_atoms.get(n, str(sub.atomtype[i]))
                for i, n in enumerate(new_names)
            ]
            sub.resname[:] = new_resname
        else:
            new_resname = residue_id.resname
            sub.resname[:] = new_resname
            # NCAA backbone: rewrite the standard backbone GAFF2 types
            # to ff14SB so the peptide bonds at the chain neighbours
            # type correctly in tLeap. Track each original GAFF2 type so
            # we can duplicate the frcmod entries below.
            if model.spec.is_chain_resident[cidx]:
                for name, ff_type in backbone_at.items():
                    sel = sub.name == name
                    for orig in sub.atomtype[sel]:
                        backbone_type_renames.setdefault(str(orig), set()).add(ff_type)
                    sub.atomtype[sel] = ff_type

        sub.resid[:] = 1
        sub.segid[:] = "A"
        sub.chain[:] = "A"
        sub.insertion[:] = ""

        # Chain-resident residues need a prepi (with HEAD / TAIL / connect
        # atoms) so tLeap can stitch them into their flanking residues
        # via standard peptide bonds. Free / scaffold residues are
        # written as CIF and converted to mol2 internally by amber.build.
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
                )
                out.topo_paths.append(topo_path)
        else:
            topo_path = os.path.join(outdir, f"{new_resname}.cif")
            if not os.path.isfile(topo_path):
                sub.write(topo_path)
                out.topo_paths.append(topo_path)

    # Copy the parmchk2 frcmod, then for canonical-anchored clusters add
    # cross-FF junction terms (bond/angle/dihedral entries spanning a
    # canonical-residue atom and a non-canonical one, with the canonical-
    # side atom types rewritten from antechamber's GAFF2 to ff14SB so
    # tLeap can resolve them when the canonical residue is loaded with
    # ff14SB types). Pure-NCAA clusters get the antechamber frcmod
    # verbatim (no FF crossing).
    #
    # Write the merged frcmod under the basename of EACH per-residue CIF
    # rather than a single shared name. amber.build's auto-handler for
    # known NCAAs / cofactors / PTMs (``_detect_cofactors_ncaa_ptm``)
    # checks ``param`` basenames against the bundled-resname registry to
    # decide whether to layer in its own prepi - if a cluster residue's
    # resname is in the registry (NLE, MSE, NAG via PTM, ...) but no
    # ``<RESNAME>.frcmod`` is in ``param``, the auto-handler kicks in and
    # collides with our cluster prepi. Naming the frcmod ``<RESNAME>.frcmod``
    # for each cluster residue suppresses the auto-add.
    needs_merged = bool(model.canonical_renames) or bool(backbone_type_renames)
    if needs_merged:
        pset = AmberParameterSet(frcmod_path)
        if model.canonical_renames:
            _add_cluster_junction_terms(pset, typed_mol, model, ff14sb_amino_lib)
        if backbone_type_renames:
            # Duplicate every frcmod entry that mentions an original
            # GAFF2 backbone type under its ff14SB rename so torsions
            # spanning the backbone-sidechain boundary still resolve
            # at build time.
            from htmd.builder.noncanonical import _duplicate_parameters

            _duplicate_parameters(
                pset,
                {orig: list(news) for orig, news in backbone_type_renames.items()},
            )
        canonical_frcmod = os.path.join(outdir, "_cluster_merged.frcmod")
        pset.write(
            canonical_frcmod,
            title="cluster: scaffold + junction parameters",
            style="frcmod",
        )
        frcmod_source = canonical_frcmod
    else:
        frcmod_source = frcmod_path
    seen_basenames = set()
    for cif_path in out.topo_paths:
        basename = os.path.splitext(os.path.basename(cif_path))[0]
        if basename in seen_basenames:
            continue
        seen_basenames.add(basename)
        per_res_frcmod = os.path.join(outdir, f"{basename}.frcmod")
        with open(frcmod_source) as fr, open(per_res_frcmod, "w") as fw:
            fw.write(fr.read())
        out.frcmod_paths.append(per_res_frcmod)

    # custombonds: one pair per cluster bond.
    for bond in model.spec.bonds:
        out.custombonds.append(
            (_atom_sel_for_unique(bond.atom_a), _atom_sel_for_unique(bond.atom_b))
        )

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


def _build_internal_cluster_spec(member_indices, cluster_bonds, mol, groups, spec_by_res_idx):
    """Construct a ClusterSpec for the cluster pipeline from the per-residue
    specs that fall inside one connected component of non-peptide bonds.
    A singleton component containing just one chain-resident NCAA is also
    valid - it gets its own 1-residue cluster with peptide-neighbour caps
    drawn from the live mol."""
    from moleculekit.tools.nonstandard_residues import (
        NCAASpec, CrosslinkedNCAASpec, ScaffoldSpec, CovalentLigandSpec,
        CanonicalRenamedSpec,
    )

    residues = []
    is_canonical = []
    is_chain_resident = []
    roles = []
    canonical_original_resnames = []
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
        if isinstance(spec, CanonicalRenamedSpec):
            is_canonical.append(True)
            is_chain_resident.append(True)
            roles.append("anchor")
            canonical_original_resnames.append(spec.residue.resname)
            # Pick the terminal-variant template only if the residue in
            # the prepared mol actually carries the terminal atoms
            # (OXT for C-term, H1/H2/H3 for N-term). Without those
            # atoms ``amber.build`` will add an ACE/NME cap, which makes
            # the residue effectively mid-chain at build time and the
            # mid-chain ff14SB template (CYX, NLN, ...) is the right
            # match instead of the terminal variant (CCYX, NCYX, ...).
            res_atoms = set(
                str(n) for n in mol.name[groups[r_idx]["atom_idx"]]
            )
            n_term_eff = spec.is_n_term and bool({"H1", "H2", "H3"} & res_atoms)
            c_term_eff = spec.is_c_term and "OXT" in res_atoms
            canonical_terminus.append(
                "n" if n_term_eff else ("c" if c_term_eff else "")
            )
            is_n_term.append(n_term_eff)
            is_c_term.append(c_term_eff)
            has_canonical = True
        elif isinstance(spec, CrosslinkedNCAASpec):
            is_canonical.append(False)
            is_chain_resident.append(True)
            roles.append("stapled_ncaa")
            canonical_original_resnames.append("")
            canonical_terminus.append("")
            is_n_term.append(spec.is_n_term)
            is_c_term.append(spec.is_c_term)
            n_chain_resident_nc += 1
            n_nc += 1
        elif isinstance(spec, NCAASpec):
            is_canonical.append(False)
            is_chain_resident.append(True)
            roles.append("ncaa")
            canonical_original_resnames.append("")
            canonical_terminus.append("")
            is_n_term.append(spec.is_n_term)
            is_c_term.append(spec.is_c_term)
            n_chain_resident_nc += 1
            n_nc += 1
        elif isinstance(spec, ScaffoldSpec):
            is_canonical.append(False)
            is_chain_resident.append(False)
            roles.append("scaffold")
            canonical_original_resnames.append("")
            canonical_terminus.append("")
            is_n_term.append(False)
            is_c_term.append(False)
            has_scaffold = True
            n_nc += 1
        elif isinstance(spec, CovalentLigandSpec):
            is_canonical.append(False)
            is_chain_resident.append(False)
            roles.append("covalent_ligand")
            canonical_original_resnames.append("")
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
        canonical_original_resnames=canonical_original_resnames,
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
    "H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1, "SI": 4, "P": 3,
    "S": 2, "CL": 1, "AS": 3, "SE": 2, "BR": 1, "I": 1,
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
    specs, mol, outdir, charge_method="gasteiger", use_pyodide=None
):
    """Parameterize every non-canonical residue in ``specs`` and return
    paths plus custombonds ready to feed :func:`htmd.builder.amber.build`.

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
    mol : :class:`moleculekit.molecule.Molecule`
        The molecule the specs describe. Must already carry covalent
        bonds (typically the post-``systemPrepare`` molecule).
    outdir : str
        Output directory for all generated CIF / frcmod files.
    charge_method : str, optional
        Charge method passed to antechamber's ``-c`` flag (``"gasteiger"``
        by default; ``"am1-bcc"`` available for higher accuracy at
        significant runtime cost).
    use_pyodide : bool or None, optional
        Force the AmberTools dispatch path (``True`` -> dispatch via
        ``antechamber_pyodide.run``; ``False`` -> native subprocess).
        ``None`` (default) auto-detects Pyodide via ``sys.platform``.

    Returns
    -------
    :class:`ClusterOutputs`
        Aggregated topology files, frcmod files, and custombonds for the
        whole system. Feed ``out.topo_paths`` to ``amber.build(topo=...)``,
        ``out.frcmod_paths`` to ``param=``, and ``out.custombonds`` to
        ``custombonds=``.

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
        pmol = systemPrepare(mol, detect_specs=specs)

        # 4. Run antechamber per cluster and split per-residue.
        out = parameterizeFromSpecs(specs, pmol, outdir="./params")

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
        NCAASpec,
        CrosslinkedNCAASpec,
        ScaffoldSpec,
        CovalentLigandSpec,
        LigandSpec,
        CanonicalRenamedSpec,
    )
    from htmd.builder.noncanonical import (
        _fftype_antechamber,
        parameterizeNonCanonicalResidues,
    )

    if use_pyodide is None:
        use_pyodide = _in_pyodide()

    os.makedirs(outdir, exist_ok=True)
    a2r, groups = _residue_groups_with_index(mol)

    # Map each spec to its residue index in mol.
    res_key_to_idx = {
        (g["segid"], g["chain"], int(g["resid"]), g["insertion"]): ri
        for ri, g in enumerate(groups)
    }
    spec_by_res_idx = {}
    for s in specs:
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
    spec_indices = set(spec_by_res_idx)
    inter_bonds = []  # (a1, a2, r1, r2) for non-peptide bonds inside the spec'd set
    for a1, a2 in mol.bonds:
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
    cluster_membership = []  # list of (sorted member_indices, bonds)
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
        if isinstance(spec, NCAASpec):
            cluster_membership.append(([r_idx], []))
            in_cluster_set.add(r_idx)

    out = ClusterOutputs()

    # Singleton NCAA clusters (one chain-resident NCAA, no cluster bonds)
    # that share resname + terminal flags are chemically identical model
    # compounds, so antechamber output is identical too. tLeap loads only
    # one unit per resname anyway. Dedup by (resname, is_n_term, is_c_term)
    # to skip the redundant antechamber runs and avoid emitting N copies
    # of the same prepi.
    seen_singleton_keys = set()

    # Run the cluster pipeline once per connected component.
    for ci, (members, cbonds) in enumerate(cluster_membership):
        if len(members) == 1 and not cbonds:
            spec = spec_by_res_idx[members[0]]
            if isinstance(spec, NCAASpec):
                g = groups[members[0]]
                key = (
                    g["resname"],
                    bool(spec.is_n_term),
                    bool(spec.is_c_term),
                )
                if key in seen_singleton_keys:
                    continue
                seen_singleton_keys.add(key)

        cluster_spec = _build_internal_cluster_spec(
            members, cbonds, mol, groups, spec_by_res_idx
        )
        cluster_outdir = os.path.join(outdir, f"cluster_{ci:03d}")
        os.makedirs(cluster_outdir, exist_ok=True)

        cluster_model = buildClusterModel(mol, cluster_spec, cluster_outdir)
        cluster_mol = Molecule(cluster_model.cif_path)
        netcharge = int(round(float(np.sum(cluster_mol.formalcharge))))
        typed_mol, typed_path, frcmod_path = _fftype_antechamber(
            cluster_mol,
            tmpdir=cluster_outdir,
            method="gaff2",
            netcharge=netcharge,
            charge_method=charge_method,
            use_pyodide=use_pyodide,
        )
        del typed_mol  # we re-load it inside prepareClusterResidues
        cluster_out = prepareClusterResidues(
            typed_path,
            frcmod_path,
            cluster_model,
            outdir=cluster_outdir,
            use_pyodide=use_pyodide,
        )
        out.topo_paths.extend(cluster_out.topo_paths)
        out.frcmod_paths.extend(cluster_out.frcmod_paths)
        out.custombonds.extend(cluster_out.custombonds)

    # Non-chain residues that ended up outside any cluster: run
    # antechamber standalone and write a per-resname CIF + frcmod pair.
    # ``LigandSpec`` is the natural case (no covalent bonds at all);
    # ``ScaffoldSpec`` / ``CovalentLigandSpec`` can also land here when
    # their only inter-residue bonds go to canonical residues that
    # aren't spec'd themselves (e.g. a metal-coordinating ligand whose
    # only neighbours are zinc ions). Chain-resident specs without a
    # cluster shouldn't happen because those are promoted to singleton
    # clusters earlier; if one slips through it's an input bug.
    standalone_types = (LigandSpec, ScaffoldSpec, CovalentLigandSpec)
    seen_ligand_resnames = set()
    for r_idx, spec in spec_by_res_idx.items():
        if r_idx in in_cluster_set:
            continue
        g = groups[r_idx]
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
        ligand_cif = os.path.join(ligand_dir, f"{g['resname']}.cif")
        sub = _write_free_residue_cif(mol, g, ligand_cif)
        netcharge = int(round(float(np.sum(sub.formalcharge))))
        typed_mol, _, frcmod_path = _fftype_antechamber(
            sub,
            tmpdir=ligand_dir,
            method="gaff2",
            netcharge=netcharge,
            charge_method=charge_method,
            use_pyodide=use_pyodide,
        )
        out_cif = os.path.join(outdir, f"{g['resname']}.cif")
        out_frcmod = os.path.join(outdir, f"{g['resname']}.frcmod")
        typed_mol.resname[:] = g["resname"]
        typed_mol.write(out_cif)
        with open(frcmod_path) as fr, open(out_frcmod, "w") as fw:
            fw.write(fr.read())
        out.topo_paths.append(out_cif)
        out.frcmod_paths.append(out_frcmod)

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

    # Re-guess intra-residue bonds for cap residues (always sliced down to
    # a few atoms) and for cluster residues that have any atom missing
    # bonds (e.g. an OXT that systemPrepare added without re-attaching the
    # C-OXT bond). Scoped one residue at a time so we don't over-bond near
    # residue boundaries.
    existing = {tuple(sorted((int(a), int(b)))) for a, b in model.bonds}
    new_bonds, new_bondtypes = [], []
    residues_to_guess = dict(cap_atoms)
    cluster_atoms_by_residue = {}
    for orig in orig_to_cluster_idx:
        live_res = int(a2r[orig])
        cluster_atoms_by_residue.setdefault(live_res, {})[orig] = None
    has_neighbor = np.zeros(model.numAtoms, dtype=bool)
    for a, b in model.bonds:
        has_neighbor[int(a)] = True
        has_neighbor[int(b)] = True
    for live_res, atoms_dict in cluster_atoms_by_residue.items():
        if any(not has_neighbor[orig_to_new[o]] for o in atoms_dict):
            residues_to_guess.setdefault(live_res, {}).update(atoms_dict)
    for residx_in_mol, atoms_dict in residues_to_guess.items():
        cap_orig_idxs_local = list(atoms_dict.keys())
        sub = mol.copy(sel=np.asarray(cap_orig_idxs_local, dtype=np.int64))
        sorted_cap_orig = sorted(cap_orig_idxs_local)
        local_to_new = {
            local: orig_to_new[orig] for local, orig in enumerate(sorted_cap_orig)
        }
        guessed = sub._guessBonds(rdkit=False)
        for la, lb in guessed:
            ga, gb = local_to_new[int(la)], local_to_new[int(lb)]
            key = tuple(sorted((ga, gb)))
            if key not in existing:
                new_bonds.append([ga, gb])
                new_bondtypes.append("1")
                existing.add(key)
    if new_bonds:
        model.bonds = np.vstack(
            [model.bonds, np.asarray(new_bonds, dtype=np.uint32)]
        )
        model.bondtype = np.hstack(
            [model.bondtype, np.asarray(new_bondtypes)]
        )

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
            allowed = {("C", "O"), ("O", "C"), ("C", "CH3"), ("CH3", "C"),
                       ("N", "H"), ("H", "N"), ("N", "CH3"), ("CH3", "N")}
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
                int(nb) for nb in model.getNeighbors(i)
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


def buildClusterModel(mol, spec, outdir):
    """Build the combined model compound for a :class:`ClusterSpec`: full
    residues + ACE/NME-style backbone caps derived from the live mol's
    chain neighbours, written as a CIF ready to feed to antechamber."""
    os.makedirs(outdir, exist_ok=True)
    nc_resnames = "_".join(
        sorted(
            r.resname
            for r, c in zip(spec.residues, spec.is_canonical)
            if not c
        )
    )
    cif_path = os.path.join(outdir, f"cluster_{spec.subtype}_{nc_resnames}.cif")
    model, atom_map, atom_to_residue, atom_to_orig_name = (
        _build_cluster_model_accurate(mol, spec, cif_path)
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


