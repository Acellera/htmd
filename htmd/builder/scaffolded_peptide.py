"""Helpers for scaffolded-peptide system building under AMBER.

The detector :func:`moleculekit.tools.nonstandard_residues.detectNonStandardResidues`
writes a model compound CIF for each scaffolded-peptide it finds: the scaffold
geometry plus a small canonical-sidechain stub (anchor atom + methyl cap) at
each anchor. The user runs :func:`htmd.builder.noncanonical._fftype_antechamber`
on that CIF to get a typed mol2 + frcmod.

:func:`prepareScaffoldedResidue` post-processes those outputs into a
matched ``<RESNAME>.cif`` + ``<RESNAME>.frcmod`` pair ready to feed
:func:`htmd.builder.amber.build`. The two files share the residue's name as
their basename so AMBER's per-residue atom-type collision-detection
(``_fix_parameterize_atomtype_collisions``) can pair them.

The CIF holds the scaffold residue alone (model-compound stubs removed). The
frcmod combines two pieces in one file:

  - antechamber's scaffold-internal terms verbatim (covers everything inside
    the residue);
  - junction terms: bonds/angles/dihedrals that span the FF boundary in the
    live system (>=1 stub atom and >=1 scaffold atom in the model compound),
    with stub-atom GAFF2 types rewritten to the canonical-FF type per
    :data:`moleculekit.tools._anchor_variants.ANCHOR_VARIANTS`.

All parameter manipulation (loading/inspecting/writing the frcmod) is done
through :class:`parmed.amber.AmberParameterSet`. There is no custom frcmod
parser in this module.
"""

from __future__ import annotations

import copy
import logging
import os
import tempfile

import numpy as np
from parmed.amber import AmberParameterSet

from moleculekit.molecule import Molecule
from moleculekit.tools._anchor_variants import ANCHOR_VARIANTS  # noqa: F401  (re-exported)
from moleculekit.util import calculateAnglesAndDihedrals

logger = logging.getLogger(__name__)


def prepareScaffoldedResidue(typed_path, frcmod_path, spec, outdir=None):
    """Post-process antechamber output for a scaffolded-peptide model compound
    into a matched ``<RESNAME>.cif`` + ``<RESNAME>.frcmod`` pair.

    Reads the typed mol2/cif and the parmchk2 frcmod for the *model compound*
    (scaffold + canonical-sidechain stubs) and writes two files into
    ``outdir``, both basenamed after the scaffold residue so AMBER's
    per-residue atom-type collision dedup pairs them:

      - ``<RESNAME>.cif``: the scaffold residue alone (stubs removed),
        ready to pass to :func:`htmd.builder.amber.build` via ``topo=``;
      - ``<RESNAME>.frcmod``: antechamber's frcmod with the cross-FF
        junction terms appended. Stub-atom GAFF2 types are rewritten to the
        canonical-FF types per :data:`ANCHOR_VARIANTS`. Pass to ``param=``.

    All frcmod manipulation goes through :class:`parmed.amber.AmberParameterSet`.

    Parameters
    ----------
    typed_path : str
        Path to a typed file produced by antechamber on the model compound
        (typically the ``typed.mol2`` file written next to the frcmod). Must
        load via :class:`moleculekit.molecule.Molecule` with ``mol.atomtype``
        populated and atom names matching ``spec.model_atom_map``.
    frcmod_path : str
        Path to the frcmod produced by parmchk2 on the model compound.
    spec : dict
        A scaffolded-peptide spec from
        :func:`moleculekit.tools.nonstandard_residues.detectNonStandardResidues`
        (must carry ``model_atom_map``).
    outdir : str or None, optional
        Directory in which to write the two output files. If ``None``, a fresh
        temporary directory is created via :func:`tempfile.mkdtemp` and used;
        its path is encoded in the returned paths.

    Returns
    -------
    cif_path : str
        Path to ``<outdir>/<RESNAME>.cif``.
    frcmod_path : str
        Path to ``<outdir>/<RESNAME>.frcmod``.
    """
    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="scaffold_params_")
    else:
        os.makedirs(outdir, exist_ok=True)

    typed_mol = Molecule(typed_path)
    resname = spec.resname
    out_cif = os.path.join(outdir, f"{resname}.cif")
    out_frcmod = os.path.join(outdir, f"{resname}.frcmod")

    _write_scaffold_only_cif(typed_mol, spec, out_cif)

    pset = AmberParameterSet(frcmod_path)
    _add_junction_terms(pset, typed_mol, spec)
    _drop_unused_terms(pset, typed_mol, spec)
    pset.write(
        out_frcmod,
        title=f"{resname}: scaffold + junction parameters",
        style="frcmod",
    )

    return out_cif, out_frcmod


def _write_scaffold_only_cif(typed_mol, spec, out_path):
    """Drop stub atoms from ``typed_mol`` and write the result as a CIF named
    after the scaffold's resname (so tLeap registers it as the right unit)."""
    atom_map = spec.model_atom_map
    if atom_map is None:
        raise ValueError(
            "spec.model_atom_map is None; was the spec produced with write_models=True?"
        )
    scaffold_only = typed_mol.copy()
    stub_mask = np.array(
        [atom_map[str(n)].role == "stub" for n in scaffold_only.name],
        dtype=bool,
    )
    if stub_mask.any():
        scaffold_only.remove(stub_mask, _logger=False)
    scaffold_only.resname[:] = spec.resname
    scaffold_only.resid[:] = 1
    scaffold_only.write(out_path)


def _add_junction_terms(pset, typed_mol, spec):
    """For every bond / angle / dihedral *instance* in ``typed_mol`` whose
    atoms span both the scaffold and a stub, add a new entry to the matching
    ``pset.<param>_types`` dict keyed on the stub-atoms-rewritten signature.

    Parameter values are looked up in ``pset`` itself by the instance's GAFF2
    type signature (antechamber/parmchk2 emitted them) and copied; parmed
    deduplicates types by ``id()`` when writing, so the junction entry must
    be a distinct object.  Impropers are left as-is (antechamber's
    empirical impropers apply by atom type and don't need rewriting at the
    junction).

    Angles and dihedrals come from
    :func:`moleculekit.util.calculateAnglesAndDihedrals`, the same helper
    used by :mod:`htmd.builder.noncanonical`."""
    atom_map = spec.model_atom_map
    if atom_map is None:
        raise ValueError(
            "spec.model_atom_map is None; was the spec produced with write_models=True?"
        )

    atom_types = np.asarray(typed_mol.atomtype)
    is_stub = np.asarray(
        [atom_map[str(n)].role == "stub" for n in typed_mol.name],
        dtype=bool,
    )
    ff_types = np.asarray(
        [atom_map[str(n)].ff_type or "" for n in typed_mol.name],
        dtype=object,
    )

    angles, dihedrals = calculateAnglesAndDihedrals(typed_mol.bonds)

    _splice_into(
        pset.bond_types, typed_mol.bonds, atom_types, is_stub, ff_types
    )
    _splice_into(
        pset.angle_types, angles, atom_types, is_stub, ff_types
    )
    _splice_into(
        pset.dihedral_types, dihedrals, atom_types, is_stub, ff_types
    )


def _splice_into(param_dict, instances, atom_types, is_stub, ff_types):
    """For each instance (an array of atom indices) that spans stub+scaffold,
    look up the parameter value at the GAFF2 key, derive the rewritten
    canonical-FF key, and write the value (deep-copied) under the new key.

    ``param_dict`` is mutated in place. ``instances`` is an ``(N, k)``
    integer array of atom indices (k=2 for bonds, 3 for angles, 4 for
    dihedrals)."""
    for idxs in instances:
        idxs = tuple(int(i) for i in idxs)
        roles = is_stub[list(idxs)]
        if not (roles.any() and (~roles).any()):
            continue  # not a mixed junction instance

        gaff_key = tuple(atom_types[list(idxs)])
        if gaff_key not in param_dict:
            continue
        new_key = _live_key(idxs, roles, atom_types, ff_types)
        if new_key is None or new_key in param_dict:
            continue
        # parmed.write deduplicates types by id(); copy so the junction entry
        # is written as a distinct line and not skipped as a "duplicate".
        val_copy = copy.copy(param_dict[gaff_key])
        param_dict[new_key] = val_copy
        param_dict[new_key[::-1]] = val_copy


def _live_key(idxs, roles, atom_types, ff_types):
    """Build the atom-type tuple a model-compound instance maps to in the
    *live* AMBER system: scaffold atoms keep their GAFF2 types; stub atoms
    are replaced by their canonical-FF type from ANCHOR_VARIANTS. Returns
    ``None`` if any stub atom is missing a canonical-FF mapping."""
    out = []
    for i, stub in zip(idxs, roles):
        if stub:
            ff = ff_types[i]
            if not ff:
                return None
            out.append(str(ff))
        else:
            out.append(str(atom_types[i]))
    return tuple(out)


def _drop_unused_terms(pset, typed_mol, spec):
    """Remove parameters and atom types that won't be referenced when the
    scaffold residue is loaded into the live AMBER system.

    A model-compound bond/angle/dihedral instance maps to a live-system key
    only if it has at least one scaffold atom; stub-only instances are
    dropped (those atoms don't exist in the live system; their canonical-FF
    counterparts get their parameters from ff14SB instead). For atom types
    in MASS/NONB sections, only GAFF2 types worn by scaffold atoms are kept.
    Canonical-FF types like ``S`` or ``2C`` are referenced by junction keys
    but aren't redefined here, since ff14SB already provides them.

    Mirrors the cleanup pattern in :func:`htmd.builder.noncanonical._clean_prm`.
    """
    atom_map = spec.model_atom_map
    is_stub = np.asarray(
        [atom_map[str(n)].role == "stub" for n in typed_mol.name],
        dtype=bool,
    )
    atom_types = np.asarray(typed_mol.atomtype)
    ff_types = np.asarray(
        [atom_map[str(n)].ff_type or "" for n in typed_mol.name],
        dtype=object,
    )

    angles, dihedrals = calculateAnglesAndDihedrals(typed_mol.bonds)

    def _live_keys(instances):
        keys = set()
        for idxs in instances:
            idxs = tuple(int(i) for i in idxs)
            roles = is_stub[list(idxs)]
            if roles.all():
                continue  # stub-only instance; absent from live system
            k = _live_key(idxs, roles, atom_types, ff_types)
            if k is None:
                continue
            keys.add(k)
            keys.add(k[::-1])
        return keys

    keep_bonds = _live_keys(typed_mol.bonds)
    keep_angles = _live_keys(angles)
    keep_dihedrals = _live_keys(dihedrals)

    # Atom types kept in MASS/NONB: GAFF2 types worn by scaffold atoms only.
    allowed_atom_types = set(str(t) for t in atom_types[~is_stub].tolist())

    for at in list(pset.atom_types):
        if at not in allowed_atom_types:
            del pset.atom_types[at]

    for d, keep in (
        (pset.bond_types, keep_bonds),
        (pset.angle_types, keep_angles),
        (pset.dihedral_types, keep_dihedrals),
    ):
        for key in list(d):
            if key not in keep:
                del d[key]

    # Impropers: drop entries that reference any deleted atom type.
    # Wildcard 'X' in impropers is always permitted.
    for d_name in ("improper_periodic_types", "improper_types"):
        d = getattr(pset, d_name, None)
        if d is None:
            continue
        for key in list(d):
            if not all(t == "X" or t in allowed_atom_types for t in key):
                del d[key]
