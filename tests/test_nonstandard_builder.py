"""End-to-end tests for non-standard residue building under AMBER.

Each end-to-end test follows the recommended five-step flow:

  1. ``Molecule(...)`` - load the input.
  2. ``detectNonStandardResidues(mol)`` - inspect ``mol`` and return one
     spec per non-standard residue. Pure analysis, no mol mutation.
  3. ``mol.templateResidueFromSmiles(...)`` - set bond orders and add
     hydrogens on each non-canonical residue from a SMILES template.
  4. ``systemPrepare(mol, detect_specs=specs)`` - protonate the canonical
     part with PDB2PQR, then apply the per-spec rename + displaced-H
     drops. Inter-residue bonds are preserved through the prep.
  5. ``parameterizeFromSpecs(specs, pmol, outdir)`` - run antechamber on
     each cluster + free residue; return CIF / frcmod paths and
     custombonds ready to feed ``amber.build``.
  6. ``amber.build(...)`` with the cluster outputs.

Two short sanity tests at the top exercise the disulfide path and
``amber.build``'s defense-in-depth rename for users that hand-build
custombonds without going through detect.
"""

import os
import shutil
from collections import Counter

import numpy as np
import pytest

from moleculekit.molecule import Molecule, UniqueAtomID, mol_equal
from moleculekit.tools.nonstandard_residues import (
    detectNonStandardResidues,
    ChainResidueSpec,
    ScaffoldSpec,
    CovalentLigandSpec,
    LigandSpec,
    PROTEIN_RESNAMES,
)
from htmd.builder.amber import _findTeLeap, _prepareMolecule
from htmd.builder.nonstandard import (
    parameterizeFromSpecs,
    _check_specs_protonated,
    _clean_frcmod_params,
    _residue_groups_with_index,
    _normalize_residue_charges,
    _backbone_charge_map,
    _load_ff14sb_amino_lib_amber,
    _load_ff14sb_amino_lib_openmm,
    _engine_for_forcefield,
    _resolve_forcefield,
    _warn_if_openff_mismatched_charges,
    _FF14SB_BACKBONE_CHARGES_BY_CLASS,
)

try:
    import openff.nagl  # noqa: F401
    import openff.nagl_models  # noqa: F401

    _nagl = True
except (ImportError, SyntaxError):
    _nagl = False


# Lazy-loaded once for the _backbone_charge_map unit tests below.
# Loader needs tleap on PATH, so wrap in the same _tleap skip guard the
# individual tests use.
def _amber_lib_or_none():
    try:
        return _load_ff14sb_amino_lib_amber()
    except (RuntimeError, FileNotFoundError):
        return None

curr_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(curr_dir, "test_nonstandard_builder")
QFZ_B_CIF = os.path.join(DATA_DIR, "8QFZ_B.cif")
QU4_A_CIF = os.path.join(DATA_DIR, "8QU4_A.cif")
VBL_PDB = os.path.join(DATA_DIR, "5vbl.pdb")
TOT_E_PDB = os.path.join(DATA_DIR, "4TOT_E.pdb")
R1J_PDB = os.path.join(DATA_DIR, "1r1j.pdb")

# Pre-reaction LFI warhead: 1,3,5-triazinane core with three bromoacetyl arms.
# After reaction with the three Cys SG atoms, the Br leaving groups are
# displaced and replaced by S-Cys connections.
LFI_SMILES = "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr"

NLE_SMILES = "CCCC[C@@H](C(=O)O)N"
MK8_SMILES = "CCCC[C@](C)(C(=O)O)N"
# Protonated benzamidine, the 3PTB ligand (from moleculekit's rdkittools docs).
BEN_SMILES = "[NH2+]=C(N)c1ccccc1"

_antechamber = shutil.which("antechamber") is not None
_tleap = _findTeLeap() is not None


# Maximum number of bonded neighbours a heavy atom of a given element
# should ever have in a built ff14SB / GAFF2 system. Bigger counts are
# almost always over-protonation bugs (e.g. a CH3 atom that ended up as
# CH4 with 5 bonds because the displaced-H drop didn't fire). Formal
# charges are not preserved through the prmtop, so this is the cheap
# sanity check we have on the built mol.
_MAX_BONDS_BY_ELEMENT = {
    "H": 1, "C": 4, "N": 4, "O": 2, "S": 4, "P": 5,
    "F": 1, "CL": 1, "BR": 1, "I": 1,
}


def _check_no_overvalent_atoms(built):
    counts = np.zeros(built.numAtoms, dtype=np.int32)
    for a, b in built.bonds:
        counts[int(a)] += 1
        counts[int(b)] += 1
    # AMBER's TIP3P water topology carries an H-H SHAKE bond, so every
    # water H legitimately has 2 bonds. Skip waters in the valence check.
    is_water = np.isin(built.resname, ("WAT", "HOH", "TIP3"))
    bad = []
    for i in range(built.numAtoms):
        if is_water[i]:
            continue
        elem = str(built.element[i]).upper()
        cap = _MAX_BONDS_BY_ELEMENT.get(elem)
        if cap is not None and int(counts[i]) > cap:
            bad.append(
                f"  atom {i} resname={str(built.resname[i])!r} "
                f"resid={int(built.resid[i])} name={str(built.name[i])!r} "
                f"element={elem!r}: {int(counts[i])} bonds (max {cap})"
            )
    assert not bad, "over-valent atoms in built mol:\n" + "\n".join(bad)


def _backbone_ring_closures(built):
    """Return the list of head-to-tail cyclic backbone closures in ``built``
    as ``(resid_a, resid_b)`` tuples. A closure is a backbone ``C``-``N`` bond
    between two residues that are not sequential neighbours (a linear peptide
    only ever bonds consecutive residues; a macrocycle additionally bonds its
    first residue back to its last)."""
    name = built.name
    resid = built.resid
    closures = []
    for a, b in built.bonds:
        a, b = int(a), int(b)
        if {str(name[a]), str(name[b])} != {"C", "N"}:
            continue
        if abs(int(resid[a]) - int(resid[b])) <= 1:
            continue
        closures.append((int(resid[a]), int(resid[b])))
    return closures


def _beta_residue_mol():
    """A minimal beta-amino-acid residue (backbone N-CA-C18-C, shortest N->C
    path = 4), like microcystin's Adda."""
    mol = Molecule().empty(5)
    mol.name[:] = ["N", "CA", "C18", "C", "O"]
    mol.element[:] = ["N", "C", "C", "C", "O"]
    mol.resname[:] = "BZA"
    mol.resid[:] = 1
    mol.chain[:] = "A"
    mol.segid[:] = "P0"
    mol.record[:] = "ATOM"
    mol.coords = np.array(
        [[0, 0, 0], [1.5, 0, 0], [2.5, 1, 0], [3.5, 0.5, 0], [3.5, -0.7, 0]],
        dtype=np.float32,
    ).reshape(5, 3, 1)
    mol.bonds = np.array([[0, 1], [1, 2], [2, 3], [3, 4]], dtype=np.uint32)
    mol.bondtype = np.array(["1"] * 4, dtype=object)
    return mol


def test_warn_if_nonalpha_backbone(caplog):
    """A chain-resident residue with a non-alpha (beta/gamma) backbone triggers
    a loud warning that its non-alpha backbone atoms get GAFF/OpenFF charges +
    analogy torsions, which are known-approximate for non-alpha amino acids."""
    import logging
    from moleculekit.molecule import UniqueResidueID
    from htmd.builder.nonstandard import _warn_if_nonalpha_backbone

    mol = _beta_residue_mol()
    spec = ChainResidueSpec(
        resname="BZA",
        residue=UniqueResidueID.fromMolecule(mol, "resid 1"),
        is_n_term=True,
        is_c_term=True,
    )
    with caplog.at_level(logging.WARNING):
        _warn_if_nonalpha_backbone(mol, [spec])
    msgs = " ".join(r.message for r in caplog.records)
    assert "BZA" in msgs, msgs
    assert "non-alpha" in msgs.lower() or "backbone" in msgs.lower(), msgs


def test_warn_if_nonalpha_backbone_silent_for_alpha(caplog):
    """A standard alpha residue (N-CA-C, path 3) triggers no such warning."""
    import logging
    from moleculekit.molecule import UniqueResidueID
    from htmd.builder.nonstandard import _warn_if_nonalpha_backbone

    mol = _beta_residue_mol()
    # collapse to an alpha backbone: bond CA directly to C (drop C18 from path)
    mol.remove("name C18", _logger=False)
    mol.bonds = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.uint32)  # N-CA-C-O
    mol.bondtype = np.array(["1"] * 3, dtype=object)
    spec = ChainResidueSpec(
        resname="BZA",
        residue=UniqueResidueID.fromMolecule(mol, "resid 1"),
        is_n_term=True,
        is_c_term=True,
    )
    with caplog.at_level(logging.WARNING):
        _warn_if_nonalpha_backbone(mol, [spec])
    assert not caplog.records, [r.message for r in caplog.records]


# N-terminal NH3+ first-hydrogen naming convention differs between
# builders: amber.build emits ``H1`` / ``H2`` / ``H3``; openff.build
# emits ``H`` / ``H2`` / ``H3``. Topologically identical. The alias is
# resname-scoped: caps (ACE / NME) carry methyl hydrogens named
# ``H1`` / ``H2`` / ``H3`` that must NOT be rewritten, otherwise the
# NME methyl ``H1`` collides with the amide ``H`` already on the
# residue.
_CAP_RESNAMES = {"ACE", "NME", "NHE", "NMA"}


def _normalised_atom_name(name, resname):
    if str(resname) in _CAP_RESNAMES:
        return str(name)
    return "H" if str(name) == "H1" else str(name)


def _assert_builds_equivalent(
    amber_built, openmm_built, amber_prmtop, openmm_prmtop, energy_tol=0.05,
):
    """Assert that the two prmtops produce the same potential energy on
    the same physical coordinates (within ``energy_tol`` kcal/mol).

    The builds may differ in atom ORDER, so we pair atoms by
    ``(sequenceID, normalised name)`` and feed each system its own
    permutation of the shared coord set. ``constraints=None`` ensures
    every X-H bond contributes to the potential (with ``HBonds``,
    OpenMM rewrites them to rigid constraints and drops them from the
    HarmonicBondForce, which would mask real divergence).
    """
    import openmm
    import openmm.app as omm_app
    import openmm.unit as omm_unit
    from moleculekit.util import sequenceID

    assert openmm_built.numAtoms == amber_built.numAtoms, (
        f"atom count differs: amber={amber_built.numAtoms}, "
        f"openmm={openmm_built.numAtoms}"
    )

    def _atom_ids(mol):
        seq = sequenceID((mol.resid, mol.insertion, mol.segid))
        return [
            (int(seq[i]), _normalised_atom_name(mol.name[i], mol.resname[i]))
            for i in range(mol.numAtoms)
        ]

    am_ids = _atom_ids(amber_built)
    om_ids = _atom_ids(openmm_built)
    assert set(am_ids) == set(om_ids), (
        f"atom identity sets differ: "
        f"in amber not openmm = {len(set(am_ids) - set(om_ids))}; "
        f"in openmm not amber = {len(set(om_ids) - set(am_ids))}"
    )

    am_idx = {k: i for i, k in enumerate(am_ids)}
    om_idx = {k: i for i, k in enumerate(om_ids)}
    shared = amber_built.coords[:, :, 0].astype(np.float64)
    om_coords = np.empty_like(shared)
    for k, j in om_idx.items():
        om_coords[j] = shared[am_idx[k]]

    def _energy(prmtop, coords):
        system = omm_app.AmberPrmtopFile(prmtop).createSystem(
            nonbondedMethod=omm_app.NoCutoff, constraints=None, rigidWater=False
        )
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
                force.setUseDispersionCorrection(False)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001))
        ctx.setPositions(coords * 0.1 * omm_unit.nanometers)
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            omm_unit.kilocalorie_per_mole
        )

    e_amber = _energy(amber_prmtop, shared)
    e_openff = _energy(openmm_prmtop, om_coords)

    assert abs(e_amber - e_openff) < energy_tol, (
        f"amber and openff builds disagree by more than {energy_tol} kcal/mol\n"
        f"  amber : {e_amber:.4f} kcal/mol\n"
        f"  openff: {e_openff:.4f} kcal/mol"
    )


def _run_pipeline(
    mol, smiles, tmp_path, build_kwargs=None
):
    """Drive the full autoSegment -> detect -> template -> systemPrepare ->
    parameterize -> amber.build pipeline and return the built Molecule.
    ``autoSegment2`` runs first so each protein chain and each
    non-protein residue (NCAAs, scaffolds, ions, ligands) ends up in
    its own segment - that's what stops amber.build from extending the
    protein chain through a HETATM ligand or auto-capping in the
    middle of a covalent group."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build

    mol = autoSegment(mol, fields=("segid", "chain"), _logger=False)
    # Strip input Hs so each residue gets a clean re-protonation from the
    # SMILES template / PDB2PQR pass, matching the canonical flow used by
    # moleculekit's systemPrepare tests.
    mol.remove("element H", _logger=False)

    specs = detectNonStandardResidues(mol)
    for resname, smi in smiles.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol, _ = systemPrepare(
        mol,
        detect_specs=specs,
    )
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )
    return amber_build(
        pmol,
        outdir=str(tmp_path / "build"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        **(build_kwargs or {}),
    )


def test_disulfide_path_unchanged(tmp_path):
    """An isolated Cys-Cys disulfide pair still gets renamed CYX with HG
    dropped through the existing disulfide rename machinery in
    ``_prepareMolecule``."""
    rows = [
        ("N", "N", 1, 0.0, 0.0, 0.0),
        ("CA", "C", 1, 1.5, 0.0, 0.0),
        ("C", "C", 1, 2.5, 1.0, 0.0),
        ("O", "O", 1, 2.0, 2.0, 0.0),
        ("CB", "C", 1, 1.5, -1.5, 0.0),
        ("SG", "S", 1, 1.5, -3.0, 0.0),
        ("HG", "H", 1, 1.5, -3.5, 0.8),
        ("N", "N", 2, 5.0, 0.0, 0.0),
        ("CA", "C", 2, 6.5, 0.0, 0.0),
        ("C", "C", 2, 7.5, 1.0, 0.0),
        ("O", "O", 2, 7.0, 2.0, 0.0),
        ("CB", "C", 2, 6.5, -1.5, 0.0),
        ("SG", "S", 2, 4.0, -3.0, 0.0),
        ("HG", "H", 2, 4.0, -3.5, 0.8),
    ]
    mol = Molecule().empty(len(rows))
    mol.name[:] = [r[0] for r in rows]
    mol.element[:] = [r[1] for r in rows]
    mol.resid[:] = [r[2] for r in rows]
    mol.resname[:] = "CYS"
    mol.segid[:] = "P"
    mol.chain[:] = "A"
    mol.record[:] = "ATOM"
    mol.coords = np.array([r[3:] for r in rows], dtype=np.float32)[:, :, np.newaxis]

    disulfide, _, _, _ = _prepareMolecule(
        mol,
        caps={},
        disulfide=[["resid 1 and name SG", "resid 2 and name SG"]],
        custombonds=None,
        remove=None,
    )
    assert len(disulfide) == 1
    assert set(mol.resname[mol.resid == 1].tolist()) == {"CYX"}
    assert set(mol.resname[mol.resid == 2].tolist()) == {"CYX"}
    assert (mol.name == "HG").sum() == 0


def test_amber_anchor_rename_via_custombonds(tmp_path):
    """Defense-in-depth: when a user passes a custombond from CYS SG to a
    scaffold atom directly to ``_prepareMolecule`` (skipping detect), the
    rename to CYX still happens and HG is dropped."""
    mol = Molecule(QFZ_B_CIF)

    # Add a dummy HG to one anchored Cys so the H-drop step has work to do.
    sg_idx = int(mol.atomselect("resid 11 and name SG", indexes=True)[0])
    sg_pos = mol.coords[sg_idx, :, mol.frame]
    new = Molecule().empty(1)
    new.name[:] = "HG"
    new.element[:] = "H"
    new.resname[:] = "CYS"
    new.resid[:] = 11
    new.chain[:] = mol.chain[sg_idx]
    new.segid[:] = mol.segid[sg_idx]
    new.coords = (sg_pos + np.array([0.0, 0.0, 1.34])).reshape(1, 3, 1).astype(np.float32)
    new.record[:] = "ATOM"
    mol.insert(new, index=sg_idx + 1, collisions=False)

    mol.segid[:] = "P"

    # Build custombonds by hand (simulating a user who skipped detect).
    cb = []
    for resid in (11, 17, 22):
        sg_uid = UniqueAtomID.fromMolecule(mol, f"resid {resid} and name SG")
        # find the LFI atom bonded to this SG
        sg_i = int(mol.atomselect(f"resid {resid} and name SG", indexes=True)[0])
        for a, b in mol.bonds:
            if int(a) == sg_i and str(mol.resname[int(b)]) == "LFI":
                lfi_uid = UniqueAtomID.fromMolecule(mol, idx=int(b))
                break
            if int(b) == sg_i and str(mol.resname[int(a)]) == "LFI":
                lfi_uid = UniqueAtomID.fromMolecule(mol, idx=int(a))
                break
        cb.append(
            (
                f'segid "{sg_uid.segid}" and resid {sg_uid.resid} and name "SG"',
                f'segid "{lfi_uid.segid}" and resid {lfi_uid.resid} and name "{lfi_uid.name}"',
            )
        )
    assert len(cb) == 3

    sg_uids_before = [
        UniqueAtomID.fromMolecule(mol, f"resid {r} and name SG") for r in (11, 17, 22)
    ]

    _prepareMolecule(mol, caps={}, disulfide=None, custombonds=cb, remove=None)

    for uid in sg_uids_before:
        sel = (
            (mol.chain == uid.chain)
            & (mol.segid == uid.segid)
            & (mol.insertion == uid.insertion)
            & (mol.name == "SG")
        )
        sg_idx = np.where(sel)[0]
        assert len(sg_idx) >= 1
        cyx_sg = sg_idx[mol.resname[sg_idx] == "CYX"]
        assert len(cyx_sg) >= 1, (uid.resid, mol.resname[sg_idx])
        resid_val = mol.resid[cyx_sg[0]]
        chain_val = mol.chain[cyx_sg[0]]
        segid_val = mol.segid[cyx_sg[0]]
        res_mask = (
            (mol.resid == resid_val)
            & (mol.chain == chain_val)
            & (mol.segid == segid_val)
        )
        assert ((mol.name == "HG") & res_mask).sum() == 0


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_8qu4_stapled_peptide_end_to_end(tmp_path):
    """8QU4 chain A: 13-mer NF-Y-derived stapled peptide. NLE272.CE -
    MK8276.CE single bond (the staple after RCM closure) joins two NCAAs
    with no canonical anchor.

    Asserts:
      - detect emits two CrosslinkedNCAASpec.
      - The build's prmtop carries exactly one NLE.CE - MK8.CE bond.
      - Each staple CE is a CH2 (exactly 2 H neighbours), not a CH3.
    """
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build

    mol = Molecule(QU4_A_CIF)
    mol.segid[:] = "P"

    # 1. Detect (no canonical anchors here, so no mutation).
    specs = detectNonStandardResidues(mol)
    assert sum(isinstance(s, ChainResidueSpec) and s.resname not in PROTEIN_RESNAMES and s.anchor_atom is not None for s in specs) == 2
    assert sum(isinstance(s, ChainResidueSpec) and s.resname in PROTEIN_RESNAMES for s in specs) == 0

    # 2. Template the two NCAAs from SMILES.
    mol.templateResidueFromSmiles("resname NLE", NLE_SMILES, addHs=True)
    mol.templateResidueFromSmiles("resname MK8", MK8_SMILES, addHs=True)

    # 3. systemPrepare (preserves bonds across PDB2PQR via bond-capture).
    pmol, _ = systemPrepare(mol, detect_specs=[])

    # 4. Parameterize.
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )
    topo_basenames = {os.path.basename(p) for p in out.topo_paths}
    frcmod_basenames = {os.path.basename(p) for p in out.frcmod_paths}
    assert topo_basenames == {"NLE.prepi", "MK8.prepi"}
    assert frcmod_basenames == {"NLE.frcmod", "MK8.frcmod"}
    assert len(out.custombonds) == 1

    # 5. Build.
    builddir = str(tmp_path / "build")
    built = amber_build(
        pmol,
        outdir=builddir,
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        caps={"P": ("none", "none")},
    )
    assert built is not None

    # Each staple CE must be a CH2 (2 H neighbours + 2 heavy = CD + staple).
    for resname in ("NLE", "MK8"):
        ce_idxs = np.where((built.resname == resname) & (built.name == "CE"))[0]
        assert len(ce_idxs) == 1
        ce = int(ce_idxs[0])
        n_h = sum(
            1 for nb in built.getNeighbors(ce)
            if str(built.element[int(nb)]) == "H"
        )
        n_heavy = sum(
            1 for nb in built.getNeighbors(ce)
            if str(built.element[int(nb)]) != "H"
        )
        assert n_h == 2, f"{resname}.CE expected 2 H neighbours, got {n_h}"
        assert n_heavy == 2

    # The staple bond must be in the built mol.
    nle_ce = set(
        np.where((built.resname == "NLE") & (built.name == "CE"))[0].tolist()
    )
    mk8_ce = set(
        np.where((built.resname == "MK8") & (built.name == "CE"))[0].tolist()
    )
    n_staple = sum(
        1 for a, b in built.bonds
        if (int(a) in nle_ce and int(b) in mk8_ce)
        or (int(b) in nle_ce and int(a) in mk8_ce)
    )
    assert n_staple == 1
    _check_no_overvalent_atoms(built)


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_8qfz_scaffolded_peptide_end_to_end(tmp_path):
    """8QFZ chain B: 12-mer cyclic peptide whose CYS11/17/22 sidechains
    are each thio-ether bonded to one of the LFI scaffold's three anchor
    carbons (a tricyclic peptide). CYS 11 is N-terminal, CYS 22 is
    C-terminal, CYS 17 is mid-chain - so they end up in three distinct
    chain-position buckets and get three distinct custom resnames.

    Asserts:
      - detect emits one ScaffoldSpec(LFI) plus three CanonicalRenamedSpec
        with three distinct new resnames (one per chain-position bucket).
      - parameterizeFromSpecs writes one CIF + frcmod per unique resname
        (3 CYS-bucket CIFs + 1 LFI.cif) and three custombonds.
      - Each CYS-bucket CIF carries ff14SB atom types on the standard
        canonical sidechain atoms.
      - The junction frcmod contains cross-FF (ff14SB / GAFF2) bonds.
    """
    from moleculekit.tools.preparation import systemPrepare
    from parmed.amber import AmberParameterSet

    mol = Molecule(QFZ_B_CIF)

    # 1. Detect: each CYS is in its own chain-position bucket, so all
    # three CYS get distinct rename targets. HG drops where present.
    specs = detectNonStandardResidues(mol)
    scaffolds = [s for s in specs if isinstance(s, ScaffoldSpec)]
    renames = [s for s in specs if isinstance(s, ChainResidueSpec) and s.resname in PROTEIN_RESNAMES]
    assert len(scaffolds) == 1 and scaffolds[0].resname == "LFI"
    assert len(renames) == 3
    cys_rename_set = {r.new_resname for r in renames}
    assert len(cys_rename_set) == 3, (
        f"expected three distinct CYS rename targets, got {cys_rename_set}"
    )
    for cys_new in cys_rename_set:
        assert len(cys_new) == 3 and cys_new.startswith("XX")

    # 2. Template LFI from SMILES.
    mol.templateResidueFromSmiles("resname LFI", LFI_SMILES, addHs=True)

    # 3. systemPrepare runs PDB2PQR on the canonical naming, then applies
    # the rename + displaced-H drop from the spec list at the end.
    pmol, _ = systemPrepare(
        mol,
        detect_specs=specs,
    )

    # 4. Parameterize.
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )
    topo_basenames = sorted(os.path.basename(p) for p in out.topo_paths)
    expected = sorted([f"{n}.prepi" for n in cys_rename_set] + ["LFI.cif"])
    assert topo_basenames == expected
    assert len(out.custombonds) == 3

    # LFI CIF retains GAFF2 (lowercase) types since it's non-canonical.
    lfi_cif = next(p for p in out.topo_paths if "LFI" in os.path.basename(p))
    lfi_mol = Molecule(lfi_cif)
    assert any(str(t).islower() for t in lfi_mol.atomtype)

    # Junction frcmod must contain at least one cross-FF bond entry
    # for the backbone-sidechain boundary: backbone atoms (N, CA, C,
    # O, H, HA) carry ff14SB types (``N``, ``CT``, ``C``, ``O``,
    # ``H``, ``H1``); the sidechain (including SG) stays GAFF2
    # (lowercase). The cluster pipeline must emit bond entries that
    # bridge the two type families so tLeap can resolve the
    # backbone-sidechain bonds.
    pset = AmberParameterSet(out.frcmod_paths[0])
    ff14sb_backbone_types = {"N", "H", "CT", "H1", "C", "O"}
    cross_ff_bonds = [
        k for k in pset.bond_types
        if any(t in ff14sb_backbone_types for t in k)
        and any(str(t).islower() for t in k)
    ]
    assert cross_ff_bonds, "junction frcmod missing cross-FF bond entries"
    # The N-CA and CA-CB bonds in particular must span the FF boundary
    # since CA is backbone (CT) and CB is sidechain (lowercase GAFF2).
    assert any(
        "CT" in k and any(str(t).islower() for t in k) for k in cross_ff_bonds
    ), "expected at least one CT-<gaff2> bond from the CA-CB junction"

    # 5. Build. Both ends of the protein chain are real terminal
    # residues (the N-terminal CYS is bonded to LFI, the C-terminal
    # CYS already carries OXT) so disable amber.build's auto-capping
    # there - otherwise it would add an ACE / NME that conflicts with
    # the existing terminal atoms.
    from htmd.builder.amber import build as amber_build

    protein_segid = sorted(set(pmol.segid[pmol.atomselect("protein")]))[0]
    builddir = str(tmp_path / "build")
    built = amber_build(
        pmol,
        outdir=builddir,
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        caps={protein_segid: ("none", "none")},
    )
    assert built is not None

    # Three SG-Cn thioether bonds in the built mol.
    sg_idxs = set(np.where(built.name == "SG")[0].tolist())
    lfi_c_idxs = set(
        np.where(
            (built.resname == "LFI")
            & np.isin(built.name, ["C10", "C11", "C12"])
        )[0].tolist()
    )
    n_thioether = sum(
        1 for a, b in built.bonds
        if (int(a) in sg_idxs and int(b) in lfi_c_idxs)
        or (int(b) in sg_idxs and int(a) in lfi_c_idxs)
    )
    assert n_thioether == 3, f"expected 3 SG-Cn bonds, got {n_thioether}"
    _check_no_overvalent_atoms(built)


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_full_pipeline_5vbl(tmp_path):
    """5VBL: chain A peptide inhibitor with five chain-resident NCAAs
    (200, ALC, HRG, NLE, OIC) and chain B protein with a free OLC ligand.
    No covalent crosslinks. The build must succeed and contain no
    over-valent atoms (no over-protonated CH3 -> CH4 etc.)."""
    mol = Molecule(VBL_PDB)

    # Chain A NCAAs are peptide-bonded on at least one side; templating
    # uses the mid-chain SMILES (carbonyl as ``C=O``) used by
    # ``moleculekit.tests.test_systemprepare`` for the same fixture.
    # OLC is a free ligand on chain B with a free hydroxyl.
    smiles = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
        "OLC": "CCCCCCCC(=O)OC[C@H](O)CO",
    }
    # Disable amber.build's auto-capping for the inhibitor segment - its
    # C-terminal NCAA (200) is parameterized with its own prepi (which
    # already includes the OXT atom), so an NME cap on top would clash.
    # ```` lets systemPrepare build a PDB2PQR template for
    # residue 200 from the spec, which stops PDB2PQR from C-terminal-capping
    # PRO 16 with a phantom OXT.
    built = _run_pipeline(
        mol,
        smiles,
        tmp_path,
        build_kwargs={"caps": {"P0": ("none", "none")}},
    )
    assert built is not None
    _check_no_overvalent_atoms(built)


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_full_pipeline_1m63(tmp_path):
    """1M63: calcineurin / cyclophilin in complex with cyclosporin. Covers
    three things the simpler fixtures do not:

      * a free iron ion (resname ``FE``) coordinating Asp/His in the
        calcineurin Fe/Zn binuclear centre. ``FE`` is not in the ion-resname
        list, so detection has to recognise it as a free metal ion from its
        element-symbol resname - otherwise the metal-coordination bond is
        mistaken for a covalent crosslink. tleap then parameterizes it (and
        the Zn / Ca ions) from AMBER's own ion library.
      * the cyclic, fully N-methylated cyclosporin undecapeptide. Its six
        NCAAs include sarcosine (``SAR``), an N-methyl-glycine whose alpha
        carbon is a CH2 carrying two hydrogens (named HA2 / HA3).
      * a head-to-tail cyclic peptide bond, which amber.build closes with no
        custombond.
    """
    mol = Molecule("1M63")
    mol.remove("water", _logger=False)

    # Mid-chain SMILES (carbonyl as C=O, backbone N as bare N; N-methylated
    # residues carry the extra methyl as "NC").
    smiles = {
        "DAL": "C[C@H](C=O)N",                     # D-alanine
        "ABA": "CC[C@@H](C=O)N",                   # L-2-aminobutyrate
        "SAR": "O=CCNC",                           # sarcosine (N-methyl-Gly)
        "MLE": "CC(C)C[C@@H](C=O)NC",              # N-methyl-L-leucine
        "MVA": "CC(C)[C@@H](C=O)NC",               # N-methyl-L-valine
        "BMT": "O=C[C@@H](NC)C(O)C(C)C/C=C/C",     # (4R)-MeBmt
    }
    built = _run_pipeline(mol, smiles, tmp_path)
    assert built is not None
    _check_no_overvalent_atoms(built)

    # Metals survive as correctly-typed free ions (FE is the key one).
    for resname, count, charge in (("FE", 2, 3.0), ("ZN", 2, 2.0), ("CA", 8, 2.0)):
        sel = built.resname == resname
        assert int(sel.sum()) == count, (
            f"{resname}: {int(sel.sum())} atoms, expected {count}"
        )
        assert np.allclose(built.charge[sel], charge), (
            f"{resname} charge {np.unique(built.charge[sel])}, expected {charge}"
        )

    # The two cyclosporin copies are each closed head-to-tail (and nothing
    # else in the system is): exactly two backbone ring closures, each
    # spanning the 11-residue macrocycle.
    closures = _backbone_ring_closures(built)
    assert len(closures) == 2, f"expected 2 cyclic closures, got {closures}"
    for ra, rb in closures:
        assert abs(ra - rb) == 10, (
            f"cyclosporin ring should span 11 residues, closure {ra}-{rb}"
        )


try:
    import openmm  # noqa: F401

    _openmm = True
except ImportError:
    _openmm = False


@pytest.mark.skipif(
    not (_antechamber and _tleap and _openmm),
    reason="OpenMM build comparison needs antechamber + teLeap + openmm",
)
def test_full_pipeline_5vbl_openmm_vs_amber(tmp_path):
    """5VBL: build the prepared mol twice, once via amber.build (tLeap
    consuming the prepi/frcmod outputs) and once via openff.build (OpenMM
    consuming the emitted ForceField XML), with solvation and ionization
    off so the two systems should carry the same atoms. Compares numAtoms
    and resname distributions to verify the XML emitter produces
    something openff.build can actually consume."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build
    from htmd.builder.openmm import build as openff_build
    from collections import Counter

    mol = Molecule(VBL_PDB)
    mol = autoSegment(mol, fields=("segid", "chain"), _logger=False)
    mol.remove("element H", _logger=False)
    mol.remove("water", _logger=False)
    # TODO: drop ZN until the openff builder can match the Zn2+ template
    # against a topology that still carries the metal-coordination bonds
    # from the input PDB's CONECT records.
    mol.remove("resname ZN", _logger=False)
    specs = detectNonStandardResidues(mol)
    smiles = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
        "OLC": "CCCCCCCC(=O)OC[C@H](O)CO",
    }
    for resname, smi in smiles.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol, _ = systemPrepare(
        mol,
        detect_specs=specs,
    )
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )

    # Same cap override as _test_full_pipeline_5vbl: the C-terminal NCAA
    # (residue 200) carries its own OXT via the prepi/XML template, so an
    # NME cap on top would clash. P0 is the inhibitor segment.
    caps = {"P0": ("none", "none")}

    amber_built = amber_build(
        pmol.copy(),
        outdir=str(tmp_path / "amber"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        caps=caps,
    )
    openmm_built, _ = openff_build(
        pmol.copy(),
        outdir=str(tmp_path / "openmm"),
        extra_xml=list(out.xml_paths),
        custombonds=out.custombonds,
        ionize=False,
        solvate=False,
        caps=caps,
    )

    # NCAA + Zn + 4 disulfides; observed gap ~0.024 kcal/mol on
    # 6000 atoms. 0.05 leaves ~2x headroom.
    _assert_builds_equivalent(
        amber_built,
        openmm_built,
        amber_prmtop=str(tmp_path / "amber" / "structure.prmtop"),
        openmm_prmtop=str(tmp_path / "openmm" / "structure.prmtop"),
        energy_tol=0.05,
    )


@pytest.mark.skipif(
    not (_tleap and _openmm),
    reason="OpenMM-vs-amber build comparison needs teLeap + openmm",
)
def test_full_pipeline_6a5j_openmm_vs_amber(tmp_path):
    """6A5J: small canonical peptide (no NCAAs, no ligands). Builds via
    amber.build (tLeap) and openff.build (OpenMM ForceField XML) and
    asserts the two systems carry identical force-field data via
    :func:`_assert_builds_equivalent`. Acts as a control for the NCAA
    builds - if 6A5J disagrees, the divergence is in the basic builder
    path rather than the parameterizeFromSpecs cluster pipeline."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build
    from htmd.builder.openmm import build as openff_build

    mol = Molecule("6A5J")
    mol = autoSegment(mol, fields=("segid", "chain"), _logger=False)
    mol.remove("element H", _logger=False)
    mol.remove("water", _logger=False)
    pmol, _ = systemPrepare(mol, verbose=False)

    amber_built = amber_build(
        pmol.copy(),
        outdir=str(tmp_path / "amber"),
        ionize=False,
    )
    openmm_built, _ = openff_build(
        pmol.copy(),
        outdir=str(tmp_path / "openmm"),
        ionize=False,
        solvate=False,
    )

    # Canonical-only ~270-atom peptide; observed gap is ~5e-5 kcal/mol.
    # 0.001 still gives ~20x headroom and catches any real regression.
    _assert_builds_equivalent(
        amber_built,
        openmm_built,
        amber_prmtop=str(tmp_path / "amber" / "structure.prmtop"),
        openmm_prmtop=str(tmp_path / "openmm" / "structure.prmtop"),
        energy_tol=0.001,
    )


@pytest.mark.skipif(
    not (_antechamber and _openmm),
    reason="OpenMM XML test needs antechamber + openmm",
)
def test_parameterize_from_specs_emits_openmm_xml(tmp_path):
    """parameterizeFromSpecs should also emit a single OpenMM ForceField
    XML covering every residue in the run. Loading it alongside ff14SB
    must produce a ForceField whose residue templates include every
    NCAA / free ligand we parameterized."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare
    import openmm.app as app

    mol = Molecule(VBL_PDB)
    mol = autoSegment(mol, fields=("segid", "chain"), _logger=False)
    mol.remove("element H", _logger=False)
    mol.remove("water", _logger=False)
    specs = detectNonStandardResidues(mol)
    smiles = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
        "OLC": "CCCCCCCC(=O)OC[C@H](O)CO",
    }
    for resname, smi in smiles.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol, _ = systemPrepare(
        mol,
        detect_specs=specs,
    )
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )

    # GAFF path emits exactly one combined OpenMM XML appended to
    # ``xml_paths`` (named ``gaff_combined.xml``); locate it.
    gaff_xml = next(
        (p for p in out.xml_paths if p.endswith("gaff_combined.xml")), None
    )
    assert gaff_xml is not None, (
        f"GAFF run should have appended a gaff_combined.xml to xml_paths, "
        f"got {out.xml_paths}"
    )
    assert os.path.isfile(gaff_xml)
    assert os.path.getsize(gaff_xml) > 0

    ff = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/tip3p.xml",
        gaff_xml,
    )

    # Each per-residue topo file should correspond to a residue template
    # in the loaded ForceField - check the resnames match.
    expected_resnames = {
        os.path.splitext(os.path.basename(p))[0] for p in out.topo_paths
    }
    loaded_resnames = set(ff._templates)
    missing = expected_resnames - loaded_resnames
    assert not missing, (
        f"Residue templates missing from the loaded ForceField: {missing}. "
        f"Loaded: {sorted(loaded_resnames & expected_resnames)}"
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_full_pipeline_4tot_e(tmp_path):
    """4TOT_E: cyclosporin-like cyclic peptide with seven N-methylated /
    backbone-modified NCAAs. The build must succeed and contain no
    over-valent atoms."""
    mol = Molecule(TOT_E_PDB)

    # Mid-chain SMILES (carbonyl as ``C=O`` rather than free-acid
    # ``C(=O)O``) since every residue in this cyclic peptide is
    # peptide-bonded on both sides. These match the SMILES strings used
    # by ``moleculekit.tests.test_systemprepare`` for the same fixture.
    smiles = {
        "33X": "CC(C=O)NC",
        "34E": "CN[C@@H]([C@H](C)CN1CCN(CCOC)CC1)C=O",
        "ABA": "CC[C@H](C=O)N",
        "BMT": "C/C=C/C[C@@H](C)[C@H]([C@@H](C=O)NC)O",
        "DAL": "C[C@H](C=O)N",
        "MLE": "CC(C)C[C@@H](C=O)NC",
        "MVA": "CC(C)[C@@H](C=O)NC",
    }
    built = _run_pipeline(mol, smiles, tmp_path)
    assert built is not None
    _check_no_overvalent_atoms(built)


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def test_full_pipeline_1r1j(tmp_path):
    """1R1J: glycoprotein with three NAG-Asn N-glycosylation sites. Each
    Asn ND2 - NAG C1 bond produces a CovalentLigandSpec(NAG) plus a
    CanonicalRenamedSpec(ASN). All three sites share one (ASN, ND2, NAG,
    mid-chain) bucket so they collapse to a single per-bucket prepi."""
    mol = Molecule(R1J_PDB)

    smiles = {
        "NAG": "CC(=O)N[C@@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",
        # OIR: N-(3-phenyl-2-sulfanylpropanoyl)phenylalanylalanine, the
        # zinc-bound peptidic inhibitor (RCSB chem-comp OIR).
        "OIR": "CC(C(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(Cc2ccccc2)S",
    }
    built = _run_pipeline(mol, smiles, tmp_path)
    assert built is not None
    _check_no_overvalent_atoms(built)


@pytest.mark.skipif(
    not _antechamber,
    reason="dedup test calls parameterizeFromSpecs which needs antechamber",
)
def test_parameterize_from_specs_dedup_4tot(tmp_path):
    """4TOT carries cyclosporin chains with three MLE residues per
    chain across multiple chains. parameterizeFromSpecs should run
    antechamber once per unique ``(resname, is_n_term, is_c_term)`` and
    emit one prepi/frcmod per unique NCAA, not one per occurrence:
    tLeap loads a single unit per resname anyway, and a duplicate
    ``loadamberprep MLE.prepi`` would silently overwrite the previous."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.nonstandard_residues import detectNonStandardResidues
    from moleculekit.tools.preparation import systemPrepare

    smiles = {
        "33X": "CC(C=O)NC",
        "34E": "CN[C@@H]([C@H](C)CN1CCN(CCOC)CC1)C=O",
        "ABA": "CC[C@H](C=O)N",
        "BMT": "C/C=C/C[C@@H](C)[C@H]([C@@H](C=O)NC)O",
        "DAL": "C[C@H](C=O)N",
        "MLE": "CC(C)C[C@@H](C=O)NC",
        "MVA": "CC(C)[C@@H](C=O)NC",
        # P6G: hexaethylene glycol, a free crystallisation-additive ligand
        # (RCSB chem-comp P6G). SO4: free sulfate dianion (also a
        # crystallisation additive, not coordinating anything in 4TOT);
        # templated with explicit -1 on each terminal O so the residue's
        # net charge is the correct -2. Without the template the PDB
        # records the four S-O bonds with no formal charges set, leaving
        # RDKit Gasteiger summing ~-0.5 and tripping the netcharge sanity
        # check.
        "P6G": "OCCOCCOCCOCCOCCOCCO",
        "SO4": "[O-]S(=O)(=O)[O-]",
    }
    mol = Molecule("4TOT")
    mol = autoSegment(mol, fields=("segid", "chain"), _logger=False)
    mol.remove("element H", _logger=False)

    specs = detectNonStandardResidues(mol)
    spec_counts = Counter(s.residue.resname for s in specs)
    # 4TOT carries 4 cyclosporin chains (each with 3 MLE plus one each
    # of 33X / 34E / ABA / BMT / DAL / MVA), 4 free P6G ligands, and
    # 2 free SO4 ligands. Without dedup that's 42 antechamber runs.
    expected_spec_counts = {
        "33X": 4, "34E": 4, "ABA": 4, "BMT": 4, "DAL": 4,
        "MLE": 12, "MVA": 4, "P6G": 4, "SO4": 2,
    }
    assert dict(spec_counts) == expected_spec_counts, (
        f"unexpected spec counts: {dict(spec_counts)}"
    )

    for resname, smi in smiles.items():
        if (mol.resname == resname).any():
            mol.templateResidueFromSmiles(
                f'resname "{resname}"', smi, addHs=True, _logger=False
            )
    pmol, _ = systemPrepare(
        mol,
        detect_specs=specs,
    )
    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"), charge_method="gasteiger"
    )

    topo_basenames = sorted(
        os.path.splitext(os.path.basename(p))[0] for p in out.topo_paths
    )
    frcmod_basenames = sorted(
        os.path.splitext(os.path.basename(p))[0] for p in out.frcmod_paths
    )
    # 42 specs collapse to 9 unique-resname prepi/frcmod files: 7 NCAA
    # singletons (one per unique resname) plus 2 free ligands (P6G, SO4).
    expected_basenames = [
        "33X", "34E", "ABA", "BMT", "DAL", "MLE", "MVA", "P6G", "SO4"
    ]
    assert topo_basenames == expected_basenames, (
        f"topo basenames {topo_basenames} != {expected_basenames}"
    )
    assert frcmod_basenames == expected_basenames, (
        f"frcmod basenames {frcmod_basenames} != {expected_basenames}"
    )


# ---------------------------------------------------------------------------
# Unit tests for the protonation sanity check (_check_specs_protonated).
# These build tiny synthetic molecules; no antechamber / tleap needed.
# ---------------------------------------------------------------------------


def _toy_mol(residues):
    """Build a multi-residue Molecule from a compact description and return
    ``(mol, groups)``.

    ``residues`` is a list of ``(resname, atoms, bonds, formalcharges)``:
      - ``atoms``    : list of element strings (one per atom in the residue)
      - ``bonds``    : list of ``(i, j)`` index pairs *local* to the residue
      - ``formalcharges`` : optional list of ints (defaults to all zero)
    """
    elements, names, resnames, resids, fcs = [], [], [], [], []
    bonds = []
    offset = 0
    for resid, (resname, atoms, res_bonds, *rest) in enumerate(residues, start=1):
        res_fc = rest[0] if rest else [0] * len(atoms)
        for k, el in enumerate(atoms):
            elements.append(el)
            names.append(f"{el}{offset + k}")
            resnames.append(resname)
            resids.append(resid)
            fcs.append(int(res_fc[k]))
        for i, j in res_bonds:
            bonds.append((offset + i, offset + j))
        offset += len(atoms)

    n = len(elements)
    mol = Molecule().empty(n)
    mol.element[:] = elements
    mol.name[:] = names
    mol.resname[:] = resnames
    mol.resid[:] = resids
    mol.formalcharge[:] = fcs
    mol.coords = np.zeros((n, 3, 1), dtype=np.float32)
    mol.bonds = (
        np.array(bonds, dtype=np.uint32).reshape(-1, 2)
        if bonds
        else np.zeros((0, 2), dtype=np.uint32)
    )
    _, groups = _residue_groups_with_index(mol)
    return mol, groups


def _alkane(resname, n_carbon):
    """Fully protonated straight-chain alkane CnH(2n+2) as a residue tuple."""
    atoms = ["C"] * n_carbon
    bonds = [(i, i + 1) for i in range(n_carbon - 1)]
    h = n_carbon
    for c in range(n_carbon):
        n_h = 3 if c in (0, n_carbon - 1) else 2
        for _ in range(n_h):
            atoms.append("H")
            bonds.append((c, h))
            h += 1
    return (resname, atoms, bonds)


# Reusable residue fragments.
_HEXANE = _alkane("HEX", 6)          # fully protonated, exercises the H-count path
_BARE_C4 = ("BC4", ["C"] * 4, [(0, 1), (1, 2), (2, 3)])  # 4-carbon chain, no H
_BENZENE = ("BNZ", ["C"] * 6 + ["H"] * 6,
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
            + [(i, 6 + i) for i in range(6)])
_ZINC = ("ZN", ["Zn"], [])
_SULFATE = ("SO4", ["S", "O", "O", "O", "O"], [(0, 1), (0, 2), (0, 3), (0, 4)],
            [0, 0, 0, -1, -1])
_SMALL_CO2 = ("CO2", ["C", "O", "O"], [(0, 1), (0, 2)])  # H-free, below the gate
# Lysine-ish heavy skeleton with only the backbone N-H / CA-HA present.
_LYS_STRIPPED = (
    "LYX",
    ["N", "C", "C", "O", "C", "C", "C", "C", "N", "H", "H"],
    [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 6), (6, 7), (7, 8),
     (0, 9), (1, 10)],
)


def test_check_specs_protonated_passes_aliphatic():
    mol, groups = _toy_mol([_HEXANE])
    _check_specs_protonated(mol, {0: None}, groups)  # must not raise


def test_check_specs_protonated_passes_aromatic():
    # The estimator treats the ring as all single bonds, so benzene's 6 H sit
    # well above the 0.2 * 12 threshold.
    mol, groups = _toy_mol([_BENZENE])
    _check_specs_protonated(mol, {0: None}, groups)  # must not raise


def test_check_specs_protonated_zero_hydrogens_raises():
    mol, groups = _toy_mol([_BARE_C4, _HEXANE])
    with pytest.raises(RuntimeError, match=r"under-protonated"):
        _check_specs_protonated(mol, {0: None}, groups)


def test_check_specs_protonated_stripped_sidechain_raises():
    # Backbone H present, every sidechain H missing -> still flagged.
    mol, groups = _toy_mol([_LYS_STRIPPED])
    with pytest.raises(RuntimeError) as exc:
        _check_specs_protonated(mol, {0: None}, groups)
    assert "LYX" in str(exc.value)
    assert "2 hydrogen(s) present" in str(exc.value)


def test_check_specs_protonated_ignores_carbonless_residues():
    # A metal ion and an oxo-anion carry no expected H even with 0 H present.
    mol, groups = _toy_mol([_HEXANE, _ZINC, _SULFATE])
    _check_specs_protonated(mol, {1: None, 2: None}, groups)  # must not raise


def test_check_specs_protonated_ignores_tiny_carbon_species():
    # CO2 estimates only ~4 expected H -> below _MIN_EXPECTED_HYDROGENS, so
    # being H-free does not trip the check.
    mol, groups = _toy_mol([_HEXANE, _SMALL_CO2])
    _check_specs_protonated(mol, {1: None}, groups)  # must not raise


def test_check_specs_protonated_reports_every_bad_residue():
    mol, groups = _toy_mol([_BARE_C4, _HEXANE, _LYS_STRIPPED])
    with pytest.raises(RuntimeError) as exc:
        _check_specs_protonated(mol, {0: None, 2: None}, groups)
    msg = str(exc.value)
    assert "BC4" in msg and "LYX" in msg
    # The fully protonated residue (index 1) is not in the spec set and must
    # not appear regardless.
    assert "HEX" not in msg


def test_check_specs_protonated_no_bonds_is_noop():
    # Defensive guard: without connectivity there is nothing to check.
    mol, groups = _toy_mol([("LIG", ["C", "C", "C"], [])])
    _check_specs_protonated(mol, {0: None}, groups)  # must not raise


# ---------------------------------------------------------------------------
# Unit tests for the frcmod parameter cleanup (_clean_frcmod_params).
# ---------------------------------------------------------------------------

# A parmchk2-style frcmod with deliberate cruft: an atom type used by nothing
# (xx), a type only present on a "cap" atom (cc), a backbone-only bond/angle
# (the protein FF provides those), and bond/angle/dihedral entries the test
# molecule's connectivity never realizes (c3-c3-c3 angle, c3-c3-c3-c3 dihe).
_DIRTY_FRCMOD = """test frcmod
MASS
c3 12.010
ca 12.010
hc 1.008
cc 12.010
N  14.010
CT 12.010
xx 99.000

BOND
c3-c3   300.0   1.526
c3-ca   300.0   1.510
ca-hc   340.0   1.090
c3-cc   300.0   1.500
N-CT    300.0   1.460
N-c3    300.0   1.470
xx-xx   100.0   1.500

ANGLE
c3-c3-ca   50.0   109.5
ca-ca-hc   50.0   120.0
c3-c3-c3   50.0   109.5
N-CT-c3    50.0   110.0
xx-xx-xx   50.0   109.5

DIHE
c3-c3-ca-ca   1   0.156   0.0   3.0
c3-c3-c3-c3   1   0.180   0.0   3.0
N-CT-c3-c3    1   0.100   0.0   3.0
xx-xx-xx-xx   1   0.000   0.0   2.0

IMPROPER
ca-ca-ca-ca   1.1   180.0   2.0
N-CT-c3-cc    1.0   180.0   2.0
xx-xx-xx-xx   1.0   180.0   2.0

NONBON
c3   1.9080   0.1094
ca   1.9080   0.0860
hc   1.4870   0.0157
cc   1.9080   0.0860
N    1.8240   0.1700
CT   1.9080   0.1094
xx   1.0000   0.0000
"""


def _typed_mol(elements, atomtypes, bonds):
    """Single-residue Molecule with explicit atom types + bonds."""
    n = len(elements)
    mol = Molecule().empty(n)
    mol.element[:] = elements
    mol.name[:] = [f"{e}{i}" for i, e in enumerate(elements)]
    mol.atomtype[:] = atomtypes
    mol.resname[:] = "LIG"
    mol.resid[:] = 1
    mol.coords = np.zeros((n, 3, 1), dtype=np.float32)
    mol.bonds = (
        np.array(bonds, dtype=np.uint32).reshape(-1, 2)
        if bonds
        else np.zeros((0, 2), dtype=np.uint32)
    )
    return mol


def _parse_frcmod(tmp_path, text=_DIRTY_FRCMOD):
    from parmed.amber import AmberParameterSet

    fp = tmp_path / "in.frcmod"
    fp.write_text(text)
    return AmberParameterSet(str(fp))


def test_clean_frcmod_free_residue(tmp_path):
    # ca-ca-ca-ca-... straight chain: c3(0)-c3(1)-ca(2)-ca(3)-hc(4).
    mol = _typed_mol(
        ["C", "C", "C", "C", "H"],
        ["c3", "c3", "ca", "ca", "hc"],
        [(0, 1), (1, 2), (2, 3), (3, 4)],
    )
    pset = _parse_frcmod(tmp_path)
    _clean_frcmod_params(pset, mol, [], np.zeros(mol.numAtoms, dtype=bool))

    # Only types the molecule actually carries survive.
    assert set(pset.atom_types) == {"c3", "ca", "hc"}
    # Realized bonds kept; unused-type and unrealized ones gone.
    assert ("c3", "c3") in pset.bond_types
    assert ("c3", "ca") in pset.bond_types
    assert ("ca", "hc") in pset.bond_types
    assert ("c3", "cc") not in pset.bond_types  # cc absent from mol
    assert ("xx", "xx") not in pset.bond_types
    assert ("N", "c3") not in pset.bond_types  # N absent from mol
    # Realized angles/dihedrals only.
    assert ("c3", "c3", "ca") in pset.angle_types
    assert ("c3", "c3", "c3") not in pset.angle_types  # not realized
    assert ("c3", "c3", "ca", "ca") in pset.dihedral_types
    assert ("c3", "c3", "c3", "c3") not in pset.dihedral_types
    # Impropers: kept if all types present, dropped otherwise.
    assert ("ca", "ca", "ca", "ca") in pset.improper_periodic_types
    assert ("xx", "xx", "xx", "xx") not in pset.improper_periodic_types
    assert ("N", "CT", "c3", "cc") not in pset.improper_periodic_types

    # Each kept dihedral key has no reversed-duplicate sibling.
    for kind in ("bond_types", "angle_types", "dihedral_types"):
        keys = list(getattr(pset, kind))
        for k in keys:
            assert tuple(reversed(k)) == k or tuple(reversed(k)) not in keys

    # The result must still serialize as a valid frcmod.
    from parmed.amber import AmberParameterSet

    out = tmp_path / "out.frcmod"
    pset.write(str(out), title="cleaned", style="frcmod")
    AmberParameterSet(str(out))  # re-parse, must not raise


def test_clean_frcmod_backbone_terms_dropped(tmp_path):
    # A chain-resident-NCAA-ish skeleton: N(0)-CT(1)-c3(2)-c3(3). N-CT is
    # purely backbone (the protein FF provides it) and must go; N-CT-c3 /
    # N-CT-c3-c3 span the backbone-sidechain boundary, are realized here, and
    # must stay; N-c3 is not realized in this fragment but touches a backbone
    # type, so it is kept conservatively (in the built system a proline-like
    # NCAA does bond N to a sidechain c3, and that torsion is realized over a
    # cap atom in the model compound).
    mol = _typed_mol(
        ["N", "C", "C", "C"],
        ["N", "CT", "c3", "c3"],
        [(0, 1), (1, 2), (2, 3)],
    )
    pset = _parse_frcmod(tmp_path)
    backbone = ["N", "H", "CT", "H1", "C", "O"]
    _clean_frcmod_params(pset, mol, backbone, np.zeros(mol.numAtoms, dtype=bool))

    assert ("N", "CT") not in pset.bond_types  # all-backbone -> dropped
    assert ("N", "c3") in pset.bond_types  # backbone-touching -> kept
    assert ("N", "CT", "c3") in pset.angle_types  # realized, mixed -> kept
    assert ("N", "CT", "c3", "c3") in pset.dihedral_types
    assert ("c3", "c3", "c3") not in pset.angle_types  # not realized, no bb
    assert ("c3", "c3", "c3", "c3") not in pset.dihedral_types
    # Backbone atom types themselves stay (referenced by mixed terms / atoms);
    # types absent from the fragment and not backbone go.
    assert "CT" in pset.atom_types and "N" in pset.atom_types and "c3" in pset.atom_types
    assert "ca" not in pset.atom_types and "xx" not in pset.atom_types


def test_clean_frcmod_cap_atoms_excluded(tmp_path):
    # cc sits only on a "cap"/padding atom -> every cc parameter must go.
    mol = _typed_mol(
        ["C", "C", "C"],
        ["c3", "c3", "cc"],
        [(0, 1), (1, 2)],
    )
    padding = np.array([False, False, True])
    pset = _parse_frcmod(tmp_path)
    _clean_frcmod_params(pset, mol, [], padding)

    assert "cc" not in pset.atom_types
    assert ("c3", "cc") not in pset.bond_types
    assert ("c3", "c3") in pset.bond_types


# ---------------------------------------------------------------------------
# Unit tests for the backbone charge pin (_normalize_residue_charges) and the
# ff14SB backbone charge lookup (_backbone_charge_map).
# ---------------------------------------------------------------------------


def _charge_mol(names, charges):
    """Single-residue Molecule with explicit atom names and partial
    charges - the inputs _normalize_residue_charges reads."""
    n = len(names)
    mol = Molecule().empty(n)
    mol.name[:] = names
    mol.element[:] = [str(nm)[0] for nm in names]
    mol.charge[:] = charges
    mol.resname[:] = "XX1"
    mol.resid[:] = 1
    mol.coords = np.zeros((n, 3, 1), dtype=np.float32)
    return mol


# Mid-chain backbone (N, H, CA, HA, C, O) plus a short sidechain.
_BACKBONE_PLUS_SIDECHAIN = ["N", "H", "CA", "HA", "C", "O", "CB", "HB2", "HB3", "SG"]
# ff14SB neutral-class backbone amide charges, used as a sample map.
_NEUTRAL_BACKBONE = _FF14SB_BACKBONE_CHARGES_BY_CLASS[0]


def test_normalize_residue_charges_pins_map_and_rebalances():
    # The charge_map atoms carry the pinned values and the residue sums
    # to its integer formal charge.
    sub = _charge_mol(_BACKBONE_PLUS_SIDECHAIN, [0.3] * len(_BACKBONE_PLUS_SIDECHAIN))
    _normalize_residue_charges(sub, net_charge=0, charge_map=_NEUTRAL_BACKBONE)
    for name, q in _NEUTRAL_BACKBONE.items():
        assert float(sub.charge[sub.name == name][0]) == pytest.approx(q, abs=1e-4)
    assert float(np.sum(sub.charge)) == pytest.approx(0.0, abs=1e-4)


def test_gasteiger_charges_renormalized_to_integer_formal_charge():
    """RDKit's Gasteiger (PEOE) does not always converge to an integer total.
    For 1FJM's microcystin DAM/FGA cluster (formal charge 0) it lands ~0.11 e
    short, even though every atom is parameterized and the per-atom charges are
    otherwise reasonable (it conserves exactly for simpler molecules). The helper
    must spread that residual so the charges sum to the integer formal charge -
    MD requires an integer total, and the downstream net-charge check expects it.
    """
    from htmd.builder._charge_helpers import _assign_rdkit_gasteiger_charges

    mol = Molecule(os.path.join(DATA_DIR, "1FJM_microcystin_cluster.cif"))
    target = int(round(float(np.sum(mol.formalcharge))))
    _assign_rdkit_gasteiger_charges(mol)
    assert float(np.sum(mol.charge)) == pytest.approx(target, abs=1e-4)


def test_normalize_residue_charges_targets_integer_formal_charge():
    # A residue with a non-zero formal charge is rebalanced to that
    # integer, not to zero.
    sub = _charge_mol(_BACKBONE_PLUS_SIDECHAIN, [-0.05] * len(_BACKBONE_PLUS_SIDECHAIN))
    _normalize_residue_charges(sub, net_charge=-1, charge_map=_NEUTRAL_BACKBONE)
    for name, q in _NEUTRAL_BACKBONE.items():
        assert float(sub.charge[sub.name == name][0]) == pytest.approx(q, abs=1e-4)
    assert float(np.sum(sub.charge)) == pytest.approx(-1.0, abs=1e-4)


def test_normalize_residue_charges_residual_spares_pinned_atoms():
    # The residual lands only on the non-pinned atoms, equally (the
    # minimum-L2 redistribution): pinned atoms keep exactly their value.
    names = _BACKBONE_PLUS_SIDECHAIN
    sub = _charge_mol(names, [0.7] * len(names))
    _normalize_residue_charges(sub, net_charge=0, charge_map=_NEUTRAL_BACKBONE)
    for name, q in _NEUTRAL_BACKBONE.items():
        assert float(sub.charge[sub.name == name][0]) == pytest.approx(q, abs=1e-4)
    free = ~np.isin(np.array(names), list(_NEUTRAL_BACKBONE))
    shifts = np.array(sub.charge)[free] - 0.7
    assert np.allclose(shifts, shifts[0], atol=1e-5)
    assert float(np.sum(sub.charge)) == pytest.approx(0.0, abs=1e-4)


def test_normalize_residue_charges_empty_map_normalises_total():
    # No atoms to pin -> the residual relative to net_charge is spread
    # equally across all atoms so the residue sums to net_charge.
    # This is the scaffold path: no backbone to pin, but the per-residue
    # total still has to be its integer formal charge.
    names = _BACKBONE_PLUS_SIDECHAIN
    n = len(names)
    sub = _charge_mol(names, [0.05] * n)  # sums to 0.5
    _normalize_residue_charges(sub, net_charge=0, charge_map={})
    # 0.5 of residual spread across n=10 atoms -> -0.05 per atom -> 0.0
    assert float(np.sum(sub.charge)) == pytest.approx(0.0, abs=1e-5)
    assert np.allclose(np.array(sub.charge), 0.0, atol=1e-5)


def test_normalize_residue_charges_empty_map_targets_integer():
    # Empty map + non-zero target: same equal-share redistribution,
    # rebalances to the residue's integer formal charge.
    names = _BACKBONE_PLUS_SIDECHAIN
    n = len(names)
    sub = _charge_mol(names, [0.0] * n)
    _normalize_residue_charges(sub, net_charge=-1, charge_map={})
    assert float(np.sum(sub.charge)) == pytest.approx(-1.0, abs=1e-5)
    assert np.allclose(np.array(sub.charge), -1.0 / n, atol=1e-5)


def test_normalize_residue_charges_pin_only_skips_normalization():
    # net_charge=None applies the pin (backbone atoms get the charge_map
    # values) but does NOT redistribute the residual onto the free
    # atoms. Used by the cluster-level normalize path, which combines
    # per-residue pinning with one cluster-wide shift at the end.
    names = _BACKBONE_PLUS_SIDECHAIN
    sub = _charge_mol(names, [0.3] * len(names))
    pinned = _normalize_residue_charges(sub, net_charge=None, charge_map=_NEUTRAL_BACKBONE)
    # Pinned atoms got the map values.
    for name, q in _NEUTRAL_BACKBONE.items():
        assert float(sub.charge[sub.name == name][0]) == pytest.approx(q, abs=1e-4)
    # Free atoms kept their original Gasteiger charges (no residual spread).
    free_mask = ~np.isin(np.array(names), list(_NEUTRAL_BACKBONE))
    assert np.allclose(np.array(sub.charge)[free_mask], 0.3, atol=1e-5)
    # Returned pinned-mask matches the names in _NEUTRAL_BACKBONE.
    expected_pinned = np.isin(np.array(names), list(_NEUTRAL_BACKBONE))
    assert np.array_equal(pinned, expected_pinned)


def test_backbone_charge_map_ncaa_charge_classes():
    # An NCAA is absent from the libraries, so the amide charges come
    # from the charge class the residue's net charge selects.
    # The fallback path doesn't actually touch ff14sb_lib for NCAAs
    # (it consults _FF14SB_BACKBONE_CHARGES_BY_CLASS), so an empty
    # dict is fine for the lib argument here.
    present = {"N", "H", "CA", "HA", "C", "O", "CB", "CG"}
    assert _backbone_charge_map("", "", present, 0, ff14sb_lib={}) == _FF14SB_BACKBONE_CHARGES_BY_CLASS[0]
    assert _backbone_charge_map("", "", present, 1, ff14sb_lib={}) == _FF14SB_BACKBONE_CHARGES_BY_CLASS[1]
    assert _backbone_charge_map("", "", present, -1, ff14sb_lib={}) == _FF14SB_BACKBONE_CHARGES_BY_CLASS[-1]
    # A net charge with no standard backbone class pins nothing.
    assert _backbone_charge_map("", "", present, 2, ff14sb_lib={}) == {}


def test_backbone_charge_map_terminal_ncaa_is_empty():
    # A terminal NCAA has no library entry and no universal terminal
    # charges, so nothing is pinned.
    present = {"N", "H", "CA", "HA", "C", "O", "OXT"}
    assert _backbone_charge_map("", "c", present, -1, ff14sb_lib={}) == {}


@pytest.mark.skipif(not _tleap, reason="ff14SB amino lib lookup needs teLeap")
def test_backbone_charge_map_proline_like_ncaa_uses_pro():
    # A proline-like NCAA (N with no bonded H) gets its backbone
    # pinned to ff14SB PRO charges - the only canonical residue with
    # this topology - per the Cornell/Cieplak + Betz/Ramos
    # convention. CD is left to refit because ring chemistry varies
    # (OIC vs. pipecolic vs. true PRO).
    present = {"N", "CA", "HA", "C", "O", "CB", "CG", "CD"}
    cmap = _backbone_charge_map(
        "", "", present, 0, n_has_bonded_h=False,
        ff14sb_lib=_amber_lib_or_none(),
    )
    assert set(cmap) == {"N", "CA", "C", "O", "HA"}
    # ff14SB PRO library values (Maier 2015 / Cornell 1995).
    assert cmap["N"] == pytest.approx(-0.2548, abs=1e-4)
    assert cmap["CA"] == pytest.approx(-0.0266, abs=1e-4)
    assert cmap["C"] == pytest.approx(0.5896, abs=1e-4)
    assert cmap["O"] == pytest.approx(-0.5748, abs=1e-4)
    assert cmap["HA"] == pytest.approx(0.0641, abs=1e-4)


def test_backbone_charge_map_proline_like_terminal_is_empty():
    # A proline-like NCAA at a chain terminus is left untouched -
    # ff14SB has no terminal-PRO entry that strips H, and applying
    # the mid-chain PRO charges to a terminal residue would mismatch
    # the terminal ionisation. Falls through to the safe empty map.
    present = {"N", "CA", "HA", "C", "O", "OXT", "CB", "CG", "CD"}
    assert (
        _backbone_charge_map(
            "", "c", present, -1, n_has_bonded_h=False, ff14sb_lib={}
        )
        == {}
    )


def test_backbone_charge_map_nterm_with_nh3_not_proline_like():
    # An N-terminal residue with NH3+ atoms named H1/H2/H3 has bonded
    # H atoms but no atom literally named "H" in present_names. The
    # name-based "H absent" check would mis-classify it as proline-
    # like; the bond-based ``n_has_bonded_h`` flag does not.
    present = {"N", "CA", "HA", "C", "O", "H1", "H2", "H3", "CB", "CG", "CD"}
    cmap = _backbone_charge_map(
        "", "", present, 0, n_has_bonded_h=True, ff14sb_lib={}
    )
    # Without the n_has_bonded_h flag, this could have wrongly fallen
    # into the proline-like PRO map; with it, the residue is treated
    # as a normal amide NCAA and the universal charge-class fallback
    # decides the map (empty here because "H" is not in present_names
    # so the {N, H, C, O} set isn't satisfied).
    assert cmap == {}


@pytest.mark.skipif(not _tleap, reason="ff14SB amino lib lookup needs teLeap")
def test_backbone_charge_map_canonical_midchain():
    # A mid-chain canonical CYS gets its whole backbone - residue-specific
    # CA / HA included - from the ff14SB CYS library entry; the sidechain
    # is left to antechamber.
    present = {"N", "H", "CA", "HA", "C", "O", "CB", "HB2", "HB3", "SG"}
    cmap = _backbone_charge_map(
        "CYS", "", present, 0, ff14sb_lib=_amber_lib_or_none()
    )
    assert set(cmap) == {"N", "H", "CA", "HA", "C", "O"}
    for name, q in _NEUTRAL_BACKBONE.items():
        assert cmap[name] == pytest.approx(q, abs=1e-4)


@pytest.mark.skipif(not _tleap, reason="ff14SB amino lib lookup needs teLeap")
def test_backbone_charge_map_canonical_terminal_variants():
    # The N-terminal variant pins the NH3+ hydrogens; the C-terminal
    # variant pins OXT. Sidechain atoms are never pinned.
    lib = _amber_lib_or_none()
    n_present = {"N", "H1", "H2", "H3", "CA", "HA", "C", "O", "CB"}
    n_map = _backbone_charge_map("ALA", "n", n_present, 1, ff14sb_lib=lib)
    assert {"H1", "H2", "H3"} <= set(n_map)
    assert "CB" not in n_map

    c_present = {"N", "H", "CA", "HA", "C", "O", "OXT", "CB"}
    c_map = _backbone_charge_map("ALA", "c", c_present, -1, ff14sb_lib=lib)
    assert "OXT" in c_map
    assert "CB" not in c_map


@pytest.mark.skipif(not _tleap, reason="ff14SB amino lib lookup needs teLeap")
def test_backbone_charge_map_his_name_falls_through_to_charge_class():
    # "HIS" is not a library key (the libraries use HID / HIE / HIP), so
    # a His anchor falls through to the charge-class fallback: a neutral
    # His gets the neutral amide, a protonated (+1) His the cationic one.
    lib = _amber_lib_or_none()
    present = {"N", "H", "CA", "HA", "C", "O", "CB", "CG", "ND1", "CE1", "NE2"}
    assert (
        _backbone_charge_map("HIS", "", present, 0, ff14sb_lib=lib)
        == _FF14SB_BACKBONE_CHARGES_BY_CLASS[0]
    )
    assert (
        _backbone_charge_map("HIS", "", present, 1, ff14sb_lib=lib)
        == _FF14SB_BACKBONE_CHARGES_BY_CLASS[1]
    )


# ---------------------------------------------------------------------------
# Forcefield dispatch helpers (engine selection, per-resname routing).
# ---------------------------------------------------------------------------


def test_engine_for_forcefield_gaff_variants_dispatch_to_antechamber():
    # The rule is "starts with 'gaff' (case-insensitive) -> antechamber,
    # everything else -> openff". Lock it so a future change doesn't
    # silently route GAFF FFs through the SMIRNOFF engine or vice-versa.
    assert _engine_for_forcefield("gaff") == "antechamber"
    assert _engine_for_forcefield("gaff2") == "antechamber"
    assert _engine_for_forcefield("GAFF2") == "antechamber"
    assert _engine_for_forcefield("Gaff") == "antechamber"


def test_engine_for_forcefield_offxml_dispatches_to_openff():
    assert _engine_for_forcefield("openff-2.3.0.offxml") == "openff"
    assert _engine_for_forcefield("openff_unconstrained-2.3.0.offxml") == "openff"
    # An offxml filename that doesn't start with "openff" still routes to
    # the SMIRNOFF engine - the rule is "anything not GAFF goes to openff".
    assert _engine_for_forcefield("custom_sage.offxml") == "openff"


def test_engine_for_forcefield_non_string_raises():
    with pytest.raises(TypeError):
        _engine_for_forcefield(None)
    with pytest.raises(TypeError):
        _engine_for_forcefield(2.3)


def test_resolve_forcefield_single_string_applies_everywhere():
    lookup = _resolve_forcefield("gaff2")
    assert lookup("BEN") == ("antechamber", "gaff2")
    assert lookup("ALC") == ("antechamber", "gaff2")
    lookup = _resolve_forcefield("openff_unconstrained-2.3.0.offxml")
    assert lookup("BEN") == ("openff", "openff_unconstrained-2.3.0.offxml")


def test_resolve_forcefield_dict_with_default_routes_per_resname():
    lookup = _resolve_forcefield(
        {
            "NLE": "openff_unconstrained-2.3.0.offxml",
            "default": "gaff2",
        }
    )
    assert lookup("NLE") == ("openff", "openff_unconstrained-2.3.0.offxml")
    # Unlisted resnames fall back to the explicit default.
    assert lookup("ALC") == ("antechamber", "gaff2")


def test_resolve_forcefield_dict_without_default_falls_back_to_gaff2():
    lookup = _resolve_forcefield({"NLE": "openff_unconstrained-2.3.0.offxml"})
    # The default is gaff2 when not specified.
    assert lookup("ALC") == ("antechamber", "gaff2")
    assert lookup("NLE") == ("openff", "openff_unconstrained-2.3.0.offxml")


def test_resolve_forcefield_rejects_non_string_values():
    with pytest.raises(TypeError):
        _resolve_forcefield({"NLE": 2.3})  # validated up front
    with pytest.raises(TypeError):
        _resolve_forcefield(2.3)  # not str or dict


@pytest.fixture
def htmd_logs_to_caplog():
    """htmd's logger has ``propagate=0`` in ``logging.ini`` so its
    records never reach pytest's caplog root handler. Re-enable
    propagation for the duration of the test and restore it after."""
    import logging
    htmd_logger = logging.getLogger("htmd")
    original = htmd_logger.propagate
    htmd_logger.propagate = True
    try:
        yield
    finally:
        htmd_logger.propagate = original


def test_warn_if_openff_mismatched_charges_warns_on_openff_plus_non_nagl(
    caplog, htmd_logs_to_caplog
):
    # SMIRNOFF + non-NAGL fires the warning.
    with caplog.at_level("WARNING"):
        _warn_if_openff_mismatched_charges(
            "openff_unconstrained-2.3.0.offxml", "gasteiger"
        )
    assert any("NAGL reproduces" in r.message for r in caplog.records)


def test_warn_if_openff_mismatched_charges_silent_on_openff_plus_nagl(
    caplog, htmd_logs_to_caplog
):
    # SMIRNOFF + NAGL is the recommended combination - no warning.
    with caplog.at_level("WARNING"):
        _warn_if_openff_mismatched_charges(
            "openff_unconstrained-2.3.0.offxml", "nagl"
        )
    assert not any("NAGL reproduces" in r.message for r in caplog.records)


def test_warn_if_openff_mismatched_charges_silent_on_gaff(
    caplog, htmd_logs_to_caplog
):
    # GAFF-only call - no SMIRNOFF in play, no warning regardless of
    # charge method. The charge-FF mismatch concern only applies to Sage.
    with caplog.at_level("WARNING"):
        _warn_if_openff_mismatched_charges("gaff2", "gasteiger")
        _warn_if_openff_mismatched_charges("gaff2", "resp")
    assert not any("NAGL reproduces" in r.message for r in caplog.records)


def test_warn_if_openff_mismatched_charges_warns_on_mixed_dict(
    caplog, htmd_logs_to_caplog
):
    # A per-resname dict containing at least one SMIRNOFF entry triggers
    # the warning when paired with a non-NAGL charge model.
    with caplog.at_level("WARNING"):
        _warn_if_openff_mismatched_charges(
            {"NLE": "openff_unconstrained-2.3.0.offxml", "default": "gaff2"},
            "gasteiger",
        )
    assert any("NAGL reproduces" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# NAGL charges with the GAFF2 backend. The Gasteiger pattern in
# _fftype_antechamber (pre-compute, set ac_charge=None, antechamber only
# types) was extended to nagl/resp/resp-multiconf. NAGL is the only one
# we can exercise without parameterize / Psi4.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (_antechamber and _nagl),
    reason="NAGL + GAFF2 needs antechamber on PATH and the openff-nagl stack",
)
def test_fftype_antechamber_with_nagl_charges(tmp_path):
    from moleculekit.molecule import Molecule
    from htmd.builder._ambertools import _fftype_antechamber

    # Protonated benzamidinium, the 3PTB ligand. Templating from the
    # SMILES locks bond orders + explicit hydrogens so antechamber can
    # type cleanly, and locks the +1 formal charge so NAGL respects it.
    mol = Molecule("3PTB")
    mol.filter("resname BEN", _logger=False)
    mol.templateResidueFromSmiles("all", BEN_SMILES, addHs=True, _logger=False)
    assert int(mol.formalcharge.sum()) == 1
    # Pre-condition: partial charges start at zero. If a future moleculekit
    # change populates mol.charge from formal charges or some other source,
    # the post-typing "charges differ from input" check below would no
    # longer detect a silent no-op in _assign_nagl_charges.
    assert np.all(mol.charge == 0.0), (
        "templateResidueFromSmiles should leave mol.charge at zeros; "
        "this test relies on the zero start state to detect a NAGL no-op"
    )

    typed, _, _ = _fftype_antechamber(
        mol,
        tmpdir=str(tmp_path),
        forcefield="GAFF2",
        netcharge=1,
        charge_method="nagl",
    )

    # GAFF2 atom typing should still run despite the externally-computed
    # charges (antechamber called with -c absent, only -at gaff2).
    # Spot-check: protonated benzamidinium should pick up the GAFF2 types
    # for sp2 aromatic carbons (ca), the central sp2 amidinium C (ce),
    # aromatic ring H (ha), amide H on N (hn), and the cationic
    # amidinium nitrogen (nv).
    assert set(typed.atomtype) >= {"ca", "ce", "ha", "hn", "nv"}

    # NAGL honours the formal charge, so the per-atom charges must sum
    # to the +1 amidinium net charge.
    assert abs(typed.charge.sum() - 1.0) < 1e-3
    # Charges actually got computed (not a NAGL no-op silently returning
    # zeros that then get rebalanced to sum to formal charge somehow).
    # Every atom in protonated benzamidinium sits in a polar / aromatic
    # / amidinium environment - no atom should have charge exactly 0.
    assert np.all(typed.charge != 0.0), (
        f"At least one atom has charge exactly 0, suggesting NAGL "
        f"was a no-op: {typed.charge}"
    )
    # Sanity bound on individual charges (avoids regressing to garbage).
    assert -1.0 < typed.charge.min() < typed.charge.max() < 1.0


# ---------------------------------------------------------------------------
# ff14SB amino-lib parity between parmed (AmberTools amino*.lib) and
# OpenMM (amber14/protein.ff14SB.xml). The two backends use these
# independently now, but parity is the chemistry-correctness invariant
# - if either upstream drifts, this drift detector trips.
# ---------------------------------------------------------------------------


# ACE/NME methyl atom-name aliases between parmed and OpenMM. Both
# libraries describe the same physical atoms; only the labels differ.
_FF14SB_CAP_NAME_ALIASES = {
    ("ACE", "H1"): "HH31",
    ("ACE", "H2"): "HH32",
    ("ACE", "H3"): "HH33",
    ("NME", "C"): "CH3",
    ("NME", "H1"): "HH31",
    ("NME", "H2"): "HH32",
    ("NME", "H3"): "HH33",
}


@pytest.mark.skipif(not _tleap, reason="parmed amino*.lib loader needs teLeap")
def test_ff14sb_amino_lib_parmed_vs_openmm_parity():
    """Assert the parmed-loaded (amino*.lib) and OpenMM-loaded
    (protein.ff14SB.xml) ff14SB amino libraries produce byte-identical
    ``(type, charge)`` for every common atom. Drift detector for either
    upstream package - if antechamber and openff build paths ever
    silently disagree on backbone-pin charges, this trips first."""
    parmed_lib = _load_ff14sb_amino_lib_amber()
    omm_lib = _load_ff14sb_amino_lib_openmm()

    # Same set of residues (78 mid-chain + N/C-terminal variants).
    assert set(parmed_lib) == set(omm_lib)

    # OpenMM types are prefixed with "protein-" in the XML; the loader
    # strips that prefix so the returned shape matches parmed's bare-
    # class type names. Verify both bare-class types and bit-identical
    # charges for every common atom.
    n_compared = 0
    for resname, p_atoms in parmed_lib.items():
        o_atoms = omm_lib[resname]
        for atom_name, (p_type, p_charge) in p_atoms.items():
            # Resolve the ACE/NME cap-atom alias on the parmed side
            # so the lookup finds the OpenMM-side atom.
            lookup_name = _FF14SB_CAP_NAME_ALIASES.get(
                (resname, atom_name), atom_name
            )
            if lookup_name not in o_atoms:
                continue  # tolerate a missing alias gracefully
            o_type, o_charge = o_atoms[lookup_name]
            assert p_type == o_type, (
                f"{resname}.{atom_name}: parmed type {p_type!r} != "
                f"openmm type {o_type!r}"
            )
            assert p_charge == o_charge, (
                f"{resname}.{atom_name}: parmed charge {p_charge!r} != "
                f"openmm charge {o_charge!r}"
            )
            n_compared += 1
    # Sanity: we should have actually compared a non-trivial number of
    # atoms. ff14SB has ~1275 atoms across the 78 residues.
    assert n_compared > 1000


# ---------------------------------------------------------------------------
# Golden test: parameterize fixed inputs and compare the per-residue
# parameterization files against committed reference copies.
# ---------------------------------------------------------------------------

_CUSTOM_PARAM_DIR = os.path.join(curr_dir, "test-custom-residue-param")


def _check_against_reference(generated_paths, ref_dir, label, regenerate):
    """Compare each generated parameterization file against a committed
    copy in ``ref_dir``, matched by basename. CIF topologies (free
    residues and scaffolds) are compared as Molecules with moleculekit's
    ``mol_equal``; prepi topologies (chain-resident residues) and frcmod
    parameter files are compared verbatim. With ``regenerate`` true,
    (re)write the references from the generated files instead."""
    generated = {os.path.basename(p): p for p in generated_paths}

    if regenerate:
        if os.path.isdir(ref_dir):
            shutil.rmtree(ref_dir)
        os.makedirs(ref_dir)
        for name, path in generated.items():
            shutil.copyfile(path, os.path.join(ref_dir, name))
        return

    assert os.path.isdir(ref_dir), f"{label}: missing reference dir {ref_dir}"
    assert set(generated) == set(os.listdir(ref_dir)), (
        f"{label}: produced files {sorted(generated)} != "
        f"reference files {sorted(os.listdir(ref_dir))}"
    )
    for name, path in sorted(generated.items()):
        ref_path = os.path.join(ref_dir, name)
        if name.endswith(".cif"):
            assert mol_equal(
                Molecule(path), Molecule(ref_path)
            ), f"{label}: {name} differs from its reference"
        else:
            with open(path) as fh:
                got = fh.read()
            with open(ref_path) as fh:
                ref = fh.read()
            assert got == ref, f"{label}: {name} differs from its reference"


def _run_8qfz_bicycle_param_reference(
    tmp_path, normalize, ref_subdir, pin_backbone_charges=True,
):
    """Shared driver for the 8QFZ_B_bicycle scaffolded cyclic peptide
    golden test under each ``normalize`` / ``pin_backbone_charges``
    combination. Templates LFI, runs systemPrepare +
    parameterizeFromSpecs with the given modes, and compares the
    emitted CIF/prepi/frcmod files against the committed references in
    ``ref_subdir``."""
    from moleculekit.tools.preparation import systemPrepare

    regenerate = bool(os.environ.get("HTMD_REGEN_REFERENCES"))
    mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, "8QFZ_B_bicycle.cif"))
    if not np.any(mol.chain != ""):
        mol.chain[:] = "A"
        mol.segid[:] = "A"
    specs = detectNonStandardResidues(mol)
    mol.templateResidueFromSmiles('resname "LFI"', LFI_SMILES, addHs=True)
    pmol = systemPrepare(mol, detect_specs=specs)[0]
    out = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / ref_subdir),
        charge_method="gasteiger",
        normalize=normalize,
        pin_backbone_charges=pin_backbone_charges,
    )
    _check_against_reference(
        out.topo_paths + out.frcmod_paths,
        os.path.join(_CUSTOM_PARAM_DIR, "reference", ref_subdir),
        ref_subdir,
        regenerate,
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference_8qfz_per_residue(tmp_path):
    """8QFZ_B_bicycle under ``normalize='per_residue'``: every emitted
    per-residue unit (XX1/XX2/XX3 prepi, LFI cif) sums to its integer
    formal charge - the AMBER tLeap convention (Betz / R.E.D. / Ramos).
    """
    _run_8qfz_bicycle_param_reference(
        tmp_path,
        normalize="per_residue",
        ref_subdir="8QFZ_B_bicycle_per_residue",
    )


def _run_5vbl_param_reference(
    tmp_path, normalize, ref_subdir, pin_backbone_charges=True,
):
    """Shared driver for the 5VBL_A peptide reference tests. The
    fixture covers five chain-resident NCAAs (HRG, ALC, OIC, NLE, 200)
    plus an isopeptide GLU(CD)-LYS(NZ) crosslink that exercises the
    canonical-anchor rename path. 200 is C-terminal."""
    from moleculekit.tools.preparation import systemPrepare

    regenerate = bool(os.environ.get("HTMD_REGEN_REFERENCES"))
    mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    # SMILES are the mid-chain "C=O" form for the NCAAs that sit
    # between peptide bonds, and the explicit "C(=O)O" carboxyl for
    # residue 200 at the C-terminus. Same set as the existing 5VBL
    # full-pipeline tests but minus OLC (a chain-B free ligand absent
    # from the chain-A fixture).
    smiles = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles.items():
        mol.templateResidueFromSmiles(
            f'resname "{resname}"', smi, addHs=True, _logger=False
        )
    # The chain-A peptide inhibitor in 5VBL_A.cif is missing the
    # sidechain atoms of LYS1 / PHE2 / ARG3 / ARG4 / LYS12 (only
    # backbone + CB), so PDB2PQR's 10%-missing-atom threshold would
    # reject the structure. ``restore_missing_sidechains=True`` runs
    # moleculekit's Dunbrack-rotamer mutator pre-PDB2PQR to fill them
    # in at canonical coordinates.
    pmol = systemPrepare(
        mol, detect_specs=specs, restore_missing_sidechains=True,
    )[0]
    out = parameterizeFromSpecs(
        specs,
        pmol,
        outdir=str(tmp_path / ref_subdir),
        charge_method="gasteiger",
        normalize=normalize,
        pin_backbone_charges=pin_backbone_charges,
    )
    _check_against_reference(
        out.topo_paths + out.frcmod_paths,
        os.path.join(_CUSTOM_PARAM_DIR, "reference", ref_subdir),
        ref_subdir,
        regenerate,
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference_5vbl_cluster(tmp_path):
    """5VBL_A under the default ``normalize='cluster'``: the cluster
    total is normalised to integer; per-residue totals stay at their
    natural Gasteiger values modulo a uniform shift."""
    _run_5vbl_param_reference(
        tmp_path, normalize="cluster", ref_subdir="5VBL_A_cluster",
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference_5vbl_per_residue(tmp_path):
    """5VBL_A under ``normalize='per_residue'``: every emitted unit
    sums to its integer formal charge (AMBER tLeap convention)."""
    _run_5vbl_param_reference(
        tmp_path, normalize="per_residue", ref_subdir="5VBL_A_per_residue",
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference_5vbl_no_pin_no_normalize(tmp_path):
    """5VBL_A under ``pin_backbone_charges=False`` and
    ``normalize=None``: no backbone pin, no rebalance. The per-atom
    charges are exactly the RDKit Gasteiger output and the cluster
    total should land near the formal-charge sum modulo antechamber's
    4-decimal mol2 rounding."""
    _run_5vbl_param_reference(
        tmp_path,
        normalize=None,
        pin_backbone_charges=False,
        ref_subdir="5VBL_A_no_pin_no_normalize",
    )
    refdir = os.path.join(
        _CUSTOM_PARAM_DIR, "reference", "5VBL_A_no_pin_no_normalize"
    )

    def _prepi_total(p):
        s = 0.0
        ina = False
        for line in open(p):
            ls = line.strip()
            if ls.startswith("CORR"):
                ina = True
                continue
            if ls.startswith(("LOOP", "IMPROPER", "DONE", "CHARGE")):
                ina = False
            if ina and len(line.split()) >= 10:
                s += float(line.split()[-1])
        return s

    # Sum every emitted unit (.prepi for chain-resident; .cif for any
    # free residues - 5VBL_A has no scaffolds, all units are .prepi).
    cluster_total = 0.0
    for fname in os.listdir(refdir):
        path = os.path.join(refdir, fname)
        if fname.endswith(".prepi"):
            cluster_total += _prepi_total(path)
        elif fname.endswith(".cif"):
            cluster_total += float(Molecule(path).charge.sum())
    # 5VBL has 7 chain-resident specs which split into 6 separate
    # antechamber runs (one XX1+XX2 multi-cluster + 5 singletons), so
    # the natural Gasteiger smear into the (discarded) ACE/NME-style
    # cap atoms compounds across clusters. With no pin and no
    # normalize, expect drift up to ~0.1; the reference-file diff
    # locks each per-atom charge, so this side-check just guards
    # against gross regressions (e.g. a future change silently
    # turning normalize back on).
    assert abs(cluster_total) < 0.15, (
        f"5VBL_A unpinned + unnormalised cluster total drifted by "
        f"{cluster_total:+.4f} (expected within rounding + cap smear)"
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference_8qfz_no_pin_no_normalize(tmp_path):
    """8QFZ_B_bicycle under ``pin_backbone_charges=False`` and
    ``normalize=None``: no backbone pin, no rebalance. The per-atom
    charges are exactly the RDKit Gasteiger output and the cluster
    total should land near the integer formal-charge sum (~0) modulo
    antechamber's 4-decimal mol2 rounding. Verifies that the rounding
    drift is small (< 0.05) and locks the raw charges down as a
    reference."""
    _run_8qfz_bicycle_param_reference(
        tmp_path,
        normalize=None,
        pin_backbone_charges=False,
        ref_subdir="8QFZ_B_bicycle_no_pin_no_normalize",
    )
    # Side-check: with neither backbone pin nor normalization, the
    # cluster total drift is only antechamber's mol2 rounding (a few
    # thousandths) - well below the looser bound we assert here.
    refdir = os.path.join(
        _CUSTOM_PARAM_DIR, "reference", "8QFZ_B_bicycle_no_pin_no_normalize"
    )

    def _prepi_total(p):
        s = 0.0
        ina = False
        for line in open(p):
            ls = line.strip()
            if ls.startswith("CORR"):
                ina = True
                continue
            if ls.startswith(("LOOP", "IMPROPER", "DONE", "CHARGE")):
                ina = False
            if ina and len(line.split()) >= 10:
                s += float(line.split()[-1])
        return s

    cluster_total = sum(
        _prepi_total(os.path.join(refdir, f"{r}.prepi")) for r in ("XX1", "XX2", "XX3")
    )
    cluster_total += float(Molecule(os.path.join(refdir, "LFI.cif")).charge.sum())
    assert abs(cluster_total) < 0.05, (
        f"unpinned + unnormalised cluster total drifted by {cluster_total:+.4f} "
        f"(expected near 0 from mol2 rounding only)"
    )


def _read_residue_template_atoms(path):
    """Parse ``path`` (a ``.prepi`` or ``.cif`` residue template emitted
    by parameterizeFromSpecs) and return ``{atom_name: (atom_type, charge)}``
    for every real atom (dummies excluded)."""
    out = {}
    if path.endswith(".cif"):
        mol = Molecule(path)
        for i in range(mol.numAtoms):
            out[str(mol.name[i])] = (
                str(mol.atomtype[i]),
                float(mol.charge[i]),
            )
        return out
    # prepi atom block lives between ``CORR`` and the next section header
    in_atoms = False
    for line in open(path):
        ls = line.strip()
        if ls.startswith("CORR"):
            in_atoms = True
            continue
        if ls.startswith(("LOOP", "IMPROPER", "DONE", "CHARGE")):
            in_atoms = False
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        name, atype = parts[1], parts[2]
        if name == "DUMM":
            continue
        out[name] = (atype, float(parts[-1]))
    return out


def _read_prepi_intra_residue_impropers(path):
    """Parse the ``IMPROPER`` section of a prepi: return a set of
    4-atom-name frozensets, restricted to **intra-residue** impropers.
    AMBER prepi uses ``-M`` / ``+M`` placeholders for main-chain atoms
    of the previous / next residue - those are cross-residue
    impropers (backbone amide / carbonyl planarity at the peptide
    bond) and naturally don't show up in our intra-residue scan of
    the built prmtop, so we drop them here. Frozenset is the right
    comparator because AMBER impropers (a,b,c,d) and (d,c,b,a) refer
    to the same torsion around the central atom."""
    out = set()
    section = None
    for line in open(path):
        ls = line.strip()
        if ls.startswith("IMPROPER"):
            section = "improper"
            continue
        if ls.startswith(("CORR", "LOOP", "DONE", "CHARGE")) or not ls:
            section = None
            continue
        if section == "improper":
            parts = line.split()
            if len(parts) == 4 and not any(p in ("-M", "+M") for p in parts):
                out.add(frozenset(parts))
    return out


def _intra_residue_impropers(built, resname):
    """All impropers in ``built`` whose four atoms lie inside the same
    residue with name ``resname``, as a set of 4-atom-name frozensets.
    Multi-instance resnames are pooled (every copy uses the same
    template, so the union over copies must equal the reference set).
    Returns ``None`` if the resname is not present."""
    if built.impropers is None or len(built.impropers) == 0:
        return None
    sel = built.resname == resname
    if not sel.any():
        return None
    out = set()
    for row in built.impropers:
        idxs = [int(x) for x in row]
        if not all(sel[i] for i in idxs):
            continue
        # Require all four atoms in the SAME residue (the prepi defines
        # one residue; cross-residue impropers are tLeap's business).
        keys = {
            (
                str(built.segid[i]),
                str(built.chain[i]),
                int(built.resid[i]),
                str(built.insertion[i]),
            )
            for i in idxs
        }
        if len(keys) != 1:
            continue
        out.add(frozenset(str(built.name[i]) for i in idxs))
    return out


def _assert_built_matches_references(built, refdir, charge_tol=1e-3):
    """For every per-residue reference file in ``refdir``, check that
    each atom carries the same ``atomtype`` and ``charge`` in the
    AMBER-built Molecule, and (for prepi templates) that the
    intra-residue impropers in the built prmtop match the prepi's
    ``IMPROPER`` block. Bonds / angles / dihedrals are not compared
    explicitly because tLeap derives them deterministically from the
    bond graph the prepi defines - if any bond were wrong, the build
    would fail or the angle/dihedral counts would diverge, but the
    direct round-trip check is on impropers (explicit in prepi)."""
    bad = []
    for fname in sorted(os.listdir(refdir)):
        name, ext = os.path.splitext(fname)
        if ext not in (".prepi", ".cif"):
            continue
        path = os.path.join(refdir, fname)
        ref = _read_residue_template_atoms(path)
        sel = built.resname == name
        if not sel.any():
            bad.append(f"{name}: not present in built mol")
            continue
        # Per-atom atomtype / charge check.
        for atom_name, (ref_type, ref_charge) in ref.items():
            mask = sel & (built.name == atom_name)
            if not mask.any():
                bad.append(f"{name}.{atom_name}: not found in built mol")
                continue
            for idx in np.where(mask)[0]:
                bt = str(built.atomtype[idx])
                bq = float(built.charge[idx])
                if bt != ref_type:
                    bad.append(
                        f"{name}.{atom_name}: atomtype {bt!r} != ref {ref_type!r}"
                    )
                if abs(bq - ref_charge) > charge_tol:
                    bad.append(
                        f"{name}.{atom_name}: charge {bq:+.4f} != ref "
                        f"{ref_charge:+.4f} (|Δ|>{charge_tol})"
                    )
        # Improper check (prepi only): every prepi-listed
        # intra-residue improper must appear in the built mol. The
        # converse is not required - tLeap adds extra impropers from
        # the ff14SB / GAFF frcmods (notably around terminal NH3+ Ns
        # and sp2 ring atoms) that the prepi doesn't list.
        if ext == ".prepi":
            ref_imp = _read_prepi_intra_residue_impropers(path)
            built_imp = _intra_residue_impropers(built, name) or set()
            for fs in sorted(ref_imp - built_imp, key=sorted):
                bad.append(
                    f"{name}: improper {sorted(fs)} in prepi but not in built"
                )
    assert not bad, "built mol disagrees with references:\n  " + "\n  ".join(bad)


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build verification needs antechamber + teLeap",
)
def test_amber_build_8qfz_matches_reference(tmp_path):
    """8QFZ_B_bicycle: run the full detect -> template -> systemPrepare ->
    parameterizeFromSpecs (default ``cluster`` mode + ``pin=True``) ->
    amber.build pipeline, then check that every atom of every spec'd
    residue in the built Molecule (XX1/XX2/XX3 prepi + LFI cif)
    carries the same ``atomtype`` and ``charge`` as our committed
    reference files. Locks the end-to-end charge / type round-trip
    through tLeap."""
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build

    mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, "8QFZ_B_bicycle.cif"))
    mol.chain[:] = "A"
    mol.segid[:] = "A"
    specs = detectNonStandardResidues(mol)
    mol.templateResidueFromSmiles('resname "LFI"', LFI_SMILES, addHs=True, _logger=False)
    pmol, _ = systemPrepare(mol, detect_specs=specs)

    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"),
        charge_method="gasteiger",
    )
    # Disable auto-capping on the protein segment: the N-terminal CYS
    # is covalently bonded to LFI and the C-terminal CYS already
    # carries an OXT from systemPrepare, so an ACE/NME cap would clash.
    protein_seg = sorted(set(pmol.segid[pmol.atomselect("protein")]))[0]
    built = amber_build(
        pmol,
        outdir=str(tmp_path / "build"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        caps={protein_seg: ("none", "none")},
    )

    _assert_built_matches_references(
        built,
        os.path.join(_CUSTOM_PARAM_DIR, "reference", "8QFZ_B_bicycle"),
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build verification needs antechamber + teLeap",
)
def test_amber_build_5vbl_matches_reference(tmp_path):
    """5VBL_A: same as the 8QFZ end-to-end build check but on the
    chain-A peptide inhibitor with 5 chain-resident NCAAs (HRG, ALC,
    OIC, NLE, 200) + GLU(CD)-LYS(NZ) isopeptide (XX1/XX2 renames).
    Exercises the prolinoid backbone pin (OIC) and the
    canonical-anchor rename path end to end."""
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build

    mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, "5VBL_A.cif"))
    specs = detectNonStandardResidues(mol)
    smiles = {
        "200": "c1cc(ccc1C[C@@H](C(=O)O)N)Cl",
        "ALC": "C1CCC(CC1)C[C@@H](C=O)N",
        "HRG": "C(CCNC(=N)N)C[C@@H](C=O)N",
        "NLE": "CCCC[C@@H](C=O)N",
        "OIC": "C1CC[C@H]2[C@@H](C1)C[C@H](N2)C=O",
    }
    for resname, smi in smiles.items():
        mol.templateResidueFromSmiles(f'resname "{resname}"', smi, addHs=True, _logger=False)
    pmol = systemPrepare(
        mol, detect_specs=specs, restore_missing_sidechains=True,
    )[0]

    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"),
        charge_method="gasteiger",
    )
    # The C-terminal NCAA (residue 200) carries its own OXT via the
    # emitted prepi template; disable amber.build auto-capping so it
    # does not stack an NME on top.
    protein_seg = sorted(set(pmol.segid[pmol.atomselect("protein")]))[0]
    built = amber_build(
        pmol,
        outdir=str(tmp_path / "build"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        caps={protein_seg: ("none", "none")},
    )

    _assert_built_matches_references(
        built,
        os.path.join(_CUSTOM_PARAM_DIR, "reference", "5VBL_A_cluster"),
    )


# Disabled until the parameterization pipeline can handle transition-metal
# ligands. RDKit's PEOE Gasteiger has no parameters for Fe (or any 3d
# metal), so a HEM-bearing residue runs the whole molecule to NaN partial
# charges and the gasteiger pathway aborts. Re-enable (remove the skip and
# run with HTMD_REGEN_REFERENCES=1 to populate the
# reference/1U5U_A_cluster directory) once the pipeline either (a) detaches
# the metal for the PEOE step and restores it with its formal charge, or
# (b) plugs in a metal-aware charge method.
@pytest.mark.skip(
    reason="HEM contains Fe; RDKit Gasteiger PEOE has no parameters for "
    "transition metals. Re-enable once parameterizeFromSpecs supports "
    "metal-containing ligands."
)
@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build verification needs antechamber + teLeap",
)
def test_amber_build_1u5u_hem_matches_reference(tmp_path):
    """1U5U_A: human catalase chain A with a heme-b (HEM) cofactor.
    The catalytic Fe is in the resting Fe(III) state, axially coordinated
    on the proximal side by a deprotonated tyrosinate (TYR-O-, -1 on the
    sidechain OH) - the canonical catalase proximal ligand. The Fe(III)
    charge in the HEM SMILES is set to +3 to match this resting state.
    Exercises a metal-bearing standalone ligand through the
    detect -> template -> systemPrepare -> parameterizeFromSpecs ->
    amber.build pipeline.

    HEM SMILES carries explicit formal charges: -1 on each of the two
    propionate oxygens, -1 on two of the four porphyrin nitrogens,
    +3 on Fe (square-planar), and dative ``->``/``<-`` bonds on the
    Fe-N coordination so RDKit perceives the porphyrin as a single
    macrocycle without inferring spurious aromatic bonds across Fe. Net
    charge of HEM as templated: -1. The active-site total (HEM -1 plus
    proximal TYR-O- -1) is -2.

    Builds against the reference in
    ``test-custom-residue-param/reference/1U5U_A_cluster`` (populated by
    running this test with ``HTMD_REGEN_REFERENCES=1``).
    """
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build

    HEM_SMILES = (
        "C=CC1=C(C)c2cc3c(C)c(CCC(=O)[O-])c4cc5[n]6->[Fe@SP2+3]7"
        "(<-[n]2c1cc1c(C)c(C=C)c(cc6C(C)=C5CCC(=O)[O-])[n-]->71)<-[n-]34"
    )
    mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, "1U5U_A.cif"))
    specs = detectNonStandardResidues(mol)
    mol.templateResidueFromSmiles(
        'resname "HEM"', HEM_SMILES, addHs=True, _logger=False
    )
    pmol = systemPrepare(mol, detect_specs=specs)[0]

    out = parameterizeFromSpecs(
        specs, pmol, outdir=str(tmp_path / "params"),
        charge_method="gasteiger",
    )
    built = amber_build(
        pmol,
        outdir=str(tmp_path / "build"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
    )

    _assert_built_matches_references(
        built,
        os.path.join(_CUSTOM_PARAM_DIR, "reference", "1U5U_A_cluster"),
    )


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="custom-residue parameterization needs antechamber + teLeap",
)
def test_custom_residue_param_reference(tmp_path):
    """Golden test for the non-canonical parameterization pipeline.

    Runs detect -> template -> (systemPrepare) -> parameterizeFromSpecs on
    two fixed inputs and compares every per-residue file it produces (CIF
    or prepi topologies and the frcmod parameter files) against committed
    reference copies, locking down the cluster pipeline and the ff14SB
    backbone-charge pinning end to end. CIF topologies are compared with
    mol_equal; prepi and frcmod files are compared verbatim.

      - 3PTB_BEN.cif: a free, +1 benzamidine ligand, parameterized with
        AM1-BCC. systemPrepare is skipped - it wraps PDB2PQR, which needs
        a biomolecule and rejects a lone ligand, and
        templateResidueFromSmiles has already added hydrogens and the
        formal charge.
      - 8QFZ_B_bicycle.cif: a scaffolded cyclic peptide whose three Cys
        anchors bond an LFI scaffold (one mid-chain, one at each chain
        terminus), exercising the canonical-anchor backbone pin and its
        N- and C-terminal variants.

    The reference files depend on the installed antechamber / AmberTools
    / PDB2PQR versions. After an intentional parameterization change or a
    tool update, regenerate them by running this test with the
    environment variable HTMD_REGEN_REFERENCES=1 set.
    """
    from moleculekit.tools.preparation import systemPrepare

    regenerate = bool(os.environ.get("HTMD_REGEN_REFERENCES"))
    cases = [
        # Exercise both charge methods: AM1-BCC for BEN, Gasteiger
        # (RDKit-computed, so it also honours the net charge) for 8QFZ.
        ("3PTB_BEN.cif", {"BEN": BEN_SMILES}, False, "am1-bcc"),
        ("8QFZ_B_bicycle.cif", {"LFI": LFI_SMILES}, True, "gasteiger"),
    ]
    for cif_name, smiles, run_prepare, charge_method in cases:
        stem = os.path.splitext(cif_name)[0]
        mol = Molecule(os.path.join(_CUSTOM_PARAM_DIR, cif_name))
        # systemPrepare and the bonded bookkeeping need a chain/segid; a
        # bare-ligand CIF may carry neither.
        if not np.any(mol.chain != ""):
            mol.chain[:] = "A"
            mol.segid[:] = "A"
        specs = detectNonStandardResidues(mol)
        for resname, smi in smiles.items():
            mol.templateResidueFromSmiles(f'resname "{resname}"', smi, addHs=True)
        pmol = systemPrepare(mol, detect_specs=specs)[0] if run_prepare else mol
        out = parameterizeFromSpecs(
            specs, pmol, outdir=str(tmp_path / stem), charge_method=charge_method
        )
        _check_against_reference(
            out.topo_paths + out.frcmod_paths,
            os.path.join(_CUSTOM_PARAM_DIR, "reference", stem),
            stem,
            regenerate,
        )
