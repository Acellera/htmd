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

from moleculekit.molecule import Molecule, UniqueAtomID
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
)

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
    out = parameterizeFromSpecs(specs, pmol, outdir=str(tmp_path / "params"))
    return amber_build(
        pmol,
        outdir=str(tmp_path / "build"),
        ionize=False,
        custombonds=out.custombonds,
        topo=out.topo_paths,
        param=out.frcmod_paths,
        **(build_kwargs or {}),
    )


def _test_disulfide_path_unchanged(tmp_path):
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

    disulfide, _, _ = _prepareMolecule(
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


def _test_amber_anchor_rename_via_custombonds(tmp_path):
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
def _test_8qu4_stapled_peptide_end_to_end(tmp_path):
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
        specs, pmol, outdir=str(tmp_path / "params")
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
def _test_8qfz_scaffolded_peptide_end_to_end(tmp_path):
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
        specs, pmol, outdir=str(tmp_path / "params")
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
def _test_full_pipeline_5vbl(tmp_path):
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
        "OLC": "CCCCCCCC(O)OC[C@H](O)CO",
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


try:
    import openmm  # noqa: F401

    _openmm = True
except ImportError:
    _openmm = False


@pytest.mark.skipif(
    not (_antechamber and _tleap and _openmm),
    reason="OpenMM build comparison needs antechamber + teLeap + openmm",
)
def _test_full_pipeline_5vbl_openmm_vs_amber(tmp_path):
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
    from htmd.builder.openff import build as openff_build
    from collections import Counter

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
        "OLC": "CCCCCCCC(O)OC[C@H](O)CO",
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
    out = parameterizeFromSpecs(specs, pmol, outdir=str(tmp_path / "params"))

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
        extra_xml=[out.xml_path],
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
def _test_full_pipeline_6a5j_openmm_vs_amber(tmp_path):
    """6A5J: small canonical peptide (no NCAAs, no ligands). Builds via
    amber.build (tLeap) and openff.build (OpenMM ForceField XML) and
    asserts the two systems carry identical force-field data via
    :func:`_assert_builds_equivalent`. Acts as a control for the NCAA
    builds - if 6A5J disagrees, the divergence is in the basic builder
    path rather than the parameterizeFromSpecs cluster pipeline."""
    from moleculekit.tools.autosegment import autoSegment
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.amber import build as amber_build
    from htmd.builder.openff import build as openff_build

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
def _test_parameterize_from_specs_emits_openmm_xml(tmp_path):
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
        "OLC": "CCCCCCCC(O)OC[C@H](O)CO",
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
    out = parameterizeFromSpecs(specs, pmol, outdir=str(tmp_path / "params"))

    assert out.xml_path is not None
    assert os.path.isfile(out.xml_path)
    assert os.path.getsize(out.xml_path) > 0

    ff = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/tip3p.xml",
        out.xml_path,
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
def _test_full_pipeline_4tot_e(tmp_path):
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
def _test_full_pipeline_1r1j(tmp_path):
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
def _test_parameterize_from_specs_dedup_4tot(tmp_path):
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
        # (RCSB chem-comp P6G). SO4 has no carbon, so the protonation guard
        # in parameterizeFromSpecs skips it.
        "P6G": "OCCOCCOCCOCCOCCOCCO",
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
    out = parameterizeFromSpecs(specs, pmol, outdir=str(tmp_path / "params"))

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


def _test_check_specs_protonated_passes_aliphatic():
    mol, groups = _toy_mol([_HEXANE])
    _check_specs_protonated(mol, {0: None}, groups)  # must not raise


def _test_check_specs_protonated_passes_aromatic():
    # The estimator treats the ring as all single bonds, so benzene's 6 H sit
    # well above the 0.2 * 12 threshold.
    mol, groups = _toy_mol([_BENZENE])
    _check_specs_protonated(mol, {0: None}, groups)  # must not raise


def _test_check_specs_protonated_zero_hydrogens_raises():
    mol, groups = _toy_mol([_BARE_C4, _HEXANE])
    with pytest.raises(RuntimeError, match=r"under-protonated"):
        _check_specs_protonated(mol, {0: None}, groups)


def _test_check_specs_protonated_stripped_sidechain_raises():
    # Backbone H present, every sidechain H missing -> still flagged.
    mol, groups = _toy_mol([_LYS_STRIPPED])
    with pytest.raises(RuntimeError) as exc:
        _check_specs_protonated(mol, {0: None}, groups)
    assert "LYX" in str(exc.value)
    assert "2 hydrogen(s) present" in str(exc.value)


def _test_check_specs_protonated_ignores_carbonless_residues():
    # A metal ion and an oxo-anion carry no expected H even with 0 H present.
    mol, groups = _toy_mol([_HEXANE, _ZINC, _SULFATE])
    _check_specs_protonated(mol, {1: None, 2: None}, groups)  # must not raise


def _test_check_specs_protonated_ignores_tiny_carbon_species():
    # CO2 estimates only ~4 expected H -> below _MIN_EXPECTED_HYDROGENS, so
    # being H-free does not trip the check.
    mol, groups = _toy_mol([_HEXANE, _SMALL_CO2])
    _check_specs_protonated(mol, {1: None}, groups)  # must not raise


def _test_check_specs_protonated_reports_every_bad_residue():
    mol, groups = _toy_mol([_BARE_C4, _HEXANE, _LYS_STRIPPED])
    with pytest.raises(RuntimeError) as exc:
        _check_specs_protonated(mol, {0: None, 2: None}, groups)
    msg = str(exc.value)
    assert "BC4" in msg and "LYX" in msg
    # The fully protonated residue (index 1) is not in the spec set and must
    # not appear regardless.
    assert "HEX" not in msg


def _test_check_specs_protonated_no_bonds_is_noop():
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


def _test_clean_frcmod_free_residue(tmp_path):
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


def _test_clean_frcmod_backbone_terms_dropped(tmp_path):
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


def _test_clean_frcmod_cap_atoms_excluded(tmp_path):
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
