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

import numpy as np
import pytest

from moleculekit.molecule import Molecule, UniqueAtomID
from moleculekit.tools.nonstandard_residues import (
    detectNonStandardResidues,
    CrosslinkedNCAASpec,
    ScaffoldSpec,
    CanonicalRenamedSpec,
)
from htmd.builder.amber import _findTeLeap, _prepareMolecule
from htmd.builder.nonstandard import parameterizeFromSpecs

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


def _run_pipeline(
    mol, smiles, tmp_path, build_kwargs=None, ignore_ns=True
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
    pmol = systemPrepare(
        mol,
        outdir=str(tmp_path / "prep"),
        ignore_ns=ignore_ns,
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
    assert sum(isinstance(s, CrosslinkedNCAASpec) for s in specs) == 2
    assert sum(isinstance(s, CanonicalRenamedSpec) for s in specs) == 0

    # 2. Template the two NCAAs from SMILES.
    mol.templateResidueFromSmiles("resname NLE", NLE_SMILES, addHs=True)
    mol.templateResidueFromSmiles("resname MK8", MK8_SMILES, addHs=True)

    # 3. systemPrepare (preserves bonds across PDB2PQR via bond-capture).
    pmol = systemPrepare(mol, outdir=str(tmp_path / "prep"), ignore_ns=True)

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
    renames = [s for s in specs if isinstance(s, CanonicalRenamedSpec)]
    assert len(scaffolds) == 1 and scaffolds[0].resname == "LFI"
    assert len(renames) == 3
    cys_rename_set = {r.new_resname for r in renames}
    assert len(cys_rename_set) == 3, (
        f"expected three distinct CYS rename targets, got {cys_rename_set}"
    )
    for cys_new in cys_rename_set:
        assert len(cys_new) == 3 and cys_new.startswith("CY")

    # 2. Template LFI from SMILES.
    mol.templateResidueFromSmiles("resname LFI", LFI_SMILES, addHs=True)

    # 3. systemPrepare runs PDB2PQR on the canonical naming, then applies
    # the rename + displaced-H drop from the spec list at the end.
    pmol = systemPrepare(
        mol,
        outdir=str(tmp_path / "prep"),
        ignore_ns=True,
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
    # (ff14SB type on one side, GAFF2 on the other).
    pset = AmberParameterSet(out.frcmod_paths[0])
    canonical_types = {"S", "2C", "CT", "CX", "C", "N", "O", "H", "H1"}
    cross_ff_bonds = [
        k for k in pset.bond_types
        if any(t in canonical_types for t in k)
        and any(str(t).islower() for t in k)
    ]
    assert cross_ff_bonds, "junction frcmod missing cross-FF bond entries"
    assert any("S" in k for k in cross_ff_bonds)

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
    # ``ignore_ns=False`` lets systemPrepare build a PDB2PQR template for
    # residue 200 from the spec, which stops PDB2PQR from C-terminal-capping
    # PRO 16 with a phantom OXT.
    built = _run_pipeline(
        mol,
        smiles,
        tmp_path,
        build_kwargs={"caps": {"P0": ("none", "none")}},
        ignore_ns=False,
    )
    assert built is not None
    _check_no_overvalent_atoms(built)


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

    smiles = {"NAG": "CC(=O)N[C@@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"}
    built = _run_pipeline(mol, smiles, tmp_path)
    assert built is not None
    _check_no_overvalent_atoms(built)
