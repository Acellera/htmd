import os
import shutil

import numpy as np
import pytest

from moleculekit.molecule import Molecule, UniqueAtomID
from moleculekit.tools.nonstandard_residues import (
    detectNonStandardResidues,
    custombondsFromSpecs,
    forceProtonationFromSpecs,
)
from htmd.builder.amber import _findTeLeap, _prepareMolecule
from htmd.builder.scaffolded_peptide import prepareScaffoldedResidue

curr_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(curr_dir, "test_scaffolded_peptide_builder")
QFZ_B_CIF = os.path.join(DATA_DIR, "8QFZ_B.cif")

# Pre-reaction LFI warhead: 1,3,5-triazinane core with three bromoacetyl arms.
# After reaction with the three Cys SG atoms, the Br leaving groups are
# displaced and replaced by S-Cys connections.
LFI_SMILES = "C1N(CN(CN1C(=O)CCBr)C(=O)CCBr)C(=O)CCBr"

_antechamber = shutil.which("antechamber") is not None
_tleap = _findTeLeap() is not None


def _load_8qfz_chainB():
    return Molecule(QFZ_B_CIF)


@pytest.mark.skipif(not _antechamber, reason="antechamber not available")
def _test_prepare_scaffolded_residue(tmp_path):
    """End-to-end: detect, antechamber, prepareScaffoldedResidue.
    The function should produce a matched <RESNAME>.cif + <RESNAME>.frcmod
    pair, with the frcmod containing both antechamber's scaffold-internal
    terms and the cross-FF junction terms appended."""
    from htmd.builder.noncanonical import _fftype_antechamber

    mol = _load_8qfz_chainB()
    specs = detectNonStandardResidues(
        mol, outdir=str(tmp_path), write_models=True, include_known=True
    )
    assert len(specs) == 1
    spec = specs[0]
    resname = spec.resname

    ac_dir = tmp_path / "ac"
    ac_dir.mkdir()
    _, typed_path, frcmod_path = _fftype_antechamber(
        Molecule(spec.model_compound_cif),
        tmpdir=str(ac_dir),
        method="gaff2",
        netcharge=0,
    )

    out = tmp_path / "params"
    cif_path, frcmod_out = prepareScaffoldedResidue(
        typed_path, frcmod_path, spec, outdir=str(out)
    )

    # Names match the resname so amber.build's resname-keyed atom-type
    # collision dedup pairs them.
    assert os.path.basename(cif_path) == f"{resname}.cif"
    assert os.path.basename(frcmod_out) == f"{resname}.frcmod"
    assert os.path.isfile(cif_path) and os.path.isfile(frcmod_out)

    # Scaffold-only CIF: stubs gone, scaffold atoms preserved, resname intact.
    scaffold_mol = Molecule(cif_path)
    n_scaffold = sum(
        1 for v in spec.model_atom_map.values() if v.role == "scaffold"
    )
    assert scaffold_mol.numAtoms == n_scaffold
    assert all(rn == resname for rn in scaffold_mol.resname)

    # Combined frcmod: parmed should round-trip-load and find at least one
    # junction BOND term carrying a canonical-FF (uppercase) atom-type name.
    from parmed.amber import AmberParameterSet

    pset = AmberParameterSet(frcmod_out)
    canonical_types = {"S", "N3", "OH", "NA", "CT", "2C", "3C", "CA", "CC", "H1"}
    has_junction_bond = any(
        any(t in canonical_types for t in key) for key in pset.bond_types
    )
    assert has_junction_bond, "junction BOND with a canonical-FF type missing"
    # The original antechamber scaffold-internal terms must still be present.
    assert ("c3", "c3") in pset.bond_types or ("c", "n") in pset.bond_types


@pytest.mark.skipif(not _antechamber, reason="antechamber not available")
def _test_prepare_scaffolded_residue_default_outdir(tmp_path):
    """When outdir=None, prepareScaffoldedResidue should fall back to
    a tempdir and still return matched <RESNAME>.cif + <RESNAME>.frcmod
    paths."""
    from htmd.builder.noncanonical import _fftype_antechamber

    mol = _load_8qfz_chainB()
    specs = detectNonStandardResidues(
        mol, outdir=str(tmp_path), write_models=True, include_known=True
    )
    spec = specs[0]
    ac_dir = tmp_path / "ac"
    ac_dir.mkdir()
    _, typed_path, frcmod_path = _fftype_antechamber(
        Molecule(spec.model_compound_cif),
        tmpdir=str(ac_dir),
        method="gaff2",
        netcharge=0,
    )

    cif_path, frcmod_out = prepareScaffoldedResidue(
        typed_path, frcmod_path, spec
    )
    assert os.path.basename(cif_path) == f"{spec.resname}.cif"
    assert os.path.basename(frcmod_out) == f"{spec.resname}.frcmod"
    assert os.path.isfile(cif_path) and os.path.isfile(frcmod_out)
    # Outputs land in the same auto-tempdir.
    assert os.path.dirname(cif_path) == os.path.dirname(frcmod_out)


def _test_amber_anchor_rename_via_custombonds(tmp_path):
    """Defense-in-depth: amber._prepareMolecule should auto-rename anchored
    Cys to CYX and drop HG when the user passes a custombond from CYS SG to
    a scaffold atom, even without going through systemPrepare first."""
    mol = _load_8qfz_chainB()

    # Make sure the Cys residues have an HG to drop. Add a dummy HG to one of
    # the anchored cysteines so we can verify the H-drop step.
    sg_idx = int(mol.atomselect("resid 11 and name SG", indexes=True)[0])
    sg_pos = mol.coords[sg_idx, :, mol.frame]
    hg_pos = sg_pos + np.array([0.0, 0.0, 1.34])
    new = Molecule().empty(1)
    new.name[:] = "HG"
    new.element[:] = "H"
    new.resname[:] = "CYS"
    new.resid[:] = 11
    new.chain[:] = mol.chain[sg_idx]
    new.segid[:] = mol.segid[sg_idx]
    new.coords = hg_pos.reshape(1, 3, 1).astype(np.float32)
    new.record[:] = "ATOM"
    mol.insert(new, index=sg_idx + 1, collisions=False)

    # Add segid (amber._prepareMolecule requires segids).
    mask = mol.segid == ""
    if mask.any():
        mol.segid[mask] = "P"
    # Dedupe segid uniformly per chain so the molecule passes _missingSegID.
    mol.segid[:] = "P"

    # Build custombonds: each Cys SG -> LFI anchor C.
    specs = detectNonStandardResidues(mol, write_models=False, include_known=True)
    cb = custombondsFromSpecs(specs)
    assert len(cb) == 3

    # Track the SG atoms by unique ID so we can find them after _prepareMolecule
    # renumbers the residues.
    sg_uids_before = [
        UniqueAtomID.fromMolecule(mol, f"resid {r} and name SG") for r in (11, 17, 22)
    ]

    # Run _prepareMolecule with empty caps/disulfide/remove and our custombonds.
    _prepareMolecule(mol, caps={}, disulfide=None, custombonds=cb, remove=None)

    # The SG atoms still exist, but the residue containing each should now
    # be CYX, with no HG hydrogen.
    for uid in sg_uids_before:
        # SG atoms keep all UniqueAtomID fields except resname (which we changed
        # to CYX) and resid (renumbered).  Match on chain/segid/insertion/name.
        sel = (
            (mol.chain == uid.chain)
            & (mol.segid == uid.segid)
            & (mol.insertion == uid.insertion)
            & (mol.name == "SG")
        )
        # Find which residue the SG is in
        sg_idx = np.where(sel)[0]
        assert len(sg_idx) >= 1
        # Among the SG hits, find the one whose resname is CYX (the renamed one).
        cyx_sg = sg_idx[mol.resname[sg_idx] == "CYX"]
        assert len(cyx_sg) >= 1, (uid.resid, mol.resname[sg_idx])
        # In the residue containing this SG there should be no HG.
        resid_val = mol.resid[cyx_sg[0]]
        chain_val = mol.chain[cyx_sg[0]]
        segid_val = mol.segid[cyx_sg[0]]
        res_mask = (
            (mol.resid == resid_val)
            & (mol.chain == chain_val)
            & (mol.segid == segid_val)
        )
        assert ((mol.name == "HG") & res_mask).sum() == 0


def _test_disulfide_path_unchanged(tmp_path):
    """Sanity check: an isolated Cys-Cys disulfide pair still gets renamed
    CYX with HG dropped, i.e. the rename logic is unchanged for the disulfide
    path even after generalisation."""
    rows = [
        # (name, element, resid, x, y, z)
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

    # The two SG atoms are too far apart for distance-based disulfide
    # detection, so we hand it explicitly via the disulfide kw.
    disulfide, _, _ = _prepareMolecule(
        mol,
        caps={},
        disulfide=[["resid 1 and name SG", "resid 2 and name SG"]],
        custombonds=None,
        remove=None,
    )

    assert len(disulfide) == 1
    # Both cysteines are now CYX, both HG are gone.
    assert set(mol.resname[mol.resid == 1].tolist()) == {"CYX"}
    assert set(mol.resname[mol.resid == 2].tolist()) == {"CYX"}
    assert (mol.name == "HG").sum() == 0


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def _test_8qfz_chainB_end_to_end(tmp_path):
    """End-to-end build of 8QFZ chain B (no solvent, no ions). Exercises the
    full workflow: SMILES templating of LFI, detector, antechamber on the
    scaffolded peptide model compound, prepareScaffoldedResidue,
    systemPrepare with force_protonation, amber.build with custombonds.

    Verifies that the resulting prmtop contains a bond between every Cys SG
    and its LFI anchor carbon."""
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.noncanonical import _fftype_antechamber
    from htmd.builder.amber import build as amber_build

    mol = _load_8qfz_chainB()
    mol.filter("not resname HOH", _logger=False)

    # User templating step (sets bond orders, formal charges, adds Hs).
    mol.templateResidueFromSmiles("resname LFI", LFI_SMILES, addHs=True)

    # 1. Detect non-standard residues + write model CIFs.
    params_dir = tmp_path / "params"
    params_dir.mkdir()
    specs = detectNonStandardResidues(mol, outdir=str(params_dir))

    scaffold_specs = [s for s in specs if s.category == "scaffolded_peptide"]
    assert len(scaffold_specs) == 1
    spec = scaffold_specs[0]
    assert spec.resname == "LFI"
    assert len(spec.anchors) == 3

    # 2. Parameterize the model compound and produce a matched <RESNAME>.cif
    #    + <RESNAME>.frcmod pair (junction terms folded into the frcmod).
    ac_dir = tmp_path / f"ac_{spec.resname}"
    ac_dir.mkdir()
    _, typed_path, frcmod_path = _fftype_antechamber(
        Molecule(spec.model_compound_cif),
        tmpdir=str(ac_dir),
        method="gaff2",
        charge_method="bcc",
    )
    cif_path, frcmod_out = prepareScaffoldedResidue(
        typed_path, frcmod_path, spec, outdir=str(params_dir)
    )

    # 3. Derive systemPrepare / amber.build args from the specs.
    fp = forceProtonationFromSpecs(specs)
    cb = custombondsFromSpecs(specs)
    assert len(fp) == 3
    assert len(cb) == 3

    # 4. Prepare protonation states (auto-renames anchored Cys -> CYX).
    prepared = systemPrepare(mol, force_protonation=fp)

    # 5. Build (no solvation, no ionization).
    build_dir = tmp_path / "build"
    built = amber_build(
        prepared,
        outdir=str(build_dir),
        ionize=False,
        custombonds=cb,
        param=[frcmod_out],
        topo=[cif_path],
    )
    assert built is not None
    assert os.path.isfile(os.path.join(str(build_dir), "structure.prmtop"))
    assert os.path.isfile(os.path.join(str(build_dir), "structure.pdb"))

    # 6. Verify the prmtop has the three Cys-SG to LFI-C bonds.
    import parmed

    s = parmed.load_file(os.path.join(str(build_dir), "structure.prmtop"))
    cys_sgs = [a for a in s.atoms if a.name == "SG" and a.residue.name == "CYX"]
    lfi_anchors = [
        a
        for a in s.atoms
        if a.name in ("C10", "C11", "C12") and a.residue.name == "LFI"
    ]
    assert len(cys_sgs) == 3
    assert len(lfi_anchors) == 3
    n_cys_lfi_bonds = 0
    for b in s.bonds:
        a1, a2 = b.atom1, b.atom2
        if (a1 in cys_sgs and a2 in lfi_anchors) or (
            a2 in cys_sgs and a1 in lfi_anchors
        ):
            n_cys_lfi_bonds += 1
    assert (
        n_cys_lfi_bonds == 3
    ), f"expected 3 Cys-SG <-> LFI-C bonds, found {n_cys_lfi_bonds}"
