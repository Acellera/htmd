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


QU4_A_CIF = os.path.join(DATA_DIR, "8QU4_A.cif")
R1J_GLYCO_CIF = os.path.join(DATA_DIR, "1R1J_glyco.cif")


@pytest.mark.skipif(
    not (_antechamber and _tleap),
    reason="end-to-end build needs antechamber + teLeap",
)
def _test_8qu4_stapled_peptide_end_to_end(tmp_path):
    """Pattern-B stapled peptide end-to-end: PDB 8QU4 chain A is a 13-mer
    NF-Y-derived peptide stapled by olefin metathesis between NLE272 and
    MK8276 (the NLE.CE - MK8.CE single bond closes the ring after RCM).

    Pipeline:
      1. SMILES-template both NCAAs to set bond orders.
      2. systemPrepare with residue_smiles to add AMBER-named hydrogens
         (PDB2PQR uses canonical H names; RDKit's addHs alone produces
         generic H1/H2/... that parameterizeNonCanonicalResidues rejects).
      3. parameterizeNonCanonicalResidues on each NCAA CIF that
         systemPrepare wrote.
      4. detectNonStandardResidues on the prepared molecule -> 1
         CrosslinkSpec.
      5. amber.build with the staple custombond + the antechamber prepi/
         frcmod for NLE and MK8.
      6. Verify the NLE.CE - MK8.CE bond is in the built Molecule's bonds.
    """
    from glob import glob
    from moleculekit.tools.preparation import systemPrepare
    from htmd.builder.noncanonical import parameterizeNonCanonicalResidues
    from htmd.builder.amber import build as amber_build
    from moleculekit.tools.nonstandard_residues import (
        detectNonStandardResidues,
        custombondsFromSpecs,
        CrosslinkSpec,
    )

    mol = Molecule(QU4_A_CIF)
    mol.segid[:] = "P"

    residue_smiles = {
        "NLE": "CCCC[C@@H](C(=O)O)N",
        "MK8": "CCCC[C@](C)(C(=O)O)N",
    }
    for resname, smi in residue_smiles.items():
        mol.templateResidueFromSmiles(f"resname {resname}", smi, addHs=True)

    # Detect on the templated molecule (PDB-derived bonds intact). systemPrepare
    # downstream rebuilds the bond list and would force the detector into
    # distance-based bond guessing, which over-bonds atoms that are close in
    # the stapled region.
    specs = detectNonStandardResidues(mol, write_models=False, include_known=True)
    cb = custombondsFromSpecs(specs)
    assert len(cb) == 1
    assert any(isinstance(s, CrosslinkSpec) for s in specs)

    # systemPrepare populates AMBER-conventional H names (HA, HB2, HB3, H,
    # ...) and writes one CIF per non-canonical residue under outdir.
    prepdir = str(tmp_path / "prepared")
    pmol = systemPrepare(
        mol, residue_smiles=residue_smiles, outdir=prepdir, _molkit_ff=False
    )

    cifs = glob(os.path.join(prepdir, "*.cif"))
    assert len(cifs) >= 2, cifs

    # Parameterize each non-canonical residue with antechamber + parmchk2.
    paramdir = str(tmp_path / "params")
    for cif in cifs:
        parameterizeNonCanonicalResidues(
            cif, paramdir, forcefield="GAFF2", charge_model="Gasteiger"
        )
    topos = glob(os.path.join(paramdir, "*.prepi"))
    params = glob(os.path.join(paramdir, "*.frcmod"))
    assert topos and params

    # Build (no solvent, no ions). The peptide is ACE/NH2-capped already so
    # disable amber's automatic terminal capping.
    builddir = str(tmp_path / "build")
    built = amber_build(
        pmol,
        outdir=builddir,
        ionize=False,
        custombonds=cb,
        topo=topos,
        param=params,
        caps={"P": ("none", "none")},
    )
    assert built is not None

    # The staple bond NLE.CE <-> MK8.CE must survive into the built mol.
    nle_ce = set(
        np.where((built.resname == "NLE") & (built.name == "CE"))[0].tolist()
    )
    mk8_ce = set(
        np.where((built.resname == "MK8") & (built.name == "CE"))[0].tolist()
    )
    assert nle_ce and mk8_ce
    n_staple = sum(
        1
        for a, b in built.bonds
        if (int(a) in nle_ce and int(b) in mk8_ce)
        or (int(b) in nle_ce and int(a) in mk8_ce)
    )
    assert n_staple == 1, (
        f"expected one NLE.CE <-> MK8.CE bond in the built prmtop, found {n_staple}"
    )


def _pdb_nag_to_glycam_0yb(mol):
    """Convert a PDB-named NAG residue to its GLYCAM 0YB equivalent in
    place. Renames the residue and the three acetyl-group atom names that
    differ between PDB and GLYCAM (C7 -> C2N, C8 -> CME, O7 -> O2N)."""
    nag_mask = mol.resname == "NAG"
    if not nag_mask.any():
        return
    mol.resname[nag_mask] = "0YB"
    pdb_to_glycam = {"C7": "C2N", "C8": "CME", "O7": "O2N"}
    for pdb_name, glycam_name in pdb_to_glycam.items():
        sel = nag_mask & (mol.name == pdb_name)
        if sel.any():
            mol.name[sel] = glycam_name


@pytest.mark.skipif(
    not _tleap, reason="end-to-end build needs teLeap"
)
def _test_1r1j_glycoprotein_end_to_end(tmp_path):
    """Glycoprotein covalent ligand end-to-end: a slim slice of PDB 1R1J
    around the ASN144 N-glycosylation site (5 protein residues + 1 NAG).

    Pipeline:
      1. Load the slim slice; set segids (P for protein, G for the glycan).
      2. detectNonStandardResidues -> CovalentLigandSpec(ASN.ND2 -> NAG.C1).
      3. forceProtonationFromSpecs -> ASN -> NLN via the new ANCHOR_VARIANTS
         entry; custombondsFromSpecs -> the glycosidic ND2-C1 bond.
      4. systemPrepare with force_protonation: renames ASN144 to NLN and
         adds AMBER-named Hs (NLN's template has no HD22, so the displaced
         amide hydrogen is correctly absent).
      5. Convert PDB NAG to GLYCAM 0YB (residue + 3 atom-name renames).
      6. amber.build with leaprc.GLYCAM_06j-1 added to the FF list.
      7. Verify NLN.ND2 - 0YB.C1 bond is in the built Molecule.

    No antechamber is needed: GLYCAM provides the 0YB sugar parameters and
    the NLN amino-acid template (Asn with HD22 removed and the C1-ND2
    junction parameters)."""
    from htmd.builder.amber import build as amber_build, defaultFf
    from moleculekit.tools.preparation import systemPrepare
    from moleculekit.tools.nonstandard_residues import (
        detectNonStandardResidues,
        custombondsFromSpecs,
        forceProtonationFromSpecs,
        CovalentLigandSpec,
    )

    mol = Molecule(R1J_GLYCO_CIF)
    mol.segid[mol.resname == "NAG"] = "G"  # glycan segment
    mol.segid[mol.segid != "G"] = "P"  # protein segment

    specs = detectNonStandardResidues(mol, write_models=False, include_known=True)
    nag_specs = [s for s in specs if isinstance(s, CovalentLigandSpec)]
    assert len(nag_specs) == 1
    cb = custombondsFromSpecs(specs)
    assert len(cb) == 1
    fp = forceProtonationFromSpecs(specs)
    assert any(variant == "NLN" for _, variant in fp), fp

    # systemPrepare adds AMBER-named Hs and applies the ASN -> NLN rename
    # via force_protonation; PDB2PQR places Hs from the NLN template, so
    # HD22 is absent in the output.
    mol = systemPrepare(mol, force_protonation=fp)

    # Convert PDB NAG naming to GLYCAM 0YB. The custombond's atomselect
    # strings target atoms by segid+chain+resid+name (no resname filter),
    # so the rename doesn't break the bond resolution.
    _pdb_nag_to_glycam_0yb(mol)

    builddir = str(tmp_path / "build")
    built = amber_build(
        mol,
        outdir=builddir,
        ionize=False,
        custombonds=cb,
        ff=defaultFf() + ["leaprc.GLYCAM_06j-1"],
    )
    assert built is not None

    # The glycosidic NLN.ND2 <-> 0YB.C1 bond must be in the built mol.
    nln_nd2 = set(
        np.where((built.resname == "NLN") & (built.name == "ND2"))[0].tolist()
    )
    yb_c1 = set(
        np.where((built.resname == "0YB") & (built.name == "C1"))[0].tolist()
    )
    assert nln_nd2 and yb_c1
    n_glyco = sum(
        1
        for a, b in built.bonds
        if (int(a) in nln_nd2 and int(b) in yb_c1)
        or (int(b) in nln_nd2 and int(a) in yb_c1)
    )
    assert n_glyco == 1, (
        f"expected one NLN.ND2 <-> 0YB.C1 bond in the prmtop, found {n_glyco}"
    )


def _test_amber_handles_stapled_peptide_custombond(tmp_path):
    """Pattern-B stapled peptide: PDB 8QU4 chain A has an i, i+4
    NLE.CE - MK8.CE hydrocarbon staple. detectNonStandardResidues
    returns a CrosslinkSpec; custombondsFromSpecs emits the staple
    bond; amber._prepareMolecule processes it without applying any anchor
    rename (CE is not in ANCHOR_VARIANTS).

    The bond emit itself is identical to the scaffolded-peptide path that
    the 8QFZ end-to-end test covers, so we verify only the
    moleculekit/htmd integration here, not the full prmtop."""
    from htmd.builder.amber import _prepareMolecule
    from moleculekit.tools.nonstandard_residues import (
        detectNonStandardResidues,
        custombondsFromSpecs,
        CrosslinkSpec,
    )

    mol = Molecule(QU4_A_CIF)
    mol.segid[:] = "P"

    specs = detectNonStandardResidues(mol, write_models=False, include_known=True)
    crosslinks = [s for s in specs if isinstance(s, CrosslinkSpec)]
    assert len(crosslinks) == 1
    cb = custombondsFromSpecs(specs)
    assert len(cb) == 1

    _prepareMolecule(mol, caps={}, disulfide=None, custombonds=cb, remove=None)

    # CE is not in ANCHOR_VARIANTS so neither stapled NCAA was renamed.
    # (Plain LYS also has a CE atom; restrict the check to the NLE/MK8
    # residues only.)
    stapled_resnames = set(
        mol.resname[
            (mol.name == "CE") & np.isin(mol.resname, ["NLE", "MK8"])
        ].tolist()
    )
    assert stapled_resnames == {"NLE", "MK8"}
    # No CYX rename since neither endpoint is a Cys SG.
    assert "CYX" not in set(mol.resname.tolist())


def _test_amber_handles_glycoprotein_asn_to_nln_rename(tmp_path):
    """Verify the new ASN -> NLN entry in ANCHOR_VARIANTS fires: a slim
    1R1J slice has Asn144.ND2 covalently bonded to NAG752.C1
    (N-glycosylation). detectNonStandardResidues returns a
    CovalentLigandSpec; forceProtonationFromSpecs emits an Asn -> NLN
    rename; amber._prepareMolecule's defense-in-depth path applies the
    rename and drops HD22 (the displaced amide hydrogen).

    Note: a full tleap build of glycoproteins additionally needs GLYCAM
    loaded with the right NAG-residue naming (0YB / 4YB / ...). Here we
    only verify the moleculekit-side spec plumbing and the htmd-side
    rename machinery - the bond emit itself is identical to the
    scaffolded-peptide path that the 8QFZ end-to-end test covers."""
    from htmd.builder.amber import _prepareMolecule
    from moleculekit.molecule import UniqueAtomID
    from moleculekit.tools.preparation import systemPrepare
    from moleculekit.tools.nonstandard_residues import (
        detectNonStandardResidues,
        custombondsFromSpecs,
        forceProtonationFromSpecs,
        CovalentLigandSpec,
    )

    mol = Molecule(R1J_GLYCO_CIF)
    mol.segid[:] = "P"

    # Detect on the raw mol (PDB-derived bonds intact). systemPrepare will
    # rebuild the bond list and would force the detector into distance-based
    # guessing.
    specs = detectNonStandardResidues(mol, write_models=False, include_known=True)
    nag_specs = [s for s in specs if isinstance(s, CovalentLigandSpec)]
    assert len(nag_specs) == 1
    cb = custombondsFromSpecs(specs)
    assert len(cb) == 1
    fp = forceProtonationFromSpecs(specs)
    assert any(variant == "NLN" for _, variant in fp), fp

    # systemPrepare adds AMBER-named Hs (HD21 / HD22 on Asn144's amide).
    # We deliberately do *not* pass force_protonation here, so the residue
    # stays ASN with HD22 intact and lets us verify amber._prepareMolecule's
    # defense-in-depth path applies the ASN -> NLN rename and drops HD22.
    # hold_nonpeptidic_bonds=False so PDB2PQR doesn't freeze Asn144 (which
    # would skip H-addition on the anchored residue).
    mol = systemPrepare(mol, hold_nonpeptidic_bonds=False)
    # Sanity: HD22 is present before the defense-in-depth rename runs.
    assert mol.atomselect("resid 144 and name HD22").sum() == 1

    # Track ND2 by UniqueAtomID so we can find it after _prepareMolecule
    # renumbers residues.
    nd2_uid = UniqueAtomID.fromMolecule(mol, "resid 144 and name ND2")

    _prepareMolecule(mol, caps={}, disulfide=None, custombonds=cb, remove=None)

    # The Asn144 residue should now be NLN, and HD22 should be gone.
    sel = (
        (mol.chain == nd2_uid.chain)
        & (mol.segid == nd2_uid.segid)
        & (mol.insertion == nd2_uid.insertion)
        & (mol.name == "ND2")
    )
    nd2_after = np.where(sel)[0]
    assert len(nd2_after) >= 1
    nln_nd2 = nd2_after[mol.resname[nd2_after] == "NLN"]
    assert len(nln_nd2) == 1, (
        "expected Asn144 to be renamed to NLN, "
        f"found resnames {set(mol.resname[nd2_after])}"
    )

    # HD22 must be absent in the renamed residue.
    resid_val = mol.resid[nln_nd2[0]]
    res_mask = (
        (mol.resid == resid_val)
        & (mol.chain == nd2_uid.chain)
        & (mol.segid == nd2_uid.segid)
    )
    assert ((mol.name == "HD22") & res_mask).sum() == 0
