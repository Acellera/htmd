# System building

System building is the work of turning a structure - a PDB, a homology model, a docked complex - into a simulation-ready system: protonated, segmented, solvated, ionised, parameterized against a force field, and written to a topology + coordinates pair that an MD engine accepts. HTMD treats this as a first-class workflow, not a script you string together by hand, and supports the full range of biomolecular systems through a single API.

## The canonical pipeline

Every HTMD build goes through the same five stages:

```{mermaid}
flowchart LR
    A[Structure input<br/>PDB / CIF / model] --> B[Preparation<br/>protonation, gaps, mutations]
    B --> C[Segmentation<br/>protein / lipid / ligand /<br/>water / ion segments]
    C --> D[Solvation &<br/>ionisation]
    D --> E[Force-field build<br/>AMBER / CHARMM / OpenFF]
    E --> F[Topology + coords<br/>ready for MD]
```

- **Segmentation** ({py:func}`~moleculekit.tools.autosegment.autoSegment`) labels each contiguous chain - and each non-protein residue - in its own segment, so the builder treats the components independently and doesn't extend protein chains through HETATM ligands.
- **Non-standard residue detection** ({py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues`) inspects the molecule and returns one spec per residue the force field doesn't natively know (ligands, NCAAs, modified residues, crosslinkers).
- **SMILES templating** ({py:meth}`~moleculekit.molecule.Molecule.templateResidueFromSmiles`) sets correct bond orders, formal charges, and hydrogen counts on each non-standard residue from its reference SMILES - a prerequisite for clean parameterization.
- **Preparation** ({py:func}`~moleculekit.tools.preparation.systemPrepare`) does pKa-aware protonation, missing-sidechain repair and partial backbone repair, and mutation; passed the `detect_specs=` list it preserves the templated non-standard residues across the PDB2PQR pass.
- **Parameterization** ({py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs`) runs antechamber on each non-standard residue, emits per-residue topology files (`topo_paths`), parameter files (`frcmod_paths`), and a `custombonds` list the builder feeds straight to tLeap.
- **Force-field build** ({py:func}`htmd.builder.amber.build` - the most fully featured backend - or {py:func}`htmd.builder.charmm.build` / {py:func}`htmd.builder.openmm.build`) consumes the prepared molecule plus the parameterization outputs and writes a simulation-ready topology + coordinates pair.

Solvation ({py:func}`~htmd.builder.solvate.solvate`) and ionisation ({py:func}`~htmd.builder.ionize.ionize`) can be called explicitly before the build step, or left to the builder's own `ionize=True` default.

The four [system-building tutorials](../tutorials/system-prep/index.md) walk this pipeline for protein-only, protein-ligand, protein-protein, and protein-in-membrane systems. The same pipeline handles every case below.

## Beyond canonical proteins

HTMD's building stack handles cases that "build a protein in CHARMM" tools usually don't:

### Non-canonical amino acids (NCAAs)

NCAAs - phosphorylated residues, fluorinated analogues, modified backbones, peptide-drug warheads - flow through the standard three-step pipeline:

1. {py:func}`~moleculekit.tools.nonstandard_residues.detectNonStandardResidues` inspects `mol` and returns one spec per non-standard residue, grouping connected non-canonicals (and their canonical anchors) into clusters automatically.
2. {py:meth}`~moleculekit.molecule.Molecule.templateResidueFromSmiles` is called once per unique non-standard residue with its reference SMILES, fixing bond orders and adding hydrogens.
3. {py:func}`~moleculekit.tools.preparation.systemPrepare` is called with `detect_specs=specs` so PDB2PQR protonates the canonical residues while preserving the templated non-standard ones, then {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` produces the GAFF2 topologies and AMBER ff14SB-anchored parameters that `amber.build` consumes.

Default charge model for NCAA parameterization is **AM1-BCC** (standard GAFF-quality charges, slower); pass `charge_method="gasteiger"` to {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` for the much faster Gasteiger model when high-quality electrostatics aren't critical.

### Stapled peptides, isopeptides, and cyclic peptides

Stapled peptides (hydrocarbon, lactam, thioether staples), isopeptide bonds (N-to-side-chain crosslinks, e.g. ubiquitin-conjugation chemistry), and head-to-tail or side-chain cyclic peptides are all crosslinked NCAA constructs. They use the same three-step flow as standalone NCAAs - `detectNonStandardResidues` groups the bridging residue together with its anchor residues (or, for cyclic peptides, the N- and C-terminal anchors) into a single cluster, and the resulting `custombonds` entries carry the crosslink across the chain.

The same {py:func}`~htmd.builder.nonstandard.parameterizeFromSpecs` call handles staples, isopeptides, and cycles; you don't have to hand-edit a tLeap script or a CHARMM PSF.

### Disulfide bonds

Disulfide handling is automatic by default: {py:func}`~htmd.builder.builder.detectDisulfideBonds` scans the prepared structure for sulphur pairs within ~3 Å on cysteine residues and returns the inferred bridges. Both {py:func}`htmd.builder.amber.build` and {py:func}`htmd.builder.charmm.build` accept a `disulfide=` argument; pass `None` to auto-detect, or a list of `(sel1, sel2)` atom-selection pairs to override.

For pathological cases - ambiguous pairings, designed mutants where the disulfide must be at a specific Cys pair - {py:func}`~htmd.builder.builder.convertDisulfide` lets you preview what the auto-detection would produce before building.

### Membrane systems

The membrane stack lives in {py:mod}`htmd.membranebuilder`:

- {py:func}`~htmd.membranebuilder.build_membrane.buildMembrane` assembles a bilayer at a given XY size from one or more lipid types, listed via `listLipids()`.
- A short LJ-fluid relaxation step ({py:mod}`~htmd.membranebuilder.ljfluid`) packs the lipids without ring penetrations ({py:mod}`~htmd.membranebuilder.ringpenetration`).
- For protein-in-membrane systems, the membrane is generated separately, then the protein is embedded via {py:func}`~htmd.builder.builder.embed` after stripping clashing lipids ({py:func}`~htmd.builder.builder.removeLipidsInProtein`).

See the {doc}`Membrane Builder how-to <../how-to/MembraneBuilder>` for the recipe.

## Force-field choice

The same prepared molecule can be built under multiple force fields without redoing preparation:

- **AMBER** ({py:func}`htmd.builder.amber.build`) via tLeap - the most fully featured backend, the one the NCAA flow is parameterized against (ff14SB anchors, GAFF2 atom types), and the recommended default for protein, protein-ligand, and crosslinked-NCAA systems.
- **CHARMM** ({py:func}`htmd.builder.charmm.build`) via VMD/psfgen - for projects standardised on CHARMM force fields and CHARMM-GUI-style membrane setups.
- **OpenFF** ({py:func}`htmd.builder.openmm.build`) for systems where you want OpenFF / SMIRNOFF small-molecule parameters and protein parameters from `openmmforcefields`. The build is driven through OpenMM but the differentiator versus the AMBER / CHARMM backends is the OpenFF small-molecule force field, not the simulation engine - all three builders produce systems that any common engine (ACEMD, OpenMM, GROMACS) can run.

Force-field-specific tooling lives in the same module as the build entry point; common topology surgery (insertions, mutations, cis-peptide detection, disulfide detection) is in {py:mod}`htmd.builder.builder`.

## What to read next

- The four [system-building tutorials](../tutorials/system-prep/index.md) for guided walk-throughs.
- {doc}`Membrane Builder how-to <../how-to/MembraneBuilder>` for the membrane-specific recipe.
- {py:mod}`htmd.builder` and {py:mod}`htmd.membranebuilder` in the [API reference](../reference/index.md) for full signatures.
