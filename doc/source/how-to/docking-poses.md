# How to generate docking poses for downstream MD

## Goal

Use HTMD's {py:func}`~htmd.dock.dock` wrapper around AutoDock Vina to generate candidate poses of a ligand inside (or around) a protein binding pocket, then write each pose out as a starting structure for a separate MD trajectory.

## Minimal example

```python
from moleculekit.molecule import Molecule
from htmd.dock import dock

protein = Molecule("protein.pdb")
ligand  = Molecule("ligand.mol2")            # or .sdf / .pdbqt

poses, scores = dock(protein, ligand,
                     center=protein.getCenter("resname POCKET_RES"),
                     extent=[20, 20, 20],    # XYZ box size in Å
                     numposes=20)

for i, pose in enumerate(poses):
    pose.write(f"./poses/pose_{i:02d}.pdb")
    print(f"pose {i:2d}: affinity {scores[i, 0]:6.2f} kcal/mol, "
          f"RMSD lb {scores[i, 1]:.2f} ub {scores[i, 2]:.2f}")
```

{py:func}`~htmd.dock.dock` calls AutoDock Vina under the hood and returns a **2-tuple**: `poses` is a list of {py:class}`~moleculekit.molecule.Molecule` objects (ligand-only, in their docked positions), and `scores` is an `(N, 3)` numpy array of Vina's per-pose `(affinity kcal/mol, RMSD lower bound, RMSD upper bound)`. Merge each pose with the protein if you need a complex.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `protein`, `ligand` | The receptor and ligand molecules. Hydrogens should already be present and correctly assigned at your target pH (use {py:func}`~moleculekit.tools.preparation.systemPrepare` upstream). |
| `center` | XYZ centre of the search box in Å. Pass the centroid of a known binding pocket residue, or a literal `[x, y, z]`. |
| `extent` | Box edge lengths in Å (`[dx, dy, dz]`). Larger boxes cover more of the protein surface but slow the search and dilute scoring. |
| `numposes` | Maximum number of poses Vina emits. Vina caps this at **20**; values above 20 silently return at most 20 poses. |
| `babelexe` | Path to `obabel` (Open Babel) for format conversion. Default `"obabel"` - must be on `$PATH`. |
| `vinaexe` | Path to the Vina binary. Default is `f"{platform.system()}-vina"` (e.g. `Linux-vina`, `Darwin-vina`), resolved via `shutil.which` - the platform-prefixed binary must be on `$PATH`. |

## Common variations

### Build an MD system from each pose

```python
from htmd.builder import amber
from htmd.builder.solvate import solvate
from moleculekit.tools.preparation import systemPrepare

for i, pose in enumerate(poses):
    complex_mol = protein.copy()
    complex_mol.append(pose)
    prepared, _ = systemPrepare(complex_mol, pH=7.4)
    solvated = solvate(prepared, pad=12)
    amber.build(solvated, outdir=f"./build_pose_{i:02d}",
                ionize=True, saltconc=0.15)
```

This is the docking → MD pipeline: generate poses, build a complete simulation system from each, then run all of them as parallel trajectories (the seed pool for an adaptive-binding study).

### Blind docking (whole-protein search box)

```python
poses, scores = dock(protein, ligand, numposes=20)
```

Omitting `center` and `extent` triggers the whole-protein auto-bounding-box: `dock` sets the box centre to the protein's geometric centre and the extent to the protein's bounding-box dimensions plus a 10 Å buffer on each side. Useful when the binding site is unknown - no need to compute the box yourself. Expect noisier scores and more time per dock.

### Targeted docking with multiple ligands

```python
from glob import glob
all_poses = {}
all_scores = {}
for ligand_path in glob("./compound_library/*.mol2"):
    lig = Molecule(ligand_path)
    poses, scores = dock(protein, lig, center=center, extent=extent, numposes=10)
    all_poses[ligand_path] = poses
    all_scores[ligand_path] = scores
```

For library screening you'd want a more featureful Vina wrapper (parallel docking, per-compound scoring) - this loop pattern handles small libraries (≤100 compounds) on a single workstation.

### Re-score Vina poses with FFEvaluate

```python
from ffevaluation.ffevaluate import FFEvaluate, loadParameters

# Build the complex once (any pose), then plug each pose's coordinates back in
complex0 = protein.copy(); complex0.append(poses[0])
prepared, _ = systemPrepare(complex0, pH=7.4)
amber.build(prepared, outdir="./scoring_build", ionize=True, saltconc=0.15)

# Load the built system: prmtop carries the topology (bonds, atom types, charges);
# the PDB provides the coordinates ParameterSet doesn't store.
prm  = loadParameters("./scoring_build/structure.prmtop")
mol  = Molecule("./scoring_build/structure.prmtop")
mol.read("./scoring_build/structure.pdb")
ffev = FFEvaluate(mol, prm, betweensets=("protein", "resname LIG"))

# Boolean mask for the ligand atoms in `mol`. The pose Molecules from dock()
# are ligand-only and share the same atom order, so we can slot pose coords
# straight into the ligand block of `mol.coords`.
lig_mask = mol.resname == "LIG"
assert lig_mask.sum() == poses[0].numAtoms, "ligand atom count must match"

for i, pose in enumerate(poses):
    mol.coords[lig_mask] = pose.coords
    score = ffev.calculateEnergies(mol.coords)["total"]
    print(f"pose {i}: ff-score {score:.2f} kcal/mol")
```

Useful sanity check: Vina's score isn't a real energy - cross-checking with a force-field re-score (or MM-GBSA later) can flag obviously-wrong poses.

## Gotchas

- **Inputs must be single-frame.** Both `protein` and `ligand` must have exactly one frame in `mol.coords` (third axis = 1). `dock` raises `NameError` otherwise - drop extra frames with `mol.dropFrames(keep=0)` before calling.
- **Vina expects PDBQT format.** HTMD's `dock` calls Open Babel internally to convert your `.pdb` / `.mol2` to PDBQT and back. If `obabel` is missing or returns a malformed PDBQT, `dock` raises a confusing Vina error - check `obabel --version` first.
- Hydrogens **must** be assigned before docking. Vina assumes the ligand and protein are at their correct protonation states; running `dock` on a heavy-atom-only PDB gives docked-but-wrong poses. Use `systemPrepare` upstream.
- A 20 Å × 20 Å × 20 Å search box is the usual default for a single pocket. For metal-binding sites, peptide-binding grooves, or large allosteric cavities, increase the box.
- The poses returned are **ligand-only** - you need to append the protein (or build the system from the complex) for any MD downstream.
- Vina scores are not free energies. Use them for ranking, not for thermodynamic claims.

## See also

- {doc}`How to use a custom force field with amber.build <system-build-custom-forcefield>` - feeding the docked complex into the build pipeline.
- {doc}`How to evaluate force-field energies on a frame <ffevaluate-energies>` - re-scoring poses with a real force field.
- [AutoDock Vina docs](https://vina.scripps.edu/) - flags / scoring details not exposed by HTMD's wrapper.
