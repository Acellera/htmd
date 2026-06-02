# How to evaluate force-field energies on a frame

## Goal

Compute bonded + non-bonded force-field energies (and optionally per-atom forces) for any coordinate snapshot, without going through an MD engine. Useful for ranking docking poses, validating a build, computing per-frame interaction energies, or scripting custom energy-based analyses.

## Minimal example

```python
from moleculekit.molecule import Molecule
from ffevaluation.ffevaluate import FFEvaluate, loadParameters

mol = Molecule("./build/structure.prmtop")
mol.read("./build/structure.pdb")

prm = loadParameters("./build/structure.prmtop")          # parmed ParameterSet
ffev = FFEvaluate(mol, prm)
energies, forces, per_atom = ffev.calculate(mol.coords)

print("total energy:", energies.sum())                    # kcal/mol
```

`FFEvaluate.calculate(coords)` returns three arrays:

- `energies` - `(6, n_frames)` per-term system totals in the order `(bond, vdw, elec, angle, dihedral, improper)`.
- `forces` - `(n_atoms, 3, n_frames)` per-atom force vector in kcal/mol/Å.
- `atmnrg` - `(n_atoms, 6, n_frames)` per-atom energy contributions in the same term order.

`FFEvaluate.calculateEnergies(coords)` is a thin wrapper that calls `calculate` and discards both `forces` and `atmnrg`, returning only the energies (as a dict when `formatted=True`).

Both `calculate` and `calculateEnergies` accept an optional `box=mol.box` argument (`(3, n_frames)`). Pass it for periodic systems so non-bonded distances get wrapped through the box; omit it for vacuum / non-periodic snapshots.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `mol` | Molecule whose topology defines the bond / angle / dihedral graph. Frames in `mol.coords` define how many coordinate snapshots get evaluated. |
| `prm` | A parmed `ParameterSet` - load from prmtop with `loadParameters("...")` or build from a CHARMM `.prm` via `loadParameters(prm_path)`. |
| `betweensets` | A tuple of two atom-selection strings - restricts non-bonded energy / force calc to interactions **between** the two sets. Computes only LJ + electrostatics, no intramolecular bonded terms. |
| `cutoff` | Non-bonded cutoff in Å. `0` (default) means no cutoff - exact electrostatics + LJ. Set to e.g. `9.0` to match a simulation's cutoff. |
| `rfa` | Use with `cutoff` to enable the reaction-field approximation for electrostatics beyond the cutoff. |
| `solventDielectric` | Solvent dielectric used when `rfa=True`. Default 78.5 (water). |

## Common variations

### Per-frame energies on a trajectory

```python
mol.read("./output.xtc")                                  # multi-frame
all_energies, all_forces, _ = ffev.calculate(mol.coords)
# all_energies shape: (n_terms, n_frames)
```

### Interaction energy between two selections

```python
ffev = FFEvaluate(mol, prm, betweensets=("resname LIG", "protein"))
e, _, _ = ffev.calculate(mol.coords)
print("protein-ligand interaction energy:", e.sum(axis=0))
```

Restricting to a pair of sets is much cheaper than the full system - no bonded terms are computed and the non-bonded loop is restricted to inter-set pairs.

### Use a CHARMM parameter file instead of prmtop

```python
prm = loadParameters("./params.prm")                      # CHARMM .prm
ffev = FFEvaluate(mol, prm)
```

`loadParameters` dispatches by extension - only `.prmtop`, `.prm`, and `.frcmod` are recognised. Any other extension raises `RuntimeError`. For CHARMM systems built from a PSF, point `loadParameters` at the matching `.prm`; for AMBER pass the `.prmtop` directly.

### Decomposed energy report

```python
e_by_term = ffev.calculateEnergies(mol.coords, formatted=True)
# Dict keyed by term name. Keys: bond, angle, dihedral, improper, vdw, elec, total.
for k, v in e_by_term.items():
    print(f"{k:>10}: {v:.3f} kcal/mol")
```

`calculateEnergies(..., formatted=True)` is the readable variant of `calculate(...)`.

## Gotchas

- `FFEvaluate` expects the same atom order in `mol.coords` and in the parameter set. If your `mol` has been reordered / filtered after the parmtop was loaded, the bond graph won't match. Always load both from the same source.
- With `betweensets`, the result excludes intra-set bonded **and** intra-set non-bonded energies. The reported number is the cross-set interaction only.
- Default `cutoff=0` is exact but O(N²) - on >100k-atom systems this is slow. Setting `cutoff > 0` only short-circuits the per-pair maths *inside* the loop; the pair iteration itself is still O(N²) (there's no neighbour list). Match your simulation's cutoff (typically 9-12 Å) for production-scale evaluations, but don't expect linear-scaling.
- The output unit is kcal/mol regardless of the parameter file format.

## See also

- {doc}`How to use a custom force field with amber.build <system-build-custom-forcefield>` - producing the prmtop that `loadParameters` consumes.
- [ffevaluation on GitHub](https://github.com/Acellera/ffevaluation) - the full API + benchmarks.
