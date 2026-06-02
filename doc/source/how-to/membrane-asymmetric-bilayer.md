# How to set up an asymmetric bilayer

## Goal

Build a membrane where the upper and lower leaflets have **different lipid compositions** - e.g. plasma membrane mimics where PS / PE / cholesterol are enriched on the inner leaflet, or experimental rafts with cholesterol concentrated on one side only.

## Minimal example

```python
from htmd.membranebuilder.build_membrane import buildMembrane

memb = buildMembrane(
    [80, 80],
    ratioupper={"popc": 0.7, "chl1": 0.2, "popg": 0.1},   # outer leaflet
    ratiolower={"pope": 0.45, "pops": 0.15,               # inner leaflet
                "popc": 0.20, "chl1": 0.20},
    minimize=300, equilibrate=1, platform="CUDA",
    outdir="./membrane",
)
```

`ratioupper` and `ratiolower` are independent dicts of **relative ratios** (auto-normalised by the builder via `ratiosAPL / ratiosAPL.sum()`); you can pass any positive numbers - `{"popc": 7, "chl1": 3}` is equivalent to `{"popc": 0.7, "chl1": 0.3}`. The builder places lipids per-leaflet according to those fractions, then runs a short LJ-fluid relaxation + (optional) OpenMM minimisation / equilibration to remove packing artefacts.

## Parameters that matter

| Parameter | What it does |
| --- | --- |
| `ratioupper`, `ratiolower` | Per-leaflet `{lipidname: mole_fraction}` dicts. The keys must be names returned by {py:func}`~htmd.membranebuilder.build_membrane.listLipids` (or by your custom `lipidf=` library). |
| `xysize` | First **positional** argument: a 2-element list `[Lx, Ly]` of XY box dimensions in Å. Make each ~10 Å wider than the protein in that direction if you're embedding one. |
| `minimize`, `equilibrate` | Steps and ns of OpenMM relaxation after the initial LJ-fluid placement. For asymmetric bilayers, set `equilibrate >= 1` (ns) - asymmetry creates initial lateral pressure that needs to relax. |
| `head_restraint_k` | Restraint on the head-group `z` position during the builder's **internal** OpenMM minimise / equilibrate step. Use a small value (~0.5 kJ/mol/Å²) to keep the leaflets from spontaneously flipping while packing relaxes. The restraint is dropped before the function returns - it never reaches production. |

## Common variations

### Outer leaflet enriched in cholesterol (raft-like)

```python
buildMembrane(
    [100, 100],
    ratioupper={"popc": 0.4, "dppc": 0.2, "chl1": 0.4},   # cholesterol-rich raft
    ratiolower={"popc": 0.7, "pope": 0.2, "chl1": 0.1},
    minimize=500, equilibrate=2, platform="CUDA",
    outdir="./raft-membrane",
)
```

### Inner leaflet only carrying the anionic lipid

```python
buildMembrane(
    [80, 80],
    ratioupper={"popc": 1.0},                             # zwitterionic only
    ratiolower={"popc": 0.7, "pops": 0.3},                # anionic enriched
    minimize=300, equilibrate=1, platform="CUDA",
    outdir="./asymm-charge",
)
```

The net `PS` charge on the inner leaflet is balanced by counter-ions added in the subsequent {py:func}`~htmd.builder.amber.build` call (with `ionize=True, saltconc=0.15`).

### With a protein embedded

```python
memb = buildMembrane(
    [100, 100],
    ratioupper={"popc": 0.7, "chl1": 0.3},
    ratiolower={"pops": 0.3, "popc": 0.5, "chl1": 0.2},
    solute=prepared,                                       # OPM-aligned protein
    minimize=300, equilibrate=1, platform="CUDA",
    outdir="./membrane",
)
```

The `solute=` argument carves the protein footprint out of both leaflets - asymmetric carve depths are handled automatically based on the protein's actual `z` range in each leaflet.

## Gotchas

- The leaflet ratios are auto-normalised - the builder does `ratios / ratios.sum()` internally, so `{"popc": 0.7, "chl1": 0.3}`, `{"popc": 7, "chl1": 3}`, and `{"popc": 0.4, "chl1": 0.17}` (incomplete) all behave identically. The "sum to 1.0" convention is for *your* readability only.
- Asymmetric bilayers have a small net lateral pressure that relaxes during equilibration. `buildMembrane`'s own `equilibrate > 0` step **already** uses OpenMM's `MonteCarloMembraneBarostat(XYIsotropic, ZFree)` internally — the membrane gets properly relaxed before the function returns. The relevant gotcha is for the **downstream ACEMD production** stage: there, use `barostatconstratio=True` on {py:func}`acemd.protocols.setup_equilibration` so X and Y scale together (as a fixed ratio) while Z relaxes independently, matching the constraint that lipid area-per-head is isotropic in the XY plane while bilayer thickness settles separately.
- Cholesterol (`chl1`) tends to flip-flop on simulation timescales (µs), even though it's placed in the leaflet you specify. If maintaining asymmetric cholesterol composition matters for your science, add a positional restraint on the cholesterol oxygens for the first few ns of production.
- Spontaneous lipid flip-flop of phospholipids is rare (>ms timescales), so asymmetric phospholipid composition is preserved naturally.

## See also

- {doc}`Build a membrane-embedded protein <../tutorials/system-prep/07-membrane>` - canonical membrane build (symmetric, but the syntax is identical).
- {doc}`How to embed a protein in a pre-equilibrated membrane <membrane-embed-preequilibrated>` - alternative when you want to drop a protein into a membrane that already exists.
