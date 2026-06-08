# How to choose a featurization with a VAMP-2 score

## Goal

Decide *which* features to build your MSM on - and how many TICA dimensions and what TICA lag to keep - objectively rather than by intuition. The cross-validated VAMP-2 score ({py:func}`~htmd.projections.vamp.vampScore`) rates how much of the slow dynamics a featurization resolves; you score your candidates and keep the best one. This is the variational feature-selection idea of McGibbon and Pande.

## Minimal example

Starting from your `sims` (a {py:func}`~htmd.simlist.simlist`), score two candidate featurizations of the same trajectories and compare. Higher is better:

```python
from htmd.ui import Metric, MetricSelfDistance, MetricDihedral, vampScore

# Candidate A: Ca-Ca contact map
metrA = Metric(sims)
metrA.set(MetricSelfDistance("protein and name CA", metric="contacts", periodic="chains"))
dataA = metrA.project()
dataA.fstep = 0.1

# Candidate B: backbone dihedrals
metrB = Metric(sims)
metrB.set(MetricDihedral())
dataB = metrB.project()
dataB.fstep = 0.1

for label, data in [("contacts", dataA), ("dihedrals", dataB)]:
    mean, std = vampScore(data, lag=2, dim=4, units="ns")
    print(f"{label}: VAMP-2 = {mean:.2f} +/- {std:.2f}")
```

Keep the featurization with the highest cross-validated score.

## How to read the score

- **Higher is better** - a higher VAMP-2 score means the features resolve more of the slow dynamics.
- **About 1.0 is the floor** - a featurization that captures no slow process (pure noise) scores around 1.0, the trivial stationary process.
- **Compare like with like** - scores are only comparable when computed at the **same `lag` and `dim`**. Changing either changes the absolute score.
- **It is cross-validated** - the score is averaged over `nfolds` folds, so a featurization that looks good only by overfitting the training data is penalised on the held-out folds. Use the mean and watch the spread (`std`).

## Common variations

Choose how many dimensions to keep. The score keeps rising as you add dimensions (it sums the squared singular values), so look for **diminishing returns** rather than a maximum:

```python
for d in (2, 3, 4, 6, 8):
    mean, std = vampScore(data, lag=2, dim=d, units="ns")
    print(f"dim={d}: VAMP-2 = {mean:.2f} +/- {std:.2f}")
```

Check robustness to the TICA lag:

```python
for lag in (0.5, 1, 2, 5):
    mean, std = vampScore(data, lag=lag, dim=4, units="ns")
    print(f"lag={lag} ns: VAMP-2 = {mean:.2f} +/- {std:.2f}")
```

Get the per-fold scores instead of the mean and standard deviation (for example to plot your own error bars):

```python
scores = vampScore(data, lag=2, dim=4, units="ns", return_scores=True)
```

If you have only a few long trajectories, cross-validating over whole trajectories is not possible; split each trajectory into blocks of `blocksize` frames instead:

```python
mean, std = vampScore(data, lag=2, dim=4, units="ns", blocksize=1000)
```

## Gotchas

- **Only compare at matched `lag` and `dim`.** A featurization scored at `dim=6` will out-score one at `dim=3` purely because of the extra dimensions, not because it is better.
- **Do not simply maximise over `dim`.** Because VAMP-2 sums squared singular values, it rises monotonically with `dim`. Read the *shape* of the curve (where the gains flatten), not the peak.
- **A tie is a real and useful answer.** Featurizations that both capture the slow process well score about the same - the score is telling you the choice does not matter for the kinetics, so decide on other grounds such as clustering robustness or cost. For example, on the villin headpiece a Cα contact map and raw Cα distances score about equally; the contact map is preferred there for its clustering robustness, not because it resolves more dynamics.
- **Needs at least 2 trajectories** for the default whole-trajectory cross-validation. With a single long trajectory, pass `blocksize` so it can be split into independent blocks.
- **Scoring is not validation.** A high VAMP-2 score says the features are promising; it does not replace building the MSM and checking the implied timescales and the Chapman-Kolmogorov test.

## See also

- {py:func}`htmd.projections.vamp.vampScore` - the API reference.
- {doc}`How to read an implied-timescales plot <msm-interpret-its-plot>` - the validation step after picking features.
- {doc}`Protein folding MSM <../tutorials/analysis/villin-folding>` and {doc}`Ligand binding MSM <../tutorials/analysis/trypsin-benzamidine-binding>` - the full pipelines these features feed into.
