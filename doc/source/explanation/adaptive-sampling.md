# Adaptive sampling

Exploring conformational space with brute-force MD wastes simulation time on regions you've already sampled. Adaptive sampling addresses this by running simulations in **sequential batches called epochs**, using the data from previous epochs to decide where the next epoch should start. The selection logic favors **under-sampled conformational regions**, which is what lets adaptive runs cross energetic barriers and reach the configurations that brute force can take orders of magnitude longer to find.

The mechanism: after each epoch HTMD discretises the conformational space sampled so far with a Markov state model, and samples the next epoch's starting frames from a distribution related to the population of each metastable state. Rare states get over-weighted; well-sampled ones get under-weighted.

The two papers below set out the algorithm and its first applications:

- S. Doerr and G. De Fabritiis. *On-the-fly learning and sampling of ligand binding by high-throughput molecular simulations.* J. Chem. Theory Comput. **2014**, 10 (5), 2064-2069. <https://pubs.acs.org/doi/abs/10.1021/ct400919u>
- S. Doerr, M. J. Harvey, F. Noé, and G. De Fabritiis. *HTMD: High-throughput molecular dynamics for molecular discovery.* J. Chem. Theory Comput. **2016**, 12 (4), 1845-1852. <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049>

## Unit of execution

Each simulation in an adaptive run is associated with a single directory containing everything needed to launch it (topology, input coordinates, MD config). To start an adaptive project you provide one or more **generators**: subdirectories of a `generators/` folder, one per starting conformation:

```text
└── generators/
    ├── gen1/
    │   ├── structure.pdb
    │   ├── input
    │   └── ...
    ├── gen2/
    │   ├── structure.pdb
    │   ├── input
    │   └── ...
```

As the run progresses HTMD creates three more folders:

```text
├── data/        # Completed simulations (auto-created)
├── filtered/    # Completed simulations with waters stripped (auto-created when filter=True,
│                # the default; selection defaults to "not water" via filtersel=)
├── generators/  # Initial generators you provided
└── input/       # Per-epoch starting directories (auto-created). Each new sim's input dir
                 # carries the coorname (default input.coor) and the boxname
                 # (default input.xsc) written from the chosen respawn frame.
```

## Naming scheme

Simulations are named with the pattern `e4s3_e2s1p0f45`. Parsed:

- `e4s3` - generated in **epoch 4**, the **3rd simulation** of that batch.
- `e2s1p0f45` - the starting conformation came from **epoch 2**, **simulation 1**, **piece 0**, **frame 45**.

Some MD engines split long simulations into pieces; the piece index is usually 0 and can be ignored.

## Simulation length

Per-simulation length is system-dependent. As a rule of thumb use **about twice the expected lag time** for the molecular process you're studying (e.g. 30-100 ns per simulation for ligand binding). Each frame seeds only coordinates - velocities are re-initialised from the Maxwell-Boltzmann distribution at the configured temperature, so velocity files are not transferred between epochs.

## Sync vs async execution

{py:class}`~htmd.adaptive.adaptiverun.AdaptiveMD` can run two ways:

- **Asynchronous** (default): the script launches the next epoch when prerequisites are met, then exits. You re-launch it on a schedule (e.g. cron).
- **Synchronous**: set `updateperiod` to a non-zero value. The script blocks, sleeping `updateperiod` seconds between adaptive iterations, and only exits once the run's stop condition fires (`nepochs` reached, `nframes` reached, or an error).

For interactive notebook work, synchronous mode is usually easier. For long campaigns on a queue, async + cron is more robust to interruptions.

## What to do next

- {doc}`Adaptive sampling tutorial <../tutorials/adaptive/adaptive-sampling>` walks through a full project.
- {doc}`Adaptive bandit tutorial <../tutorials/adaptive/adaptive-bandit>` covers the bandit variant.
- {py:class}`~htmd.adaptive.adaptiverun.AdaptiveMD` is the main class; {py:class}`~htmd.adaptive.adaptivebandit.AdaptiveBandit` and {py:class}`~htmd.adaptive.adaptivegoal.AdaptiveGoal` are the other entry points.
