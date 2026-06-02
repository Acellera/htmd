---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Build a protein-RNA complex

**You will learn:** how to build a solvated, ionised protein-RNA complex - same flow as a canonical protein build, with the RNA force field already covered by HTMD's defaults.

**Prerequisites:**
- HTMD installed.
- {doc}`Build a protein <01-protein>` covers the canonical build steps used here.

## The system

PDB `1AUD` is the **NMR solution structure of the human U1A protein bound to its own 3'-UTR polyadenylation inhibition element (PIE)**. U1A is one of the workhorse models for studying RNP recognition: the protein binds the RNA loop nucleotides through a classic RRM β-sheet, and the RNA backbone wraps around the protein surface.

The only things that change versus {doc}`tutorial 01 <01-protein>` are (a) dropping the NMR ensemble down to a single frame, and (b) the AMBER RNA force field handles the nucleotides automatically because `leaprc.RNA.OL3` is already in {py:func}`~htmd.builder.amber.defaultFf`.

## Setup

```{code-cell} python
from moleculekit.molecule import Molecule
from moleculekit.tools.autosegment import autoSegment
from moleculekit.tools.preparation import systemPrepare
from htmd.builder import amber
from htmd.builder.solvate import solvate
```

```{code-cell} python
:tags: [remove-input]
from acellera_docs_theme.molstar import show3d
```

## Step 1 - Load and pick a single frame

```{code-cell} python
mol = Molecule("1AUD")
print(f"loaded {mol.numFrames} NMR conformers, {mol.numAtoms} atoms")

mol.dropFrames(keep=0)
```

```{code-cell} python
:tags: [remove-input]
show3d(mol)
```

NMR depositions ship as an **ensemble**: 1AUD has 31 conformers. The downstream prep / parameterization / build pipeline works on a single frame at a time, so we keep frame 0. Pick a different frame, or run the prep on each frame separately, depending on your study design.

## Step 2 - Segment

```{code-cell} python
mol = autoSegment(mol, fields=("segid", "chain"))
```

{py:func}`~moleculekit.tools.autosegment.autoSegment` walks the connectivity and assigns one segment per chain. For 1AUD the result is two segments: the 30-nt RNA construct as one strand (the residue-number jump from 30 to 33 is not a chain break - the phosphodiester bond from O3'(30) to P(33) is intact at 1.61 Å, and autoSegment's spatial validation recognises this) and the protein as the other.

## Step 3 - Prepare

```{code-cell} python
prepared, specs = systemPrepare(mol, pH=7.4)
print(f"prepared system has {prepared.numAtoms} atoms")
```

{py:func}`~moleculekit.tools.preparation.systemPrepare` protonates the protein side chains at pH 7.4 the same way it does for a soluble protein. The RNA contributes nothing to the protonation problem (RNA backbones are uniformly deprotonated above pH ~2), but `systemPrepare` may emit a *"dubious protonation state"* warning for histidines whose pKa is close to the chosen pH - in 1AUD that's HIS 9, sitting near the RNA backbone where the local electrostatic field shifts its pKa. Inspect those manually if your study cares about the exact charge state.

## Step 4 - Solvate

```{code-cell} python
solvated = solvate(prepared, pad=12)
```

`pad=12` extends the simulation box by 12 Å on each side of the solute's bounding box (waters that would clash with the solute are removed, so the *actual* water shell can be a bit thinner than 12 Å around bulges). RNA backbones are flexible and need a bit more padding than a globular protein - `pad=10` would also work for a tutorial-sized system but `pad=12` is safer for the typical RNA breathing motions.

## Step 5 - Build under AMBER

```{code-cell} python
amber.build(solvated, outdir="./build", ionize=True, saltconc=0.15)
```

That single call:

- Loads `leaprc.protein.ff14SB` + `leaprc.RNA.OL3` + `leaprc.water.tip3p` (all part of {py:func}`~htmd.builder.amber.defaultFf`), so the protein, the RNA, and the water are each typed by their dedicated force field.
- Adds **Na⁺ counter-ions** to neutralise the RNA backbone. RNA carries -1 per nucleotide, so a 30-nt construct contributes -30 in total - expect ~20-25 neutralising Na⁺ (the protein also carries some net charge that offsets a few). Then enough NaCl is added on top to reach 0.15 M.
- No disulfide bonds in 1AUD; auto-detection finds nothing and moves on.

The output `./build/structure.prmtop` + `structure.pdb` (and `structure.crd`, the AMBER inpcrd that `setup_equilibration` actually reads) is now ready for an equilibration protocol (e.g. {py:func}`acemd.protocols.setup_equilibration`).

## Gotchas

- For an NMR ensemble, decide upfront whether you want one frame (default) or all of them. Each frame is a slightly different conformer, so running prep + build on each gives you an ensemble of starting structures for replica MD.
- RNA has **roughly one negative backbone charge per nucleotide**, so even a small protein-RNA system needs a substantial number of neutralising Na⁺ ions on top of the 0.15 M NaCl. Don't be surprised if `Adding N cations` is in the dozens.
- For protein-DNA complexes the flow is identical - `leaprc.DNA.bsc1` is also in {py:func}`~htmd.builder.amber.defaultFf`, so a DNA-only or DNA-bound build works the same way (`Molecule("1BNA")` is the Dickerson dodecamer fixture used elsewhere as a reference).
- If you have a modified nucleotide (5mC, m6A, ...) treat it as a non-canonical residue: detect → SMILES template → `parameterizeFromSpecs` → `amber.build`, exactly as in {doc}`tutorial 02 <02-protein-ligand>`. The cluster pipeline handles nucleic-acid backbones the same way it handles peptide backbones.

## See also

- {doc}`Build a protein <01-protein>` - the canonical build this tutorial mirrors.
- {doc}`Build a protein with a ligand <02-protein-ligand>` - what to add when the nucleic acid carries non-canonical modifications.
- {doc}`System-building overview <../../explanation/system-building>` - where the nucleic-acid force fields fit in the wider stack.
