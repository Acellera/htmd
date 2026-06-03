# Adaptive sampling

In this tutorial we showcase how to run adaptive sampling simulations on a
molecular system. The sample system here is the NTL9 protein.

```{note}
This tutorial drives real MD simulations through a GPU queue, so it is not
executed when the documentation is built — the code blocks below are a
read-through walkthrough. Run them on a machine with a GPU and a configured
queue to reproduce the campaign.
```

Let's import HTMD and set up some definitions:

```python
from htmd.ui import *
```

## Get the generators folder structure

Get the data for this tutorial [here](https://ndownloader.figshare.com/files/65181267).
Alternatively, download it with `wget`:

```python
import os, glob
assert os.system('wget -q https://ndownloader.figshare.com/files/65181267 -O adaptive_generators.zip && unzip -q -o adaptive_generators.zip && rm adaptive_generators.zip') == 0
for file in glob.glob('./generators/*/run.sh'):
    os.chmod(file, 0o755)
```

The `generators` folder holds one subfolder per starting structure, each a
ready-to-run simulation directory:

```text
generators
├── ntl9_1ns_0
│   ├── input
│   ├── input.coor
│   ├── input.xsc
│   ├── parameters
│   ├── run.sh
│   ├── structure.pdb
│   └── structure.psf
├── ntl9_1ns_1
│   └── ...
└── ntl9_1ns_2
    └── ...
```

## Adaptive classes

HTMD has two types of adaptive sampling:

* {py:class}`~htmd.adaptive.adaptiverun.AdaptiveMD` (free exploration)
* {py:class}`~htmd.adaptive.adaptivegoal.AdaptiveGoal` (exploration + exploitation)

Create a directory for each type of adaptive run and copy the generators into them:

```python
os.makedirs('./adaptivemd', exist_ok=True)
os.makedirs('./adaptivegoal', exist_ok=True)
shutil.copytree('./generators', './adaptivemd/generators')
shutil.copytree('./generators', './adaptivegoal/generators')
```

## AdaptiveMD

Let's change into the `adaptivemd` directory and work there:

```python
os.chdir('./adaptivemd')
```

* Set up the queue that will be used for the simulations.
* Tell it to store completed trajectories in the `data` folder, as this is
  where `AdaptiveMD` expects them to be by default.

```python
queue = LocalGPUQueue()
queue.datadir = './data'
```

```python
ad = AdaptiveMD()
ad.app = queue
```

* Set the `nmin`, `nmax` and `nepochs`:

```python
ad.nmin = 1
ad.nmax = 3
ad.nepochs = 3
```

* Choose which projection to use for the construction of the Markov model:

```python
protsel = 'protein and name CA'
ad.projection = MetricSelfDistance(protsel)
```

* Set the `updateperiod` of the adaptive run to define how often it polls for
  completed simulations and redoes the analysis:

```python
ad.updateperiod = 120  # execute every 2 minutes
```

Launch the `AdaptiveMD` run:

```python
ad.run()
```

The run processes one epoch at a time: it submits simulations, sleeps for
`updateperiod` seconds, retrieves completed trajectories, builds an MSM, picks
the next starting frames, and repeats until `nepochs` is reached.

## AdaptiveGoal

Now let's change into the `adaptivegoal` directory and work there instead:

```python
os.chdir('../adaptivegoal')
```

* Most of the class arguments are identical to `AdaptiveMD`:

```python
adg = AdaptiveGoal()
adg.app = queue
adg.nmin = 1
adg.nmax = 3
adg.nepochs = 2
adg.generatorspath = './generators'
adg.projection = MetricSelfDistance('protein and name CA')
adg.updateperiod = 120  # execute every 2 minutes
adg.goalfunction = None  # set to None just as an example
```

* It also requires the `goalfunction` argument, which defines a goal.
* We can define a variety of different goal functions.

## The goal function

The goal function will:

* take as input a {py:class}`~moleculekit.molecule.Molecule` object of a simulation, and
* produce as output a score for each frame of that simulation.
* The higher the score, the more desirable that simulation frame is for being respawned.

### RMSD goal function

For this goal function we use a crystal structure of NTL9. You can download the
structure from the [NTL9 crystal structure](http://pub.htmd.org/tutorials/adaptive-sampling/ntl9_crystal.pdb)
link and save it in the `adaptivegoal` directory, or use `wget`:

```python
assert os.system('wget -q http://pub.htmd.org/tutorials/adaptive-sampling/ntl9_crystal.pdb') == 0
```

We can define a simple goal function that uses the RMSD between the sampled
conformation and a reference (here, the crystal structure), and returns a score
to be evaluated by the `AdaptiveGoal` algorithm:

```python
ref = Molecule('./ntl9_crystal.pdb')

def mygoalfunction(mol):
    rmsd = MetricRmsd(ref, 'protein and name CA').project(mol)
    return -rmsd  # or even 1/rmsd

adg.goalfunction = mygoalfunction
```

`AdaptiveGoal` ranks conformations from high to low score. For RMSD, since we
want a lower RMSD to give a higher score, the symmetric value is returned
instead (the inverse would also work).

Launch the `AdaptiveGoal` run:

```python
adg.run()
```

### Functions with multiple arguments

The goal function can also take multiple arguments. This allows flexibility and
on-the-fly comparisons to non-static conformations (i.e. comparing against
different references as the run progresses). Here we redefine the previous goal
function with multiple arguments:

```python
def newgoalfunction(mol, crystal):
    rmsd = MetricRmsd(crystal, 'protein and name CA').project(mol)
    return -rmsd  # or even 1/rmsd
```

Now we clean the previous `AdaptiveGoal` run and start a new one with the new
goal function:

```python
# clean previous run
shutil.rmtree('./input')
shutil.rmtree('./data')
shutil.rmtree('./filtered')

# run with new goal
ref = Molecule('./ntl9_crystal.pdb')
adg.goalfunction = (newgoalfunction, (ref,))
adg.run()
```

### Other goal function examples

HTMD includes two other goal functions: the secondary-structure goal function
and the contacts goal function.

#### Secondary structure goal function

```python
ref = Molecule('./ntl9_crystal.pdb')

def ssGoal(mol, crystal):
    crystalSS = MetricSecondaryStructure().project(crystal)[0]
    proj = MetricSecondaryStructure().project(mol)
    # How many crystal SS elements match the simulation SS
    ss_score = np.sum(proj == crystalSS, axis=1) / proj.shape[1]
    return ss_score

adg.goalfunction = (ssGoal, (ref,))
```

#### Contacts goal function

```python
ref = Molecule('./ntl9_crystal.pdb')

def contactGoal(mol, crystal):
    crystalCO = MetricSelfDistance('protein and name CA', periodic=None,
                                   metric='contacts',
                                   threshold=10).project(crystal)
    proj = MetricSelfDistance('protein and name CA',
                              metric='contacts',
                              threshold=10).project(mol)
    # How many crystal contacts are seen?
    co_score = np.sum(proj[:, crystalCO] == 1, axis=1)
    co_score /= np.sum(crystalCO)
    return co_score

adg.goalfunction = (contactGoal, (ref,))
```

Many more goal functions can be devised.