# Adaptive Bandit

```{note}
This tutorial drives real MD simulations through a GPU queue, so it is not
executed when the documentation is built — the code blocks below are a
read-through walkthrough. Run them on a machine with a GPU and a configured
queue to reproduce the campaign.
```

## Concept

* Adaptive sampling algorithms usually employ empirical policies, and they are
  not based on any mathematical decision process.
* We describe adaptive sampling in terms of a multi-armed bandit problem to
  develop a novel adaptive sampling algorithm, Adaptive Bandit, providing strong
  fundamentals to tackle the exploration-exploitation dilemma faced in adaptive
  sampling.
* Adaptive Bandit is framed into a reinforcement-learning based framework, using
  an action-value function and an upper confidence bound selection algorithm,
  improving adaptive sampling's performance and versatility when faced against
  different free energy landscapes.
* Discretized conformational states are defined as actions, and each action has
  an associated reward distribution. When an action is picked, the algorithm
  computes the associated reward for that action, based on MSM free-energy
  estimations, and applies a policy to select the next action.
* `AdaptiveBandit` relies on the UCB1 algorithm to optimize the action-picking
  policy, defining an upper confidence bound for each action based on the number
  of times the agent has picked that action and the total number of actions
  taken:

$$
a_t = \arg\max_{a\in\mathcal{A}}\left[ Q_t(a) + c\sqrt{\frac{\ln t}{N_t(a)}} \right]
$$

A. Pérez, P. Herrera-Nieto, S. Doerr and G. De Fabritiis,
[AdaptiveBandit: A multi-armed bandit framework for adaptive sampling in molecular simulations](https://arxiv.org/abs/2002.12582).
arXiv preprint 2020; arXiv:2002.12582.

Auer P. [Using confidence bounds for exploitation-exploration trade-offs](http://www.jmlr.org/papers/volume3/auer02a/auer02a.pdf).
Journal of Machine Learning Research. 2002; 3(Nov):397-422.

## Getting started

This tutorial shows you how to properly set up an `AdaptiveBandit` project,
highlighting the main differences with respect to standard adaptive sampling. As
an example, we will perform some folding simulations using the chicken villin
headpiece (PDB: 2F4K).

Let's start by importing HTMD and the {py:class}`~htmd.adaptive.adaptivebandit.AdaptiveBandit` class:

```python
from htmd.ui import *
from htmd.adaptive.adaptivebandit import AdaptiveBandit
```

`AdaptiveBandit` uses the same project structure as adaptive sampling, with each
simulation associated to a single directory containing all the files needed to
run it.

To begin, get the starting generators [here](https://ndownloader.figshare.com/files/65181204).
You can also download the data with `wget -O gen.tar.gz https://ndownloader.figshare.com/files/65181204`.
Uncompress that `tar.gz` file and allow execution on all `run.sh` files:

```python
for file in glob('./generators/*/run.sh'):
    os.chmod(file, 0o755)
```

These generators contain prepared unfolded structures of villin, which we want
to simulate long enough to reach the folded native structure.

## AdaptiveBandit

We start our `AdaptiveBandit` project the same way as with adaptive sampling, by
defining the queue used for simulations.

```python
queue = LocalGPUQueue()
queue.datadir = './data'
```

```python
ab = AdaptiveBandit()
ab.app = queue
```

Then we define `nmin`, `nmax` and `nframes` to set the maximum amount of
simulated frames:

```python
ab.nmin = 0
ab.nmax = 2
ab.nframes = 1000000
```

And we choose the projection and clustering method used to construct a Markov
model at each epoch:

```python
ab.clustmethod = MiniBatchKMeans
ab.projection = MetricSelfDistance('protein and name CA')
```

Up until now, the setup is exactly the same as with `AdaptiveMD`. However,
`AdaptiveBandit` has an additional parameter, which sets the $c$ parameter from
the UCB1 equation:

```python
ab.exploration = 0.01
```

Additionally, `AdaptiveBandit` accepts a goal function as input that is used to
initialize the action-value estimates. In this example we use the contacts goal
function defined in the previous tutorial to initialize the $Q(a)$ values. The
`goal_init` parameter sets an $N_t(a)$ initial value proportional to the maximum
frames per cluster at the end of the run, which represents the statistical
certainty we give to the goal function.

```python
ref = Molecule('2F4K')

def contactGoal(mol, crystal):
    crystalCO = MetricSelfDistance('protein and name CA', periodic=None,
                                   metric='contacts',
                                   threshold=10).project(crystal).squeeze()
    proj = MetricSelfDistance('protein and name CA',
                              metric='contacts',
                              threshold=10).project(mol)
    # How many crystal contacts are seen?
    co_score = np.sum(proj[:, crystalCO] == 1, axis=1).astype(float)
    co_score = co_score / np.sum(crystalCO)
    return co_score

ab.goalfunction = (contactGoal, (ref,))
ab.goal_init = 0.3
```

And now we just need to launch our `AdaptiveBandit` run:

```python
ab.run()
```
