Adaptive sampling
=================

Concept
-------

-  Exploration of the conformational space using MD simulations can
   waste lots of simulation time sampling the same conformational
   regions which does not provide any new knowledge.

-  Instead it would be desirable to explore more under-sampled regions
   of the conformational space, to overcome energetic barriers and
   eventually reach the desired conformation (folded protein / bound
   ligand etc.)

-  In adaptive sampling, instead of launching thousands of simulations
   at once from a small set of initial structures as in naive
   high-throughput MD, simulations are launched in sequential batches
   called epochs utilizing knowledge of the conformational space
   obtained from all previous epochs.

-  The starting points of the simulations in each epoch are chosen based
   on some criteria; in this case, it selects conformations from the
   most under-sampled conformational regions detected in all previous
   simulations. This is done by using Markov state models which
   discretize the conformational space into a set of most metastable
   states. Then the starting conformations are sampled based on a
   distribution related to the population of each state.

S. Doerr and G. De Fabritiis, `On-the-fly learning and sampling of
ligand binding by high-throughput molecular
simulations <http://pubs.acs.org/doi/abs/10.1021/ct400919u>`__, J. Chem.
Theory Comput. 10 (5), pp 2064–2069(2014).

S. Doerr, M. J. Harvey, Frank Noé, and G. De Fabritiis, `HTMD:
High-Throughput Molecular Dynamics for Molecular
Discovery <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049>`__ J.
Chem. Theory Comput. 2016 12 (4), 1845-1852

Unit of execution
-----------------

Each simulation in adaptive is associated to a single directory which
contains all files to run it. To run a project it is therefore necessary
to provide one or more initial simulation directories, called
generators.

How to start
------------

Adaptive creates multiple directories which it uses to organize the
simulations. The user only needs to provide a generators folder
containing one sub-folder for each starting conformation containing all
files needed to execute that simulation. For example:

::

    └── generators/
        ├── gen1/
        │   ├── structure.pdb
        │   ├── input
        │   └── ...
        ├── gen2/
        │   ├── structure.pdb
        │   ├── input
        │   └── ...

Then the adaptive will generate an ``input``, ``data`` and later a
``filtered`` folder as well, looking like this:

::

    ├── data/          # Contains the completed simulations (automatically created)
    ├── filtered/      # Contains the completed simulations without water atoms (automatically created)
    ├── generators/    # Contains the initial generators provided by the user
    └── input/         # Contains the files needed to start all simulations of all epochs (automatically created)

Adaptive uses a naming scheme for simulations which follows the pattern:
``e4s3_e2s1p0f45``. This name tells us that this simulation was
generated in epoch 4 as the 3rd simulation of the batch. The starting
conformation was taken from simulation 1 of epoch 2 from the first piece
of the simulation [*]_ and from frame 45 of that simulation piece.

.. [*] some MD software might fragment simulations into pieces. Usually
       though this number will be 0 and can be ignored.

Simulation length
-----------------

1. The length of each simulation is really system dependent.
2. It could be anything like tens of nanoseconds to hundred of
   nanoseconds.
3. As a rule of thumb use twice the expected lag time for your molecular
   process (e.g. for binding anything between 30 and 100 ns).

Simulation details
------------------

As only the coordinates files are seeded for every new epoch,
simulations cannot use a velocity file. Velocities are therefore
reinitialized to the Maxwell Boltzmann distribution at the given
temperature.

E.g. if setting up the simulations with the :class:`~htmd.protocols.production_v5.Production` class:

.. code:: python

    from htmd.protocols.production_v5 import Production
    md = Production()
    md.adaptive = True
    [...]

or directly modifying the ACEMD ``input`` file of the simulations and
removing the binvelocities line.

Adaptive script example
-----------------------

The power of adaptive sampling is accessible on HTMD through the :class:`~htmd.adaptive.adaptiverun.AdaptiveMD` class:

.. code:: python

    from htmd.ui import *
    app = LocalGPUQueue()
    app.datadir = './data'
    md = AdaptiveMD()
    md.nmin=5
    md.nmax=10
    md.nepochs = 30
    md.app = app
    md.projection = MetricDistance('name CA', '(resname BEN) and ((name C7) or (name C6))', metric='contacts')
    md.ticadim = 3
    md.updateperiod = 14400 # execute every 4 hours
    md.run()

Execution in a notebook
-----------------------

1. It is possible to run the adaptive scheme syncronosly or
   asyncrounsly.
2. The option ``updateperiod`` controls this behaviour.
3. The default is to run and exit, so ``updateperiod`` needs to be specified
   if adaptive should be run synchronously

Setting a simple cron job
-------------------------

1. This is useful for having the script execute automatically every x
   hours.
2. Do not set ``updateperiod`` then, or set it to zero such that the
   program will execute and exit

.. code:: bash

    #!/bin/bash -login
    # cron.sh file
    # use crontab -e to add this line:
    # 0 */4 * * * cd /pathtomydir/; ./cron.sh
    #
    python conf.py

Visualizing the starting conformations
--------------------------------------

If we want to look at what structures were chosen as the starting
conformations of a given epoch we can use a code snippet like the
following:

.. code:: python

    for s in glob('input/e28s*'):  # Visualize all starting conf of epoch 28
       mol = Molecule(s+'/structure.pdb')
       mol.read(s+'/input.coor')
       mol.view()
