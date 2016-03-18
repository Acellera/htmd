
Adaptive sampling
=================

by G. De Fabritiis

Concept
-------

-  AdaptiveRun executes adaptive runs on a given resource. It starts
   from a set of initial input directories containing all input files
   and generates new simulations based on these generators by building a
   (Markov) model of the data on-the-fly and deciding using a criteria
   where sampling is more needed.

S. Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand
binding by high-throughput molecular simulations, J. Chem. Theory
Comput. 10 (5), pp 2064–2069(2014).

Unit of execution
-----------------

Each simulation is associated to a single directory which contains all
files to run it. To run a project it is therefore necessary to provide
one or more initial simulation directories, called generators.

How it works
------------

1. Adaptive takes one or more generator simulations and construct
   on-the-fly new input coordinate files for new simulations based on
   these generators.
2. Generator simulations consist individual subdirectories in the
   ''genereratorpath'' directory.
3. Each directory must corresponds to a single simulation and names of
   directories are not important.

How it works
------------

1. Instead of launching many simulations at once, this scheme launches
   simulations in sequential batches called epochs.

2. At every epoch new input coordinatates are created while velocities
   are reinitialized to the Maxwell Boltzmann distribution at the given
   temperature.

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
simulations cannot use a velocity file. E.g. using Acemd

``python acemd.binvelocities = None # remove binvelocities for respawning acemd.temperature = 300 # set velocities generation automatically``

Adaptive script example
-----------------------

.. code:: python

    md = AdaptiveRun()
    md.nmin=5
    md.nmax=10
    md.nepochs = 30
    md.app = AcemdLocal()
    md.dryrun = True  # creates everything but does not submit anything
    md.metricsel1 = 'name CA'
    md.metricsel2 = '(resname BEN) and ((name C7) or (name C6))'
    md.metrictype = 'contacts'
    md.metricticadim = 3
    md.updateperiod = 14400 # execute every 4 hours
    md.run()

Execution in a notebook
-----------------------

1. It is possible to run the adaptive scheme syncronosly or
   asyncrounsly.
2. The command ''updateperiod'' controls this behaviour.
3. The default is to run and exist, so updateperiod needs to be
   specified if adaptive is run from the notebook

Setting a simple cron job
-------------------------

1. This is useful for having the script to execute automatically every x
   hours.
2. Do not set updateperiod then, or set it to zero such that the program
   will execute and exit

``bash #!/bin/bash -login # cron.sh file # use crontab -e to add this line: # 0 */4 * * * cd /pathtomydir/; ./cron.sh # python conf.py``

