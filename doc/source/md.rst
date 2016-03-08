HTMD Basic Molecular Dynamics
=============================


Introduction
------------

HTMD optionally includes simple Molecular Dynamics capability. This MD engine is less performant and has fewer capabilities than ACEMD but is adequate for simple, short simulations.

Installation
------------

The ligand parameterisation tool is installed using the command:

.. parsed-literal::
  
    conda install mdy -y

MD requires a GPU for best performance. In the absense of a GPU, a multicore CPU can be used.

Execution
---------

The program can be run from the command-line with 'mdy':

.. parsed-literal::

 $ mdy

 HTMD Molecular Dynamics 2016 (c) Acellera

  Syntax: cli [args]:
    --input [inputfile]      : Specify the input file. Default value 'input'
    --device [device-number] : Set the GPU to run on 
    --platform [name]        : OpenMM platform. Defaults to fastest available
    --command ([command])    : Show help for an input-file command
    --verbose                : Show full simulation configuration before run
    --ff                     : Show supported force-field parameter files
    --license                : Show the software license
    --help

 HTMD Molecular Dynamics incorporates OpenMM. See http://wiki.simtk.org/openmm/License
 No re-distribution in whole or part

As with ACEMD, simulations are configured in an input file. MD's input file grammar is a simplified version of ACEMD's. 
The full command set, or help on an individual command, can be obtained with the '--command' option, eg:

.. parsed-literal::

 $ mdy --command run

   run [ number (ps|ns|us) ]

 The length of simulation to run.

   Default: 1000 


The GPU to be used by the simulation can be selected with '--device'. Likewise, the exeuction platform can be selected explicitly ( '--platform CUDA | OpenCL | CPU' ). The defualt behaviour is to use the fastest mode available.


Where commands differ between ACEMD and MD, a warning will be printed out, eg:

.. parsed-literal::

 Command 'hydrogenscale'   is deprecated and is no longer required
 Command 'langevin'  is deprecated and replaced by 'thermostat'

Deprecated commands can be used in preference to their replacements, but will be removed in a future version.

MD from Python
--------------

It's also possible to invoke an MD simulation  directly from Python code:

.. parsed-literal::
 from mdy import Configuration, Simulation

.. parsed-literal::
 # Set up a Configuration object.
 # The fields match the commands available in the input file

 config = Configuration()
 config.coordinates = "system.pdb"
 config.structure   = "system.psf"
 config.parameters  = "all_charmm36_prot.inp"
 config.minimize    = 500
 config.run         = 1000000


.. parsed-literal::
 # Start a simulation

 sim = Simulation( config, device=0, platform="CUDA"  )


The simulation will run to completion before returning.
