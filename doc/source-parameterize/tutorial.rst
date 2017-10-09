Parameterization of benzamidine
===============================

This tutorial showcases the features of ``parameterize`` using benzamidine as an example.

Contents:

.. toctree::
    :maxdepth: 1

Prepare molecule
----------------

The structure of benzamidine can be obtained from PBD database. The structure comes in SDF format, but for
``parameterize`` it has to be converted to MOL2 format.

.. code:: bash

    wget https://files.rcsb.org/ligands/view/BEN_ideal.sdf
    obabel BEN_ideal.sdf -O benzamidine.mol2

``benzamidine.mol2`` file containes a molecule ready for parameterization.

Get command-line options
-------------------------

``parameterize`` command line-options can be shown using:

.. code:: bash

    parameterize -h

.. argparse::
    :ref: htmd.parameterization.cli.cli_parser
    :prog: parameterize

Set total charge
----------------

Verifying if the small molecule has the correct protonation state and charge is still a user-dependent task. Some
combinations of a total charge with a given set of atoms may cause the underlying QM software to error due to
incompatibility of the multiplicity with the total amount of electrons of the system. For example, when taking our
``ben.mol2``, and leaving the charge as the default one, deduced from the MOL2 charges (a total charge of +1), the tool
errors as it should:

.. code:: bash

    parameterize benzamidine.mol2 --no-esp --no-dihed --output benzamidine-fail

::

    RuntimeError: sanity check failed! A multiplicity of 1 with 63 electrons is impossible.
    Please check your input
    Running QM Calculations: 100% (1/1) [################################################################################################] eta --:-- -
    Traceback (most recent call last):
      File "/home/joao/miniconda3/bin/parameterize", line 6, in <module>
        sys.exit(htmd.parameterization.cli.main_parameterize())
      File "/home/joao/maindisk/software/repos/Acellera/htmd/htmd/parameterization/cli.py", line 180, in main_parameterize
        mol.minimize()
      File "/home/joao/maindisk/software/repos/Acellera/htmd/htmd/parameterization/ffmolecule.py", line 163, in minimize
        raise RuntimeError("QM Optimization failed")
    RuntimeError: QM Optimization failed

This is because that combination is wrong. For that set of atoms, with that protonation of the nitrogens, the total charge
should be set to 0 (zero). We can use the ``-c`` flag for this:

.. code:: bash

    parameterize -m benzamidine.mol2 --charge 1 --no-esp --no-dihed --output benzamidine-pass

This last command just performs the minization of the small molecule and outputs that minimized structure with the
CGENFF/GAFF2 parameters. Please note that a combination of total charge and a given set of atoms working on
``parameterize`` does not mean that protonation state is the most common or relevant in normal pH conditions.

Choose force field
------------------

The ``parameterize`` tool has an option to use a given FF as initial guess for the parameters, set by the flag
``--forcefield``. By default, it outputs parameters for both CHARMM (through CGENFF) and AMBER (through GAFF2). If one
wants parameters in the original GAFF, one may use ``--forcefield GAFF``.

The most simple use of ``parameterize``, which does not go through any optimization procedure, is to provide the guessed
initial parameters. This can be done by setting the flags ``--no-min``, ``--no-esp`` and ``--no-torsions``:

.. code:: bash

    parameterize benzamidine.mol2 --no-min --no-esp --no-dihed --output benzamidine-noqm

Inside the output directory ``benzamidine-noqm`` (specified in the flag ``--output``), one can find parameters for each both
CHARMM and AMBER.

List parameterizable dihedral angles
------------------------------------

Before doing any parameterization, one can list the soft torsions that the molecule has. This can be easily done by
using the ``--list`` flag:

.. code:: bash

    parameterize benzamidine.mol2 --list

::

    C2-C1-C7-N1
    C1-C7-N1-H6
    C1-C7-N2-H7

that are 3 sets of 4 atoms describing the 3 detected soft torsions.

Choose QM code
--------------

By default, ``parameterize`` uses an open-source QM code, PSI4 for doing the QM calculations. If one has access to
Gaussian, ``parameterize`` supports it and one may change the QM code using the flag ``--qmcode``:

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --code Gaussian --output benzamidine-qm-gaussian

Choose QM level
---------------

By default, ``parameterize`` uses settings that account for the most accurate QM settings available in the tool. This
means using a higher level of theory (B3LYP, can be set through flag ``--theory``), a larger basis set (cc-pVDZ, can be
set through the flag ``--basis``), and a solvation model for the QM calculations (turned on by default, can be
turned off by using the flag ``--vacuum``). By keeping all of the settings on default, the parameterization of the small
molecule will, however, be more computationally demanding (also depending on its number of atoms, number of soft
torsions, and on the resources available for the parameterization). One can do one of these parameterization using:

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --output benzamidine-qm


One may one to sacrifice on accuracy to increase the speed of parameterization. This can be done by using the flags
related with accuracy mentioned in the last section.

One may start by turning off the solvation model, and doing the calculations in vacuum. This can be one of the first
options if one wants to increase the speed of parameterization:

.. code:: bash

    parameterize ben.mol2 --charge 1 --vacuum --output benzamidine-vacuum

In conjunction with this option (or alone, if one wants to keep the solvation model), one may also choose a smaller
basis set, 6-31g-star:

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --vacuum --basis 6-31g-star -o ./param-ben-vacuum-lowbasis

Finally, one may parameterize with the least accuracy possible, but at the highest speed, by lowering the level of theory
to HF:

.. code:: bash

    parameterize -m ben.mol2 --charge 1 --vacuum --basis 6-31g-star --theory HF -o ./param-ben-fast

Note that the output directories are changed for simplicity. In fact, when changing these three options, the outputs are
always automatically sent to differently named directories, with the format ``<theory>-<basis-set>-<vacuum/water>``.

Control QM job execution
------------------------

Furthermore, by default, the tool runs the QM calculations on the local machine by (``-e`` flag set to inline by
default), guessing the number of CPUs to use from the maximum available in the local machine (``-n`` flag).

One may want to use less resources for a given parameterization, or maybe run two parameterizations in parallel on the
same machine. For this, one may override the number of CPUs to be used by:

.. code:: bash

    parameterize -m benzamidine.mol2 --charge 1 -ncpus 8 --output benzamidine-qm

for using 8 CPU nodes. Obviously, the less amount of CPUs used, the slower the parameterization will be.

Alternatively, one may have access for remote resources that are ready to run the QM calculations. The ``parameterize``
allows the interface with several execution back-ends, such as LSF, Slurm, and AceCloud. If one has a cluster available
running on Slurm and properly set up, one may use the ``-e`` flag to send the QM calculations to that cluster
automatically:

.. code:: bash

    parameterize -m benzamindine.mol2 --charge 1 --queue Slurm --output benzamidine-fullqm

The tool will run locally, which is very computationally inexpensive, and all the computationally expensive QM jobs
will be sent for the user to the Slurm queue system, in this case.

Restart QM jobs
---------------

TODO