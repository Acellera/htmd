Parameterization of benzamidine
===============================

This tutorial showcases the features of ``parameterize`` using benzamidine as an example.

Contents:

.. toctree::
    :maxdepth: 1

Prepare molecule
----------------

The structure of a molecule has to be in MOL2 format.

.. note::

    The protonation state and charge of the molecule has to be correct.

For benzamidine, it can be downloaded from *HTMD* repository:

.. code:: bash

    wget https://raw.githubusercontent.com/Acellera/htmd/master/htmd/data/test-param/benzamidine.mol2

``benzamidine.mol2`` file contains a molecule ready for parameterization.

Get command-line options
-------------------------

``parameterize`` command-line options can be shown using:

.. code:: bash

    parameterize -h

.. argparse::
    :ref: htmd.parameterization.cli.cli_parser
    :prog: parameterize

Set total charge
----------------

The total charge of the molecule is set to a sum of the atomic partial charges from the structure file. If necessary,
it can be overridden with ``--charge`` flags, i.e. ``--charge -2`` set the total charge of a molecule to -2.

.. warning::

    An incorrect combination of the total charge and a protonation state may result into QM code failure.

In case of benzamidine, it is in a protonate state with the total charge of +1. ``benzamidine.mol2`` does not have
partial atomic charges, so ``--charge 1`` has to be set.

Choose force field
------------------

The initial forcefield guess for AMBER

TODO: finish!
``parameterize`` support


``parameterize`` has an option to use a given FF as initial guess for the parameters, set by the flag
``--forcefield``. By default, it outputs parameters for both CHARMM (through CGENFF) and AMBER (through GAFF2). If one
wants parameters in the original GAFF, one may use ``--forcefield GAFF``.


The tool can be used to obtain just GAFF/GAFF2 or CGenFF parameters, i.e. no QM calculation are performed. This can be
done by setting ``--no-min``, ``--no-esp`` and ``--no-torsions`` flags.

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --forcefield GAFF2 CGENFF --no-min --no-esp --no-dihed --outdir initial

The parameters are written to ``intial`` directory (specified with ``--outdir`` flag):::



List parameterizable dihedral angles
------------------------------------

Parametrizable dihedral angles for a given molecule can be listed using ``--list `` flag.

.. code:: bash

    parameterize benzamidine.mol2 --list

::
 === Parameterizable dihedral angles ===

  C1-C7-N1-H8
  C2-C1-C7-N1


.. note::

    Symmetry equivalent dihedral angles are taken into account and are not shown in the list.

Choose QM code
--------------

By default, *Psi4* is used for all QM calculations. QM code can be changed with ``--code`` flag, i.e.
``--code Gaussian`` switches *Psi4* to *Gaussian 09*.

.. note::

    *Gaussian 09* is not distributed with *HTMD*. It has to be installed separately.

Choose QM level
---------------

The default QM level is the denisty functional theory (DFT) with B3LYP exchange-correlation functional and DFT-D3
dispersion correction. The level of theory can be changed with ``--theory`` flag, i.e. ``--theory HF`` switches to
Hartree-Fock metohd.

The default basis sets are ``cc-pVZD``, though for a negatively charged molecule, more diffuse ``aug-cc-pVZD`` are used.
The basis sets can be changed with ``--basis`` flag, i.e. ``--basis 6-31G*``.

The default QM environment (solvation model) is vacuum. It can be changed with ``--environment`` flag, i.e.
``--environment PCM`` switches to the polarizable continuum model (PCM).

TODO: finish!

By default, ``parameterize`` uses settings that account for the most accurate QM settings available in the tool.

By keeping all of the settings on default, the parameterization of the small
molecule will, however, be more computationally demanding (also depending on its number of atoms, number of soft
torsions, and on the resources available for the parameterization). One can do one of these parameterization using:

Note that the output directories are changed for simplicity. In fact, when changing these three options, the outputs are
always automatically sent to differently named directories, with the format ``<theory>-<basis-set>-<vacuum/water>``.

Control QM job execution
------------------------

TODO: finish!

By default, QM calculations are performed on the local machine by using all available CPUs/cores.

(``--queue`` flag set to inline by
default), guessing the number of CPUs to use from the maximum available in the local machine (``--ncpus`` flag).

One may want to use less resources for a given parameterization, or maybe run two parameterizations in parallel on the
same machine. For this, one may override the number of CPUs to be used by:

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --ncpus 1 --outdir local

.. note::

    Benzamidine parametrization can takes 2-4 hours, depending on your machine.

for using 8 CPU nodes. Obviously, the less amount of CPUs used, the slower the parameterization will be.

The QM calculation can be parallized using one of several queue systems (Slurm, LSF, and PBS).

Alternatively, one may have access for remote resources that are ready to run the QM calculations. The ``parameterize``
allows the interface with several execution back-ends, such as LSF, Slurm, and AceCloud. If one has a cluster available
running on Slurm and properly set up, one may use the ``-e`` flag to send the QM calculations to that cluster
automatically:

.. code:: bash

    parameterize benzamidine.mol2 --charge 1 --queue Slurm --outdir slurm

The tool will run locally, which is very computationally inexpensive, and all the computationally expensive QM jobs
will be sent for the user to the Slurm queue system, in this case.

Reuse QM and reparametrize
--------------------------

TODO

``--seed``

Find and validate parameters
----------------------------

TODO
