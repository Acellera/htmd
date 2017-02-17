Parameterization
================

Parameterization of small molecules is an important aspect of molecular simulations, especially when dealing with
molecules which are not readily available by default on CHARMM or AMBER. Even when parameters can be obtained from
other sources such as CGENFF or GAFF, they are not guaranteed to be fully optimized.

The ``parameterize`` command-line tool comes with HTMD. It takes a MOL2 file, parameterizes it, and outputs force-field
files in common CHARMM and AMBER formats. There are many options that are shown in detail below. The structure inside
the output directory (controlled by the ``-o`` flag) is the following::

    .
    ├── dihedral-opt
    ├── esp
    ├── minimize
    └── parameters

The output inside ``dihedral-opt/``, ``esp/``, and ``minimize/`` is related with the QM calculations performed during
the parameterization of the molecule and work as check-points for the different steps of the parameterization process.
The ``parameters`` is where the relevant outputs are written with following format:
``parameters/<force-field>/<theory>-<basis-set>-<vacuum/water>/``. Inside this directory there are the force-field
parameter files, as well as a ``plots/`` directory, where one can check the optimization carried on by the tool.

The method followed by the ``parameterize`` tool is inspired by the GAAMP method [#]_. A detailed description of our
method will soon be published [#]_.

.. [#]  L. Huang and B. Roux, Automated Force Field Parameterization for Nonpolarizable and Polarizable
        Atomic Models Based on Ab Initio Target Data, J. Chem. Theory Comput., 2013, 9 (8), pp 3543–3556.
        doi: `10.1021/ct4003477 <http://dx.doi.org/10.1021/ct4003477>`_
.. [#]  M. J. Harvey, J. M. Damas, G. Martínez and G. De Fabritiis, Small Molecule Forcefield Parameterization with
        HTMD, (manuscript).

------------

.. argparse::
    :ref: htmd.parameterization.cli.cli_parser
    :prog: parameterize