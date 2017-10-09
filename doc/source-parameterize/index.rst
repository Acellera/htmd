Parameterize
============

``parameterize`` is a tool to obtain force field parameters for new molecules easily and efficiently.

Commonly used AMBER and CHARMM force fields contain parameters for biomolecules (proteins, nucleotides, saccharides,
lipids, etc), but lack parameters for other biologically relevant molecules (co-factors, drugs, etc).

GAFF and CGenFF address this by using empirical rules and pre-computed data sets to derive parameters for arbitrary
organic molecules. However, these parameters are not guaranteed to be transferable to all possible chemical
environments.

``parameterize`` improves the quality of the parameters by using QM data, i.e. refitting ESP charges and rotatable
dihedral angle parameters. It fundamentally solves the problem of transferability.

Contents:

.. toctree::
    :maxdepth: 1

    Method <method>
    Tutorial <tutorial>

Capabilities
------------

* Supported force fields:
    * AMBER
    * CHARMM
* Atom typing and initial parameters:
    * GAFF and GAFF2 (for AMBER)
    * CGenFF (for CHARMM)
* Supported QM codes:
    * Psi4
    * Gaussian (experimental)
* Atomic charges:
    * Fitted to reproduce the electrostatic potential (ESP) of QM
    * Selected charges can be fixed during the fitting
* Dihedral angle parameters:
    * Automatic rotatable dihedral detection
    * Automatic identical dihedral detection
    * Dihedral scanning with and without QM structure optimization
    * Parameter fitting to reproduce QM energies
* QM job execution:
    * Local machine
    * LSF/PBS/Slum queuing systems
    * AceCloud (experimental)

.. comment:
``parameterize`` comes with HTMD. It takes a MOL2 file, parameterizes it, and outputs force field files in CHARMM and
AMBER formats. There are many options that are shown in detail below. The structure inside the output directory
 (controlled by the ``--output`` flag) is the following::
    .
    ├── dihedral-opt
    ├── esp
    ├── minimize
    └── parameters

.. comment:
The output inside ``dihedral-opt/``, ``esp/``, and ``minimize/`` is related with the QM calculations performed during
the parameterization of the molecule and work as check-points for the different steps of the parameterization process.
The ``parameters`` is where the relevant outputs are written with following format:
``parameters/<force-field>/<theory>-<basis-set>-<vacuum/water>/``. Inside this directory there are the force-field
parameter files, as well as a ``plots/`` directory, where one can check the optimization carried on by the tool.