Parameterization of Benzamidine
===============================

Let's imagine we want to work in the binding of benzamidine to trypsin. Most force-fields will not have parameters
available for benzamidine, and it will need to be parameterized in order to be simulated.

HTMD comes with a ``parameterize`` tool, a command-line interface to generate parameters for small molecules and output
them in CHARMM or AMBER formats.

Get Benzamidine
---------------

First of all, we would need to get benzamidine. An ideal structure in SDF format can be taken from the PDB, and then
babel can be used to convert it to MOL2 format

.. code:: bash

    wget https://files.rcsb.org/ligands/view/BEN_ideal.sdf
    babel -isdf BEN_ideal.sdf -omol2 ben.mol2

These commands will get you a ``ben.mol2`` file ready to be parameterized.

Get Help for parameterize
-------------------------

All options available for ``parameterize`` can be shown using:

.. code:: bash

    parameterize -h

::

    Please cite HTMD: Doerr et al.(2016)JCTC,12,1845. https://dx.doi.org/10.1021/acs.jctc.6b00049

    HTMD Documentation at: https://www.htmd.org/docs/latest/

    You are on the latest HTMD version (unpackaged : /home/joao/maindisk/software/repos/Acellera/htmd/htmd).

    usage: parameterize [-h] -m <input.mol2> [-l] [-c CHARGE] [--rtf RTF]
                        [--prm PRM] [-o OUTDIR] [-t A1-A2-A3-A4] [-n NCPUS]
                        [-f {GAFF,GAFF2,CGENFF,all}] [-b {6-31g-star,cc-pVDZ}]
                        [--theory {RHF,B3LYP}] [--vacuum] [--no-min] [--no-esp]
                        [--no-torsions] [-e {inline,LSF,PBS,Slurm,AceCloud}]
                        [--qmcode {Gaussian,PSI4,TeraChem}] [--freeze-charge A1]
                        [--no-geomopt]

    Acellera Small Molecule Parameterization Tool

    optional arguments:
      -h, --help            show this help message and exit
      -m <input.mol2>, --mol2 <input.mol2>
                            Molecule to parameterise, in mol2 format
      -l, --list, --list-torsions
                            List parameterisable torsions
      -c CHARGE, --charge CHARGE
                            Net charge on molecule (default: sum of the partial
                            charges on the .mol2 file)
      --rtf RTF             Inital RTF parameters (req --prm)
      --prm PRM             Inital PRM parameters (req --rtf)
      -o OUTDIR, --outdir OUTDIR
                            Output directory (default: ./)
      -t A1-A2-A3-A4, --torsion A1-A2-A3-A4
                            Torsion to parameterise (default: all)
      -n NCPUS, --ncpus NCPUS
                            Number of CPUs to use (default: 12)
      -f {GAFF,GAFF2,CGENFF,all}, --forcefield {GAFF,GAFF2,CGENFF,all}
                            Inital FF guess to use (default: all)
      -b {6-31g-star,cc-pVDZ}, --basis {6-31g-star,cc-pVDZ}
                            QM Basis Set (default: cc-pVDZ)
      --theory {RHF,B3LYP}  QM Theory (default: B3LYP)
      --vacuum              Perform QM calculations in vacuum (default: False)
      --no-min              Do not perform QM minimisation (default: False)
      --no-esp              Do not perform QM charge fitting (default: False)
      --no-torsions         Do not perform torsion fitting (default: False)
      -e {inline,LSF,PBS,Slurm,AceCloud}, --exec {inline,LSF,PBS,Slurm,AceCloud}
                            Mode of execution for the QM calculations (default:
                            inline)
      --qmcode {Gaussian,PSI4,TeraChem}
                            QM code (default: PSI4)
      --freeze-charge A1    Freeze the charge of the named atom (default: None)
      --no-geomopt          Do not perform QM geometry optimisation when fitting
                            torsions (default: True)

This tutorial will aim to showcase some of the most important features of the tool using benzamidine.

Get the CGENFF/GAFF parameters for benzamidine
----------------------------------------------

The ``parameterize`` tool has an option to use a given FF as initial guess for the parameters, set by the flag ``-f`` or
``--forcefield``. By default, it outputs parameters for both CHARMM (through CGENFF) and AMBER (through GAFF2). If one
wants parameters in the original GAFF, one may use ``-f GAFF``.

The most simple use of ``parameterize``, which does not go through any optimization procedure, is to provide the guessed
initial parameters. This can be done by setting the flags ``--no-min``, ``--no-esp`` and ``--no-torsions``:

.. code:: bash

    parameterize -m ben.mol2 --no-min --no-esp --no-torsions -o ./param-ben-noqm

Inside the output directory ``./param-ben-noqm`` (specified in the flag ``-o``), one can find parameters for each both
CHARMM and AMBER.

List torsions of benzamidine
----------------------------

Before doing any parameterization, one can list the soft torsions that the molecule has. This can be easily done by
using the ``-l`` flag:

.. code:: bash

    parameterize -m ben.mol2 -l

which for benzamidine should output:

::

    C2-C1-C7-N1
    C1-C7-N1-H6
    C1-C7-N2-H7

that are 3 sets of 4 atoms describing the 3 detected soft torsions.

Setting up the charge of benzamidine
------------------------------------

Verifying if the small molecule has the correct protonation state and charge is still a user-dependent task. Some
combinations of a total charge with a given set of atoms may cause the underlying QM software to error due to
incompatibility of the multiplicity with the total amount of electrons of the system. For example, when taking our
``ben.mol2``, and leaving the charge as the default one, deduced from the MOL2 charges (a total charge of +1), the tool
errors as it should:

.. code:: bash

    parameterize -m ben.mol2 --no-esp --no-torsions -o ./param-ben-fail

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

    parameterize -m ben.mol2 -c 0 --no-esp --no-torsions -o ./param-ben-pass

This last command just performs the minization of the small molecule and outputs that minimized structure with the
CGENFF/GAFF2 parameters. Please note that a combination of total charge and a given set of atoms working on
``parameterize`` does not mean that protonation state is the most common or relevant in normal pH conditions.

Default, most accurate parameterization of benzamidine
------------------------------------------------------

By default, ``parameterize`` uses settings that account for the most accurate QM settings available in the tool. This
means using a higher level of theory (B3LYP, can be set through flag ``--theory``), a larger basis set (cc-pVDZ, can be
set through the flag ``-b`` or ``--basis``), and a solvation model for the QM calculations (turned on by default, can be
turned off by using the flag ``--vacuum``). By keeping all of the settings on default, the parameterization of the small
molecule will, however, be more computationally demanding (also depending on its number of atoms, number of soft
torsions, and on the resources available for the parameterization). One can do one of these parameterization using:

.. code:: bash

    parameterize -m ben.mol2 -c 0 -o ./param-ben-fullqm


Get faster but less accurate parameterizations of benzamidine
-------------------------------------------------------------

One may one to sacrifice on accuracy to increase the speed of parameterization. This can be done by using the flags
related with accuracy mentioned in the last section.

One may start by turning off the solvation model, and doing the calculations in vacuum. This can be one of the first
options if one wants to increase the speed of parameterization:

.. code:: bash

    parameterize -m ben.mol2 -c 0 --vacuum -o ./param-ben-vacuum

In conjunction with this option (or alone, if one wants to keep the solvation model), one may also choose a smaller
basis set, 6-31g-star:

.. code:: bash

    parameterize -m ben.mol2 -c 0 --vacuum -b 6-31g-star -o ./param-ben-vacuum-lowbasis

Finally, one may parameterize with the least accuracy possible, but at the highest speed, by lowering the level of theory
to RHF:

.. code:: bash

    parameterize -m ben.mol2 -c 0 --vacuum -b 6-31g-star --theory RHF -o ./param-ben-fast

Note that the output directories are changed for simplicity. In fact, when changing these three options, the outputs are
always automatically sent to differently named directories, with the format ``<theory>-<basis-set>-<vacuum/water>``.

Control the execution of the parameterization
---------------------------------------------

By default, ``parameterize`` uses an open-source QM code, PSI4 for doing the QM calculations. If one has access to
Gaussian, ``parameterize`` supports it and one may change the QM code using the flag ``--qmcode``:

.. code:: bash

    parameterize -m ben.mol2 -c 0 --qmcode Gaussian -o ./param-ben-gaussian

Furthermore, by default, the tool runs the QM calculations on the local machine by (``-e`` flag set to inline by
default), guessing the number of CPUs to use from the maximum available in the local machine (``-n`` flag).

One may want to use less resources for a given parameterization, or maybe run two parameterizations in parallel on the
same machine. For this, one may override the number of CPUs to be used by:

.. code:: bash

    parameterize -m ben.mol2 -c 0 -n 8 -o ./param-ben-fullqm

for using 8 CPU nodes. Obviously, the less amount of CPUs used, the slower the parameterization will be.

Alternatively, one may have access for remote resources that are ready to run the QM calculations. The ``parameterize``
allows the interface with several execution back-ends, such as LSF, Slurm, and AceCloud. If one has a cluster available
running on Slurm and properly set up, one may use the ``-e`` flag to send the QM calculations to that cluster
automatically:

.. code:: bash

    parameterize -m ben.mol2 -c 0 -e Slurm -o ./param-ben-fullqm

The tool will run locally, which is very computationally inexpensive, and all the computationally expensive QM jobs
will be sent for the user to the Slurm queue system, in this case.