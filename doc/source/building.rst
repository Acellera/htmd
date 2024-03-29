Building
========

Building is the most important module for preparing Molecules for simulation. HTMD includes tools to solvate molecular
systems in water, prepare proteins and nucleic acids for simulation by assigning appropriate protonation states, and can build all the
necessary files to simulate the system in either two of the most widely used force-fields, CHARMM and AMBER.

The workflow of the building process is normally:

#. Obtain structures
#. Clean structures
#. Define segments
#. Combine structures
#. Solvate
#. Choose force-field
#. Build and ionize

Contents:

.. toctree::
    :maxdepth: 1

    System preparation <../moleculekit/moleculekit.tools.preparation.rst>
    Solvating <htmd.builder.solvate>
    CHARMM builder <htmd.builder.charmm>
    AMBER builder <htmd.builder.amber>
