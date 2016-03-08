Ligand Parameterisation
=======================


Introduction
------------

The parameterise package produces optimised force field parameters for small modules according to the `GAAMP <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3819940/>`__ method.

Starting from an initial geometry it will optimise with respect to quantum mechanical calculations:

 - atomic partial charges
 - bond and angle reference lengths
 - torsional (rotameric) potentials

Installation
------------

The ligand parameterisation tool is installed using the command:

.. parsed-literal::
  
    conda install parameterisation -y

Parameterise will perform all QC calculations with Gaussian if the command g09 is found at runtime. If Gaussian isn't available, Parameterise will fall back to using its internal copy of the QC code Psi4.

Parameterisation
----------------

The Parameterise tool requires two input files: Firstly, a ligand starting geometry in mol2 or PDB format. Secondly, an input file containing, at a minimum:

.. parsed-literal::
    JobName    [some-name]
    FileName   [ligand.pdb|mol2]
    NetCharge  [charge]

Where "some-name" is a unique identifier for referencing the parameterisation.
Additional input file commands are available, and can be seen with "parameterise --command", but are generally won't need to be changed from their default values.

Run the parameterisation with the command:

.. parsed-literal::
    parameterise --input [input-file-name]

The output will look like this:

.. parsed-literal::
some stuff

Controlling QC execution
------------------------

The amount of memory and CPUs used by the QC calculations will be automatically determined from the specification of the machine running parameterise. These settings can be override with the input file commands

.. parsed-literal::
    NCORES   [num-cpus]
    MEMORY   [memory-in-GB]


The default behaviour when running the QC calculations is to run them sequentially directly on the machine on which "parameterise" was invoked. Instead, the jobs can be submitted to a batch system using the command:

.. parsed-literal::
    ExectionMode  Inline|PBS|LSF

Note: 
 *  this assumes a common file system between interactive session and cluster exection resources. Only PBS Pro and LSF resource management systems are supported at the moment. 
 * NCORES and MEMORY settings will be used to qualify the resource specification of the submitted jobs.
 * Presently the water-fitting step always runs its QC calculations inline. This means the machine running "parameterise" must have sufficient free resources for the QC.


Output
------

A successful parameterisation will produce the following output files:

 * mol.pdb -- A QC-minimised geometry
 * mol.rtf, mol.prm -- CHARMM format force-field parameters 
 * mol-dihedral-N.svg -- A plot of the potentials of each fitted torsion including orginal MM potential, QC potential and fitted MM potential

Parameterise from Python
------------------------

It's also possible to invoke a parameterisation directly from Python code:

.. parsed-literal::
 from parameterise import Configuration, Parameterisation

.. parsed-literal::
 # Set up a Configuration object.
 # The fields match the commands available in the input file

 config = Configuration()
 config.FileName = "ligand.pdb"
 config.JobName  = "parameterise-test"
 config.NCORES   = 4
 config.MEMORY   = 4


.. parsed-literal::
 # Start a parameterisation

 param = Parameterisation( config = config )


.. parsed-literal::
 # Plot each of the dihedral fits
 param.plotDihedrals()

.. parsed-literal::
 # Plot dihedrals to files
 fileset = param.plotDihedrals( show=False)

.. parsed-literal::
 # Returns a dictionary of temporary files 
 # with keys 'RTF', 'PRM' and 'PDB'

 fileset=param.getParameters()

