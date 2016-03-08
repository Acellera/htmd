
System ready for equilibration
==============================

.. code:: python

    from htmd.console import *

.. code:: python

    path=mkdtemp() + '/equilibrate'
    copytree(home()+'/data/equilibrate', path )
    chdir(path)

These steps show how to prepare your system for equilibration.

Check for the acemd protocols available within your HTMD distribution,

.. code:: python

    Acemd.protocols()


.. parsed-literal::

    ['equilibration.py']


.. code:: python

    acemd = Acemd('equilibration.py')
    acemd.show()


.. parsed-literal::

    
    set numsteps 2500000
    set temperature 300
    proc calcforces_init {} {
      berendsenpressure  off}
    proc calcforces {} {	
      global numsteps
      set step [ getstep ]
      if { $step > 500 } {
        berendsenpressure  on
      } else {
        berendsenpressure  off}
      if { $step > [expr $numsteps/2] } {
        constraintscaling 0
      } else {
        constraintscaling [expr 1 + $step*(0.1-1)*2/$numsteps]}
    }
    proc calcforces_endstep { } { }
    1-4scaling	1.0
    cutoff	9
    hydrogenscale	4
    consref	structure.pdb
    tclforces	on
    constraints	on
    langevindamping	1
    switching	on
    minimize	500
    structure	structure.psf
    restartfreq	5000
    temperature	$temperature
    exclude	scaled1-4
    berendsenpressuretarget	1.01325
    constraintscaling	1.0
    pme	on
    switchdist	7.5
    xtcfile	output.xtc
    langevintemp	$temperature
    rigidbonds	all
    fullelectfrequency	2
    berendsenpressure	on
    outputname	output
    langevin	on
    restart	on
    xtcfreq	25000
    energyfreq	1000
    timestep	4
    berendsenpressurerelaxationtime	800
    run	$numsteps
    pmegridspacing	1.0
    parameters	parameters
    coordinates	structure.pdb
    


For use of constraints during the equilibration set occupancy and beta
columns of protein and ligand to 1

.. code:: python

    s = Molecule('structure.pdb')
    s.set('occupancy',0)
    s.set('beta',0)
    s.set('beta',1,sel='segid L and noh')
    s.set('beta',1,sel='segid P and noh')
    s.write('structure.pdb')

Calculate ans set the size of the system periodic box,

.. code:: python

    box = amax(s.get('coords','water'),axis=0)-amin(s.get('coords','water'),axis=0)

.. code:: python

    acemd.celldimension = str(box[0])+' '+str(box[1])+' '+str(box[2])

Provide here the name of the directory you used to build your system

.. code:: python

    acemd.load('./')

Generate the input folder ready for running an equilibration of your
system using ACEMD

.. code:: python

    acemd.save('./equil/')

To start the equilibration of this system, just run ACEMD from the
./equil directory

.. code:: python

    loc = AcemdLocal()
    loc.submit('./equil/')
