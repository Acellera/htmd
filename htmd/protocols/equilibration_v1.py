# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.molecule.molecule import Molecule
from htmd.userinterface import UserInterface
from htmd.acemd.acemd import Acemd
import os
import numpy as np


class Equilibration(UserInterface):
    """ Equilibration protocol

        Equilibration protocol for globular and membrane proteins
        It includes a flatbottom potential box to retrain a ligand
        for example within this box.

        Parameters
        ----------
        numsteps: int
            Number of steps to run the simulations in units of 4fs
        temperature: float
            Temperature of the thermostat in Kelvin
        useconstantratio: bool
            For membrane protein simulations set it to true so that the barostat
            does not modify the xy aspect ratio. Default 0
        k: float
            Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5
        reference: str
            Reference selection to use as dynamic center of the flatbottom box.
        selection: str
            Selection of atoms to apply the flatbottom potential
        box: str
            Position of the flatbottom box in term of the reference center given as
            [xmin, xmax, ymin, ymax, zmin, zmax]
        inputdir: str
            Input directory where to find the simulation files
        outputdir: str
            Output directory to create with the simulation ready to run

        Example
        -------
        >>> from htmd.protocols.equilibration_v1 import Equilibration
        >>> md = Equilibration()
        >>> md.numsteps = 10000000
        >>> md.temperature = 300
        >>> md.useconstantratio = True  # only for membrane sims
        >>> # this is only needed for setting the flatbottom potential, otherwise remove it
        >>> md.reference = 'protein and resid 293'
        >>> md.selection = 'segname L and noh'
        >>> md.box = [-25, 25, -25, 25, 43, 45]
        >>> md.k = 5
        >>> md.write('./build','./equil')
    """
    def __init__(self):
        self._commands = {'acemd': None, 'numsteps': 0, 'temperature': 300, 'k': 0, 'reference': 'none',
                          'selection': 'none', 'box': [0,0,0,0,0,0], 'useconstantratio': False, 'constraints': None}
        for k in self._commands:
            self.__dict__[k] = self._commands[k]

        self.constraints = {'protein and noh and not name CA': 0.1, 'protein and name CA': 1}

        self.acemd = Acemd()
        self.acemd.coordinates = 'structure.pdb'
        self.acemd.structure = 'structure.psf'
        self.acemd.parameters = 'parameters'
        self.acemd.temperature = '$temperature'
        self.acemd.restart = 'on'
        self.acemd.restartfreq = '5000'
        self.acemd.outputname = 'output'
        self.acemd.xtcfile = 'output.xtc'
        self.acemd.xtcfreq = '25000'
        self.acemd.timestep = '4'
        self.acemd.rigidbonds = 'all'
        self.acemd.hydrogenscale = '4'
        self.acemd.switching = 'on'
        self.acemd.switchdist = '7.5'
        self.acemd.cutoff = '9'
        self.acemd.exclude = 'scaled1-4'
        self.acemd.scaling14 = '1.0'
        self.acemd.langevin = 'on'
        self.acemd.langevintemp = '$temperature'
        self.acemd.langevindamping = '1'
        self.acemd.pme = 'on'
        self.acemd.pmegridspacing = '1.0'
        self.acemd.fullelectfrequency = '2'
        self.acemd.energyfreq = '1000'
        self.acemd.constraints = 'on'
        self.acemd.consref = 'structure.pdb'
        self.acemd.constraintscaling = '1.0'
        self.acemd.berendsenpressure = 'on'
        self.acemd.berendsenpressuretarget = '1.01325'
        self.acemd.berendsenpressurerelaxationtime = '800'
        self.acemd.tclforces = 'on'
        self.acemd.minimize = '500'
        self.acemd.run = '$numsteps'
        self.acemd.TCL='''
set numsteps NUMSTEPS
set temperature TEMPERATURE
set refindex { REFINDEX }
set selindex { SELINDEX }
set box { BOX }
set K KCONST
#
proc flatbot1d {x xm xM K} {
  set f 0
  if {$x < $xm} {
    set f [expr $K*[expr $xm-$x]]
  }
  if {$x > $xM} {
    set f [expr $K*[expr $xM-$x]]
  }
  return $f
}
proc calcforces_init {} {
  global ref sel refindex selindex
  berendsenpressure  off
  set ref [addgroup  $refindex]
  set sel [addgroup  $selindex]
}
proc calcforces {} {
  global ref sel numsteps K box 
  loadcoords coords
##FLATBOTTOM
  if {$K>0} {
    set r0 $coords($ref)
    set r1 $coords($sel)
    set dr  [vecsub $r1 $r0]
    set fx [flatbot1d [lindex $dr 0] [lindex $box 0] [lindex $box 1] $K]
    set fy [flatbot1d [lindex $dr 1] [lindex $box 2] [lindex $box 3] $K]
    set fz [flatbot1d [lindex $dr 2] [lindex $box 4] [lindex $box 5] $K]
    #print "dr: $dr  fx: $fx fy: $fy fz: $fz"
    addforce $sel [list $fx $fy $fz]
  }
##EQUIL
  set step [ getstep ]
  if { $step > 500 } {
    berendsenpressure  on
  } else {
    berendsenpressure  off}
  if { $step > [expr $numsteps/2] } {
    constraintscaling 0
  } else {
    constraintscaling [expr 1 + $step*(0.05-1)*2/$numsteps]}
}
proc calcforces_endstep { } { }
'''

    def write(self, inputdir=None, outputdir=None):
        """ Write the equilibration protocol

        Writes the equilibration protocol and files into a folder for execution
        using files inside the inputdir directory

        Parameters
        ----------
        inputdir : str
            Path to a directory containing the files produced by a build process.
        outputdir : str
            Directory where to write the equilibration setup files.

        Examples
        --------
        >>> md = Equilibration()
        >>> md.write('./build','./equil')
        """
        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        self.acemd.TCL = self.acemd.TCL.replace('NUMSTEPS', str(self.numsteps))
        self.acemd.TCL = self.acemd.TCL.replace('TEMPERATURE', str(self.temperature))
        self.acemd.TCL = self.acemd.TCL.replace('KCONST', str(self.k))
        self.acemd.TCL = self.acemd.TCL.replace('REFINDEX', ' '.join(map(str, inmol.get('index', self.reference))))
        self.acemd.TCL = self.acemd.TCL.replace('SELINDEX', ' '.join(map(str, inmol.get('index', self.selection))))
        self.acemd.TCL = self.acemd.TCL.replace('BOX', ' '.join(map(str, self.box)))
        if 'celldimension' not in self.acemd.__dict__ and 'extendedsystem' not in self.acemd.__dict__:
            coords = inmol.get('coords', sel='water')
            if coords.size == 0:  # It's a vacuum simulation
                coords = inmol.get('coords', sel='all')
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
                dim = dim + 12. 
            else:
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
            self.acemd.celldimension = '{} {} {}'.format(dim[0], dim[1], dim[2])
        if self.useconstantratio:
            self.acemd.useconstantratio = 'on'
        self.acemd.setup(inputdir, outputdir, overwrite=True)

        # Adding constraints
        pdbfile = os.path.join(outputdir, self.acemd.coordinates)
        mol = Molecule(pdbfile)
        mol.set('occupancy', 0)
        mol.set('beta', 0)
        for sel in self.constraints:
            mol.set('beta', self.constraints[sel], sel)
        mol.write(pdbfile)

    def getConstraints(self):
        """ Get the currently applied constraints

        Returns
        -------
        constraints : dict
            A dictionary containing as keys the atomselections of the constraints and as values the constraint scaling
            factor. 0 scaling means no constraint, 1 full constraints and in between values are used for scaling.
        """
        return self.constraints

    def setConstraints(self, constraints):
        """ Set equilibration constraints

        Overwrites all existing constraints with your own dictionary of constraints.

        Parameters
        ----------
        constraints : dict
            A dictionary containing as keys the atomselections of the constraints and as values the constraint scaling
            factor. 0 scaling means no constraint, 1 full constraints and in between values are used for scaling.

        Example
        -------
        >>> # The order with which the constraints are written is random.
        >>> # So if you do:
        >>> eq.setConstraints({'protein and noh': 0.1, 'protein and name CA': 1})
        >>> # The CA's might get factor 0.1 because of the first constraint being applied after the second.
        >>> # So make atomselects mutually exclusive to be sure you get the correct constraints:
        >>> eq.setConstraints({'protein and noh and not name CA': 0.1, 'protein and name CA': 1})
        """
        self.constraints = constraints

    def addConstraint(self, atomselect, factor=1):
        """ Add a new constraint to existing constraints

        Parameters
        ----------
        atomselect : str
            Atom selection of atoms we want to constrain
        factor : float
            The scaling factor of the constraints applied to the atoms

        Example
        -------
        >>> eq.addConstraint('chain X', 0.3)
        """
        self.constraints[atomselect] = factor

if __name__ == "__main__":
    import htmd
    eq = Equilibration()
    eq.numsteps = 1000000
    eq.temperature = 300
    eq.reference = 'protein and name CA'
    eq.selection = 'segname L and noh'
    eq.box = [-20, 20, -20, 20, 43, 45]
    eq.k = 5
    eq.write(htmd.home() + '/data/equilibrate', '/tmp/equil1')
