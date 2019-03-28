# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.molecule import Molecule
from htmd.apps.acemd import Acemd
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
import os
import numpy as np
import logging
from htmd.decorators import _Deprecated
logger = logging.getLogger(__name__)


@_Deprecated('1.13.6')
class Equilibration(ProtocolInterface):
    """ Equilibration protocol

        Equilibration protocol for globular and membrane proteins
        It includes a flatbottom potential box to retrain a ligand
        for example within this box.

        Parameters
        ----------
        numsteps : int, default=0
            Number of steps to run the simulations in units of 4fs
        temperature : float, default=300
            Temperature of the thermostat in Kelvin
        k : float, default=0
            Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5
        reference : str, default='none'
            Reference selection to use as dynamic center of the flatbottom box.
        selection : str, default='none'
            Selection of atoms to apply the flatbottom potential
        box : list, default=[0, 0, 0, 0, 0, 0]
            Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]
        useconstantratio : bool, default=False
            For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.
        constraints : dict, default={'protein and noh and not name CA': 0.1, 'protein and name CA': 1}
            A dictionary containing as keys the atomselections of the constraints and as values the constraint scaling factor. 0 factor means no constraint, 1 full constraints and in between values are used for scaling. The order with which the constraints are applied is random, so make atomselects mutually exclusive to be sure you get the correct constraints.

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
        super().__init__()
        self._cmdObject('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Acemd)
        self._cmdValue('numsteps', 'int', 'Number of steps to run the simulations in units of 4fs', 0, TYPE_INT, RANGE_0POS)
        self._cmdValue('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, TYPE_FLOAT, RANGE_ANY)
        self._cmdValue('k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0, TYPE_FLOAT, RANGE_ANY)
        self._cmdString('reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.', 'none')
        self._cmdString('selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none')
        self._cmdList('box', 'list', 'Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]', [0,0,0,0,0,0])
        self._cmdBoolean('useconstantratio', 'bool', 'For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.', False)
        self._cmdDict('constraints', 'dict', 'A dictionary containing as keys the atomselections of the constraints '
                                             'and as values the constraint scaling factor. 0 factor means no constraint'
                                             ', 1 full constraints and in between values are used for scaling.'
                                             ' The order with which the constraints are applied is random, so make '
                                             'atomselects mutually exclusive to be sure you get the correct constraints.'
                                             , {'protein and noh and not name CA': 0.1, 'protein and name CA': 1})

        self.acemd = Acemd()
        self.acemd.coordinates = None
        self.acemd.structure = None
        self.acemd.parameters = None
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
        self.acemd.consref = None
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
proc calcforces_terminate { } { }
'''

    def _findFiles(self, inputdir):
        # Tries to find default files if the given don't exist
        defaults = {'coordinates': ('structure.pdb', ),
                    'structure': ('structure.psf', 'structure.prmtop'),
                    'parameters': ('parameters', 'structure.prmtop')}

        for field in defaults:
            userval = self.acemd.__dict__[field]
            if userval is not None and not os.path.exists(os.path.join(inputdir, userval)):
                self.acemd.__dict__[field] = None

            if self.acemd.__dict__[field] is None:
                for val in defaults[field]:
                    if os.path.exists(os.path.join(inputdir, val)):
                        self.acemd.__dict__[field] = val
                        break

            if userval is not None and self.acemd.__dict__[field] is not None and self.acemd.__dict__[field] != userval:
                logger.warning('Could not locate structure file {}. Using {} instead.'.format(
                    os.path.join(inputdir, userval), os.path.join(inputdir, self.acemd.__dict__[field])
                ))
            elif self.acemd.__dict__[field] is None:
                raise RuntimeError('Could not locate any {f:} file in {i:}. '
                                   'Please set the Equilibration.acemd.{f:} property to '
                                   'point to the {f:} file'.format(f=field, i=inputdir))

        if self.acemd.consref is None:
            self.acemd.consref = self.acemd.coordinates

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.parameters.endswith('structure.prmtop'):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None
            self.acemd.scaling14 = '0.8333333'
            self.acemd.amber = 'on'

    def write(self, inputdir, outputdir):
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
        self._findFiles(inputdir)
        self._amberFixes()

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        self.acemd.TCL = self.acemd.TCL.replace('NUMSTEPS', str(self.numsteps))
        self.acemd.TCL = self.acemd.TCL.replace('TEMPERATURE', str(self.temperature))
        self.acemd.TCL = self.acemd.TCL.replace('KCONST', str(self.k))
        self.acemd.TCL = self.acemd.TCL.replace('REFINDEX', ' '.join(map(str, inmol.get('index', self.reference))))
        self.acemd.TCL = self.acemd.TCL.replace('SELINDEX', ' '.join(map(str, inmol.get('index', self.selection))))
        self.acemd.TCL = self.acemd.TCL.replace('BOX', ' '.join(map(str, self.box)))
        if self.acemd.celldimension is None and self.acemd.extendedsystem is None:
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
        inmol.set('occupancy', 0)
        inmol.set('beta', 0)
        for sel in self.constraints:
            inmol.set('beta', self.constraints[sel], sel)
        outfile = os.path.join(outputdir, self.acemd.coordinates)
        inmol.write(outfile)

    def addConstraint(self, atomselect, factor=1):
        """ Convenience function for adding a new constraint to existing constraints.

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
    import htmd.home
    eq = Equilibration()
    eq.numsteps = 1000000
    eq.temperature = 300
    eq.reference = 'protein and name CA'
    eq.selection = 'segname L and noh'
    eq.box = [-20, 20, -20, 20, 43, 45]
    eq.k = 5
    eq.write(htmd.home.home() + '/data/equilibrate', '/tmp/equil1')
