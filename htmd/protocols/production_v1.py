# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from moleculekit.molecule import Molecule
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
from htmd.apps.acemd import Acemd
import os
import logging
from htmd.decorators import _Deprecated
logger = logging.getLogger(__name__)


@_Deprecated('1.13.6')
class Production(ProtocolInterface):
    ''' Production protocol for globular and membrane proteins
        It also includes a possible flatbottom potential box

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
    '''
    def __init__(self):
        super().__init__()
        self._cmdObject('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Acemd)
        self._cmdValue('numsteps', 'int', 'Number of steps to run the simulations in units of 4fs', 0, TYPE_INT, RANGE_0POS)
        self._cmdValue('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, TYPE_FLOAT, RANGE_ANY)
        self._cmdValue('k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0, TYPE_FLOAT, RANGE_ANY)
        self._cmdString('reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.', 'none')
        self._cmdString('selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none')
        self._cmdList('box', 'list', 'Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]', [0,0,0,0,0,0])

        self.acemd = Acemd()
        #self.acemd.binindex='input.idx'
        self.acemd.extendedsystem='input.xsc'
        self.acemd.coordinates = None
        self.acemd.structure = None
        self.acemd.parameters = None
        self.acemd.temperature = '300'
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
        self.acemd.langevintemp = '300'
        self.acemd.langevindamping = '0.1'
        self.acemd.pme = 'on'
        self.acemd.pmegridspacing = '1.0'
        self.acemd.fullelectfrequency = '2'
        self.acemd.energyfreq = '5000'
        self.acemd.run = '$numsteps'
        self._TCL='''
set numsteps NUMSTEPS
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
  global ref sel K box
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
}
proc calcforces_endstep { } { }
proc calcforces_terminate { } { }
'''

    def _findFiles(self, inputdir):
        # Tries to find default files if the given don't exist
        defaults = {'coordinates': ('structure.pdb',),
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

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.parameters.endswith('structure.prmtop'):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None
            self.acemd.scaling14 = '0.8333333'
            self.acemd.amber = 'on'

    def write(self, inputdir, outputdir):
        """ Writes the production protocol and files into a folder.

        Parameters
        ----------
        inputdir : str
            Path to a directory containing the files produced by a equilibration process.
        outputdir : str
            Directory where to write the production setup files.
        """
        self._findFiles(inputdir)
        self._amberFixes()

        self.acemd.temperature = str(self.temperature)
        self.acemd.langevintemp = str(self.temperature)
        if self.k > 0: #use TCL only for flatbottom
            mol = Molecule(os.path.join(inputdir, self.acemd.coordinates))
            self.acemd.tclforces = 'on'
            TCL = self._TCL
            TCL = TCL.replace('NUMSTEPS', str(self.numsteps))
            TCL = TCL.replace('KCONST', str(self.k))
            TCL = TCL.replace('REFINDEX', ' '.join(map(str, mol.get('index', self.reference))))
            TCL = TCL.replace('SELINDEX', ' '.join(map(str, mol.get('index', self.selection))))
            TCL = TCL.replace('BOX', ' '.join(map(str, self.box)))
            self.acemd.TCL = TCL
        else:
            self.acemd.TCL = 'set numsteps {}\n'.format(self.numsteps)
        self.acemd.setup(inputdir, outputdir, overwrite=True)

if __name__ == "__main__":
    import htmd.home
    from htmd.util import tempname
    md = Production()
    md.temperature = 300
    md.reference = 'protein and name CA'
    md.selection = 'segname L and noh'
    md.acemd.extendedsystem = None  # use different data
    md.acemd.binindex = None  # use different data
    md.box = [-20, 20, -20, 20, 43, 45]
    md.k = 5
    md.write(htmd.home.home() +'/data/equilibrate', tempname())
    md.k = 0
    md.write(htmd.home.home() +'/data/equilibrate', tempname())
