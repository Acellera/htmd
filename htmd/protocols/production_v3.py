# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from moleculekit.molecule import Molecule
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
from htmd.apps.acemd import Acemd
import os
import htmd
import logging
from htmd.decorators import _Deprecated
logger = logging.getLogger(__name__)


@_Deprecated('1.13.6')
class Production(ProtocolInterface):
    ''' Production protocol v3

        Production protocol for globular and membrane proteins
        It also includes a possible flatbottom potential box
        It is also possible to define constraints for the production run

        Parameters
        ----------
        runtime : float, default=0
            Running time of the simulation.
        timeunits : str, default='steps'
            Units for runtime. Can be 'steps', 'ns' etc.
        temperature : float, default=300
            Temperature of the thermostat in Kelvin
        fb_k : float, default=0
            Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5
        fb_reference : str, default='none'
            Reference selection to use as dynamic center of the flatbottom box.
        fb_selection : str, default='none'
            Selection of atoms to apply the flatbottom potential
        fb_box : list, default=[0, 0, 0, 0, 0, 0]
            Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]
        useconstraints: bool, default=False
            Apply constraints to the production simulation, defined by the constraints parameter
        constraints : dict, default={}
            A dictionary of atomselections and values of the constraint to be applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.

    '''
    def __init__(self):
        super().__init__()
        self._cmdObject('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Acemd)
        self._cmdValue('runtime', 'float', 'Running time of the simulation.', 0, TYPE_FLOAT, RANGE_0POS)
        self._cmdString('timeunits', 'str', 'Units for runtime. Can be \'steps\', \'ns\' etc.', 'steps')
        self._cmdValue('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, TYPE_FLOAT, RANGE_ANY)
        self._cmdValue('fb_k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0, TYPE_FLOAT, RANGE_ANY)
        self._cmdString('fb_reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.', 'none')
        self._cmdString('fb_selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none')
        self._cmdList('fb_box', 'list', 'Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]', [0,0,0,0,0,0])
        self._cmdList('useconstraints', 'bool', 'Apply constraints to the production simulation, defined by the constraints parameter', False)
        self._cmdDict('constraints', 'dict', 'A dictionary of atomselections and values of the constraint to be applied '
                                             '(in kcal/mol/A^2). Atomselects must be mutually exclusive.'
                                             , {})

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
        self.acemd.TCL = ('''
set numsteps {NUMSTEPS}
set fb_refindex {{ {REFINDEX} }}
set fb_selindex {{ {SELINDEX} }}
set fb_box {{ {BOX} }}
set fb_K {KCONST}
#
''',
'''
proc flatbot1d {x xm xM fb_K} {
  set f 0
  if {$x < $xm} {
    set f [expr $fb_K*[expr $xm-$x]]
  }
  if {$x > $xM} {
    set f [expr $fb_K*[expr $xM-$x]]
  }
  return $f
}
proc calcforces_init {} {
  global ref sel fb_refindex fb_selindex
  berendsenpressure  off
  set ref [addgroup  $fb_refindex]
  set sel [addgroup  $fb_selindex]
}
proc calcforces {} {
  global ref sel fb_K fb_box
  loadcoords coords
##FLATBOTTOM
  if {$fb_K>0} {
    set r0 $coords($ref)
    set r1 $coords($sel)
    set dr  [vecsub $r1 $r0]
    set fx [flatbot1d [lindex $dr 0] [lindex $fb_box 0] [lindex $fb_box 1] $fb_K]
    set fy [flatbot1d [lindex $dr 1] [lindex $fb_box 2] [lindex $fb_box 3] $fb_K]
    set fz [flatbot1d [lindex $dr 2] [lindex $fb_box 4] [lindex $fb_box 5] $fb_K]
    #print "dr: $dr  fx: $fx fy: $fy fz: $fz"
    addforce $sel [list $fx $fy $fz]
  }
}
proc calcforces_endstep { } { }
''')

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

        from htmd.units import convert
        numsteps = convert(self.timeunits, 'timesteps', self.runtime, timestep=self.acemd.timestep)
        self.acemd.temperature = str(self.temperature)
        self.acemd.langevintemp = str(self.temperature)
        if self.fb_k > 0: #use TCL only for flatbottom
            mol = Molecule(os.path.join(inputdir, self.acemd.coordinates))
            self.acemd.tclforces = 'on'
            tcl = list(self.acemd.TCL)
            tcl[0] = tcl[0].format(NUMSTEPS=numsteps, KCONST=self.fb_k,
                                   REFINDEX=' '.join(map(str, mol.get('index', self.fb_reference))),
                                   SELINDEX=' '.join(map(str, mol.get('index', self.fb_selection))),
                                   BOX=' '.join(map(str, self.fb_box)))
            self.acemd.TCL = tcl[0] + tcl[1]
        else:
            self.acemd.TCL = 'set numsteps {}\n'.format(numsteps)
        if self.useconstraints:
            # Turn on constraints
            self.acemd.constraints = 'on'
            self.acemd.constraintscaling = '1.0'
        else:
            if len(self.constraints) != 0:
                logger.warning('You have setup constraints to {} but constraints are turned off. '
                               'If you want to use constraints, define useconstraints=True'.format(self.constraints))
        self.acemd.setup(inputdir, outputdir, overwrite=True)

        # Adding constraints
        if self.useconstraints:
            inmol = Molecule(os.path.join(inputdir, self.acemd.coordinates))
            inmol.set('occupancy', 0)
            inmol.set('beta', 0)
            if len(self.constraints) == 0:
                raise RuntimeError('You have set the production to use constraints (useconstraints=True), but have not '
                               'defined any constraints (constraints={}).')
            else:
                for sel in self.constraints:
                    inmol.set('beta', self.constraints[sel], sel)
                outfile = os.path.join(outputdir, self.acemd.coordinates)
                inmol.write(outfile)
                self.acemd.consref = self.acemd.coordinates

if __name__ == "__main__":
    import htmd.home
    from htmd.util import tempname
    md = Production()
    md.temperature = 300
    md.fb_reference = 'protein and name CA'
    md.fb_selection = 'segname L and noh'
    md.acemd.extendedsystem = None  # use different data
    md.acemd.binindex = None  # use different data
    md.fb_box = [-20, 20, -20, 20, 43, 45]
    md.fb_k = 5
    md.write(htmd.home.home() +'/data/equilibrate', tempname())
    md.fb_k = 0
    md.write(htmd.home.home() +'/data/equilibrate', tempname())
    md.useconstraints = True
    md.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
    md.write(htmd.home.home() +'/data/equilibrate', tempname())