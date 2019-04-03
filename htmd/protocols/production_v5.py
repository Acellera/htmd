# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from moleculekit.molecule import Molecule
from htmd.protocols.oldprotocolinterface import TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
from htmd.protocols.oldprotocolinterface import ProtocolInterface as OldProtocolInterface
from htmd.decorators import _Deprecated
from htmd.apps.acemd import Acemd
import os
import numpy as np
import logging
logger = logging.getLogger(__name__)


@_Deprecated('1.13.6')
class Production(OldProtocolInterface):
    ''' Production protocol v5

        Production protocol for globular and membrane proteins. You can optionally define a flatbottom potential box and
        atom constraints for the production run.

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
        useconstantratio : bool, default=False
            For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.
        useconstraints : bool, default=False
            Apply constraints to the production simulation, defined by the constraints parameter
        constraints : dict, default={}
            A dictionary of atomselections and values of the constraint to be applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.
        adaptive : bool, default=False
            Set to True if making production runs for adaptive sampling.
    '''
    def __init__(self):
        super().__init__()
        self._cmdObject('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Acemd)
        self._cmdValue('runtime', 'float', 'Running time of the simulation.', 0, TYPE_FLOAT, RANGE_0POS)
        self._cmdString('timeunits', 'str', 'Units for runtime. Can be \'steps\', \'ns\' etc.', 'steps')
        self._cmdValue('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, TYPE_FLOAT, RANGE_ANY)
        self._cmdValue('fb_k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0,
                       TYPE_FLOAT, RANGE_ANY)
        self._cmdString('fb_reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.',
                        'none')
        self._cmdString('fb_selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none')
        self._cmdList('fb_box', 'list', 'Position of the flatbottom box in term of the reference center given as '
                                        '[xmin, xmax, ymin, ymax, zmin, zmax]', [0, 0, 0, 0, 0, 0])
        self._cmdBoolean('useconstantratio', 'bool', 'For membrane protein simulations set it to true so that the '
                                                     'barostat does not modify the xy aspect ratio.', False)
        self._cmdList('useconstraints', 'bool', 'Apply constraints to the production simulation, defined by the '
                                                'constraints parameter', False)
        self._cmdDict('constraints', 'dict', 'A dictionary of atomselections and values of the constraint to be '
                                             'applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.', {})
        self._cmdBoolean('adaptive', 'bool', 'Set to True if making production runs for adaptive sampling.', False)

        self.acemd = Acemd()
        #self.acemd.binindex='input.idx'
        self.acemd.binindex = None
        self.acemd.binvelocities = None
        self.acemd.bincoordinates = 'output.coor'
        self.acemd.extendedsystem='output.xsc'
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
        self.acemd.langevindamping = '0.1'
        self.acemd.pme = 'on'
        self.acemd.pmegridspacing = '1.0'
        self.acemd.fullelectfrequency = '2'
        self.acemd.energyfreq = '5000'
        self.acemd.consref = None
        self.acemd.run = '$numsteps'
        self.acemd.TCL = ('''
set numsteps {NUMSTEPS}
set temperature {TEMPERATURE}
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
                raise RuntimeError('Could not locate file {} set by the user for argument '
                                   'Production.acemd.{}'.format(os.path.join(inputdir, userval), field))

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
                                   'Please set the Production.acemd.{f:} property to '
                                   'point to the {f:} file'.format(f=field, i=inputdir))

            if self.useconstraints and self.acemd.consref is None:
                self.acemd.consref = self.acemd.coordinates

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.structure.endswith('.prmtop'):
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

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        if np.any(inmol.atomselect('lipids')) and not self.useconstantratio:
            logger.warning('Lipids detected in input structure. We highly recommend setting useconstantratio=True '
                           'for membrane simulations.')

        if self.fb_k > 0:  # use TCL only for flatbottom
            self.acemd.tclforces = 'on'
            tcl = list(self.acemd.TCL)
            tcl[0] = tcl[0].format(NUMSTEPS=numsteps, TEMPERATURE=self.temperature, KCONST=self.fb_k,
                                   REFINDEX=' '.join(map(str, inmol.get('index', self.fb_reference))),
                                   SELINDEX=' '.join(map(str, inmol.get('index', self.fb_selection))),
                                   BOX=' '.join(map(str, self.fb_box)))
            self.acemd.TCL = tcl[0] + tcl[1]
        else:
            self.acemd.TCL = 'set numsteps {NUMSTEPS}\n' \
                             'set temperature {TEMPERATURE}\n'.format(NUMSTEPS=numsteps, TEMPERATURE=self.temperature)

        if self.useconstraints:
            # Turn on constraints
            self.acemd.constraints = 'on'
            self.acemd.constraintscaling = '1.0'
        else:
            if len(self.constraints) != 0:
                logger.warning('You have setup constraints to {} but constraints are turned off. '
                               'If you want to use constraints, define useconstraints=True'.format(self.constraints))

        if self.useconstantratio:
            self.acemd.useconstantratio = 'on'

        if self.adaptive:
            self.acemd.binvelocities = None

        self.acemd.setup(inputdir, outputdir, overwrite=True)

        # Adding constraints by writing them to the consref file
        if self.useconstraints:
            inconsreffile = os.path.join(inputdir, self.acemd.consref)
            consrefmol = Molecule(inconsreffile)
            consrefmol.set('occupancy', 0)
            consrefmol.set('beta', 0)
            if len(self.constraints) == 0:
                raise RuntimeError('You have set the production to use constraints (useconstraints=True), but have not '
                                   'defined any constraints (constraints={}).')
            else:
                for sel in self.constraints:
                    consrefmol.set('beta', self.constraints[sel], sel)
            outconsreffile = os.path.join(outputdir, self.acemd.consref)
            consrefmol.write(outconsreffile)

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
    from htmd.home import home
    from htmd.util import tempname
    import filecmp
    from glob import glob

    pdbid = '3PTB'
    pd = Production()
    pd.runtime = 4
    pd.timeunits = 'ns'
    pd.temperature = 350
    pd.useconstraints = True
    pd.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
    pd.fb_reference = 'protein and name CA'
    pd.fb_selection = 'segname L and noh'
    pd.fb_box = [-20, 20, -20, 20, 43, 45]
    pd.fb_k = 5
    tmpdir = tempname()
    pd.write(home(dataDir=os.path.join('test-protocols', 'equilibration', pdbid, 'postrun')), tmpdir)

    # Compare with reference
    refdir = home(dataDir=os.path.join('test-production', pdbid, 'prerun'))
    files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
    match, mismatch, error = filecmp.cmpfiles(refdir, tmpdir, files, shallow=False)

    if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError('Different results produced by Equilibration.write for '
                               'test {} between {} and {} in files {}.'.format(pdbid, refdir, tmpdir, mismatch))
