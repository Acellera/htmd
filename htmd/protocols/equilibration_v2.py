# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.molecule.molecule import Molecule
from htmd.apps.acemd import Acemd
from htmd.protocols.oldprotocolinterface import TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
from htmd.protocols.oldprotocolinterface import ProtocolInterface as OldProtocolInterface
from protocolinterface import ProtocolInterface, val
import os
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Equilibration(OldProtocolInterface):
    """ Equilibration protocol v2

        Equilibration protocol for globular and membrane proteins
        It includes a flatbottom potential box to retrain a ligand
        for example within this box.

        Parameters
        ----------
        runtime : float, default=0
            Running time of the simulation.
        timeunits : str, default='steps'
            Units for time arguments. Can be 'steps', 'ns' etc.
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
        constraints : dict, default={'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
            A dictionary of atomselections and values of the constraint to be applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.
        nvtsteps : int, default=500
            Number of initial steps to apply NVT in units of 4fs.
        constraintsteps : int, default=None
            Number of initial steps to apply constraints in units of 4fs. Defaults to half the simulation time.

        Example
        -------
        >>> from htmd.protocols.equilibration_v2 import Equilibration
        >>> md = Equilibration()
        >>> md.runtime = 4
        >>> md.timeunits = 'ns'
        >>> md.temperature = 300
        >>> md.useconstantratio = True  # only for membrane sims
        >>> # this is only needed for setting the flatbottom potential, otherwise remove it
        >>> md.fb_reference = 'protein and resid 293'
        >>> md.fb_selection = 'segname L and noh'
        >>> md.fb_box = [-25, 25, -25, 25, 43, 45]
        >>> md.fb_k = 5
        >>> md.write('./build','./equil')
    """
    def __init__(self):
        super().__init__()
        self._cmdObject('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Acemd)
        self._cmdValue('runtime', 'float', 'Running time of the simulation.', 25000, TYPE_FLOAT, RANGE_0POS)
        self._cmdString('timeunits', 'str', 'Units for time arguments. Can be \'steps\', \'ns\' etc.', 'steps')
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
        self._cmdDict('constraints', 'dict', 'A dictionary of atomselections and values of the constraint to be '
                                             'applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.'
                                             , {'protein and noh and not name CA': 0.1, 'protein and name CA': 1})
        self._cmdValue('nvtsteps', 'int', 'Number of initial steps to apply NVT in units of 4fs.', 500, TYPE_INT,
                       RANGE_ANY)
        self._cmdValue('constraintsteps', 'int', 'Number of initial steps to apply constraints in units of 4fs. '
                                                 'Defaults to half the simulation time.', None, TYPE_INT, RANGE_ANY)

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
        self.acemd.TCL = ('''
set numsteps {NUMSTEPS}
set temperature {TEMPERATURE}
set nvtsteps {NVTSTEPS}
set constraintsteps {CONSTRAINTSTEPS}
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
  global ref sel numsteps fb_K fb_box nvtsteps constraintsteps
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
##EQUIL
  set step [ getstep ]
  if { $step > $nvtsteps } {
    berendsenpressure  on
  } else {
    berendsenpressure  off
  }
  if { $step > $constraintsteps } {
    constraintscaling 0
  } else {
    constraintscaling [expr 1 - 0.95*$step/$constraintsteps]
  }
}
proc calcforces_endstep { } { }
''')

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

        from htmd.units import convert
        numsteps = convert(self.timeunits, 'timesteps', self.runtime, timestep=self.acemd.timestep)

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        if np.any(inmol.atomselect('lipids')) and not self.useconstantratio:
            logger.warning('Lipids detected in input structure. We highly recommend setting useconstantratio=True '
                           'for membrane simulations.')

        if self.constraintsteps is None:
            constrsteps = int(numsteps / 2)
        else:
            constrsteps = int(self.constraintsteps)

        if isinstance(self.acemd.TCL, tuple):
            tcl = list(self.acemd.TCL)
            tcl[0] = tcl[0].format(NUMSTEPS=numsteps, KCONST=self.fb_k,
                                   REFINDEX=' '.join(map(str, inmol.get('index', self.fb_reference))),
                                   SELINDEX=' '.join(map(str, inmol.get('index', self.fb_selection))),
                                   BOX=' '.join(map(str, self.fb_box)),
                                   NVTSTEPS=self.nvtsteps, CONSTRAINTSTEPS=constrsteps, TEMPERATURE=self.temperature)
            self.acemd.TCL = tcl[0] + tcl[1]
        else:
            logger.warning('{} default TCL was already formatted.'.format(self.__class__.__name__))

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

        # Adding constraints by writing them to the consref file
        inconsreffile = os.path.join(inputdir, self.acemd.consref)
        consrefmol = Molecule(inconsreffile)
        consrefmol.set('occupancy', 0)
        consrefmol.set('beta', 0)
        if len(self.constraints) == 0:
            raise RuntimeError('You have not defined any constraints for the Equilibration (constraints={}).')
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


from htmd.apps.acemd3 import Acemd3, _Restraint, GroupRestraint, AtomRestraint
class EquilibrationAcemd3(ProtocolInterface):
    """ Equilibration protocol v2

        Equilibration protocol for globular and membrane proteins
        It includes a flatbottom potential box to retrain a ligand
        for example within this box.

        Parameters
        ----------
        runtime : float, default=0
            Running time of the simulation.
        timeunits : str, default='steps'
            Units for time arguments. Can be 'steps', 'ns' etc.
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
        constraints : dict, default={'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
            A dictionary of atomselections and values of the constraint to be applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.
        nvtsteps : int, default=500
            Number of initial steps to apply NVT in units of 4fs.
        constraintsteps : int, default=None
            Number of initial steps to apply constraints in units of 4fs. Defaults to half the simulation time.

        Example
        -------
        >>> from htmd.protocols.equilibration_v2 import Equilibration
        >>> md = Equilibration()
        >>> md.runtime = 4
        >>> md.timeunits = 'ns'
        >>> md.temperature = 300
        >>> md.useconstantratio = True  # only for membrane sims
        >>> # this is only needed for setting the flatbottom potential, otherwise remove it
        >>> md.fb_reference = 'protein and resid 293'
        >>> md.fb_selection = 'segname L and noh'
        >>> md.fb_box = [-25, 25, -25, 25, 43, 45]
        >>> md.fb_k = 5
        >>> md.write('./build','./equil')
    """
    def __init__(self):
        super().__init__()
        self._arg('acemd', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, val.Object(Acemd3))
        self._arg('runtime', 'float', 'Running time of the simulation.', 0, val.Number(float, '0POS'))
        self._arg('timeunits', 'str', 'Units for runtime. Can be \'steps\', \'ns\' etc.', 'steps', val.String())
        self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, val.Number(float, 'ANY'))
        self._arg('restraints', 'list', 'A list of restraint objects', None, val.Object(_Restraint), nargs='*')
        self._arg('useconstantratio', 'bool', 'For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.', False, val.Boolean())
        self._arg('constraintsteps', 'int', 'Number of initial steps to apply constraints in units of 4fs. Defaults to half the simulation time.', None, val.Number(int, 'ANY'))

        self.acemd = Acemd3()
        self.acemd.coordinates = None
        self.acemd.structure = None
        self.acemd.parameters = None
        self.acemd.restart = 'on'
        self.acemd.trajectoryfile = 'output.xtc'
        self.acemd.trajectoryfreq = 25000
        self.acemd.timestep = 4
        self.acemd.switching = 'on'
        self.acemd.switchdist = 7.5
        self.acemd.cutoff = 9
        self.acemd.thermostat = 'on'
        self.acemd.thermostatdamping = 1
        self.acemd.pme = 'on'
        self.acemd.barostat = 'on'
        self.acemd.barostatpressure = 1.01325
        self.acemd.minimize = 500

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

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.parameters.endswith('structure.prmtop'):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None

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

        from htmd.units import convert
        numsteps = convert(self.timeunits, 'timesteps', self.runtime, timestep=self.acemd.timestep)
        self.acemd.temperature = self.temperature
        self.acemd.thermostattemp = self.temperature
        self.acemd.run = str(numsteps)

        if self.constraintsteps is None:
            constrsteps = int(numsteps / 2)
        else:
            constrsteps = int(self.constraintsteps)

        # Adding the default restraints of the equilibration
        restraints = list()
        restraints.append(AtomRestraint('protein and noh and not name CA', 0, [(0.1, 0), (0, constrsteps)]))
        restraints.append(AtomRestraint('protein and name CA', 0, [(1, 0), (0, constrsteps)]))
        if self.restraints is not None:
            restraints += self.restraints
        self.acemd.restraints = restraints

        if self.acemd.celldimension is None and self.acemd.extendedsystem is None:
            inmol = Molecule(os.path.join(inputdir, self.acemd.coordinates))
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


if __name__ == "__main__":
    from htmd.home import home
    from htmd.util import tempname
    import filecmp
    from glob import glob
    from os.path import join

    def _cutfirstline(infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(outfile, 'w') as fout:
            fout.writelines(data[1:])


    def _compareResultFolders(compare, tmpdir, pid):
        ignore_ftypes = ('.log', '.txt')
        files = []
        deletefiles = []
        for f in glob(join(compare, '*')):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith('prmtop'):
                _cutfirstline(f, join(compare, fname + '.mod'))
                _cutfirstline(join(tmpdir, fname), os.path.join(tmpdir, fname + '.mod'))
                files.append(os.path.basename(f) + '.mod')
                deletefiles.append(join(compare, fname + '.mod'))
            else:
                files.append(os.path.basename(f))

        match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError(
                'Different results produced by amber.build for test {} between {} and {} in files {}.'.format(pid, compare, tmpdir, mismatch))

        for f in deletefiles:
            os.remove(f)

    pdbid = '3PTB'
    eq = Equilibration()
    eq.runtime = 4
    eq.timeunits = 'ns'
    eq.temperature = 300
    eq.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
    eq.fb_reference = 'protein and name CA'
    eq.fb_selection = 'segname L and noh'
    eq.fb_box = [-20, 20, -20, 20, 43, 45]
    eq.fb_k = 5
    tmpdir = tempname()
    eq.write(home(dataDir=os.path.join('test-amber-build', pdbid)), tmpdir)

    # Compare with reference
    refdir = home(dataDir=os.path.join('test-equilibration', pdbid, 'prerun'))
    files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
    # match, mismatch, error = filecmp.cmpfiles(refdir, tmpdir, files, shallow=False)
    _compareResultFolders(refdir, tmpdir, '3PTB')

    # if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
    #         raise RuntimeError('Different results produced by Equilibration.write for '
    #                            'test {} between {} and {} in files {}.'.format(pdbid, refdir, tmpdir, mismatch))

    # from htmd.protocols.production_v5 import ProductionAcemd3, GroupRestraint, AtomRestraint
    r = list()
    r.append(GroupRestraint('segid P2', 5, [(10, '10ns'), (5, '15ns'), (0, '20ns')], axes='z'))
    r.append(AtomRestraint('segid P1 and name CA', 0.1, [(10, '10ns'), (5, '15ns'), (0, '20ns')]))

    eq = EquilibrationAcemd3()
    eq.runtime = 4
    eq.timeunits = 'ns'
    eq.temperature = 300
    eq.restraints = r
    # eq.write('/workspace5/pablo/bound_KIX_cMYB/1_build/2-struct/', '/tmp/testdir2/')


