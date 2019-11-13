# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from protocolinterface import ProtocolInterface, val
from htmd.apps.acemd import Acemd as Acemd2
from htmd.mdengine.acemd.acemd import Acemd, _Restraint, GroupRestraint, AtomRestraint
import os
import numpy as np
import logging
from htmd.config import _config
logger = logging.getLogger(__name__)


class Production(ProtocolInterface):
    """ Production protocol v6

        Production protocol for globular and membrane proteins. You can optionally define a flatbottom potential box and
        atom constraints for the production run.

        An Acemd class object is stored in the Production object which can be used to modify futher options.
        For documentation on further options see :class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>`

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
    """
    def __init__(self, _version=_config['acemdversion']):
        super().__init__()
        self._version = _version
        self._arg('acemd', ':class:`Acemd2 <htmd.apps.acemd.Acemd>` or :class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>`'
                           ' object', 'Acemd class object', None, val.Object([Acemd2, Acemd]))
        self._arg('runtime', 'float', 'Running time of the simulation.', 25000, val.Number(float, '0POS'))
        self._arg('timeunits', 'str', 'Units for runtime. Can be \'steps\', \'ns\' etc.', 'steps', val.String())
        self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, val.Number(float, 'ANY'))
        self._arg('fb_k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0,
                  val.Number(float, 'ANY'))
        self._arg('fb_reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.',
                  'none', val.String())
        self._arg('fb_selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none', val.String())
        self._arg('fb_box', 'list', 'Position of the flatbottom box in term of the reference center given as '
                                    '[xmin, xmax, ymin, ymax, zmin, zmax]',
                  [0, 0, 0, 0, 0, 0], val.Number(float, 'ANY'), nargs=6)
        self._arg('useconstantratio', 'bool', 'For membrane protein simulations set it to true so that the barostat '
                                              'does not modify the xy aspect ratio.', False, val.Boolean())
        self._arg('useconstraints', 'bool', 'Apply constraints to the production simulation, defined by the '
                                            'constraints parameter', False, val.Boolean())
        self._arg('constraints', 'dict', 'A dictionary of atomselections and values of the constraint to be '
                                         'applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.', {},
                  val.Dictionary(key_type=str))
        self._arg('adaptive', 'bool', 'Set to True if making production runs for adaptive sampling.', False,
                  val.Boolean())
        self._arg('restraints', 'list', 'A list of restraint objects. Only works with {}(_version=3),'
                                'see :class:`AtomRestraint <htmd.mdengine.acemd.acemd.AtomRestraint>` and'
                                ':class:`GroupRestraint <htmd.mdengine.acemd.acemd.GroupRestraint>`'
                                ')'.format(self.__class__.__name__), None, val.Object(_Restraint), nargs='*')

        if self._version == 2:
            self.acemd = Acemd2()
            self.acemd.binindex = None
            self.acemd.binvelocities = None
            self.acemd.bincoordinates = 'output.coor'
            self.acemd.extendedsystem = 'output.xsc'
            self.acemd.coordinates = None
            self.acemd.structure = None
            self.acemd.parameters = None
            self.acemd.temperature = None
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
            self.acemd.langevintemp = None
            self.acemd.langevindamping = '0.1'
            self.acemd.pme = 'on'
            self.acemd.pmegridspacing = '1.0'
            self.acemd.fullelectfrequency = '2'
            self.acemd.energyfreq = '5000'
            self.acemd.consref = None
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
        elif self._version == 3:
            self.acemd = Acemd()
            self.acemd.binvelocities = None
            self.acemd.bincoordinates = 'output.coor'
            self.acemd.extendedsystem = 'output.xsc'
            self.acemd.coordinates = None
            self.acemd.structure = None
            self.acemd.parameters = None
            self.acemd.restart = 'on'
            self.acemd.trajectoryfile = 'output.xtc'
            self.acemd.trajectoryperiod = 25000
            self.acemd.timestep = 4
            self.acemd.switching = 'on'
            self.acemd.switchdistance = 7.5
            self.acemd.cutoff = 9
            self.acemd.thermostat = 'on'
            self.acemd.thermostatdamping = 0.1
            self.acemd.pme = 'on'
        else:
            raise ValueError('_version can not be {}. Choose either 2 or 3.'.format(self._version))

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

            if self._version == 2:
                if self.useconstraints and self.acemd.consref is None:
                    self.acemd.consref = self.acemd.coordinates

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.structure.endswith('.prmtop'):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None
            if self._version == 2:
                self.acemd.scaling14 = '0.8333333'
                self.acemd.amber = 'on'

    def _fb_potential2restraints(self, inputdir):
        from moleculekit.molecule import Molecule
        restraints = list()

        fb_box = np.array(self.fb_box)
        # convert fb_box to width
        width = list(np.concatenate(np.diff(np.array([fb_box[::2], fb_box[1::2]]), axis=0)))

        # If fb_box is not symmetrical
        if not np.all(fb_box[::2] == -fb_box[1::2]):
            # convert fb_box and fb_reference to fbcentre and width
            mol = Molecule(os.path.join(inputdir, self.acemd.structure))
            mol.read(os.path.join(inputdir, self.acemd.coordinates))
            fb_refcentre = mol.get('coords', sel=self.fb_reference).mean(axis=0).squeeze()

            fbcentre = list(np.around(np.mean(np.array([fb_box[::2], fb_box[1::2]]), axis=0) + fb_refcentre, 3))
            restraints.append(GroupRestraint(self.fb_selection, width, [(self.fb_k, 0)], fbcentre=fbcentre))
        else:
            restraints.append(GroupRestraint(self.fb_selection, width, [(self.fb_k, 0)], fbcentresel=self.fb_reference))

        return restraints

    def _constraints2restraints(self):

        restraints = list()
        for constr in sorted(self.constraints):
            restraints.append(AtomRestraint(constr, 0, [(self.constraints[constr], 0)]))

        return restraints

    def write(self, inputdir, outputdir):
        """ Writes the production protocol and files into a folder.

        Parameters
        ----------
        inputdir : str
            Path to a directory containing the files produced by a equilibration process.
        outputdir : str
            Directory where to write the production setup files.
        """
        from moleculekit.molecule import Molecule

        # Do version consistency check
        if (self._version == 2 and not isinstance(self.acemd, Acemd2)) and \
                (self._version == 3 and not isinstance(self.acemd, Acemd)):
            raise RuntimeError('Acemd object version ({}) inconsistent with protocol version at instantiation '
                               '({})'.format(type(self.acemd), self._version))

        self._findFiles(inputdir)
        self._amberFixes()

        if self._version == 2:
            self.acemd.temperature = str(self.temperature)
            self.acemd.langevintemp = str(self.temperature)
        elif self._version == 3:
            self.acemd.temperature = self.temperature
            self.acemd.thermostattemperature = self.temperature

        from htmd.units import convert
        numsteps = convert(self.timeunits, 'timesteps', self.runtime, timestep=self.acemd.timestep)
        if self._version == 3:
            self.acemd.run = str(numsteps)

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        from htmd.builder.builder import detectCisPeptideBonds
        detectCisPeptideBonds(inmol)

        if np.any(inmol.atomselect('lipids')) and not self.useconstantratio:
            logger.warning('Lipids detected in input structure. We highly recommend setting useconstantratio=True '
                           'for membrane simulations.')

        if self._version == 2:
            if self.restraints:
                raise RuntimeWarning('restraints are only available on {}(_version=3)'.format(self.__class__.__name__))
            if self.fb_k > 0:  # use TCL only for flatbottom
                self.acemd.tclforces = 'on'
                if isinstance(self.acemd.TCL, tuple):
                    tcl = list(self.acemd.TCL)
                    tcl[0] = tcl[0].format(NUMSTEPS=numsteps, TEMPERATURE=self.temperature, KCONST=self.fb_k,
                                           REFINDEX=' '.join(map(str, inmol.get('index', self.fb_reference))),
                                           SELINDEX=' '.join(map(str, inmol.get('index', self.fb_selection))),
                                           BOX=' '.join(map(str, self.fb_box)))
                    self.acemd.TCL = tcl[0] + tcl[1]
                else:
                    logger.warning('{} default TCL was already formatted.'.format(self.__class__.__name__))
            else:
                self.acemd.TCL = 'set numsteps {NUMSTEPS}\n'.format(NUMSTEPS=numsteps)
            if self.useconstraints:
                # Turn on constraints
                self.acemd.constraints = 'on'
                self.acemd.constraintscaling = '1.0'
            else:
                if len(self.constraints) != 0:
                    logger.warning('You have setup constraints to {} but constraints are turned off. '
                                   'If you want to use constraints, define '
                                   'useconstraints=True'.format(self.constraints))
        elif self._version == 3:
            if self.restraints is not None:
                logger.info('Using user-provided restraints and ignoring constraints and fb_potential')
                self.acemd.restraints = self.restraints
            else:
                restraints = list()
                if self.fb_k > 0:
                    logger.warning('Converting fb_potential to restraints. This is a convenience '
                                   'functional conversion. We recommend start using restraints with '
                                   '{}(_version=3)'.format(self.__class__.__name__))
                    restraints += self._fb_potential2restraints(inputdir)
                if self.useconstraints:
                    logger.warning('Converting constraints to restraints. This is a convenience '
                                   'functional conversion. We recommend start using restraints with '
                                   '{}(_version=3)'.format(self.__class__.__name__))
                    restraints += self._constraints2restraints()
                else:
                    if len(self.constraints) != 0:
                        logger.warning('You have setup constraints to {} but constraints are turned off. '
                                       'If you want to use constraints, define '
                                       'useconstraints=True'.format(self.constraints))
                if len(restraints) != 0:
                    self.acemd.restraints = restraints

        if self.useconstantratio:
            self.acemd.useconstantratio = 'on'

        if self.adaptive:
            self.acemd.binvelocities = None

        self.acemd.setup(inputdir, outputdir, overwrite=True)

        if self._version == 2:
            # Adding constraints by writing them to the consref file
            if self.useconstraints:
                inconsreffile = os.path.join(inputdir, self.acemd.consref)
                consrefmol = Molecule(inconsreffile)
                consrefmol.set('occupancy', 0)
                consrefmol.set('beta', 0)
                if len(self.constraints) == 0:
                    raise RuntimeError('You have set the production to use constraints (useconstraints=True), but have '
                                       'not defined any constraints (constraints={}).')
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


import unittest
class _TestProduction(unittest.TestCase):
    def _cutfirstline(self, infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(outfile, 'w') as fout:
            fout.writelines(data[1:])

    def _compareResultFolders(self, compare, tmpdir, pid):
        from glob import glob
        import os
        import filecmp

        ignore_ftypes = ('.log', '.txt')
        files = []
        deletefiles = []
        for f in glob(os.path.join(compare, '*')):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith('prmtop'):
                self._cutfirstline(f, os.path.join(compare, fname + '.mod'))
                self._cutfirstline(os.path.join(tmpdir, fname), os.path.join(tmpdir, fname + '.mod'))
                files.append(os.path.basename(f) + '.mod')
                deletefiles.append(os.path.join(compare, fname + '.mod'))
            else:
                files.append(os.path.basename(f))

        match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError(
                'Different results produced by amber.build for test {} between {} and {} in files {}.'.format(pid, compare, tmpdir, mismatch))

        for f in deletefiles:
            os.remove(f)

    def test_acemd2(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import os

        pd = Production(_version=2)
        pd.runtime = 4
        pd.timeunits = 'ns'
        pd.temperature = 300
        pd.useconstraints = True
        pd.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
        pd.fb_reference = 'protein and name CA'
        pd.fb_selection = 'resname MOL and noh'
        pd.fb_box = [-21, 21, -19, 19, 29, 30]
        pd.fb_k = 5
        tmpdir = tempname()
        pd.write(home(dataDir=os.path.join('test-protocols', 'equilibration', 'acemd2', 'protLig', 'postrun')), tmpdir)

        # Compare with reference
        refdir = home(dataDir=os.path.join('test-protocols', 'production', 'acemd2', 'protLig', 'prerun'))
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
        self._compareResultFolders(refdir, tmpdir, 'protLig')

    def test_acemd3(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import os

        pd = Production(_version=3)
        pd.runtime = 4
        pd.timeunits = 'ns'
        pd.temperature = 300
        pd.useconstraints = True
        pd.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
        pd.fb_reference = 'protein and name CA'
        pd.fb_selection = 'resname MOL and noh'
        pd.fb_box = [-21, 21, -19, 19, 29, 30]
        pd.fb_k = 5
        tmpdir = tempname()
        pd.write(home(dataDir=os.path.join('test-protocols', 'equilibration', 'acemd3', 'protLig', 'postrun')), tmpdir)

        # Compare with reference
        refdir = home(dataDir=os.path.join('test-protocols', 'production', 'acemd3', 'protLig', 'prerun'))
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
        self._compareResultFolders(refdir, tmpdir, 'protLig')

    @unittest.skipUnless('ACE3ARG' in os.environ, 'Untrusted PR')
    def test_run_water(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        from subprocess import check_output
        import subprocess
        import shutil
        import os

        acemd3exe = shutil.which('acemd3', mode=os.X_OK)
        if not acemd3exe:
            raise NameError('Could not find acemd3, or no execute permissions are given')

        for system in ['amber-equil-completed', 'charmm-equil-completed']:
            pd = Production(_version=3)
            pd.runtime = 5
            pd.timeunits = 'steps'
            pd.temperature = 300
            pd.constraints = {}
            # Set these down for tiny box size of water
            pd.acemd.cutoff = 3
            pd.acemd.switchdistance = 2
            ######
            tmpdir = tempname()
            pd.write(home(dataDir=os.path.join('test-acemd', 'tiny-water', system)), tmpdir)
            try:
                res = check_output(['acemd3', '--platform', 'CPU', os.getenv('ACE3ARG')], cwd=tmpdir)
            except subprocess.CalledProcessError as exc:
                assert False, f'Failed to run due to error: {exc}\n\n ---> Error log:\n\n{exc.output.decode("ascii")}'
            res = res.decode('utf-8').strip()
            print(res)
            assert res.endswith('Completed simulation!'), 'Failed at system ' + system


if __name__ == "__main__":
    unittest.main(verbosity=2)


