# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.apps.acemd import Acemd as Acemd2
from htmd.mdengine.acemd.acemd import Acemd, _Restraint, GroupRestraint, AtomRestraint
from htmd.config import _config
from protocolinterface import ProtocolInterface, val
import os
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Equilibration(ProtocolInterface):
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
    def __init__(self, _version=_config['acemdversion']):
        super().__init__()
        self._version = _version
        self._arg('acemd', ':class:`Acemd2 <htmd.apps.acemd.Acemd>` or :class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>`'
                           ' object', 'Acemd class object', None, val.Object([Acemd2, Acemd]))
        self._arg('runtime', 'float', 'Running time of the simulation.', 25000, val.Number(float, '0POS'))
        self._arg('timeunits', 'str', 'Units for time arguments. Can be \'steps\', \'ns\' etc.', 'steps', val.String())
        self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin', 300, val.Number(float, 'ANY'))

        self._arg('fb_k', 'float', 'Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5', 0,
                  val.Number(float, 'ANY'))
        self._arg('fb_reference', 'str', 'Reference selection to use as dynamic center of the flatbottom box.', 'none',
                  val.String())
        self._arg('fb_selection', 'str', 'Selection of atoms to apply the flatbottom potential', 'none', val.String())
        self._arg('fb_box', 'list', 'Position of the flatbottom box in term of the reference center given as '
                                    '[xmin, xmax, ymin, ymax, zmin, zmax]',
                  [0, 0, 0, 0, 0, 0], val.Number(float, 'ANY'), nargs=6)

        self._arg('constraints', 'dict', 'A dictionary of atomselections and values of the constraint to be '
                                         'applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.',
                  {'protein and noh and not name CA': 0.1, 'protein and name CA': 1}, val.Dictionary(key_type=str))
        self._arg('useconstantratio', 'bool', 'For membrane protein simulations set it to true so that the barostat '
                                              'does not modify the xy aspect ratio.', False, val.Boolean())

        self._arg('nvtsteps', 'int', 'Number of initial steps to apply NVT in units of 4fs.', 500,
                  val.Number(int, 'ANY'))
        self._arg('constraintsteps', 'int', 'Number of initial steps to apply constraints in units of 4fs. Defaults '
                                            'to half the simulation time.', None, val.Number(int, 'ANY'))
        self._arg('restraints', 'list', 'A list of restraint objects. Only works with {}(_version=3),'
                                        'see :class:`AtomRestraint <htmd.mdengine.acemd.acemd.AtomRestraint>` and'
                                        ':class:`GroupRestraint <htmd.mdengine.acemd.acemd.GroupRestraint>`'
                                        ')'.format(self.__class__.__name__), None, val.Object(_Restraint), nargs='*')

        if self._version == 2:
            self.acemd = Acemd2()
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
        elif self._version == 3:
            self.acemd = Acemd()
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
            self.acemd.thermostatdamping = 1
            self.acemd.pme = 'on'
            self.acemd.barostat = 'on'
            self.acemd.barostatpressure = 1.01325
            self.acemd.minimize = 500
        else:
            raise ValueError('_version can not be {}. Choose either 2 or 3.'.format(self._version))

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
                                   'Please set the {name:}.acemd.{f:} property to '
                                   'point to the {f:} file'.format(f=field, i=inputdir, name=self.__class__.__name__))

        if self._version == 2:
            if self.acemd.consref is None:
                self.acemd.consref = self.acemd.coordinates

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.parameters.endswith('structure.prmtop'):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None
            if self._version == 2:
                self.acemd.scaling14 = '0.8333333'
                self.acemd.amber = 'on'

    def _constraints2restraints(self, constrsteps):

        restraints = list()
        for constr in sorted(self.constraints):
            restraints.append(AtomRestraint(constr, 0, [(self.constraints[constr], 0), (0, constrsteps)]))

        return restraints

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

        from moleculekit.molecule import Molecule

        # Do version consistency check
        if (self._version == 2 and not isinstance(self.acemd, Acemd2)) and \
                (self._version == 3 and not isinstance(self.acemd, Acemd)):
            raise RuntimeError('Acemd object version ({}) inconsistent with protocol version at instantiation '
                               '({})'.format(type(self.acemd), self._version))

        self._findFiles(inputdir)
        self._amberFixes()

        from htmd.units import convert
        numsteps = convert(self.timeunits, 'timesteps', self.runtime, timestep=self.acemd.timestep)
        if self._version == 3:
            self.acemd.temperature = self.temperature
            self.acemd.thermostattemperature = self.temperature
            self.acemd.run = str(numsteps)

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        from htmd.builder.builder import detectCisPeptideBonds
        detectCisPeptideBonds(inmol)

        if np.any(inmol.atomselect('lipids')) and not self.useconstantratio:
            logger.warning('Lipids detected in input structure. We highly recommend setting useconstantratio=True '
                           'for membrane simulations.')

        if self.constraintsteps is None:
            constrsteps = int(numsteps / 2)
        else:
            constrsteps = int(self.constraintsteps)

        if self._version == 2:
            if self.restraints:
                raise RuntimeWarning('restraints are only available on {}(_version=3)'.format(self.__class__.__name__))
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
        elif self._version == 3:
            if self.restraints is not None:
                logger.info('Using user-provided restraints and ignoring constraints and fb_potential')
                self.acemd.restraints = self.restraints
            else:
                logger.warning('Converting constraints and fb_potential to restraints. This is a convenience '
                               'functional conversion. We recommend start using restraints with '
                               '{}(_version=3)'.format(self.__class__.__name__))
                restraints = list()
                # convert constraints to restraints and add them
                if self.constraints is not None:
                    restraints += self._constraints2restraints(constrsteps)
                # convert fb_potential to restraints and add them
                if self.fb_k > 0:
                    restraints += self._fb_potential2restraints(inputdir)
                self.acemd.restraints = restraints

        if ((self._version == 2) and self.acemd.celldimension is None and self.acemd.extendedsystem is None) or \
            ((self._version == 3) and self.acemd.boxsize is None and self.acemd.extendedsystem is None):
            coords = inmol.get('coords', sel='water')
            if coords.size == 0:  # It's a vacuum simulation
                coords = inmol.get('coords', sel='all')
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
                dim += 12.
            else:
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
            if self._version == 2:
                self.acemd.celldimension = '{} {} {}'.format(dim[0], dim[1], dim[2])
            else:
                self.acemd.boxsize = '{} {} {}'.format(dim[0], dim[1], dim[2])

        if self.useconstantratio:
            self.acemd.useconstantratio = 'on'

        self.acemd.setup(inputdir, outputdir, overwrite=True)

        if self._version == 2:
            # Adding constraints by writing them to the consref file
            inconsreffile = os.path.join(inputdir, self.acemd.consref)
            consrefmol = Molecule(inconsreffile)
            consrefmol.set('occupancy', 0)
            consrefmol.set('beta', 0)
            if len(self.constraints) == 0:
                raise RuntimeError('You have not defined any constraints for the {} ('
                                   'constraints={{}}).'.format(self.__class__.__name__))
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
class _TestEquilibration(unittest.TestCase):
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

        eq = Equilibration(_version=2)
        eq.runtime = 4
        eq.timeunits = 'ns'
        eq.temperature = 300
        # keep protein on the same place during all equilibration
        from htmd.units import convert
        eq.constraintsteps = convert(eq.timeunits, 'timesteps', eq.runtime, timestep=eq.acemd.timestep)
        eq.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
        eq.fb_reference = 'protein and name CA'
        eq.fb_selection = 'resname MOL and noh'
        eq.fb_box = [-21, 21, -19, 19, 29, 30]
        eq.fb_k = 5
        tmpdir = tempname()
        eq.write(home(dataDir=os.path.join('test-protocols', 'build', 'protLig')), tmpdir)

        # Compare with reference
        refdir = home(dataDir=os.path.join('test-protocols', 'equilibration', 'acemd2', 'protLig', 'prerun'))
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
        self._compareResultFolders(refdir, tmpdir, 'protLig')

    def test_acemd3(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import os

        eq = Equilibration(_version=3)
        eq.runtime = 4
        eq.timeunits = 'ns'
        eq.temperature = 300
        # keep protein on the same place during all equilibration
        from htmd.units import convert
        eq.constraintsteps = convert(eq.timeunits, 'timesteps', eq.runtime, timestep=eq.acemd.timestep)
        eq.constraints = {'protein and name CA': 1, 'protein and noh and not name CA': 0.1}
        eq.fb_reference = 'protein and name CA'
        eq.fb_selection = 'resname MOL and noh'
        eq.fb_box = [-21, 21, -19, 19, 29, 30]
        eq.fb_k = 5
        tmpdir = tempname()
        eq.write(home(dataDir=os.path.join('test-protocols', 'build', 'protLig')), tmpdir)

        # Compare with reference
        refdir = home(dataDir=os.path.join('test-protocols', 'equilibration', 'acemd3', 'protLig', 'prerun'))
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
        self._compareResultFolders(refdir, tmpdir, 'protLig')

    @unittest.skipUnless('ACE3ARG' in os.environ, 'Untrusted PR')
    def test_run_water(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import subprocess
        from subprocess import check_output
        import shutil
        import os

        acemd3exe = shutil.which('acemd3', mode=os.X_OK)
        if not acemd3exe:
            raise NameError('Could not find acemd3, or no execute permissions are given')

        for system in ['amber-build', 'charmm-build']:
            eq = Equilibration(_version=3)
            eq.runtime = 5
            eq.timeunits = 'steps'
            eq.temperature = 300
            eq.constraints = {}
            # Set these down for tiny box size of water
            eq.acemd.cutoff = 3
            eq.acemd.switchdistance = 2
            eq.acemd.minimize = 50  # Do few steps to finish fast
            ######
            tmpdir = tempname()
            eq.write(home(dataDir=os.path.join('test-acemd', 'tiny-water', system)), tmpdir)
            try:
                res = check_output(['acemd3', '--platform', 'CPU', os.getenv('ACE3ARG')], cwd=tmpdir)
            except subprocess.CalledProcessError as exc:
                assert False, f'Failed to run due to error: {exc}\n\n ---> Error log:\n\n{exc.output.decode("ascii")}'
            res = res.decode('utf-8').strip()
            assert 'Completed minimization' in res, 'Failed at system ' + system
            assert res.endswith('Completed simulation!'), 'Failed at system ' + system

if __name__ == "__main__":
    unittest.main(verbosity=2)