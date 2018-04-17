# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from protocolinterface import ProtocolInterface, val
from htmd.util import ensurelist
import logging

logger = logging.getLogger(__name__)


class _Restraint:
    def __init__(self, contype, selection, width, constraints, axes='xyz'):
        self.type = contype
        self.selection = selection
        self.width = width
        self.constraints = constraints
        self.axes = axes

    def format(self, maxwidth=None):
        if maxwidth is None:
            res = '{type}Restraint '.format(type=self.type)
        else:
            res = '{type: <{maxwidth}}\t'.format(type=self.type+'Restraint', maxwidth=maxwidth)
        res += '"{sel}" axes {axes} width {width} setpoints'.format(sel=self.selection,
                                                                    axes=self.axes,
                                                                    width=self.width)
        for p in self.constraints:
            res += ' {con}@{time}'.format(con=p[0], time=p[1])
        return res

    def __str__(self):
        return self.format()


class GroupRestraint(_Restraint):
    def __init__(self, selection, width, constraints, axes='xyz'):
        """ A restraint applied on the center of mass of a group of atoms.

        Parameters
        ----------
        selection : str
            Atom selection string of the atoms on whose center of mass the restraints will be applied.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        width : float
            The width of the flat-bottom potential box in Angstrom
        constraints : list of tuples
            A list of contraint/time pairs indicating which value the constraints should have at each point in time.
            Constraints are given in kcal/mol units. Time can be either a timestep or time if one of the suffices
            "us", "ns", "ps", "fs" is used.
        axes : str
            The axes on which to apply the potential.

        Examples
        --------
        >>> gr = GroupRestraint('resname MOL', 10, [(20, '10ns'), (10, '20ns'), (0, '30ns')])
        """
        super().__init__('group', selection, width, constraints, axes)


class AtomRestraint(_Restraint):
    def __init__(self, selection, width, constraints, axes='xyz'):
        """ A restraint applied individually on each atom in the selection

        Parameters
        ----------
        selection : str
            Atom selection string of the atoms on whose center of mass the restraints will be applied.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        width : float
            The width of the flat-bottom potential box in Angstrom
        constraints : list of tuples
            A list of contraint/time pairs indicating which value the constraints should have at each point in time.
            Constraints are given in kcal/mol units. Time can be either a timestep or time if one of the suffices
            "us", "ns", "ps", "fs" is used.
        axes : str
            The axes on which to apply the potential.

        Examples
        --------
        >>> gr = AtomRestraint('name CA', 10, [(20, '10ns'), (10, '20ns'), (0, '30ns')])
        """
        super().__init__('atom', selection, width, constraints, axes)


class Acemd(ProtocolInterface):
    """ ACEMD class.

    Parameters
    ----------
    temperature : float, default=None
        Temperature of the thermostat in Kelvin.
    restart : ('on', 'off'), str, default=None
        Restart simulation.
    restartfreq : int, default=None
        Restart file frequency.
    outputname : str, default=None
        Output file name.
    xtcfile : str, default=None
        Output XTC file name.
    xtcfreq : int, default=None
        XTC sampling frequency in steps.
    timestep : float, default=None
        Simulation timestep.
    rigidbonds : ('none', 'all'), str, default=None
        Enable holonomic constraints on all hydrogen-heavy atom bond terms
    hydrogenscale : float, default=None
        Amount by which to scale the mass of H atoms
    switching : ('on', 'off'), str, default=None
        Enable to apply smoothing and switching functions to electrostatic and VdW forces.
    switchdist : float, default=None
        Range beyond which to begin to apply the switching functions.
    cutoff : str, default=None
        Cutoff distance for direct-space electrostatic interaction evaluation.
    exclude : ('none', '1-2', '1-3', '1-4', 'scaled1-4'), str, default=None
        Which pairs of bonded atoms should be excluded from non-bonded interactions
    scaling14 : float, default=None
        Scaling factor for 1-4 electrostatic interations.
    langevin : ('on', 'off'), str, default=None
        Enable the Langevin thermostat.
    langevintemp : float, default=None
        The set point in K for the Langevin thermostat.
    langevindamping : float, default=None
        Langevin damping constant gamma (1/ps)
    pme : ('on', 'off'), str, default=None
        Enable the use of PME for long-range electrostatics.
    pmegridspacing : float, default=None
        The spacing of the PME mesh in 1/A.
    fullelectfrequency : int, default=None
        The frequency in interations between successive calculations of long-range (PME) electrostatics.
    energyfreq : int, default=None
        The frequency with which ACEMD will calculate system energies.
    constraints : ('on', 'off'), str, default=None
        Set to enable positional constraints on specified atoms.
    consref : str, default=None
        Specify a PDB file giving reference positions for constrained atoms.
    constraintscaling : float, default=None
        The harmonic constraint energy function is multiplied by this parameter.
    berendsenpressure : ('on', 'off'), str, default=None
        Set to enable the Berendsen pressure bath barostatic control.
    berendsenpressuretarget : float, default=None
        The target pressure (Bar) for the barostat.
    berendsenpressurerelaxationtime : float, default=None
        Relaxation time for the barostat (fs).
    tclforces : ('on', 'off'), str, default=None
        Enable TCL force scripting.
    minimize : int, default=None
        Number of steps of conjugate-gradient minimisation to perform.
    run : str, default=None
        The number of simulation iterations to perform.
    celldimension : str, default=None
        Dimensions of the unit cell.
    useconstantratio : ('on', 'off'), str, default=None
        Keep the ratio of the X-Y dimensions constant while allowing Z to fluctuate independently.
    amber : ('on', 'off'), str, default=None
        Indicate that the Amber force field is to be used.
    dielectric : float, default=None
        Dielectric constant.
    pairlistdist : str, default=None
        Specify the buffer size for grid cells.
    TCL : str, default=None

    bincoordinates : str, default=None
        Filename for initial structure coordinates, in NAMD Bincoor format.
    binvelocities : str, default=None
        Initial velocity field, in NAMD Bincoor format.
    binindex : str, default=None
        Filename for index file to set initial timestep (as made by a check-point)
    structure : str, default=None
        CHARMM structure topology in PSF format
    parameters : str, default=None
        CHARMM force-field parameter file (PRM)
    extendedsystem : str, default=None
        If set, specifies an extended system .xsc file, from which a cell dimension will be read.
    coordinates : str, default=None
        Filename for initial structure coordinates, in PDB format.
    velocities : str, default=None
        Initial velocity field, in PDB format.
    parmfile : str, default=None
        The filename of the Amber PRMTOP parameter file.

    """
    _defaultfnames = {'bincoordinates': 'input.coor', 'binvelocities': 'input.vel', 'binindex': 'input.idx',
                      'structure': 'structure.*', 'parameters': 'parameters', 'extendedsystem': 'input.xsc',
                      'coordinates': 'structure.pdb', 'velocities': 'velocity.pdb', 'consref': 'structure.pdb',
                      'parmfile': 'parameters'}


    def __init__(self, version=3):
        super().__init__()
        self._version = version
        self._file_data = {}
        self._outnames = {}

        if version == 3:
            # ACEMD3 Options
            self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin.', None, val.Number(float, 'ANY'))
            self._arg('restart', 'str', 'Restart simulation.', None, val.String())
            self._arg('trajectoryfile', 'str', 'Output file name.', None, val.String())
            self._arg('trajectoryfreq', 'int', 'Trajectory sampling frequency in steps.', None, val.Number(int, 'POS'))
            self._arg('timestep', 'int', 'Simulation timestep.', None, val.Number(int, 'POS'))
            self._arg('pme', 'str', 'Particle-mesh Ewald summation.', None, val.String())
            self._arg('switching', 'str', 'Apply switching function to the van der Waals potential.', None, val.String())
            self._arg('switchdist', 'float', 'Distance in Angstrom at which to begin applying the switching function.', None, val.Number(float, '0POS'))
            self._arg('cutoff', 'float', 'Real-space cutoff in Angstroms for electrostatics and van der Waals.', None, val.Number(float, '0POS'))
            self._arg('thermostat', 'str', 'Enable thermostatic control', None, val.String())
            self._arg('thermostattemp', 'float', 'Target temperature (K) for thermostatic control', None, val.Number(float, '0POS'))
            self._arg('thermostatdamping', 'float', 'Damping constant for the Langevin thermostat in ps^-1', None, val.Number(float, '0POS'))
            self._arg('restraints', 'str', 'Restraining potentials', None, val.Object(_Restraint), nargs='*')
            self._arg('barostat', 'str', 'Enable pressure control', None, val.String())
            self._arg('barostatpressure', 'float', 'The target pressure in bar', None, val.Number(float, '0POS'))
            self._arg('useflexiblecell', 'str', 'Allow X, Y and Z unit cell dimensions to vary independently', None, val.String())
            self._arg('useconstantarea', 'str', 'Constrain the X,Y dimensions of the unit cell. Allow Z to vary independently.', None, val.String())
            self._arg('useconstantratio', 'str', 'Constrain the X:Y ratio of the unit cell dimensions. Allow Z to vary independently.', None, val.String())
            self._arg('minimize', 'int', 'The number of energy minimization steps to perform before commencing dynamics.', None, val.Number(int, '0POS'))
            self._arg('run', 'str', 'The length of simulation ro run. May be specified as a number of steps or as a time if one of the suffices "us", "ns", "ps", "fs" is used.', None, val.String())
            self._arg('celldimension', 'str', 'The dimensions of the unit cell in Angstrom. Note that the unit cell must be cuboid. Overrides any dimension given in the "coordinates" PDB.', None, val.String())
            self._arg('implicit', 'str', 'Set to True to enable implicit solvent simulations in AMBER.', None, val.String())

            # Files
            self._arg('bincoordinates', 'str', 'Optional initial system geometry in NAMD BINCOOR format. If specified, overrides "coordinates"', None, val.String())
            self._arg('binvelocities', 'str', 'Optional initial system velocity field in NAMD BINCOOR format. If specified, overrides field generated by "temperature" and "velocities"', None, val.String())
            self._arg('structure', 'str', 'The filename of a CHARMM PSF file', None, val.String())
            self._arg('parameters', 'str', 'The filename of a CHARMM PAR file', None, val.String())
            self._arg('parmfile', 'str', 'The filename of an Amber PRMTOP file', None, val.String())
            self._arg('extendedsystem', 'str', 'Filename of a NAMD XSC format file giving the periodic cell dimensions. Overrides "celldimension" and any dimensions in the "coordinates" PDB', None, val.String())
            self._arg('coordinates', 'str', 'Mandatory initial system geometry in PDB format', None, val.String())
            self._arg('velocities', 'str', 'Optional initial system velocity field in NAMD BINCOOR format. If specified, overrides field generated by "temperature"', None, val.String())


        if version == 2:
            logger.warning('The default Acemd engine for HTMD has changed to Acemd3. You can keep on using Acemd2 but will not be able to use new features introduced in Acemd3')
            # ACEMD2 Options
            self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin.', None, val.Number(float, '0POS'))
            self._arg('restart', 'str', 'Restart simulation.', None, val.String(), valid_values=('on', 'off'))
            self._arg('restartfreq', 'int', 'Restart file frequency.', None, val.Number(int, '0POS'))
            self._arg('outputname', 'str', 'Output file name.', None, val.String())
            self._arg('xtcfile', 'str', 'Output XTC file name.', None, val.String())
            self._arg('xtcfreq', 'int', 'XTC sampling frequency in steps.', None, val.Number(int, '0POS'))
            self._arg('timestep', 'float', 'Simulation timestep.', None, val.Number(float, '0POS'))
            self._arg('rigidbonds', 'str', 'Enable holonomic constraints on all hydrogen-heavy atom bond terms', None, val.String(), valid_values=('none', 'all'))
            self._arg('hydrogenscale', 'float', 'Amount by which to scale the mass of H atoms', None, val.Number(float, '0POS'))
            self._arg('switching', 'str', 'Enable to apply smoothing and switching functions to electrostatic and VdW forces.', None, val.String(), valid_values=('on', 'off'))
            self._arg('switchdist', 'float', 'Range beyond which to begin to apply the switching functions.', None, val.Number(float, '0POS'))
            self._arg('cutoff', 'str', 'Cutoff distance for direct-space electrostatic interaction evaluation.', None, val.Number(float, '0POS'))
            self._arg('exclude', 'str', 'Which pairs of bonded atoms should be excluded from non-bonded interactions', None, val.String(), valid_values=('none', '1-2', '1-3', '1-4', 'scaled1-4'))
            self._arg('scaling14', 'float', 'Scaling factor for 1-4 electrostatic interations.', None, val.Number(float, '0POS'))
            self._arg('langevin', 'str', 'Enable the Langevin thermostat.', None, val.String(), valid_values=('on', 'off'))
            self._arg('langevintemp', 'float', 'The set point in K for the Langevin thermostat.', None, val.Number(float, '0POS'))
            self._arg('langevindamping', 'float', 'Langevin damping constant gamma (1/ps)', None, val.Number(float, '0POS'))
            self._arg('pme', 'str', 'Enable the use of PME for long-range electrostatics.', None, val.String(), valid_values=('on', 'off'))
            self._arg('pmegridspacing', 'float', 'The spacing of the PME mesh in 1/A.', None, val.Number(float, '0POS'))
            self._arg('fullelectfrequency', 'int', 'The frequency in interations between successive calculations of long-range (PME) electrostatics.', None, val.Number(int, '0POS'))
            self._arg('energyfreq', 'int', 'The frequency with which ACEMD will calculate system energies.', None, val.Number(int, '0POS'))
            self._arg('constraints', 'str', 'Set to enable positional constraints on specified atoms.', None, val.String(), valid_values=('on', 'off'))
            self._arg('consref', 'str', 'Specify a PDB file giving reference positions for constrained atoms.', None, val.String())
            self._arg('constraintscaling', 'float', 'The harmonic constraint energy function is multiplied by this parameter.', None, val.Number(float, 'ANY'))
            self._arg('berendsenpressure', 'str', 'Set to enable the Berendsen pressure bath barostatic control.', None, val.String(), valid_values=('on', 'off'))
            self._arg('berendsenpressuretarget', 'float', 'The target pressure (Bar) for the barostat.', None, val.Number(float, '0POS'))
            self._arg('berendsenpressurerelaxationtime', 'float', 'Relaxation time for the barostat (fs).', None, val.Number(float, '0POS'))
            self._arg('tclforces', 'str', 'Enable TCL force scripting.', None, val.String(), valid_values=('on', 'off'))
            self._arg('minimize', 'int', 'Number of steps of conjugate-gradient minimisation to perform.', None, val.Number(int, '0POS'))
            self._arg('run', 'str', 'The number of simulation iterations to perform.', None)
            self._arg('celldimension', 'str', 'Dimensions of the unit cell.', None, val.Number(float, 'ANY'), nargs=3)
            self._arg('useconstantratio', 'str', 'Keep the ratio of the X-Y dimensions constant while allowing Z to fluctuate independently.', None, val.String(), valid_values=('on', 'off'))
            self._arg('amber', 'str', 'Indicate that the Amber force field is to be used.', None, val.String(), valid_values=('on', 'off'))
            self._arg('dielectric', 'float', 'Dielectric constant.', None, val.Number(float, 'ANY'))
            self._arg('pairlistdist', 'str', 'Specify the buffer size for grid cells.', None, val.Number(float, 'ANY'))
            self._arg('TCL', 'str', '', None, val.String(), nargs='*')

            # Files
            self._arg('bincoordinates', 'str', 'Filename for initial structure coordinates, in NAMD Bincoor format.', None, val.String())
            self._arg('binvelocities', 'str', 'Initial velocity field, in NAMD Bincoor format.', None, val.String())
            self._arg('binindex', 'str', 'Filename for index file to set initial timestep (as made by a check-point)', None, val.String())
            self._arg('structure', 'str', 'CHARMM structure topology in PSF format', None, val.String())
            self._arg('parameters', 'str', 'CHARMM force-field parameter file (PRM)', None, val.String())
            self._arg('extendedsystem', 'str', 'If set, specifies an extended system .xsc file, from which a cell dimension will be read.', None, val.String())
            self._arg('coordinates', 'str', 'Filename for initial structure coordinates, in PDB format.', None, val.String())
            self._arg('velocities', 'str', 'Initial velocity field, in PDB format.', None, val.String())
            self._arg('parmfile', 'str', 'The filename of the Amber PRMTOP parameter file.', None, val.String())


    def load(self, path='.'):
        """ Loads all files required to run a simulation and apply eventually configured protocols to it

        Parameters
        ----------
        path : str
            Working directory relative to which the configuration file is read
        """
        # load files and reset filenames
        for cmd in self._defaultfnames.keys():
            if cmd in self.__dict__ and self.__dict__[cmd] is not None:
                fname = os.path.join(path, self.__dict__[cmd])
                f = open(fname, 'rb')  # read all as binary
                data = f.read()
                f.close()
                self._file_data[cmd] = data

                defaultname = self._defaultfnames[cmd]
                if defaultname.endswith('*'):
                    defaultname = '{}.{}'.format(os.path.splitext(defaultname)[0], os.path.splitext(fname)[1][1:])
                self._outnames[cmd] = defaultname


    def save(self, path, overwrite=False):
        """ Create a directory with all necessary input to run acemd.

        Parameters
        ----------
        path : str
            New execution directory to create
        overwrite : bool
            Overwrite output directory if it exists.
        """
        if os.path.exists(path):
            if overwrite:
                import shutil
                shutil.rmtree(path)
            else:
                raise RuntimeError('Directory exists, use overwrite=True or remove directory')
        os.makedirs(path)

        # Write all files using default filenames
        for f in self._file_data:
            fo = open(os.path.join(path, self._outnames[f]), 'wb')  # write all as binary
            fo.write(self._file_data[f])
            fo.close()

        self._writeBashRun(os.path.join(path, 'run.sh'))

        self._file_data = {}  # empty file dictionary after writing

        self.writeConf(os.path.join(path, 'input'))


    def setup(self, indir='./', outdir='./run', overwrite=False):
        """ Convenience method performing load and save.

        Parameters
        ----------
        indir : str
            The directory to load from.
        outdir : str
            The directory to save to.
        overwrite : bool
            Overwrite output directory if it exists.
        """
        self.load(indir)
        self.save(outdir, overwrite)


    def show(self, quiet=False):
        """ Returns the Acemd configuration file string

        Parameters
        ----------
        quiet : bool
            If true it prints the string to stdout

        Returns
        -------
        conf : str
            The string of the configuration file
        """
        text = ''
        if 'TCL' in self.__dict__ and self.__dict__['TCL'] is not None:
            text = self.__dict__['TCL']

        text += '#\n'

        maxwidth = np.max([len(k) for k in self.__dict__.keys()])

        keys = sorted(list(self.__dict__.keys()))
        keys = keys + [keys.pop(keys.index('restraints')), keys.pop(keys.index('run'))]  # Move the run command to the end
        for cmd in keys:
            if cmd == 'restraints' and self.restraints is not None:
                for r in ensurelist(self.restraints):
                    text += '{}\n'.format(r.format(maxwidth))
            elif not cmd.startswith('_') and self.__dict__[cmd] is not None and cmd != 'TCL':
                val = self.__dict__[cmd]
                if cmd in self._outnames:
                    val = self._outnames[cmd]
                name = cmd
                if cmd == 'scaling14':  # variables cannot start with numbers. We need to rename it here for acemd
                    name = '1-4scaling'
                text += '{name: <{maxwidth}}\t{val:}\n'.format(name=name, val=val, maxwidth=maxwidth)

        if not quiet:
            print(text)
        else:
            return text


    def _writeBashRun(self, fname):
        with open(fname, 'w') as f:
            if self._version == 3:
                f.write('#!/bin/bash\nacemd3 >log.txt 2>&1')
            else:
                f.write('#!/bin/bash\nacemd >log.txt 2>&1')
        os.chmod(fname, 0o700)

    def __repr__(self):
        return self.show(quiet=True)


    def writeConf(self, fname='input'):
        """ Write an acemd configuration file

        Parameters
        ----------
        fname : output file name
        """
        text = self.show(quiet=True)
        fo = open(fname, 'w')
        fo.write(text)
        fo.close()


from unittest import TestCase
class TestAcemd(TestCase):
    def test_acemd2(self):
        from htmd.home import home
        from htmd.util import tempname
        import filecmp
        from glob import glob
        import sys

        tmpdir = tempname()

        acemd = Acemd()
        acemd.temperature = 350
        acemd.restart = 'on'
        acemd.restartfreq = 25
        acemd.outputname = 'output'
        acemd.xtcfile = 'myout.xtc'
        acemd.xtcfreq = 500
        acemd.timestep = 4
        acemd.rigidbonds = 'all'
        acemd.hydrogenscale = 0.3
        acemd.switching = 'on'
        acemd.switchdist = 9
        acemd.cutoff = 5
        acemd.exclude = 'scaled1-4'
        acemd.scaling14 = 5.6
        acemd.langevin = 'on'
        acemd.langevintemp = 300
        acemd.langevindamping = 0.8
        acemd.pme = 'on'
        acemd.pmegridspacing = 5
        acemd.fullelectfrequency = 3
        acemd.energyfreq = 10
        acemd.constraints = 'on'
        acemd.consref = '5dhfr_cube.pdb'
        acemd.constraintscaling = 3.2
        acemd.berendsenpressure = 'on'
        acemd.berendsenpressuretarget = 14
        acemd.berendsenpressurerelaxationtime = 2
        acemd.tclforces = 'on'
        acemd.minimize = 150
        acemd.run = 555
        acemd.celldimension = [3, 56, 2]
        acemd.useconstantratio = 'on'
        acemd.amber = 'off'
        acemd.dielectric = 16
        acemd.pairlistdist = 9
        acemd.TCL = 'on'
        acemd.bincoordinates = None
        acemd.binvelocities = None
        acemd.binindex = None
        acemd.structure = '5dhfr_cube.psf'
        acemd.parameters = 'par_all22_prot.inp'
        acemd.extendedsystem = None
        acemd.coordinates = '5dhfr_cube.pdb'
        acemd.velocities = None
        acemd.parmfile = None

        acemd.coordinates = '5dhfr_cube.pdb'
        acemd.setup(home(dataDir='dhfr/'), tmpdir, overwrite=True)

        # Compare with reference
        refdir = home(dataDir='test-acemd-v2')
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, '*'))]
        match, mismatch, error = filecmp.cmpfiles(refdir, tmpdir, files, shallow=False)

        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError('Different results produced by Acemd.write between {} and {} '
                               'in files {}.'.format(refdir, tmpdir, mismatch))

        import shutil
        shutil.rmtree(tmpdir)

    def test_acemd3(self):
        import htmd.home
        homedir = htmd.home.home()

        r = list()
        r.append(GroupRestraint('resname MOL', 5, [(10, '10ns'), (5, '15ns'), (0, '20ns')], axes='z'))
        r.append(AtomRestraint('name CA', 0.1, [(10, '10ns'), (5, '15ns'), (0, '20ns')]))

        acemd = Acemd(version=3)
        acemd.structure = '5dhfr_cube.psf'
        acemd.parameters = 'par_all22_prot.inp'
        acemd.coordinates = '5dhfr_cube.pdb'
        acemd.restraints = r
        acemd.setup(homedir + '/data/dhfr', '/tmp/testdir', overwrite=True)
        print(acemd)


if __name__ == "__main__":
    import unittest
    unittest.main(verbosity=2)

