# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import htmd
import shutil
import numpy as np
from htmd.protocols.oldprotocolinterface import ProtocolInterface as OldProtocolInterface
from protocolinterface import ProtocolInterface, val


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
                      'parmfile': 'structure.prmtop'}

    def __init__(self, version=2):
        super().__init__()
        self._version = version
        self._files = {}
        self._outnames = {}

        # Options
        self._arg('temperature', 'float', 'Temperature of the thermostat in Kelvin.', None, val.Number(float, '0POS'))
        self._arg('restart', 'str', 'Restart simulation.', None, val.String(), valid_values=('on', 'off'))
        self._arg('restartfreq', 'int', 'Restart file frequency.', None, val.Number(int, '0POS'))
        self._arg('outputname', 'str', 'Output file name.', None, val.String())
        self._arg('xtcfile', 'str', 'Output XTC file name.', None, val.String())
        self._arg('xtcfreq', 'int', 'XTC sampling frequency in steps.', None, val.Number(int, '0POS'))
        self._arg('timestep', 'float', 'Simulation timestep.', None, val.Number(float, '0POS'))
        self._arg('rigidbonds', 'str', 'Enable holonomic constraints on all hydrogen-heavy atom bond terms',
                  None, val.String(), valid_values=('none', 'all'))
        self._arg('hydrogenscale', 'float', 'Amount by which to scale the mass of H atoms',
                  None, val.Number(float, '0POS'))
        self._arg('switching', 'str',
                  'Enable to apply smoothing and switching functions to electrostatic and VdW forces.',
                  None, val.String(), valid_values=('on', 'off'))
        self._arg('switchdist', 'float', 'Range beyond which to begin to apply the switching functions.',
                  None, val.Number(float, '0POS'))
        self._arg('cutoff', 'str', 'Cutoff distance for direct-space electrostatic interaction evaluation.',
                  None, val.Number(float, '0POS'))
        self._arg('exclude', 'str', 'Which pairs of bonded atoms should be excluded from non-bonded interactions',
                  None, val.String(), valid_values=('none', '1-2', '1-3', '1-4', 'scaled1-4'))
        self._arg('scaling14', 'float', 'Scaling factor for 1-4 electrostatic interations.',
                  None, val.Number(float, '0POS'))
        self._arg('langevin', 'str', 'Enable the Langevin thermostat.', None, val.String(), valid_values=('on', 'off'))
        self._arg('langevintemp', 'float', 'The set point in K for the Langevin thermostat.',
                  None, val.Number(float, '0POS'))
        self._arg('langevindamping', 'float', 'Langevin damping constant gamma (1/ps)',
                  None, val.Number(float, '0POS'))
        self._arg('pme', 'str', 'Enable the use of PME for long-range electrostatics.',
                  None, val.String(), valid_values=('on', 'off'))
        self._arg('pmegridspacing', 'float', 'The spacing of the PME mesh in 1/A.', None, val.Number(float, '0POS'))
        self._arg('fullelectfrequency', 'int',
                  'The frequency in interations between successive calculations of long-range (PME) electrostatics.',
                  None, val.Number(int, '0POS'))
        self._arg('energyfreq', 'int', 'The frequency with which ACEMD will calculate system energies.',
                  None, val.Number(int, '0POS'))
        if self._version == 2:
            self._arg('constraints', 'str', 'Set to enable positional constraints on specified atoms.',
                      None, val.String(), valid_values=('on', 'off'))
            self._arg('consref', 'str', 'Specify a PDB file giving reference positions for constrained atoms.',
                      None, val.String())
            self._arg('constraintscaling', 'float',
                      'The harmonic constraint energy function is multiplied by this parameter.',
                      None, val.Number(float, 'ANY'))
        if self._version == 3:
            self._arg('atomrestraint', 'str', '', None)
            self._arg('grouprestraint', 'str', '', None)
        self._arg('berendsenpressure', 'str', 'Set to enable the Berendsen pressure bath barostatic control.',
                  None, val.String(), valid_values=('on', 'off'))
        self._arg('berendsenpressuretarget', 'float', 'The target pressure (Bar) for the barostat.',
                  None, val.Number(float, '0POS'))
        self._arg('berendsenpressurerelaxationtime', 'float', 'Relaxation time for the barostat (fs).',
                  None, val.Number(float, '0POS'))
        self._arg('tclforces', 'str', 'Enable TCL force scripting.', None, val.String(), valid_values=('on', 'off'))
        self._arg('minimize', 'int', 'Number of steps of conjugate-gradient minimisation to perform.',
                  None, val.Number(int, '0POS'))
        self._arg('run', 'str', 'The number of simulation iterations to perform.', None)
        self._arg('celldimension', 'str', 'Dimensions of the unit cell.', None, val.Number(float, 'ANY'), nargs=3)
        self._arg('useconstantratio', 'str',
                  'Keep the ratio of the X-Y dimensions constant while allowing Z to fluctuate independently.',
                  None, val.String(), valid_values=('on', 'off'))
        self._arg('amber', 'str', 'Indicate that the Amber force field is to be used.',
                  None, val.String(), valid_values=('on', 'off'))
        self._arg('dielectric', 'float', 'Dielectric constant.', None, val.Number(float, 'ANY'))
        self._arg('pairlistdist', 'str', 'Specify the buffer size for grid cells.', None, val.Number(float, 'ANY'))
        if self._version == 2:
            self._arg('TCL', 'str', '', None, val.String(), nargs='*')

        # Files
        self._arg('bincoordinates', 'str', 'Filename for initial structure coordinates, in NAMD Bincoor format.',
                  None, val.String())
        self._arg('binvelocities', 'str', 'Initial velocity field, in NAMD Bincoor format.',
                  None, val.String())
        self._arg('binindex', 'str', 'Filename for index file to set initial timestep (as made by a check-point)',
                  None, val.String())
        self._arg('structure', 'str', 'CHARMM structure topology in PSF format', None, val.String())
        self._arg('parameters', 'str', 'CHARMM force-field parameter file (PRM)', None, val.String())
        self._arg('extendedsystem', 'str',
                  'If set, specifies an extended system .xsc file, from which a cell dimension will be read.',
                  None, val.String())
        self._arg('coordinates', 'str', 'Filename for initial structure coordinates, in PDB format.',
                  None, val.String())
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
                self._files[cmd] = data

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
                shutil.rmtree(path)
            else:
                raise NameError('Directory exists, use overwrite=True or remove directory')
        os.makedirs(path)

        # Write all files using default filenames
        for f in self._files:
            fo = open(os.path.join(path, self._outnames[f]), 'wb')  # write all as binary
            fo.write(self._files[f])
            fo.close()

        self._writeBashRun(os.path.join(path, 'run.sh'))

        self._files = {}  # empty file dictionary after writing

        self.writeConf(os.path.join(path, 'input'))

    def setup(self, cwd='./', outdir='./run', overwrite=False):
        """ Convenience method performing load and save.

        Parameters
        ----------
        cwd : str
            The directory to load from.
        outdir : str
            The directory to save to.
        overwrite : bool
            Overwrite output directory if it exists.
        """
        self.load(cwd)
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
        if self.__dict__['TCL'] is not None:
            text = self.__dict__['TCL']

        text += '#\n'

        maxwidth = np.max([len(k) for k in self.__dict__.keys()])

        keys = sorted(list(self.__dict__.keys()))
        keys = keys + [keys.pop(keys.index('run'))]  # Move the run command to the end
        for cmd in keys:
            if not cmd.startswith('_') and self.__dict__[cmd] is not None and cmd != 'TCL':
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


if __name__ == "__main__":
    # l=Acemd.protocols(quiet=True)
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
