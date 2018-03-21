# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import htmd
import shutil
import numpy as np
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY


class Acemd(ProtocolInterface):
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
        self._cmdString('temperature', 'float', 'Temperature of the thermostat in Kelvin.', None)
        self._cmdString('restart', 'str', 'Restart simulation.', None)
        self._cmdString('restartfreq', 'str', 'Restart file frequency.', None)
        self._cmdString('outputname', 'str', 'Output file name.', None)
        self._cmdString('xtcfile', 'str', 'Output XTC file name.', None)
        self._cmdString('xtcfreq', 'str', 'XTC sampling frequency in steps.', None)
        self._cmdString('timestep', 'str', 'Simulation timestep.', None)
        self._cmdString('rigidbonds', 'str', '', None)
        self._cmdString('hydrogenscale', 'str', '', None)
        self._cmdString('switching', 'str', '', None)
        self._cmdString('switchdist', 'str', '', None)
        self._cmdString('cutoff', 'str', '', None)
        self._cmdString('exclude', 'str', '', None)
        self._cmdString('scaling14', 'str', '', None)
        self._cmdString('langevin', 'str', '', None)
        self._cmdString('langevintemp', 'str', '', None)
        self._cmdString('langevindamping', 'str', '', None)
        self._cmdString('pme', 'str', '', None)
        self._cmdString('pmegridspacing', 'str', '', None)
        self._cmdString('fullelectfrequency', 'str', '', None)
        self._cmdString('energyfreq', 'str', '', None)
        if self._version == 2:
            self._cmdString('constraints', 'str', '', None)
            self._cmdString('consref', 'str', '', None)
            self._cmdString('constraintscaling', 'str', '', None)
        if self._version == 3:
            self._cmdString('atomrestraint', 'str', '', None)
            self._cmdString('grouprestraint', 'str', '', None)
        self._cmdString('berendsenpressure', 'str', '', None)
        self._cmdString('berendsenpressuretarget', 'str', '', None)
        self._cmdString('berendsenpressurerelaxationtime', 'str', '', None)
        self._cmdString('tclforces', 'str', '', None)
        self._cmdString('minimize', 'str', '', None)
        self._cmdString('run', 'str', '', None)
        self._cmdString('celldimension', 'str', '', None)
        self._cmdString('useconstantratio', 'str', '', None)
        self._cmdString('amber', 'str', '', None)
        self._cmdString('dielectric', 'str', '', None)
        self._cmdString('pairlistdist', 'str', '', None)
        if self._version == 2:
            self._cmdString('TCL', 'str', '', None)

        # Files
        self._cmdString('bincoordinates', 'str', '', None)
        self._cmdString('binvelocities', 'str', '', None)
        self._cmdString('binindex', 'str', '', None)
        self._cmdString('structure', 'str', '', None)
        self._cmdString('parameters', 'str', '', None)
        self._cmdString('extendedsystem', 'str', '', None)
        self._cmdString('coordinates', 'str', '', None)
        self._cmdString('velocities', 'str', '', None)
        self._cmdString('parmfile', 'str', '', None)

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
    #l=Acemd.protocols(quiet=True)
    import htmd.home
    acemd = Acemd()
    acemd.structure = '5dhfr_cube.psf'
    acemd.parameters = 'par_all22_prot.inp'
    homedir = htmd.home.home()
    acemd.coordinates = '5dhfr_cube.pdb'
    acemd.setup(homedir + '/data/dhfr', '/tmp/testdir', overwrite=True)
    print(acemd)

