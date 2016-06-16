# Code to setup simulations with Amber
# Description of parameters is mainly being copied from the AMBER 15 manual
# http://ambermd.org/doc12/Amber15.pdf
# However, some parameter options have different names from Amber to 
# to keep protocol creation as consistent as possible with acemd

import os
import htmd
import shutil
import numpy as np
from htmd.protocols.protocolinterface import ProtocolInterface #TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY

class Amber(ProtocolInterface):
    _defaultfnames = {'bincoordinates': 'input.nc', 'structure': 'structure.*', 
                      'parameters': 'parameters', 'coordinates': 'structure.rst',
                      'parmfile': 'structure.prmtop'}

    def __init__(self):
        super().__init__()

        self._files = {}

        # Options
        self._cmdString('minimisation', 'int', 'Whether to perform energy minimisation', None)
        self._cmdString('maxcycle', 'int', 'Maximum number of cycles of minimisation to perform', None)
        self._cmdString('ncyc', 'int', 'Cycle at which to switch from steepest to conjugate descent')
        self._cmdString('cutoff', 'float', '', None)
        self._cmdString('constraints', 'str', '', None)
        self._cmdString('printoutnsteps', 'int','Number of steps between to mdout and mdinfo print out', None)
        self._cmdString('initcoordread', 'int', 'Option to read the initial coordinates, velocities and box size from the inpcrd file', None)
        self._cmdString('finalcoordwrite', 'int', 'Format of the final coordinates, velocities, and box size', None)
        self._cmdString('restart', 'int', 'Flag to restart a simulation.', None)
        self._cmdString('ncfreq', 'int', 'NetCFD sampling frequency in steps.', None)
        self._cmdString('coorfiletype', 'int', 'The format of coordinate and velocity trajectory files', None)

        self._cmdString('temperature', 'float', 'Temperature of the thermostat in Kelvin.', None)
        self._cmdString('restartfreq', 'str', 'Restart file frequency.', None)
        self._cmdString('outputname', 'str', 'Output file name.', None)
        self._cmdString('ncfile', 'str', 'Output NetCFD file name.', None)
        self._cmdString('timestep', 'str', 'Simulation timestep.', None)
        self._cmdString('rigidbonds', 'str', '', None)
        self._cmdString('hydrogenscale', 'str', '', None)
        self._cmdString('switching', 'str', '', None)
        self._cmdString('switchdist', 'str', '', None)
        self._cmdString('exclude', 'str', '', None)
        self._cmdString('langevin', 'str', '', None)
        self._cmdString('langevintemp', 'str', '', None)
        self._cmdString('langevindamping', 'str', '', None)
        self._cmdString('pme', 'str', '', None)
        self._cmdString('pmegridspacing', 'str', '', None)
        self._cmdString('fullelectfrequency', 'str', '', None)
        self._cmdString('energyfreq', 'str', '', None)
        self._cmdString('consref', 'str', '', None)
        self._cmdString('constraintscaling', 'str', '', None)
        self._cmdString('berendsenpressure', 'str', '', None)
        self._cmdString('berendsenpressuretarget', 'str', '', None)
        self._cmdString('berendsenpressurerelaxationtime', 'str', '', None)
        self._cmdString('run', 'str', '', None)
        self._cmdString('celldimension', 'str', '', None)
        self._cmdString('useconstantratio', 'str', '', None)

        # Files
        self._cmdString('bincoordinates', 'str', '', None)
        self._cmdString('structure', 'str', '', None)
        self._cmdString('parameters', 'str', '', None)
        self._cmdString('coordinates', 'str', '', None)
        self._cmdString('consref', 'str', '', None)
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
                f = open(fname, 'rb')  #read all as binary
                data = f.read()
                f.close()
                self._files[cmd] = data

                defaultname = self._defaultfnames[cmd]
                if defaultname.endswith('*'):
                    defaultname = '{}.{}'.format(os.path.splitext(defaultname)[0], os.path.splitext(fname)[1][1:])

                self.__dict__[cmd] = defaultname  # use default file names

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
                name = cmd
                if cmd == 'scaling14':  # variables cannot start with numbers. We need to rename it here for acemd
                    name = '1-4scaling'
                text += '{name: <{maxwidth}}\t{val:}\n'.format(name=name, val=str(self.__dict__[cmd]), maxwidth=maxwidth)

        if not quiet:
            print(text)
        else:
            return text

    def __str__(self):
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