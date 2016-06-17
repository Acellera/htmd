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

        #  General flags describing the calculation (Manual section 18.6.1)
        self._cmdList(key='imin', datatype='int',
                      descr='Flag to run minimization.',
                      default=0, valid_values=[0, 1, 5])
        self._cmdList(key='nmropt', datatype='int',
                      descr='Flag to apply NMR restraints',
                      default=0, valid_values=[0, 1, 2])

        #  Nature and format of the input (Manual section 18.6.2)
        self._cmdList(key='ntx', datatype='int',
                      descr="""Option to read the initial coordinates,
                            velocities and box size from the inpcrd file.""",
                      default=1, valid_values=[1, 5])
        self._cmdList(key='irest', datatype='int',
                      descr='Flag to restart a simulation.',
                      default=0, valid_values=[0, 1])

        #  Nature and format of the output (Manual section 18.6.3)
        self._cmdList(key='ntxo', datatype='int',
                      descr="""Format of the final coordinates, velocities,
                            and box size (if constant volume or pressure run)
                            written to file "restrt".""",
                      default=2, valid_values=[1, 2])
        self._cmdValue(key='ntpr', datatype='int',
                       descr="""Every ntpr steps, energy information will be
                             printed in human-readable form to files "mdout"
                             and "mdinfo".""",
                       default=50)
        self._cmdValue(key='ntave', datatype='int',
                       descr="""Every ntave steps of dynamics, running averages
                             of average energies and fluctuations over the last
                             ntave steps will be printed out.""",
                       default=0)
        self._cmdValue(key='ntwr', datatype='int',
                       descr="""Every ntwr steps during dynamics, the “restrt”
                             file will be written.""",
                       default=None)
        self._cmdList(key='iwrap', datatype='int',
                      descr="""If iwrap = 1, the coordinates written to the
                            restart and trajectory files will be "wrapped"
                            into a primary box. If iwrap = 0, no wrapping will
                            be performed.""",
                      default=0, valid_values=[0, 1])
        self._cmdValue(key='ntwx', datatype='int',
                       descr="""Every ntwx steps, the coordinates will be
                             written to the mdcrd file. If ntwx = 0, no
                             coordinate trajectory file will be written.""",
                       default=0)
        self._cmdValue(key='ntwv', datatype='int',
                       descr="""Every ntwv steps, the velocities will be written
                             to the mdvel file. If ntwv = 0, no velocity
                             trajectory file will be written. If ntwv = -1,
                             velocities will be written to mdcrd, which then
                             becomes a combined coordinate/velocity trajectory
                             file, at the interval defined by ntwx.""",
                       default=0)
        self._cmdValue(key='ntwf', datatype='int',
                       descr="""Every ntwf steps, the forces will be written to
                             the mdfrc file. If ntwf = 0, no force trajectory
                             file will be written. If ntwf = -1, forces will be
                             written to the mdcrd, which then becomes a combined
                             coordinate/force trajectory file, at the interval
                             defined by ntwx.""",
                       default=0)
        self._cmdValue(key='ntwe', datatype='int',
                       descr="""Every ntwe steps, the energies and temperatures
                             will be written to file "mden" in a compact form.
                             If ntwe = 0 then no mden file will be written.""",
                       default=0)
        self._cmdList(key='ioutfm', datatype='int',
                      descr="""The format of coordinate and velocity trajectory
                            files (mdcrd, mdvel and inptraj). 0 is ASCII and 1
                            is Binary NetCDF""",
                      default=1, valid_values=[0, 1])
        self._cmdValue(key='ntwprt', datatype='int',
                       descr="""The number of atoms to include in trajectory
                             files (mdcrd and mdvel). If ntwprt = 0, all atoms
                             will be included.""",
                       default=0)
        self._cmdList(key='idecomp', datatype='int',
                      descr="""Perform energy decomposition according to a
                            chosen scheme.""",
                      default=0, valid_values=[0, 1, 2, 3, 4])




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
        with open(fname, 'w') as f:
            f.write(text)
