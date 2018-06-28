# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Author: Juan Eiros Zamora (je714@ic.ac.uk) & Gil Ferreira Hoben (gf712@ic.ac.uk) Organisation: Imperial College London
# Modified by: Stefan Doerr

# Code to setup simulations with Amber
# Description of parameters is mainly being copied from the AMBER 15 manual
# http://ambermd.org/doc12/Amber15.pdf
# However, some parameter options have different names from Amber to
# to keep protocol creation as consistent as possible with acemd

import os
import shutil
import numpy as np
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
from htmd.decorators import _Deprecated


@_Deprecated('1.13.6')
class Pmemd(ProtocolInterface):
    _defaultfnames = {'bincoordinates': 'input.nc',
                      'parameters': 'parameters', 'coordinates': 'structure.rst',
                      'parmfile': 'structure.prmtop'}

    def __init__(self):
        super().__init__()

        self._files = {}

        # GENERAL MINIMIZATION AND DYNAMICS PARAMETERS (MANUAL SECTION 18.6)
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
                      default=2, valid_values=[1, 2, 4, 5, 6])
        self._cmdList(key='irest', datatype='int',
                      descr='Flag to restart a simulation.',
                      default=1, valid_values=[0, 1])

        #  Nature and format of the output (Manual section 18.6.3)
        self._cmdList(key='ntxo', datatype='int',
                      descr="""Format of the final coordinates, velocities,
                            and box size (if constant volume or pressure run)
                            written to file "restrt".""", default=2,
                      valid_values=[1, 2])
        self._cmdValue(key='ntpr', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntpr steps, energy information will be
                             printed in human-readable form to files "mdout"
                             and "mdinfo".""", default=50, valid_range=RANGE_POS)
        self._cmdValue(key='ntave', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntave steps of dynamics, running averages
                             of average energies and fluctuations over the last
                             ntave steps will be printed out.""", default=0,
                       valid_range=RANGE_0POS)
        self._cmdValue(key='ntwr', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntwr steps during dynamics, the “restrt”
                             file will be written.""", default=None, valid_range=RANGE_POS)
        self._cmdList(key='iwrap', datatype='int',
                      descr="""If iwrap = 1, the coordinates written to the
                            restart and trajectory files will be "wrapped"
                            into a primary box. If iwrap = 0, no wrapping will
                            be performed.""", default=0, valid_values=[0, 1])
        self._cmdValue(key='ntwx', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntwx steps, the coordinates will be
                             written to the mdcrd file. If ntwx = 0, no
                             coordinate trajectory file will be written.""",
                       default=0, valid_range=RANGE_0POS)
        self._cmdValue(key='ntwv', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntwv steps, the velocities will be written
                             to the mdvel file. If ntwv = 0, no velocity
                             trajectory file will be written. If ntwv = -1,
                             velocities will be written to mdcrd, which then
                             becomes a combined coordinate/velocity trajectory
                             file, at the interval defined by ntwx.""",
                       default=0, valid_range=RANGE_ANY)
        self._cmdValue(key='ntwf', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntwf steps, the forces will be written to
                             the mdfrc file. If ntwf = 0, no force trajectory
                             file will be written. If ntwf = -1, forces will be
                             written to the mdcrd, which then becomes a combined
                             coordinate/force trajectory file, at the interval
                             defined by ntwx.""", default=0, valid_range=RANGE_ANY)
        self._cmdValue(key='ntwe', datatype='int', realdatatype=TYPE_INT,
                       descr="""Every ntwe steps, the energies and temperatures
                             will be written to file "mden" in a compact form.
                             If ntwe = 0 then no mden file will be written.""",
                       default=0, valid_range=RANGE_0POS)
        self._cmdList(key='ioutfm', datatype='int',
                      descr="""The format of coordinate and velocity trajectory
                            files (mdcrd, mdvel and inptraj). 0 is ASCII and 1
                            is Binary NetCDF""", default=1, valid_values=[0, 1])
        self._cmdValue(key='ntwprt', datatype='int', realdatatype=TYPE_INT,
                       descr="""The number of atoms to include in trajectory
                             files (mdcrd and mdvel). If ntwprt = 0, all atoms
                             will be included.""", default=0, valid_range=RANGE_0POS)
        self._cmdList(key='idecomp', datatype='int',
                      descr="""Perform energy decomposition according to a
                            chosen scheme.""", default=0, valid_values=[0, 1, 2, 3, 4])

        #  Frozen or restrained atoms (Manual section 18.6.4)
        self._cmdList(key='ibelly', datatype='int',
                      descr="""Flag for belly type dynamics. If set to 1, a
                            subset of the atoms in the system will be allowed
                            to move, and the coordinates of the rest will be
                            frozen.""", default=0, valid_values=[0, 1])
        self._cmdValue(key='ntr', datatype='int', realdatatype=TYPE_INT,
                       descr="""Flag for restraining specified atoms in
                             Cartesian space using a harmonic potential,
                             if ntr > 0.""", default=0, valid_range=RANGE_0POS)
        self._cmdValue(key='restraint_wt', datatype='int', realdatatype=TYPE_INT,
                       descr="""The weight (in kcal/mol−A^̊2) for the positional
                             restraints.""", default=None, valid_range=RANGE_0POS)
        self._cmdString(key='restraintmask', datatype='str',
                        descr="""String that specifies the restrained atoms
                              when ntr=1.""", default=None)
        self._cmdString(key='bellymask', datatype='str',
                        descr="""String that specifies the moving atoms when
                        ibelly=1.""", default=None)

        #  Energy minimization (Manual section 18.6.5)
        self._cmdValue(key='maxcyc', datatype='int', realdatatype=TYPE_INT,
                       descr='The maximum number of cycles of minimization.',
                       default=1, valid_range=RANGE_POS)
        self._cmdValue(key='ncyc', datatype='int', realdatatype=TYPE_INT,
                       descr="""If NTMIN is 1 then the method of minimization
                             will be switched from steepest descent to
                             conjugate gradient after NCYC cycles.""",
                       default=10, valid_range=RANGE_POS)
        self._cmdList(key='ntmin', datatype='int',
                      descr='Flag for the method of minimization.',
                      default=1, valid_values=[0, 1, 2, 3, 4])
        self._cmdValue(key='dx0', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='The initial step length.', default=0.01,
                       valid_range=RANGE_POS)
        self._cmdValue(key='drms', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""The convergence criterion for the energy
                             gradient: minimization will halt when the RMS of
                             the Cartesian elements of the gradient
                             is < DRMS.""", default=1e-4, valid_range=RANGE_POS)

        #  Molecular dynamics (Manual section 18.6.6)
        self._cmdValue(key='nstlim', datatype='int', realdatatype=TYPE_INT,
                       descr='Number of MD-steps to be performed.', default=1,
                       valid_range=RANGE_POS)
        self._cmdValue(key='nscm', datatype='int', realdatatype=TYPE_INT,
                       descr="""Flag for the removal of translational and
                             rotational center-of-mass (COM) motion at regular
                             intervals.""", default=1000, valid_range=RANGE_POS)
        self._cmdValue(key='t', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""The time at the start (ps) this is for your own
                             reference and is not critical.""", default=0.0,
                       valid_range=RANGE_0POS)
        self._cmdValue(key='dt', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""The time step (ps).""", default=0.001,
                       valid_range=RANGE_POS)
        self._cmdValue(key='nrespa', datatype='int', realdatatype=TYPE_INT,
                       descr="""This variable allows the user to evaluate
                             slowly-varying terms in the force field less
                             frequently.""", default=None, valid_range=RANGE_POS)

        #  Temperature regulation (Manual section 18.6.7)
        self._cmdList(key='ntt', datatype='int',
                      descr='Switch for temperature scaling.', default=None,
                      valid_values=[0, 1, 2, 3, 9, 10])
        self._cmdValue(key='temp0', datatype='int', realdatatype=TYPE_INT,
                       descr="""Reference temperature at which the system is to
                             be kept, if ntt > 0.""", default=300, valid_range=RANGE_POS)
        self._cmdValue(key='temp0les', datatype='int', realdatatype=TYPE_INT,
                       descr="""This is the target temperature for all LES
                       particles (see Chapter 6 in the Amber Manual).""",
                       default=-1, valid_range=RANGE_ANY)
        self._cmdValue(key='tempi', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""Initial temperature. For the initial dynamics
                             run, (NTX < 3) the velocities are assigned from
                             a Maxwellian distribution at TEMPI K. If
                             TEMPI = 0.0, the velocities will be calculated
                             from the forces instead. TEMPI has no effect if
                             NTX > 3.""", default=0.0, valid_range=RANGE_0POS)
        self._cmdValue(key='ig', datatype='int', realdatatype=TYPE_INT,
                       descr='The seed for the pseudo-random numbergenerator.',
                       default=-1, valid_range=RANGE_ANY)
        self._cmdValue(key='tautp', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""Time constant (ps) for heat bath coupling for
                             the system, if ntt = 1.""", default=1.0,
                       valid_range=RANGE_0POS)
        self._cmdValue(key='gamma_ln', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='The collision frequency γ, in ps−1, when ntt = 3.',
                       default=0.0, valid_range=RANGE_0POS)
        self._cmdValue(key='vrand', datatype='int', realdatatype=TYPE_INT,
                       descr="""If vrand>0 and ntt=2, the velocities will be
                             randomized to temperature TEMP0 every vrand steps.""",
                       default=1000, valid_range=RANGE_0POS)
        self._cmdValue(key='vlimit', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""If not equal to 0.0, then any component of the
                             velocity that is greater than abs(VLIMIT) will be
                             reduced to VLIMIT (preserving the sign).""",
                       default=20.0, valid_range=RANGE_0POS)
        self._cmdValue(key='nkija', datatype='int', realdatatype=TYPE_INT,
                       descr="""For use with ntt=9 and ntt=10., For ntt=9, this
                             is the number of substeps of dt when integrating
                             the thermostat equations of motion, for greater
                             accuracy. For ntt=10, this specifies the number of
                             additional auxiliary velocity variables v1 and v2,
                             which will total nkija×v1 +nkija×v2""",
                       default=1, valid_range=RANGE_0POS)
        self._cmdValue(key='idistr', datatype='int', realdatatype=TYPE_INT,
                       descr="""For use with ntt=9, this is the frequency at
                             which the thermostat velocity distribution
                             functions are accumulated.""",
                       default=None, valid_range=RANGE_0POS)
        self._cmdValue(key='sinrtau', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""this specifies the time scale for determining
                             the masses associated with the two auxiliary
                             velocity variables v1 and v2 (e.g. thermostat
                             velocities)""", default=1.0, valid_range=RANGE_0POS)

        # Pressure regulation (Manual section 18.6.8)
        self._cmdList(key='ntp', datatype='int',
                      descr='Flag for constant pressure dynamics.',
                      default=0, valid_values=[0, 1, 2, 3])
        self._cmdList(key='barostat', datatype='int',
                      descr="""Flag used to control which barostat to use in
                            order to control the pressure.""", default=1,
                      valid_values=[1, 2])
        self._cmdValue(key='mcbarint', datatype='int', realdatatype=TYPE_INT,
                       descr="""Number of steps between volume change attempts
                             performed as part of the Monte Carlo barostat.""",
                       default=100, valid_range=RANGE_POS)
        self._cmdValue(key='pres0', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""Reference pressure (in units of bars) at which
                             the system is maintained (when NTP > 0).""",
                       default=1.0, valid_range=RANGE_POS)
        self._cmdValue(key='comp', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""compressibility of the system when NTP > 0. The
                             units are in 1.0 × 10-6 bar-1""", default=44.6,
                       valid_range=RANGE_0POS)
        self._cmdValue(key='taup', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='Pressure relaxation time (in ps), when NTP > 0.',
                       default=1.0, valid_range=RANGE_POS)
        self._cmdList(key='csurften', datatype='int',
                      descr='Flag for constant surface tension dynamics.',
                      default=0, valid_values=[0, 1, 2, 3])
        self._cmdValue(key='gamma_ten', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='Surface tension value in units of dyne/cm.',
                       default=0.0, valid_range=RANGE_0POS)
        self._cmdValue(key='ninterface', datatype='int', realdatatype=TYPE_INT,
                       descr="""Number of interfaces in the periodic box. There
                             must be at least two interfaces in the periodic box.""",
                       default=2, valid_range=RANGE_0POS)

        # SHAKE bond length constraints (Manual section 18.6.9)
        self._cmdList(key='ntc', datatype='int',
                      descr='Flag for SHAKE to perform bond length constraints.',
                      default=1, valid_values=[1, 2, 3])
        self._cmdValue(key='tol', datatype='float', realdatatype=TYPE_FLOAT,
                       descr="""Relative geometrical tolerance for coordinate
                             resetting in shake.""", default=0.00001,
                       valid_range=RANGE_0POS)
        self._cmdList(key='jfastw', datatype='int',
                      descr='Fast water definition flag.', default=0, valid_values=[0, 4])
        self._cmdString(key='noshakemask', datatype='str',
                        descr="""String that specifies atoms that are not to be
                              shaken (assuming that ntc>1).""", default='')

        # Water cap (Manual section 18.6.10)
        self._cmdList(key='ivcap', datatype='int',
                      descr="""Flag to control cap option. The "cap" refers to
                            a spherical portion of water centered on a point in
                            the solute and restrained by a soft half-harmonic
                            potential. For the best physical realism, this
                            option should be combined with igb=10, in order to
                            include the reaction field of waters that are beyond
                            the cap radius.""", default=0, valid_values=[0, 1, 2, 5])
        self._cmdValue(key='fcap', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='The force constant for the cap restraint potential.',
                       default=None, valid_range=RANGE_POS)
        self._cmdValue(key='cutcap', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='Radius of the cap, if ivcap=1 is used.',
                       default=None, valid_range=RANGE_POS)
        self._cmdValue(key='xcap', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='Location of the cap center, if ivcap=1 is used.',
                       default=None, valid_range=RANGE_ANY)

        # TODO: NMR refinement options (Manual section 18.6.11)
        # TODO: EMAP restraints (Manual section 18.6.12)
        # The pmemd cuda engine does not support these (and many other options
        # that I already have written). Still, users can use sander or the CPU
        # version of pmemd if they like

        # POTENTIAL FUNCTION PARAMETERS (MANUAL SECTION 18.7)
        # Only implementing the bold ones in the Manual (most important ones)
        # for the moment

        # Generic parameters (Manual Section 18.7.1)
        self._cmdList(key='ntf', datatype='int', descr='Force evaluation',
                      default=1, valid_values=[1, 2, 3, 4, 5, 6, 7, 8])
        self._cmdList(key='ntb', datatype='int',
                      descr="""Whether or not periodic boundaries are imposed
                            on the system during the calculation of non-bonded
                            interactions.""", default=None, valid_values=[0, 1, 2])
        self._cmdValue(key='cut', datatype='float', realdatatype=TYPE_FLOAT,
                       descr='Nonbonded cutoff, in Å',
                       default=None, valid_range=RANGE_POS)
        self._cmdList(key='ipol', datatype='int',
                      descr="""When set to 1, use a polarizable force field.
                            See Section 18.7.5 for more information.""",
                      default=0, valid_values=[0, 1])
        self._cmdList(key='ifqnt', datatype='int',
                      descr="""Flag for QM/MM run; if set to 1, you must also
                            include a &qmmm namelist.""", default=0, valid_values=[0, 1])
        self._cmdList(key='igb', datatype='int',
                      descr="""Flag for using the generalized Born or Poisson-Boltzmann
                            implicit solvent models.""",
                      default=0, valid_values=[0, 1, 2, 3, 4, 5, 6, 7, 8, 10])
        self._cmdList(key='irism', datatype='int',
                      descr="""Flag for 3D-reference interaction site model
                            molecular solvation method.""", default=0, valid_values=[0, 1])
        self._cmdList(key='ievb', datatype='int',
                      descr="""If set to 1, use the empirical valence bond method
                            to compute energies and forces.""", default=0, valid_values=[0, 1])
        self._cmdList(key='iamoeba', datatype='int',
                      descr="""Flag for using the amoeba polarizable potentials
                            of Ren and Ponder.""", default=0, valid_values=[0, 1])

        self._cmdString('outputnc', 'str', 'Name of output NetCDF file', None)

        # Assign everything to a single string
        self._cmdString('FORTRAN', 'str', 'Protocol in FORTRAN', None)
        self._cmdString('bash', 'str', '', None)
        # Files
        self._cmdString('bincoordinates', 'str', '', None)  # coordinate binary file .nc (-x)
        self._cmdString('coordinates', 'str', '', None)  # frame coordinates .rst (-r)
        self._cmdString('consref', 'str', '', None)  # constraints reference (-ref)
        self._cmdString('parmfile', 'str', '', None)  # topology file .prmtop (-p)

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

                self.__dict__[cmd] = defaultname  # use default file names

    def save(self, path, overwrite=False):
        """ Create a directory with all necessary input to run with AMBER.

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
            fo = open(os.path.join(path, self.__dict__[f]), 'wb')  # write all as binary
            fo.write(self._files[f])
            fo.close()

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
        """ Returns the AMBER configuration file string

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
        if self.__dict__['FORTRAN'] is not None:
            text = self.__dict__['FORTRAN']

        text += '#\n'

        maxwidth = np.max([len(k) for k in self.__dict__.keys()])

        keys = sorted(list(self.__dict__.keys()))
        # keys = keys + [keys.pop(keys.index('run'))]  # Move the run command to the end
        for cmd in keys:
            if not cmd.startswith('_') and self.__dict__[cmd] is not None and cmd != 'TCL':
                name = cmd
                if cmd == 'scaling14':  # variables cannot start with numbers. We need to rename it here for acemd
                    name = '1-4scaling'
                text += '{name: <{maxwidth}}\t{val:}\n'.format(name=name, val=str(self.__dict__[cmd]),
                                                               maxwidth=maxwidth)

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
