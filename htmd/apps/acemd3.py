# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import htmd
import shutil
import numpy as np
from protocolinterface import ProtocolInterface, val


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
            Atomselection containing the atoms on whose center of mass the restraints will be applied
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
            Atomselection containing the atoms on whose center of mass the restraints will be applied
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


class Acemd3(ProtocolInterface):
    _defaultfnames = {'bincoordinates': 'input.coor', 'binvelocities': 'input.vel', 'structure': 'structure.*',
                      'parameters': 'parameters', 'extendedsystem': 'input.xsc',
                      'coordinates': 'structure.pdb', 'velocities': 'velocity.pdb', 'parmfile': 'structure.prmtop'}

    def __init__(self):
        """
        temperature : float, default=None
            Temperature of the thermostat in Kelvin.
        restart : str, default=None
            Restart simulation.
        trajectoryfile : str, default=None
            Output file name.
        trajectoryfreq : int, default=None
            Trajectory sampling frequency in steps.
        timestep : int, default=None
            Simulation timestep.
        pme : str, default=None
            Particle-mesh Ewald summation.
        switching : str, default=None
            Apply switching function to the van der Waals potential.
        switchdist : float, default=None
            Distance in Angstrom at which to begin applying the switching function.
        cutoff : float, default=None
            Real-space cutoff in Angstroms for electrostatics and van der Waals.
        thermostat : str, default=None
            Enable thermostatic control
        thermostattemp : float, default=None
            Target temperature (K) for thermostatic control
        thermostatdamping : float, default=None
            Damping constant for the Langevin thermostat in ps^-1
        restraints : str, default=None
            Restraining potentials
        barostat : str, default=None
            Enable pressure control
        barostatpressure : float, default=None
            The target pressure in bar
        useflexiblecell : str, default=None
            Allow X, Y and Z unit cell dimensions to vary independently
        useconstantarea : str, default=None
            Constrain the X,Y dimensions of the unit cell. Allow Z to vary independently.
        useconstantratio : str, default=None
            Constrain the X:Y ratio of the unit cell dimensions. Allow Z to vary independently.
        minimize : int, default=None
            The number of energy minimization steps to perform before commencing dynamics.
        run : str, default=None
            The length of simulation ro run. May be specified as a number of steps or as a time if one of the suffices "us", "ns", "ps", "fs" is used.
        celldimension : str, default=None
            The dimensions of the unit cell in Angstrom. Note that the unit cell must be cuboid. Overrides any dimension given in the "coordinates" PDB.
        bincoordinates : str, default=None
            Optional initial system geometry in NAMD BINCOOR format. If specified, overrides "coordinates"
        binvelocities : str, default=None
            Optional initial system velocity field in NAMD BINCOOR format. If specified, overrides field generated by "temperature" and "velocities"
        structure : str, default=None
            The filename of a CHARMM PSF file
        parameters : str, default=None
            The filename of a CHARMM PAR file
        parmfile : str, default=None
            The filename of an Amber PRMTOP file
        extendedsystem : str, default=None
            Filename of a NAMD XSC format file giving the periodic cell dimensions. Overrides "celldimension" and any dimensions in the "coordinates" PDB
        coordinates : str, default=None
            Mandatory initial system geometry in PDB format
        velocities : str, default=None
            Optional initial system velocity field in NAMD BINCOOR format. If specified, overrides field generated by "temperature"
        """
        super().__init__()
        self._files = {}
        self._outnames = {}

        # Options
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
        from htmd.util import ensurelist
        text = ''
        maxwidth = np.max([len(k) for k in self.__dict__.keys()])

        keys = sorted(list(self.__dict__.keys()))
        keys = keys + [keys.pop(keys.index('restraints')), keys.pop(keys.index('run'))]  # Move the run command to the end
        for cmd in keys:
            if cmd == 'restraints' and self.restraints is not None:
                for r in ensurelist(self.restraints):
                    text += '{}\n'.format(r.format(maxwidth))
            elif not cmd.startswith('_') and self.__dict__[cmd] is not None:
                val = self.__dict__[cmd]
                if cmd in self._outnames:
                    val = self._outnames[cmd]
                text += '{name: <{maxwidth}}\t{val:}\n'.format(name=cmd, val=val, maxwidth=maxwidth)

        if not quiet:
            print(text)
        else:
            return text

    def _writeBashRun(self, fname):
        with open(fname, 'w') as f:
            f.write('#!/bin/bash\nacemd3 >log.txt 2>&1')
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
    homedir = htmd.home.home()

    r = list()
    r.append(GroupRestraint('resname MOL', 5, [(10, '10ns'), (5, '15ns'), (0, '20ns')], axes='z'))
    r.append(AtomRestraint('name CA', 0.1, [(10, '10ns'), (5, '15ns'), (0, '20ns')]))

    acemd = Acemd3()
    acemd.structure = '5dhfr_cube.psf'
    acemd.parameters = 'par_all22_prot.inp'
    acemd.coordinates = '5dhfr_cube.pdb'
    acemd.restraints = r
    acemd.setup(homedir + '/data/dhfr', '/tmp/testdir', overwrite=True)
    print(acemd)

