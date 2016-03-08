# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import htmd
import shutil
import numpy as np
from htmd.userinterface import UserInterface


class Acemd(UserInterface):
    _cmdfile = {'bincoordinates': 'input.coor', 'binvelocities': 'input.vel', 'binindex': 'index.idx',
               'structure': 'structure.psf', 'parameters': 'parameters', 'extendedsystem': 'input.xsc',
               'coordinates': 'structure.pdb', 'velocities': 'velocity.pdb', 'consref': 'structure.pdb'}

    _cmdopt =  {'temperature':None,'restart':'on','restartfreq':'5000','outputname':'output',
               'xtcfile':'output.xtc','xtcfreq':'25000','timestep':'4','rigidbonds':'all',
               'hydrogenscale':'4','switching':'on','switchdist':'7.5','cutoff':'9','exclude':'scaled1-4',
               'scaling14':'1.0','langevin':'on','langevintemp':'300','langevindamping':'1','pme':'on',
               'pmegridspacing':'1.0','fullelectfrequency':'2','energyfreq':'1000','constraints':'off',
               'consref':None,'constraintscaling':'1.0','berendsenpressure':'off','berendsenpressuretarget':'1.01325',
               'berendsenpressurerelaxationtime':'800','tclforces':'off','minimize':None,'run':None,'TCL':'',
                'celldimension':None,'useconstantratio':None}

    def __init__(self):
        self._commands = Acemd._cmdfile.copy()  # needs to act on dictionary to avoid going via__setattr__
        self._commands.update(Acemd._cmdopt) # merging

    def load(self, path='.'):
        """ Loads all files required to run a simulation and apply eventually configured protocols to it

        Files content are stored in _FILE attributes. Filenames are not needed and set to default names,
        e.g. FILE_bincoordinates will be present containing the input coordinate file, if bincoordinates was present.
        bincoordinates is reset to the canonical value of 'input.coor'.

        Parameters
        ----------
        path : str
            Working directory relative to which the configuration file is read
        """
        #load files and reset filenames
        for cmd in self._cmdfile.keys():
            if cmd in self.__dict__.keys():
                if self.__dict__[cmd]:
                    absfname = os.path.join(path, self.__dict__[cmd])
                    fo = open(absfname, 'rb')  #read all as binary
                    a = fo.read()
                    fo.close()
                    self.__dict__['_FILE' + cmd] = a
                    self.__dict__[cmd] = self._cmdfile[cmd]  # use default file names

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
    
        commands = self.__dict__
        for cmdorg in list(commands.keys()):  #make copy to delete while iterating for py3
            if cmdorg.startswith('_FILE'):
                cmd = cmdorg[5:]  # _FILE is long 5
                fname = commands[cmd]
                fo = open(os.path.join(path, fname), 'wb')   #write all as binary
                fo.write(commands[cmdorg])
                fo.close()
                del(commands[cmdorg])
    
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
        if 'TCL' in self.__dict__.keys():
            text = self.__dict__['TCL']

        text += '#\n'

        keys = list(self._commands.keys())
        keys = keys + [keys.pop(keys.index('run'))]  # Move the run command to the end
        for cmd in keys:
            if cmd in self.__dict__.keys():
                if cmd[0] != '_':
                    if cmd == 'scaling14':  # 1-4 is not accetable field
                        text += '1-4scaling\t' + str(self.__dict__[cmd]) + '\n'
                    elif cmd[0:5] != '_FILE' and cmd != 'TCL':
                        text += cmd + '\t' + str(self.__dict__[cmd]) + '\n'

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

if __name__ == "__main__":
    #l=Acemd.protocols(quiet=True)
    acemd = Acemd()
    acemd.structure = '5dhfr_cube.psf'
    acemd.consref = None
    acemd.constraintscaling = None
    acemd.parameters = 'par_all22_prot.inp'
    homedir = htmd.home()
    acemd.coordinates = '5dhfr_cube.pdb'
    acemd.setup(homedir + '/data/dhfr', '/tmp/testdir', overwrite=True)
    print(acemd)

