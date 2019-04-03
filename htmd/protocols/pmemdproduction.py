# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Author: Juan Eiros Zamora (je714@ic.ac.uk) & Gil Ferreira Hoben (gf712@ic.ac.uk) Organisation: Imperial College London
# Modified by: Stefan Doerr

from moleculekit.molecule import Molecule
from htmd.apps.pmemd import Pmemd
from htmd.protocols.oldprotocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
import os
import numpy as np
import logging
import shutil
import re
from htmd.decorators import _Deprecated
logger = logging.getLogger(__name__)


@_Deprecated('1.13.6')
class Production(ProtocolInterface):

    """ Production protocol for globular and membrane proteins
        Example
        -------
        >>> from htmd.protocols.equilibration_v1 import Equilibration
        >>> md = Equilibration()
        >>> md.numsteps = 10000000
        >>> md.temperature = 300
        >>> md.useconstantratio = True  # only for membrane sims
        >>> # this is only needed for setting the flatbottom potential, otherwise remove it
        >>> md.reference = 'protein and resid 293'
        >>> md.selection = 'segname L and noh'
        >>> md.box = [-25, 25, -25, 25, 43, 45]
        >>> md.k = 5
        >>> md.write('./build','./equil')
    """

    def __init__(self):
        super().__init__()
        self._cmdObject('amber', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Pmemd)

        self.amber = Pmemd()
        self.amber.imin = 0
        self.amber.ntx = 5
        self.amber.irest = 1
        self.amber.ntc = 2
        self.amber.ntf = 2
        self.amber.nstlim = 12500000
        self.amber.ntt = 3
        self.amber.gamma_ln = 1.0
        self.amber.temp0 = 300.0
        self.amber.ntpr = 5000
        self.amber.ntwr = 50000
        self.amber.ntwx = 5000
        self.amber.dt = 0.004
        self.amber.ig = -1
        self.amber.ntb = 2
        self.amber.ntp = 1
        self.amber.cut = 8.0
        self.amber.barostat = 2
        self.amber.ioutfm = 1
        self.amber.ntxo = 2
        self.amber.outputnc = 'Production.nc'
        self.amber.FORTRAN = """ Production
 &cntrl
   imin=IMIN, ntx=NTX, irest=IREST,
   ntc=NTC, ntf=NTF,
   nstlim=NSTLIM, ntt=NTT, gamma_ln=GAMMA_LN,
   temp0=TEMP0,
   ntpr=NTPR, ntwr=NTWR, ntwx=NTWX,
   dt=DT, ig=IG,
   ntb=NTB, ntp=NTP, cut=CUT, barostat=BAROSTAT, ioutfm=IOUTFM, ntxo=NTXO,
 / \n"""
        # self.amber.FORTRAN = ''' HEATING\n &cntrl\n'''
        self.amber.bash = '''ENGINE -O -i INPUT -o OUTPUT -p TOPOLOGY -c RESTART -x TRAJOUT -r RESTOUT'''  # FIXME: add -ref REFERENCE

    def _findFiles(self, inputdir):
        # Tries to find default files if the given don't exist
        defaults = {'coordinates': ('structure.rst', 'structure.inpcrd', 'structure.mdcrd'),
                    'parmfile': ('structure.prmtop',)}

        for field in defaults:
            userval = self.amber.__dict__[field]
            if userval is not None and not os.path.exists(os.path.join(inputdir, userval)):
                self.amber.__dict__[field] = None

            if self.amber.__dict__[field] is None:
                for val in defaults[field]:
                    if os.path.exists(os.path.join(inputdir, val)):
                        self.amber.__dict__[field] = val
                        break

            if userval is not None and self.amber.__dict__[field] is not None and self.amber.__dict__[field] != userval:
                logger.warning('Could not locate structure file {}. Using {} instead.'.format(
                    os.path.join(inputdir, userval), os.path.join(
                        inputdir, self.amber.__dict__[field])
                ))
            elif self.amber.__dict__[field] is None:
                raise RuntimeError('Could not locate any {f:} file in {i:}. '
                                   'Please set the Equilibration.acemd.{f:} property to '
                                   'point to the {f:} file'.format(f=field, i=inputdir))

    # FIXME: add inputdir and outputdir to write method
    # need to decide how to organize MD engine input files!!!
    # Probably will create a bash script here with all the information needed to run pmemd, e.g.:
    # pmemd.cuda_SPFP -O -i Equil.in -o Equil.out -c structure.rst -p structure \
    # -r EquilRestart.rst -x EquilRestart.nc (-ref Min.rst)

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
        # this code automatically loops through the amber attributes and assigns the given value to
        # the parameter.

        self._findFiles(inputdir)

        restartfile = os.path.join(inputdir, self.amber.coordinates)
        topology = os.path.join(inputdir, self.amber.parmfile)

        if self.amber.consref is not None:
            reference = os.path.join(inputdir, self.amber.consref)

        # FOR NOW WE WILL HAVE A FIXED STRING LIKE IN THE ACEMD PROTOCOL
        # THIS WAY WE CAN CONTROL BETTER WHAT PARAMETERS CAN BE SPECIFIED
        # IN THE PROTOCOL

        # protocol=[]
        # i=0
        # for key, value in self.amber.__dict__.items():

        #     if not key.startswith('_') and key not in ['bincoordinates', 'coordinates',
        #                                                'consref', 'parmfile','FORTRAN',
        #                                                'bash', 'outputnc']:

        #         # cleans up the .in file a bit
        #         if i % 3 == 0 and i>0:
        #             protocol.append('\n   {}={}'.format(key,value))
        #         else:
        #             protocol.append('{}={}'.format(key,value))

        #         i+=1

        # # format the FORTRAN file with MD parameters
        # self.amber.FORTRAN = self.amber.FORTRAN+'   '+(', '.join(protocol)) + "\n /"

        # preparing the FORTRAN file
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'IMIN', str(self.amber.imin))
        self.amber.FORTRAN = re.sub(r'\bNTX\b', str(
            self.amber.ntx), self.amber.FORTRAN)
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'IREST', str(self.amber.irest))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTC', str(self.amber.ntc))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTF', str(self.amber.ntf))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NSTLIM', str(self.amber.nstlim))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTT', str(self.amber.ntt))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'GAMMA_LN', str(self.amber.gamma_ln))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTR', str(self.amber.ntr))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'IG', str(self.amber.ig))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'TEMP0', str(self.amber.temp0))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTPR', str(self.amber.ntpr))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTWR', str(self.amber.ntwr))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTWX', str(self.amber.ntwx))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'DT', str(self.amber.dt))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTB', str(self.amber.ntb))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTP', str(self.amber.ntp))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'CUT', str(self.amber.cut))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'IOUTFM', str(self.amber.ioutfm))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'BAROSTAT', str(self.amber.barostat))
        self.amber.FORTRAN = self.amber.FORTRAN.replace(
            'NTXO', str(self.amber.ntxo))

        # write out the FORTRAN file, for now let's call it Production.in
        with open(os.path.join(outputdir, 'Production.in'), 'w') as text_file:
            text_file.write(self.amber.FORTRAN)

        # prepare the bash file
        # putting these strings as placeholders for now
        self.amber.bash = self.amber.bash.replace('INPUT', 'Production.in')
        self.amber.bash = self.amber.bash.replace(
            'OUTPUT', 'Production.out')
        self.amber.bash = self.amber.bash.replace(
            'TOPOLOGY', self.amber.parmfile)
        self.amber.bash = self.amber.bash.replace(
            'RESTART', self.amber.coordinates)
        self.amber.bash = self.amber.bash.replace(
            'TRAJOUT', self.amber.outputnc)
        self.amber.bash = self.amber.bash.replace(
            'RESTOUT', 'Production_new.rst')
        # FIXME: add ref
        # self.amber.bash = self.amber.bash.replace('REFERENCE', self.amber.reference)

        with open(os.path.join(outputdir, 'MD.sh'), 'w') as text_file:
            text_file.write(self.amber.bash)

        shutil.copy(restartfile, outputdir)
        shutil.copy(topology, outputdir)
        # FIXME: add ref
        # shutil.copy(consref, outputdir)
