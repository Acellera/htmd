# File to write out Amber protocols

from htmd.molecule.molecule import Molecule
from htmd.userinterface import UserInterface
from htmd.acemd.amber import Amber
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_INT, TYPE_FLOAT, RANGE_0POS, RANGE_POS, RANGE_ANY
import os
import numpy as np
import logging
import shutil
logger = logging.getLogger(__name__)

class Equilibration(ProtocolInterface):

    """ Equilibration protocol

        Equilibration protocol for globular and membrane proteins
        It includes a flatbottom potential box to retrain a ligand
        for example within this box.

        Parameters
        ----------
        numsteps : int, default=0
            Number of steps to run the simulations in units of 4fs
        temperature : float, default=300
            Temperature of the thermostat in Kelvin
        k : float, default=0
            Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5
        reference : str, default='none'
            Reference selection to use as dynamic center of the flatbottom box.
        selection : str, default='none'
            Selection of atoms to apply the flatbottom potential
        box : list, default=[0, 0, 0, 0, 0, 0]
            Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]
        useconstantratio : bool, default=False
            For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.
        constraints : dict, default={'protein and noh and not name CA': 0.1, 'protein and name CA': 1}
            A dictionary containing as keys the atomselections of the constraints and as values the constraint scaling factor. 0 factor means no constraint, 1 full constraints and in between values are used for scaling. The order with which the constraints are applied is random, so make atomselects mutually exclusive to be sure you get the correct constraints.

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
        self._cmdObject('amber', ':class:`MDEngine <htmd.apps.app.App>` object', 'MD engine', None, Amber)

        self.amber = Amber()
        self.amber.imin = 0
        self.amber.nmropt = 0
        self.amber.ntx = 1
        self.amber.irest = 0
        self.amber.ntxo = 2
        self.amber.ntpr = 50
        self.amber.ntave = 0
#        self.amber.ntwr = #TODO: WHAT GOES HERE
        self.amber.iwrap = 0
        self.amber.ntwx = 0
        self.amber.ntwv = 0
        self.amber.ntwf = 0
        self.amber.ntwe = 0
        self.amber.ioutfm = 1
        self.amber.ntwprt = 0
        self.amber.idecomp = 0
        #  Frozen or restrained atoms (Manual section 18.6.4)
        self.amber.ibelly = 0
        self.amber.ntr = 0
        # These values need conditionals
        # Check Amber class in htmd.acemd.amber
#        self.amber.restraint_wt = #TODO: WHAT GOES HERE
#        self.amber.restraintmask = #TODO: WHAT GOES HERE
#        self.amber.bellymask = #TODO: WHAT GOES HERE
        #  Energy minimization (Manual section 18.6.5)
        self.amber.maxcyc = 1
        self.amber.ncyc = 10
        self.amber.ntmin = 1
        self.amber.dx0 = 0.01
        self.amber.drms = 1e-4
        #  Molecular dynamics (Manual section 18.6.6)
        
        self.amber.outputnc = 'Equilibration.nc'
        self.amber.FORTRAN = ''' HEATING\n &cntrl\n'''
        self.amber.bash = '''ENGINE -O -i INPUT -o OUTPUT -p TOPOLOGY -r RESTART -x TRAJOUT''' #FIXME: add -ref REFERENCE



    def _findFiles(self, inputdir):
        # Tries to find default files if the given don't exist
        defaults = {'coordinates': ('structure.rst', 'structure.inpcrd', 'structure.mdcrd'),
                    'structure': ('structure.psf', 'structure.prmtop'),
                    'parameters': ('parameters.in')}

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
                    os.path.join(inputdir, userval), os.path.join(inputdir, self.amber.__dict__[field])
                ))
            elif self.amber.__dict__[field] is None:
                raise RuntimeError('Could not locate any {f:} file in {i:}. '
                                   'Please set the Equilibration.acemd.{f:} property to '
                                   'point to the {f:} file'.format(f=field, i=inputdir))

    # FIXME: add inputdir and outputdir to write method
    # need to decide how to organize MD engine input files!!!
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
        self.amber.FORTRAN = ''' HEATING\n &cntrl\n'''

        restartfile = os.path.join(inputdir, self.amber.coordinates)
        topology= os.path.join(inputdir, self.amber.parmfile)

        if self.amber.consref is not None:
            reference = os.path.join(inputdir, self.amber.consref) 

        protocol=[]
        i=0
        for key, value in self.amber.__dict__.items():
            
            if key != 'FORTRAN' and key[0] != '_' and key not in ['bincoordinates', 'coordinates', 'consref', 'parmfile']:
                
                # cleans up the .in file a bit
                if i % 3 == 0 and i>0:
                    protocol.append('\n   {}={}'.format(key,value))
                else:
                    protocol.append('{}={}'.format(key,value))
                
                i+=1

        # format the FORTRAN file with MD parameters
        self.amber.FORTRAN = self.amber.FORTRAN+'   '+(', '.join(protocol)) + "\n /"

        # write out the FORTRAN file
        with open(os.path.join(outputdir,'Equilibration.in'), 'w') as text_file:
            text_file.write(self.amber.FORTRAN)

        # prepare the bash file
        self.amber.bash = self.amber.bash.replace('INPUT', 'Equilibration.in') # putting these strings as placeholders for now
        self.amber.bash = self.amber.bash.replace('OUTPUT', 'Equilibration.out')
        self.amber.bash = self.amber.bash.replace('TOPOLOGY', self.amber.parmfile)
        self.amber.bash = self.amber.bash.replace('RESTART', self.amber.coordinates)
        self.amber.bash = self.amber.bash.replace('TRAJOUT', self.amber.outputnc)
        #FIXME: add ref
        #self.amber.bash = self.amber.bash.replace('REFERENCE', self.amber.reference)

        with open(os.path.join(outputdir,'Equilibration.sh'), 'w') as text_file:
            text_file.write(self.amber.bash)

        shutil.copy(restartfile, outputdir)
        shutil.copy(topology, outputdir)
        #FIXME: add ref
        #shutil.copy(consref, outputdir)
