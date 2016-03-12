from htmd.projections.projection import Projection
from htmd.molecule.molecule import Molecule
import logging
import numpy

import subprocess
import os
import ctypes
import tempfile

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


# The struct typedef'd in src/wrapper/Plumed.h
class c_plumed_type(ctypes.Structure):
    _fields_ = [("p", ctypes.c_void_p)]


class PlumedGroup():
    """ An atom group for use in the Plumed interface

    Parameters
    ----------
    mol: Molecule
        The molecule
    sel: str
        The atom selection defining the group
    """

    def __init__(self, label, mol, sel):
        self._atomlist=mol.atomselect(sel)
        self._label=label

    def __repr__(self):
        out="%s: GROUP ATOMS=%s" % (self._label,",".join(self._al))
        return out


class MetricPlumed2(Projection):
    """ Calculates generic collective variables through Plumed 2

    Parameters
    ----------
    plumed_inp_str: string
        The PLUMED script defining CVs - possibly a list of strings
    """

    def __init__(self, plumed_inp_str):
        # I am not sure at all about opening files here is good style
        self._enginePreinitialized = False
        self._plib = False

        self._metainp = tempfile.NamedTemporaryFile(mode="w+", suffix=".meta_inp", dir="/tmp")
        self._colvar_name = self._tempFileName(suffix=".colvar")

        if not isinstance(plumed_inp_str,str):
            plumed_inp_str="\n".join(plumed_inp_str)
        plumed_inp_str += '\nPRINT ARG=* FILE=%s\n# FLUSH STRIDE=1\n' % self._colvar_name
        self._metainp.write(plumed_inp_str)
        self._metainp.flush()


        logger.info("Plumed temporary files are %s (in) and %s (out)" %
                    (self._metainp.name, self._colvar_name))


    def _initEngine(self, mol):
        # This may throw exception, and I'm happy to stop if it does. 
        # TODO: Fallback to Matt's prepackaged library
        plumed_dir = subprocess.check_output(['plumed', 'info', '--root'])
        plumed_so = os.path.join(plumed_dir.strip(), b'src/lib/libplumed.so')

        self._plib = ctypes.cdll.LoadLibrary(plumed_so)
        logger.info("Using plumed shared library at %s" % plumed_so)

        plumed_create = self._plib.plumed_create
        plumed_create.restype = c_plumed_type
        self._pmain = plumed_create()


        # Stuff which could be needed.  http://plumed.github.io/doc-v2.1/developer-doc/html/_how_to_plumed_your_m_d.html       

        #        // Calls to pass data to plumed
        #        plumed_cmd(plumedmain,"setRealPrecision",&real_precision);     // Pass a pointer to an integer containing the size of a real number (4 or 8)
        #        plumed_cmd(plumedmain,"setMDEnergyUnits",&energyUnits);        // Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
        #        plumed_cmd(plumedmain,"setMDLengthUnits",&lengthUnits);        // Pass a pointer to the conversion factor between the length unit used in your code and nm 
        #        plumed_cmd(plumedmain,"setMDTimeUnits",&timeUnits);            // Pass a pointer to the conversion factor between the time unit used in your code and ps
        #        plumed_cmd(plumedmain,"setPlumedDat",&plumedInput);            // Pass the name of the plumed input file from the md code to plumed
        #        plumed_cmd(plumedmain,"setMPIComm",&MPI_COMM_WORLD);           // Pass a pointer to the MPI communicator to plumed
        #        // notice that from fortran the command "setMPIFComm" should be used instead
        #        plumed_cmd(plumedmain,"setNatoms",&natoms);                    // Pass a pointer to the number of atoms in the system to plumed
        #        plumed_cmd(plumedmain,"setMDEngine","gromacs");                // Pass the name of your md engine to plumed (now it is just a label) 
        #        plumed_cmd(plumedmain,"setLog",fplog);                         // Pass the file on which to write out the plumed log (if the file is already open)
        #        plumed_cmd(plumedmain,"setLogFile",fplog);                     // Pass the file  on which to write out the plumed log (to be created)
        #        plumed_cmd(plumedmain,"setTimestep",&delta_t);                 // Pass a pointer to the molecular dynamics timestep to plumed
        #        
        #        plumed_cmd(plumedmain,"setKbT",&kbT);                          // Tell to PLUMED the value of kbT - ONLY VALID IF API VERSION > 1



        # Calls to do the actual initialization (all the above commands must appear before this call)
        # plumed_cmd(plumedmain,"init",NULL);                            // Do all the initialization of plumed

    def getApiVersion(self):
        if not self._plib:
            raise Exception("Engine not initalized yet")

        version = ctypes.c_int()
        self._plib.plumed_cmd(self._pmain, b"getApiVersion", ctypes.byref(version))
        return version.value

    # Only called if single topology
    def _precalculate(self, mol):
        logger.info("In _precalculate")
        self._initEngine(mol)
        self._enginePreinitialized = True

    def _getEngine(self, mol):
        if not self._enginePreinitialized:
            self._initEngine(mol)

    def getMapping(self, mol):
        # TODO
        return

    def _tempFileName(self,prefix="",suffix=""):
        return os.path.join(tempfile._get_default_tempdir(),
                            prefix+next(tempfile._get_candidate_names())+suffix)

    # Assumptions: file begins with #! FIELDS time
    # first line is time
    # only colvars follow
    def _readColvar(self):
        data=[]
        with open(self._colvar_name,"r") as file:
            headerline=file.readline()
            cvnames=headerline.strip().replace('#! FIELDS time ','').split()

            for line in file:
                if line[0] == "#":
                    continue
                cols_str=line.split()[1:]
                cols=[float(x) for x in cols_str]
                data.append(cols)
        return data



    # Arguments are actually self, mol
    def project(self, mol):
        """ Project molecule.

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>`
            A :class:`Molecule <htmd.molecule.molecule.Molecule>` object to project.

        Returns
        -------
        data : np.ndarray
            An array containing the projected data.
        """

        self._getEngine(mol)
        cmd = self._plib.plumed_cmd
        cmd(self._pmain, b"setRealPrecision", ctypes.byref(ctypes.c_int(4)))
        cmd(self._pmain, b"setNatoms", ctypes.byref(ctypes.c_int(mol.numAtoms)))
        cmd(self._pmain, b"setMDEngine", b"HTMD")
        cmd(self._pmain, b"setPlumedDat", bytes(self._metainp.name,"ASCII"))
        cmd(self._pmain, b"setTimestep", ctypes.byref(ctypes.c_float(1.0)))
        # cmd(self._pmain, b"setLog", b"plumed.log")
        # cmd(self._pmain, b"setLogFile", b"plumed.logfile")
        cmd(self._pmain, b"init", 0)

        c_float_3natoms = ctypes.c_float * (3 * mol.numAtoms)
        c_float_9 = ctypes.c_float * 9

        c_masses= numpy.ctypeslib.as_ctypes(mol.masses)
        c_charges= numpy.ctypeslib.as_ctypes(mol.charge)
        c_forces=c_float_3natoms(0)
        c_virial=c_float_9(0)
        c_cell=c_float_9(0)

        # Not sure if something can be simplified here
        for fr in range(mol.numFrames):
            cmd(self._pmain, b"setStep", ctypes.byref(ctypes.c_int(fr)))

            coo = mol.coords[:, :, fr]
            coo = coo.reshape(-1)  # Make it a long vector (x1 y1 z1 x2 ...)
            coo = numpy.array(coo)  # Remove the slicing bc. not supported
            c_coo = numpy.ctypeslib.as_ctypes(coo)
            cmd(self._pmain, b"setPositions", c_coo)

            cmd(self._pmain, b"setMasses", c_masses)
            cmd(self._pmain, b"setCharges", c_charges)
            cmd(self._pmain, b"setForces", c_forces)

            cmd(self._pmain, b"setVirial", c_virial)
            cmd(self._pmain, b"setBox", c_cell)

            cmd(self._pmain, b"calc", 0)

        logger.info("Done passing %d frames" % fr)
        # cmd(self._pmain,b"runFinalJobs") # Undocumented?
        self._plib.plumed_finalize(self._pmain)

        data=self._readColvar()
        return data


if __name__ == "__main__":
    from htmd.home import home
    from os import path

    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))

    metric = MetricPlumed2(['d1: DISTANCE ATOMS=2,3',
                            'd2: DISTANCE ATOMS=5,6'] )
    data = metric.project(mol)

    # print("Plumed API is version %d" % pl.getApiVersion())


    #    metr = MetricDihedral(protsel='protein')
    #    data = metr.project(mol)
    #
    #    calcdata = np.array([ 0.91631763,  0.40045224,  0.90890707,  0.41699872, -0.99956623,
    #                          0.02945084,  0.52407037, -0.85167496, -0.67766999, -0.73536616,
    #                          0.53415969, -0.8453836 , -0.66133656, -0.7500893 ,  0.55669439,
    #                         -0.83071738, -0.90348715, -0.42861517,  0.5950773 , -0.80366847,
    #                         -0.5837572 , -0.81192828,  0.71012313, -0.70407751, -0.95668505,
    #                         -0.29112493,  0.53835619, -0.84271739, -0.9231271 ,  0.38449493,
    #                         -0.12417973,  0.99225974, -0.93138983,  0.36402332,  0.37667118,
    #                         -0.92634703, -0.14376672, -0.98961161,  0.37357125, -0.92760149,
    #                         -0.93655808,  0.35051244,  0.64918191, -0.76063319, -0.93758286,
    #                         -0.34776195,  0.51787137, -0.8554585 , -0.96970912,  0.2442626 ])
    #
    #    assert np.all(np.abs(calcdata - data[147, 500:550]) < 0.001), 'Diherdals calculation is broken'
