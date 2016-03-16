from htmd.projections.projection import Projection
from htmd.molecule.molecule import Molecule
import logging
import numpy
import subprocess
import shutil

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def _getTempFileName(prefix="",suffix=""):
        return os.path.join(tempfile._get_default_tempdir(),
                            prefix+next(tempfile._get_candidate_names())+suffix)

# Assumptions: file begins with #! FIELDS time
# first line is time
# only colvars follow
def _readColvar(colvar):
    data=[]
    with open(colvar,"r") as file:
        headerline=file.readline()
        cvnames=headerline.strip().replace('#! FIELDS time ','').split()

        for line in file:
            if line[0] == "#":
                continue
            cols_str=line.split()[1:]
            cols=[float(x) for x in cols_str]
            data.append(cols)
    return numpy.array(data)


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
        al=mol.get("serial",sel)
        al=list(al)
        self.label=label
        self.mol=mol
        self.sel=sel
        self.code = "%s: GROUP ATOMS=%s" % (label,",".join(map(str,al)))

    def __repr__(self):
        return self.code


class MetricPlumed2(Projection):
    """ Calculates generic collective variables through Plumed 2

    Parameters
    ----------
    plumed_inp_str: string
        The PLUMED script defining CVs - possibly a list of strings
    """

    def __init__(self, plumed_inp_str):
        # I am not sure at all about opening files here is good style
        self._precalculation_enabled = False
        self._plumed_exe = shutil.which("plumed")

        if not isinstance(plumed_inp_str,str):
            plumed_inp_str="\n".join(plumed_inp_str)
        self._plumed_inp=plumed_inp_str


    # Only called if single topology
    def _precalculate(self, mol):
        logger.info("In _precalculate")
        self._precalculation_enabled = True


    def getMapping(self, mol):
        # Useful?
        return





    # Arguments are actually self, mol
    def project(self, mol):
        """ Project molecule.

        Parameters
        ------------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>`
            A :class:`Molecule <htmd.molecule.molecule.Molecule>` object to project.

        Returns
        -------
        data : np.ndarray
            An array containing the projected data.
        """

        # --standalone-executable driver --box 100000,100000,100000 --mf_dcd /var/tmp/vmdplumed.8003/temp.dcd
        # --pdb /var/tmp/vmdplumed.8003/temp.pdb --plumed /var/tmp/vmdplumed.8003/META_INP

        # PDB
        pdb=_getTempFileName(suffix=".pdb")
        mol.write(pdb)

        # DCD
        dcd=_getTempFileName(suffix=".dcd")
        mol.write(dcd)
        logger.info("Done writing %d frames in %s" % (mol.numFrames,dcd))

        # Colvar
        colvar = _getTempFileName(suffix=".colvar")

        # Metainp
        metainp = _getTempFileName(suffix=".metainp")
        metainp_fp= open(metainp,"w+")
        metainp_fp.write("UNITS  LENGTH=A  ENERGY=kcal/mol  TIME=ps\n")
        metainp_fp.write(self._plumed_inp)
        metainp_fp.write('\nPRINT ARG=* FILE=%s\n# FLUSH STRIDE=1\n' % colvar)
        metainp_fp.close()

        logger.info("Plumed temporary files are %s (in) and %s (out)" %
                    (metainp, colvar))

        subprocess.check_output([self._plumed_exe,
                                 '--standalone-executable',
                                 'driver',
                                 '--mf_dcd', dcd,
                                 '--pdb', pdb,
                                 '--plumed',metainp])

        data = _readColvar(colvar)

        os.remove(pdb)
        os.remove(dcd)
        os.remove(colvar)
        os.remove(metainp)
        return data




if __name__ == "__main__":
    from htmd.home import home
    from os import path
    from glob import glob
    from htmd.simlist import simlist
    from htmd import *

    # One simulation
    mol = Molecule(path.join(home(), 'data', '1kdx', '1kdx_0.pdb'))
    mol.read(path.join(home(), 'data', '1kdx', '1kdx.dcd'))

    metric = MetricPlumed2(['d1: DISTANCE ATOMS=1,200',
                            'd2: DISTANCE ATOMS=5,6'] )
    data = metric.project(mol)



    # Simlist
    # datadirs=glob(path.join(home(), 'data', 'adaptive', 'data', '*' )
    # fsims=simlist(glob(path.join(home(), 'data', 'adaptive', 'data', '*', '/')),
    #              path.join(home(), 'data', 'adaptive', 'generators', '1','structure.pdb'))

    fsims=simlist(['/home/toni/work/htmd/htmd/htmd/data/adaptive/data/e1s1_1/',
             '/home/toni/work/htmd/htmd/htmd/data/adaptive/data/e1s2_1/'],
            '/home/toni/work/htmd/htmd/htmd/data/adaptive/generators/1/structure.pdb')

    metr=Metric(fsims)
    metr.projection(MetricPlumed2(['d1: DISTANCE ATOMS=2,3',
                                   'd2: DISTANCE ATOMS=5,6'] ))
    data=metr.project()
    pass
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
