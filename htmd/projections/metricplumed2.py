from htmd.projections.projection import Projection
from htmd.projections.metric import Metric
from htmd.molecule.molecule import Molecule

from abc import ABC
import logging
import numpy
import subprocess
import shutil
import os
import tempfile

logger = logging.getLogger(__name__)

cvcounter = 1

# Local utility functions --------------------------------------------------------

def _getTempDirName(prefix=""):
    return os.path.join(tempfile._get_default_tempdir(),
                        prefix + next(tempfile._get_candidate_names()))


def _getPlumedPath():
    """ Return path to plumed executable, or raise an exception if not found. """
    pr = subprocess.check_output(["plumed", "--standalone-executable", "info", "--root"])
    pr = pr.strip().decode("utf-8")
    return pr


# Plumed statement wrappers --------------------------------------------------------

# ABC should prevent instantiation but doesn't
class PlumedStatement(ABC):
    """ Abstract base class for Plumed statements. Do not use directly. """
    pass


class PlumedCV(PlumedStatement):
    """ Define a Plumed2 CV

    Parameters
    ----------
    cv : str
        The CV keyword as a string (e.g.: "DISTANCE"). Case-insensitive.
    label : str
        The label assigned to the CV
    args :
        Other named arguments to the CV (e.g. ATOMS=1, NOPBC=True). Objects of
        class PlumedGroup and PlumedCOM are expanded and prepended as appropriate.
    verbatim : str, optional
        Code which will be added as-is to the CV line

    Examples
    --------
    >>> PlumedCV("GYRATION", "rgyr", ATOMS="10-20", TYPE="RADIUS")
    rgyr: GYRATION ATOMS=10-20 TYPE=RADIUS
    >>> m=Molecule("3ptb")
    >>> grp=PlumedGroup(m,"grp","serial 10 to 20")
    >>> PlumedCV("GYRATION", "rgyr2", ATOMS=grp, TYPE="ASPHERICITY", NOPBC=True)
    grp: GROUP ATOMS=10,11,12,13,14,15,16,17,18,19,20
    rgyr2: GYRATION ATOMS=grp NOPBC TYPE=ASPHERICITY
    >>> protCA=PlumedCOM(m,"protCA","chain A and name CA")
    >>> lig=PlumedCOM(m,"lig","resname BEN and noh")
    >>> PlumedCV("DISTANCE", "dist", ATOMS=[protCA,lig])
    lig: COM ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    protCA: COM ATOMS=2,10,17,21,25,37,44,50,54,59,67,74,81,88,100,109,116,122,130,138,144,148,160,170,181,187,191,195,201,209,217,225,231,240,254,261,268,274,279,284,294,300,312,321,327,331,339,348,355,366,374,378,387,395,403,411,419,426,433,442,446,454,463,472,483,491,497,502,508,517,523,531,538,548,555,561,573,581,587,595,602,610,618,626,634,642,650,658,666,675,683,692,698,703,708,714,722,730,736,747,754,759,765,773,779,787,794,801,807,813,818,824,829,833,840,849,855,863,871,877,881,895,899,907,914,923,929,935,939,946,952,964,971,979,986,994,1003,1009,1017,1026,1031,1038,1046,1054,1060,1068,1074,1080,1086,1095,1101,1106,1118,1125,1129,1138,1146,1153,1159,1167,1175,1186,1192,1197,1201,1213,1221,1230,1234,1238,1247,1255,1261,1267,1276,1280,1288,1294,1298,1302,1309,1316,1323,1329,1335,1339,1348,1356,1365,1369,1377,1384,1390,1404,1408,1414,1418,1424,1429,1438,1447,1455,1464,1471,1475,1482,1494,1501,1510,1517,1523,1531,1543,1550,1556,1570,1578,1587,1596,1603,1611,1616,1622,1631
    dist: DISTANCE ATOMS=protCA,lig
    >>> PlumedCV("GYRATION", "rgyr3", ATOMS=m)
    """

    def __init__(self, cv, label, verbatim=None, **kw):
        self.label = label
        self.cv = cv
        self.args = kw
        self.verbatim = verbatim
        self.prereq = []

        for k in self.args:
            v = self.args[k]
            # Check all possible RHS types...
            # If the arg is a Plumed group, add it to prereq and replace by label
            if isinstance(v, PlumedGenericGroup):
                self.prereq.append(v)
                self.args[k] = v.label
            # If the args is a Molecule object, convert to a PlumedGroup
            elif isinstance(v, Molecule):
                tmpGrp=PlumedGroup(label=None,mol=v,sel="all")
                self.prereq.append(tmpGrp)
                self.args[k] = tmpGrp.label
            # Boolean, as in PBC=True, is ok
            elif isinstance(v, bool):
                pass
            # String is ok
            elif isinstance(v, str):
                pass
            elif isinstance(v, int):
                self.args[k]=str(v)
            # Ditto if it is a list-like object, plus expand with commas
            elif hasattr(v, '__iter__'):
                for l in range(len(v)):
                    le = v[l]
                    if isinstance(le, PlumedGenericGroup):
                        self.prereq.append(le)
                        v[l] = le.label
                    elif isinstance(le, Molecule):
                        tmpGrp = PlumedGroup(label=None, mol=le, sel="all")
                        self.prereq.append(tmpGrp)
                        v[l] = tmpGrp.label
                    elif isinstance(le, str):
                        pass # already a string
                    elif isinstance(le, int):
                        v[l] = str(le)
                    else:
                        raise TypeError("Unexpected type passed at argument (list): " + k)
                self.args[k] = ",".join(v)
            else:
                raise TypeError("Unexpected type passed at argument: " + k)

    def __repr__(self):
        r = ""

        # Prerequisites (uniquely)
        for p in sorted([str(i) for i in set(self.prereq)]):
            r = r + p + "\n"

        # Label
        r = r + self.label + ": " + self.cv + " "

        # Code
        if self.verbatim:
            r = r + self.verbatim + " "

        # Args
        for k, v in sorted(self.args.items()):
            if v == True:
                r = r + k + " "
            else:
                r = r + k + "=" + str(v) + " "

        return r.strip()


class PlumedGenericGroup(PlumedStatement):
    """ Abstract class from which PLUMED groups are inherited. Do not use directly. """

    def __init__(self, mol, label, sel, type=""):
        global cvcounter
        al = mol.get("serial", sel)
        al = list(al)
        if label:
            self.label = label
        else:
            self.label = "lab_" + str(cvcounter)
            cvcounter = cvcounter + 1
        self.mol = mol
        self.sel = sel
        self.code = "%s: %s ATOMS=%s" % (self.label, type, ",".join(map(str, al)))

    def __repr__(self):
        return self.code


class PlumedGroup(PlumedGenericGroup):
    """ An atom GROUP for use in the Plumed interface

    Parameters
    ----------
    mol: Molecule
        The molecule
    label: str
        The label assigned to the group.
    sel: str
        The atom selection defining the group

    Example
    -------
    >>> m=Molecule("3PTB")
    >>> PlumedGroup(m,"ben","resname BEN")
    ben: GROUP ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    """

    def __init__(self, mol, label, sel):
        return super(PlumedGroup, self).__init__(mol, label, sel, "GROUP")


class PlumedCOM(PlumedGenericGroup):
    """ An atom center-of-mass for use in the Plumed interface

    Parameters
    ----------
    mol: Molecule
        The molecule
    label: str
        The label assigned to the group
    sel: str
        The atom selection defining the group

    Example
    -------
    >>> m=Molecule("3PTB")
    >>> PlumedCOM(m,"ben_cm","resname BEN")
    ben_cm: COM ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    """

    def __init__(self, mol, label, sel):
        return super(PlumedCOM, self).__init__(mol, label, sel, "COM")


# Plumed projector --------------------------------------------------------

class MetricPlumed2(Projection):
    """ Calculates generic collective variables through Plumed 2

    The collective variables are defined in PLUMED 2's syntax. PLUMED needs be installed
    separately; see http://www.plumed.org/.

    Parameters
    ----------
    plumed_inp_str: string
        The PLUMED script defining CVs - a string or a list of strings (which are concatenated)

    Examples
    --------
    >>> dd = htmd.home(dataDir="adaptive")
    >>> fsims = htmd.simlist([dd + '/data/e1s1_1/', dd + '/data/e1s2_1/'], dd + '/generators/1/structure.pdb')
    >>> metr = Metric(fsims)
    >>> metr.projection(MetricPlumed2( ['d1: DISTANCE ATOMS=2,3', 'd2: DISTANCE ATOMS=5,6']))
    >>> data=metr.project()
    >>> data.dat
    array([ array([[ 1.68597198,  1.09485197], ...
    """

    def __init__(self, plumed_inp_str):
        # I am not sure at all about opening files here is good style
        self._precalculation_enabled = False
        self._plumed_exe = shutil.which("plumed")
        self.colvar = None
        self.cvnames = None

        try:
            pp = _getPlumedPath()
            logger.info("Plumed path is " + pp)
        except Exception as e:
            raise Exception("To use MetricPlumed2 please ensure PLUMED 2's executable is installed and in path")

        if not isinstance(plumed_inp_str, str):
            plumed_inp_str = "\n".join(plumed_inp_str)
        self._plumed_inp = plumed_inp_str

    def _readColvar(self):
        # Assumptions: file begins with #! FIELDS time
        # first line is time
        # only colvars follow
        assert self.colvar, "colvar variable not defined"
        data = []
        with open(self.colvar, "r") as file:
            headerline = file.readline()
            self.cvnames = headerline.strip().replace('#! FIELDS time ', '').split()

            for line in file:
                if line[0] == "#":
                    continue
                cols_str = line.split()[1:]
                cols = [float(x) for x in cols_str]
                data.append(cols)
        return numpy.array(data, dtype=numpy.float32)

    # Only called if single topology
    def _precalculate(self, mol):
        logger.info("In _precalculate")
        self._precalculation_enabled = True

    def getMapping(self):
        """ Return the labels of the colvars used in this projection.

        Can only be used after the projection has been executed.

        Returns
        -------
        cvnames
            A list of cv names
        """
        if self.cvnames:
            return self.cvnames
        else:
            raise Exception("MetricPlumed's getMapping can only be called after the projection")

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

        logger.debug("_precalculate was called? %d" % self._precalculation_enabled)

        # --standalone-executable driver --box 100000,100000,100000 --mf_dcd /var/tmp/vmdplumed.8003/temp.dcd
        # --pdb /var/tmp/vmdplumed.8003/temp.pdb --plumed /var/tmp/vmdplumed.8003/META_INP

        td = _getTempDirName("metricplumed2-")
        os.mkdir(td)

        # PDB
        pdb = os.path.join(td, "temp.pdb")
        mol.write(pdb)

        # DCD
        dcd = os.path.join(td, "temp.dcd")
        mol.write(dcd)
        logger.info("Done writing %d frames in %s" % (mol.numFrames, dcd))

        # Colvar
        colvar = os.path.join(td, "temp.colvar")
        self.colvar = colvar
        logger.info("Colvar file is " + colvar)

        # Metainp
        metainp = os.path.join(td, "temp.metainp")
        metainp_fp = open(metainp, "w+")
        metainp_fp.write("UNITS  LENGTH=A  ENERGY=kcal/mol  TIME=ps\n")
        metainp_fp.write(self._plumed_inp)
        metainp_fp.write('\nPRINT ARG=* FILE=%s\n# FLUSH STRIDE=1\n' % colvar)
        metainp_fp.close()

        cmd = [self._plumed_exe, '--standalone-executable',
               'driver',
               '--mf_dcd', dcd,
               '--pdb', pdb,
               '--plumed', metainp]
        logger.info("Invoking " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logger.error("Error from PLUMED (stdout): " + e.stdout.decode("utf-8"))
            logger.error("Error from PLUMED (stderr):" + e.stderr.decode("utf-8"))
            logger.error("Leaving temporary data in " + td)
            raise e

        data = self._readColvar()
        shutil.rmtree(td)

        return data


# Test --------------------------------------------------------

if __name__ == "__main__":
    import sys
    import numpy as np
    import htmd
    from htmd.projections.metricplumed2 import MetricPlumed2

    try:
        _getPlumedPath()
    except:
        print("Tests in %s skipped because plumed executable not found." % __file__)
        sys.exit()

    import doctest

    doctest.testmod()

    # One simulation
    mol = Molecule(os.path.join(htmd.home(), 'data', '1kdx', '1kdx_0.pdb'))
    mol.read(os.path.join(htmd.home(), 'data', '1kdx', '1kdx.dcd'))

    metric = MetricPlumed2(['d1: DISTANCE ATOMS=1,200',
                            'd2: DISTANCE ATOMS=5,6'])
    #    metric = MetricPlumed2([''])  # to test exceptions
    data = metric.project(mol)
    ref = np.array([0.536674, 21.722393, 22.689391, 18.402114, 23.431387, 23.13392, 19.16376, 20.393544,
                    23.665517, 22.298349, 22.659769, 22.667669, 22.484084, 20.893447, 18.791701,
                    21.833056, 19.901318])
    assert np.all(np.abs(ref - data[:, 0]) < 0.01), 'Plumed demo calculation is broken'

    # Simlist
    # datadirs=glob(os.path.join(home(), 'data', 'adaptive', 'data', '*' )
    # fsims=simlist(glob(os.path.join(home(), 'data', 'adaptive', 'data', '*', '/')),
    #              os.path.join(home(), 'data', 'adaptive', 'generators', '1','structure.pdb'))

    dd = htmd.home(dataDir="adaptive")
    fsims = htmd.simlist([dd + '/data/e1s1_1/', dd + '/data/e1s2_1/'],
                         dd + '/generators/1/structure.pdb')
    metr = Metric(fsims)
    metr.projection(MetricPlumed2(
        ['d1: DISTANCE ATOMS=2,3',
         'd2: DISTANCE ATOMS=5,6']))
    data2 = metr.project()
    pass
