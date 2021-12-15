# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from protocolinterface import ProtocolInterface, val
from htmd.mdengine.acemd.acemd import Acemd, _Restraint, GroupRestraint, AtomRestraint
import os
import numpy as np
import logging

logger = logging.getLogger(__name__)


class Production(ProtocolInterface):
    """Production protocol v6

    Production protocol for globular and membrane proteins. You can optionally define a flatbottom potential box and
    atom constraints for the production run.

    An Acemd class object is stored in the Production object which can be used to modify futher options.
    For documentation on further options see :class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>`

    Parameters
    ----------
    runtime : float, default=0
        Running time of the simulation.
    timeunits : str, default='steps'
        Units for runtime. Can be 'steps', 'ns' etc.
    temperature : float, default=300
        Temperature of the thermostat in Kelvin
    fb_k : float, default=0
        Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5
    fb_reference : str, default='none'
        Reference selection to use as dynamic center of the flatbottom box.
    fb_selection : str, default='none'
        Selection of atoms to apply the flatbottom potential
    fb_box : list, default=[0, 0, 0, 0, 0, 0]
        Position of the flatbottom box in term of the reference center given as [xmin, xmax, ymin, ymax, zmin, zmax]
    useconstantratio : bool, default=False
        For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.
    useconstraints : bool, default=False
        Apply constraints to the production simulation, defined by the constraints parameter
    constraints : dict, default={}
        A dictionary of atomselections and values of the constraint to be applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.
    adaptive : bool, default=False
        Set to True if making production runs for adaptive sampling.
    """

    def __init__(self):
        super().__init__()
        self._arg(
            "acemd",
            ":class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>`" " object",
            "Acemd class object",
            None,
            val.Object(Acemd),
        )
        self._arg(
            "runtime",
            "float",
            "Running time of the simulation.",
            25000,
            val.Number(float, "0POS"),
        )
        self._arg(
            "timeunits",
            "str",
            "Units for runtime. Can be 'steps', 'ns' etc.",
            "steps",
            val.String(),
        )
        self._arg(
            "temperature",
            "float",
            "Temperature of the thermostat in Kelvin",
            300,
            val.Number(float, "ANY"),
        )
        self._arg(
            "fb_k",
            "float",
            "Force constant of the flatbottom potential in kcal/mol/A^2. E.g. 5",
            0,
            val.Number(float, "ANY"),
        )
        self._arg(
            "fb_reference",
            "str",
            "Reference selection to use as dynamic center of the flatbottom box.",
            "none",
            val.String(),
        )
        self._arg(
            "fb_selection",
            "str",
            "Selection of atoms to apply the flatbottom potential",
            "none",
            val.String(),
        )
        self._arg(
            "fb_box",
            "list",
            "Position of the flatbottom box in term of the reference center given as "
            "[xmin, xmax, ymin, ymax, zmin, zmax]",
            [0, 0, 0, 0, 0, 0],
            val.Number(float, "ANY"),
            nargs=6,
        )
        self._arg(
            "useconstantratio",
            "bool",
            "For membrane protein simulations set it to true so that the barostat "
            "does not modify the xy aspect ratio.",
            False,
            val.Boolean(),
        )
        self._arg(
            "useconstraints",
            "bool",
            "Apply constraints to the production simulation, defined by the "
            "constraints parameter",
            False,
            val.Boolean(),
        )
        self._arg(
            "constraints",
            "dict",
            "A dictionary of atomselections and values of the constraint to be "
            "applied (in kcal/mol/A^2). Atomselects must be mutually exclusive.",
            {},
            val.Dictionary(key_type=str),
        )
        self._arg(
            "adaptive",
            "bool",
            "Set to True if making production runs for adaptive sampling.",
            False,
            val.Boolean(),
        )
        self._arg(
            "restraints",
            "list",
            "A list of restraint objects."
            "See :class:`AtomRestraint <htmd.mdengine.acemd.acemd.AtomRestraint>` and"
            ":class:`GroupRestraint <htmd.mdengine.acemd.acemd.GroupRestraint>`"
            ")",
            None,
            val.Object(_Restraint),
            nargs="*",
        )

        self.acemd = Acemd()
        self.acemd.binvelocities = None
        self.acemd.bincoordinates = "output.coor"
        self.acemd.extendedsystem = "output.xsc"
        self.acemd.coordinates = None
        self.acemd.structure = None
        self.acemd.parameters = None
        self.acemd.restart = "on"
        self.acemd.trajectoryfile = "output.xtc"
        self.acemd.trajectoryperiod = 25000
        self.acemd.timestep = 4
        self.acemd.switching = "on"
        self.acemd.switchdistance = 7.5
        self.acemd.cutoff = 9
        self.acemd.thermostat = "on"
        self.acemd.thermostatdamping = 0.1
        self.acemd.pme = "on"

    def _findFiles(self, inputdir):
        # Tries to find default files if the given don't exist
        defaults = {
            "coordinates": ("structure.pdb",),
            "structure": ("structure.psf", "structure.prmtop"),
            "parameters": ("parameters", "structure.prmtop"),
        }

        for field in defaults:
            userval = self.acemd.__dict__[field]
            if userval is not None and not os.path.exists(
                os.path.join(inputdir, userval)
            ):
                raise RuntimeError(
                    "Could not locate file {} set by the user for argument "
                    "Production.acemd.{}".format(os.path.join(inputdir, userval), field)
                )

            if self.acemd.__dict__[field] is None:
                for vv in defaults[field]:
                    if os.path.exists(os.path.join(inputdir, vv)):
                        self.acemd.__dict__[field] = vv
                        break

            if (
                userval is not None
                and self.acemd.__dict__[field] is not None
                and self.acemd.__dict__[field] != userval
            ):
                logger.warning(
                    "Could not locate structure file {}. Using {} instead.".format(
                        os.path.join(inputdir, userval),
                        os.path.join(inputdir, self.acemd.__dict__[field]),
                    )
                )
            elif self.acemd.__dict__[field] is None:
                raise RuntimeError(
                    "Could not locate any {f:} file in {i:}. "
                    "Please set the Production.acemd.{f:} property to "
                    "point to the {f:} file".format(f=field, i=inputdir)
                )

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.structure.endswith(".prmtop"):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None

    def _fb_potential2restraints(self, inputdir):
        from moleculekit.molecule import Molecule

        restraints = list()

        fb_box = np.array(self.fb_box)
        # convert fb_box to width
        width = list(
            np.concatenate(np.diff(np.array([fb_box[::2], fb_box[1::2]]), axis=0))
        )

        # If fb_box is not symmetrical
        if not np.all(fb_box[::2] == -fb_box[1::2]):
            # convert fb_box and fb_reference to fbcentre and width
            mol = Molecule(os.path.join(inputdir, self.acemd.structure))
            mol.read(os.path.join(inputdir, self.acemd.coordinates))
            fb_refcentre = (
                mol.get("coords", sel=self.fb_reference).mean(axis=0).squeeze()
            )

            fbcentre = list(
                np.around(
                    np.mean(np.array([fb_box[::2], fb_box[1::2]]), axis=0)
                    + fb_refcentre,
                    3,
                )
            )
            restraints.append(
                GroupRestraint(
                    self.fb_selection, width, [(self.fb_k, 0)], fbcentre=fbcentre
                )
            )
        else:
            restraints.append(
                GroupRestraint(
                    self.fb_selection,
                    width,
                    [(self.fb_k, 0)],
                    fbcentresel=self.fb_reference,
                )
            )

        return restraints

    def _constraints2restraints(self):

        restraints = list()
        for constr in sorted(self.constraints):
            restraints.append(AtomRestraint(constr, 0, [(self.constraints[constr], 0)]))

        return restraints

    def write(self, inputdir, outputdir):
        """Writes the production protocol and files into a folder.

        Parameters
        ----------
        inputdir : str
            Path to a directory containing the files produced by a equilibration process.
        outputdir : str
            Directory where to write the production setup files.
        """
        from moleculekit.molecule import Molecule
        from htmd.units import convert
        from htmd.builder.builder import detectCisPeptideBonds

        self._findFiles(inputdir)
        self._amberFixes()

        self.acemd.temperature = self.temperature
        self.acemd.thermostattemperature = self.temperature

        numsteps = convert(
            self.timeunits, "timesteps", self.runtime, timestep=self.acemd.timestep
        )
        self.acemd.run = str(numsteps)

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        detectCisPeptideBonds(inmol)

        if np.any(inmol.atomselect("lipids")) and not self.useconstantratio:
            logger.warning(
                "Lipids detected in input structure. We highly recommend setting useconstantratio=True "
                "for membrane simulations."
            )

        if self.restraints is not None:
            logger.info(
                "Using user-provided restraints and ignoring constraints and fb_potential"
            )
            self.acemd.restraints = self.restraints
        else:
            restraints = list()
            if self.fb_k > 0:
                logger.warning(
                    "Converting fb_potential to restraints. This is a convenience "
                    "functional conversion. We recommend to start using restraints."
                )
                restraints += self._fb_potential2restraints(inputdir)
            if self.useconstraints:
                logger.warning(
                    "Converting constraints to restraints. This is a convenience "
                    "functional conversion. We recommend to start using restraints"
                )
                restraints += self._constraints2restraints()
            else:
                if len(self.constraints) != 0:
                    logger.warning(
                        "You have setup constraints to {} but constraints are turned off. "
                        "If you want to use constraints, define "
                        "useconstraints=True".format(self.constraints)
                    )
            if len(restraints) != 0:
                self.acemd.restraints = restraints

        if self.useconstantratio:
            self.acemd.barostatconstratio = "on"

        if self.adaptive:
            self.acemd.binvelocities = None

        self.acemd.setup(inputdir, outputdir, overwrite=True)

    def addConstraint(self, atomselect, factor=1):
        """Convenience function for adding a new constraint to existing constraints.

        Parameters
        ----------
        atomselect : str
            Atom selection of atoms we want to constrain
        factor : float
            The scaling factor of the constraints applied to the atoms

        Example
        -------
        >>> eq.addConstraint('chain X', 0.3)
        """
        self.constraints[atomselect] = factor


import unittest


class _TestProduction(unittest.TestCase):
    def _cutfirstline(self, infile, outfile):
        # Cut out the first line of prmtop which has a build date in it
        with open(infile, "r") as fin:
            data = fin.read().splitlines(True)
        with open(outfile, "w") as fout:
            fout.writelines(data[1:])

    def _compareResultFolders(self, compare, tmpdir, pid):
        from glob import glob
        import os
        import filecmp

        ignore_ftypes = (".log", ".txt")
        files = []
        deletefiles = []
        for f in glob(os.path.join(compare, "*")):
            fname = os.path.basename(f)
            if os.path.splitext(f)[1] in ignore_ftypes:
                continue
            if f.endswith("prmtop"):
                self._cutfirstline(f, os.path.join(compare, fname + ".mod"))
                self._cutfirstline(
                    os.path.join(tmpdir, fname), os.path.join(tmpdir, fname + ".mod")
                )
                files.append(os.path.basename(f) + ".mod")
                deletefiles.append(os.path.join(compare, fname + ".mod"))
            else:
                files.append(os.path.basename(f))

        match, mismatch, error = filecmp.cmpfiles(tmpdir, compare, files, shallow=False)
        if len(mismatch) != 0 or len(error) != 0 or len(match) != len(files):
            raise RuntimeError(
                "Different results produced by amber.build for test {} between {} and {} in files {}.".format(
                    pid, compare, tmpdir, mismatch
                )
            )

        for f in deletefiles:
            os.remove(f)

    def test_acemd3(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import os

        pd = Production()
        pd.runtime = 4
        pd.timeunits = "ns"
        pd.temperature = 300
        pd.useconstraints = True
        pd.constraints = {
            "protein and name CA": 1,
            "protein and noh and not name CA": 0.1,
        }
        pd.fb_reference = "protein and name CA"
        pd.fb_selection = "resname MOL and noh"
        pd.fb_box = [-21, 21, -19, 19, 29, 30]
        pd.fb_k = 5
        tmpdir = tempname()
        pd.write(
            home(
                dataDir=os.path.join(
                    "test-protocols", "equilibration", "acemd3", "protLig", "postrun"
                )
            ),
            tmpdir,
        )

        # Compare with reference
        refdir = home(
            dataDir=os.path.join(
                "test-protocols", "production", "acemd3", "protLig", "prerun"
            )
        )
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, "*"))]
        self._compareResultFolders(refdir, tmpdir, "protLig")

    @unittest.skipUnless("ACE3ARG" in os.environ, "Untrusted PR")
    def test_run_water(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        from subprocess import check_output
        import subprocess
        import shutil
        import os

        acemd3exe = shutil.which("acemd3", mode=os.X_OK)
        if not acemd3exe:
            raise NameError(
                "Could not find acemd3, or no execute permissions are given"
            )

        for system in ["amber-equil-completed", "charmm-equil-completed"]:
            pd = Production()
            pd.runtime = 5
            pd.timeunits = "steps"
            pd.temperature = 300
            pd.constraints = {}
            # Set these down for tiny box size of water
            pd.acemd.cutoff = 3
            pd.acemd.switchdistance = 2
            ######
            tmpdir = tempname()
            pd.write(
                home(dataDir=os.path.join("test-acemd", "tiny-water", system)), tmpdir
            )
            try:
                res = check_output(
                    ["acemd3", "--platform", "CPU", os.getenv("ACE3ARG")], cwd=tmpdir
                )
            except subprocess.CalledProcessError as exc:
                assert (
                    False
                ), f'Failed to run due to error: {exc}\n\n ---> Error log:\n\n{exc.output.decode("ascii")}'
            res = res.decode("utf-8").strip()
            print(res)
            assert res.endswith("Completed simulation!"), "Failed at system " + system


if __name__ == "__main__":
    unittest.main(verbosity=2)
