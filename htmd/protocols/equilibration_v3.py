# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.apps.acemd import Acemd as Acemd2
from htmd.mdengine.acemd.acemd import Acemd, _Restraint, GroupRestraint, AtomRestraint
from htmd.config import _config
from protocolinterface import ProtocolInterface, val
import os
import numpy as np
import logging

logger = logging.getLogger(__name__)


class Equilibration(ProtocolInterface):
    """Equilibration protocol v3

    Equilibration protocol for globular and membrane proteins
    Supporst extra restraints like a flatbottom potential box to retain a ligand
    for example within this box.

    Parameters
    ----------
    runtime : float, default=0
        Running time of the simulation.
    timeunits : str, default='steps'
        Units for time arguments. Can be 'steps', 'ns' etc.
    temperature : float, default=300
        Temperature of the thermostat in Kelvin
    restraints : list, default=None
        A list of restraint objects. See :class:`AtomRestraint <htmd.mdengine.acemd.acemd.AtomRestraint>` and:class:`GroupRestraint<htmd.mdengine.acemd.acemd.GroupRestraint>`)
    useconstantratio : bool, default=False
        For membrane protein simulations set it to true so that the barostat does not modify the xy aspect ratio.
    restraintsteps : int, default=None
        Number of initial steps to apply restraints in units of 4fs. Defaults to half the simulation time.

    Example
    -------
    >>> from htmd.protocols.equilibration_v3 import Equilibration
    >>> from htmd.mdengine.acemd.acemd import GroupRestraint
    >>> md = Equilibration()
    >>> md.runtime = 4
    >>> md.timeunits = 'ns'
    >>> md.temperature = 300
    >>> md.useconstantratio = True  # only for membrane sims
    Use a 10A flat bottom potential to prevent the ligand from diffusing from original position during equilibration
    >>> width = np.array([10, 10, 10])
    >>> flatbot = GroupRestraint('segname L and noh', width, [(5, '0ns')])
    Add also the default restraints for protein
    >>> mol = Molecule("./build/structure.pdb")
    >>> md.restraints = [flatbot,] + md.defaultEquilRestraints('2ns', mol)
    >>> md.write('./build','./equil')
    """

    def __init__(self, _version=_config["acemdversion"]):
        if _version == 2:
            raise RuntimeError("Equilibration v3 only supports _version=3")

        super().__init__()
        self._arg(
            "acemd",
            ":class:`Acemd <htmd.mdengine.acemd.acemd.Acemd>` object",
            "Acemd class object",
            None,
            val.Object([Acemd2, Acemd]),
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
            "Units for time arguments. Can be 'steps', 'ns' etc.",
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
            "useconstantratio",
            "bool",
            "For membrane protein simulations set it to true so that the barostat "
            "does not modify the xy aspect ratio.",
            False,
            val.Boolean(),
        )
        self._arg(
            "restraintsteps",
            "int",
            "Number of initial steps to apply restraints in units of 4fs. Defaults "
            "to half the simulation time.",
            None,
            val.Number(int, "ANY"),
        )
        self._arg(
            "restraints",
            "list",
            "A list of restraint objects. If None will apply defaultEquilRestraints decaying over half the runtime. If no restraints are required set to empty list []."
            "See :class:`AtomRestraint <htmd.mdengine.acemd.acemd.AtomRestraint>` and"
            ":class:`GroupRestraint <htmd.mdengine.acemd.acemd.GroupRestraint>`)",
            None,
            val.Object(_Restraint),
            nargs="*",
        )

        self.acemd = Acemd()
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
        self.acemd.thermostatdamping = 1
        self.acemd.pme = "on"
        self.acemd.barostat = "on"
        self.acemd.barostatpressure = 1.01325
        self.acemd.minimize = 500

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
                self.acemd.__dict__[field] = None

            if self.acemd.__dict__[field] is None:
                for val in defaults[field]:
                    if os.path.exists(os.path.join(inputdir, val)):
                        self.acemd.__dict__[field] = val
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
                    "Please set the {name:}.acemd.{f:} property to "
                    "point to the {f:} file".format(
                        f=field, i=inputdir, name=self.__class__.__name__
                    )
                )

    def _amberFixes(self):
        # AMBER specific fixes
        if self.acemd.parameters.endswith("structure.prmtop"):
            self.acemd.parmfile = self.acemd.parameters
            self.acemd.parameters = None

    def defaultEquilRestraints(self, decay, mol=None):
        """Get the default equilibration restraints

        Parameters
        ----------
        decay : str
            The restrains will get scaled to 0 over this much time.

        Returns
        -------
        restraints : list
            A list of default protein restraints

        Examples
        --------
        >>> md = Equilibration()
        >>> res = md.defaultEquilRestraints('20ns')
        """
        caatoms = AtomRestraint("protein and name CA", 0, [(1, 0), (0, decay)])
        notcaatoms = AtomRestraint(
            "protein and noh and not name CA", 0, [(0.1, 0), (0, decay)]
        )
        nucleic = AtomRestraint("nucleic and backbone", 0, [(1, 0), (0, decay)])
        nucleicside = AtomRestraint(
            "nucleic and not backbone and noh", 0, [(0.1, 0), (0, decay)]
        )
        restraints = [caatoms, notcaatoms, nucleic, nucleicside]
        if mol is not None:
            restraints = [r for r in restraints if mol.atomselect(r.selection).sum()]

        return restraints

    def write(self, inputdir, outputdir):
        """Write the equilibration protocol

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

        from moleculekit.molecule import Molecule

        self._findFiles(inputdir)
        self._amberFixes()

        from htmd.units import convert

        numsteps = convert(
            self.timeunits, "timesteps", self.runtime, timestep=self.acemd.timestep
        )

        self.acemd.temperature = self.temperature
        self.acemd.thermostattemperature = self.temperature
        self.acemd.run = str(numsteps)

        pdbfile = os.path.join(inputdir, self.acemd.coordinates)
        inmol = Molecule(pdbfile)

        from htmd.builder.builder import detectCisPeptideBonds

        detectCisPeptideBonds(inmol)

        if np.any(inmol.atomselect("lipids")) and not self.useconstantratio:
            logger.warning(
                "Lipids detected in input structure. We highly recommend setting useconstantratio=True "
                "for membrane simulations."
            )

        if self.restraintsteps is None:
            constrsteps = int(numsteps / 2)
        else:
            constrsteps = int(self.restraintsteps)

        if self.restraints is not None:
            self.acemd.restraints = self.restraints
        else:
            self.acemd.restraints = self.defaultEquilRestraints(constrsteps, mol=inmol)

        if self.acemd.boxsize is None and self.acemd.extendedsystem is None:
            coords = inmol.get("coords", sel="water")
            if coords.size == 0:  # It's a vacuum simulation
                coords = inmol.get("coords", sel="all")
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
                dim += 12.0
            else:
                dim = np.max(coords, axis=0) - np.min(coords, axis=0)
            self.acemd.boxsize = "{} {} {}".format(dim[0], dim[1], dim[2])

        if self.useconstantratio:
            self.acemd.useconstantratio = "on"

        self.acemd.setup(inputdir, outputdir, overwrite=True)


import unittest


class _TestEquilibration(unittest.TestCase):
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
        from htmd.units import convert
        from moleculekit.molecule import Molecule
        from glob import glob
        import os
        from htmd.mdengine.acemd.acemd import GroupRestraint

        eq = Equilibration()
        eq.runtime = 4
        eq.timeunits = "ns"
        eq.temperature = 300
        eq.restraintsteps = convert(
            eq.timeunits, "timesteps", eq.runtime, timestep=eq.acemd.timestep
        )
        ligres = GroupRestraint(
            "resname MOL and noh",
            [42, 38, 1],
            [(5, 0)],
            fbcentre=[-0.178, -0.178, 29.195],
        )
        mol = Molecule(
            home(
                dataDir=os.path.join(
                    "test-protocols", "build", "protLig", "structure.pdb"
                )
            )
        )
        eq.restraints = eq.defaultEquilRestraints(1000000, mol=mol) + [ligres]
        tmpdir = tempname()
        eq.write(
            home(dataDir=os.path.join("test-protocols", "build", "protLig")), tmpdir
        )

        # Compare with reference
        refdir = home(
            dataDir=os.path.join(
                "test-protocols", "equilibration", "acemd3", "protLig", "prerun"
            )
        )
        files = [os.path.basename(f) for f in glob(os.path.join(refdir, "*"))]
        self._compareResultFolders(refdir, tmpdir, "protLig")

    @unittest.skipUnless("ACE3ARG" in os.environ, "Untrusted PR")
    def test_run_water(self):
        from htmd.util import tempname
        from htmd.home import home
        from glob import glob
        import subprocess
        from subprocess import check_output
        import shutil
        import os

        acemd3exe = shutil.which("acemd3", mode=os.X_OK)
        if not acemd3exe:
            raise NameError(
                "Could not find acemd3, or no execute permissions are given"
            )

        for system in ["amber-build", "charmm-build"]:
            eq = Equilibration()
            eq.runtime = 5
            eq.timeunits = "steps"
            eq.temperature = 300
            eq.restraints = []
            # Set these down for tiny box size of water
            eq.acemd.cutoff = 3
            eq.acemd.switchdistance = 2
            eq.acemd.minimize = 50  # Do few steps to finish fast
            ######
            tmpdir = tempname()
            eq.write(
                home(dataDir=os.path.join("test-acemd", "tiny-water", system)), tmpdir
            )
            print(tmpdir)
            try:
                res = check_output(
                    ["acemd3", "--platform", "CPU", os.getenv("ACE3ARG")], cwd=tmpdir
                )
            except subprocess.CalledProcessError as exc:
                assert (
                    False
                ), f'Failed to run due to error: {exc}\n\n ---> Error log:\n\n{exc.output.decode("ascii")}'
            res = res.decode("utf-8").strip()
            assert "Completed minimization" in res, "Failed at system " + system
            assert res.endswith("Completed simulation!"), "Failed at system " + system


if __name__ == "__main__":
    unittest.main(verbosity=2)
