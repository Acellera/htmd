# (c) 2015-2019 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from scipy import constants as const

from moleculekit.molecule import Molecule
from htmd.qm.base import QMBase, QMResult


class TeraChem(QMBase):
    """
    Class to set up and run QM calculations with TeraChem

    Capabilities
    ------------
    - Single-point energy calculations with HF and DFT.
    - Electronic properties: dipole and Mulliken charges
    - Geometry optimization with/without dihedral restraints

    Attributes
    ----------
    molecule : :class:`Molecule`
        Molecule
    multiplity : int, default = 1
        Multiplicity of the molecule (1 = singlet, 2 = doublet, 3 = triplet, etc.)
    theory : str, default = 'BLYP3'
        Level of theory. A full list is at Psi4.THEORIES
    basis : str, default = '6-31G*'
        Basis set. A full list is at Psi4.BASIS_SETS
    correction : str, default = 'none'
        Empirical dispersion correction. A full list is at Psi4.CORRECTIONS
    solvent : str, default = 'vacuum'
         Implicit solvent model. A full isit is at Psi4.SOLVENTS
    esp_points : nunpy.ndarray, default = None
        Point coordinates to compute ESP values. The array shape has to be (number_of_points, 3)
    optimize : bool, default = False
        Run geometery optimization
    restrained_dihedrals : nunpy.ndarray, default = None
        List of restrained dihedrals (0-based atom indices). The array shape has to be (number_of_restrained_dihedrals, 4).
    queue : :class:`SimQueue`, defualt = :class:`LocalCPUQueue`
        Queue object used to run simulations
    directory : str, default = '.'
        Working directoy

    Examples
    --------

    Create an object of H2 molecule
    >>> import os
    >>> from htmd.home import home
    >>> from moleculekit.molecule import Molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2-0.74.mol2')
    >>> mol = Molecule(molFile)

    Create a Psi4 object
    >>> from htmd.qm import TeraChem
    >>> qm = TeraChem()
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.terachem.TeraChem object at 0x...>

    Run single-point QM calculation of H2 with BLYP and cc-pVDZ
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as tmp:
    ...     qm.molecule = mol
    ...     qm.theory = 'BLYP'
    ...     qm.basis = 'cc-pVDZ'
    ...     qm.directory = tmp
    ...     result = qm.run()

    The QM results are returned as a list of htmd.qm.QMResult objects. See htmd.qm.QMResult documentation for details.
    >>> result # doctest: +ELLIPSIS
    [<htmd.qm.base.QMResult object at 0x...>]
    >>> result[0].errored
    False
    >>> result[0].energy # doctest: +ELLIPSIS
    -728.97068164...
    >>> result[0].mulliken # doctest: +ELLIPSIS
    [...0.0, ...0.0]

    Run the geometry optimization of H2 with BLYP, but change basis to 3-21G.
    NOTE: the `directory` attribut needs to be set to empty or non-existing directory, overwise the previous
    calculation results are read.
    >>> with TemporaryDirectory() as tmp:
    ...     qm.basis = '3-21G'
    ...     qm.optimize = True
    ...     qm.directory = tmp
    ...     result = qm.run()

    The initial coordinates and optimize coordinates of H2 can be compared
    >>> mol.coords
    array([[[ 0.  ],
            [ 0.  ],
            [-0.37]],
    <BLANKLINE>
           [[ 0.  ],
            [ 0.  ],
            [ 0.37]]], dtype=float32)
    >>> result[0].coords
    array([[[ 0.        ],
            [ 0.        ],
            [-0.37542403]],
    <BLANKLINE>
           [[-0.        ],
            [-0.        ],
            [ 0.37542403]]], dtype=float32)

    The QM calculations run using LocalCPUQueue by default, but this can be changed to the others.
    >>> from jobqueues.slurmqueue import SlurmQueue
    >>> qm.queue = SlurmQueue() # doctest: +SKIP
    """

    @property
    def _command(self):
        if "TeraChem" not in os.environ:
            raise RuntimeError("Environment variable TeraChem is undefined")
        if not os.path.exists(os.path.join(os.environ["TeraChem"], "SetTCVars.sh")):
            raise RuntimeError("$TeraChem/SetTCVars.sh does not exist")
        return (
            "source $TeraChem/SetTCVars.sh && " "terachem terachem.in &> terachem.out"
        )

    def _completed(self, directory):

        outFile = os.path.join(directory, "terachem.out")

        if not os.path.exists(outFile):
            return False

        with open(outFile) as fd:
            for line in fd.readlines():
                if line.startswith(" Job terminated:"):
                    return True

        return False

    def _writeInput(self, directory, frame):

        xyzFile = os.path.join(directory, "terachem.xyz")
        self.molecule.frame = frame
        self.molecule.write(xyzFile)

        with open(os.path.join(directory, "terachem.in"), "w") as f:

            f.write(
                f'method {("R" if self.multiplicity == 1 else "U") + self.theory}\n'
            )
            f.write(f"basis {self.basis}\n")
            if self.correction != "none":
                f.write(
                    f'dispersion {"D2" if self.correction == "D" else self.correction}\n'
                )
            if self.solvent == "PCM":
                f.write("pcm cosmo\n")

            f.write(f"charge {self.charge}\n")
            f.write(f"spinmult {self.multiplicity}\n")
            f.write("coordinates terachem.xyz\n")

            if self.esp_points is not None:
                raise NotImplemented("ESP is not available")

            if self.optimize:
                f.write("new_minimizer yes\n")
                if self._restrained_dihedrals is not None:
                    f.write("$constraint_freeze\n")
                    for dihedral in self._restrained_dihedrals:
                        f.write(f'  dihedral {"_".join(map(str, dihedral))}\n')
                    f.write("$end\n")
                f.write("run minimize\n")
            else:
                f.write("run energy\n")

            f.write("precision double\n")

            f.write(f"end\n")

    def _readOutput(self, directory):

        result = QMResult()
        result.completed = True

        with open(os.path.join(directory, "terachem.out")) as fd:
            for line in fd.readlines():
                if line.startswith("FINAL ENERGY:"):
                    result.energy = float(line.split()[2])
                    result.energy *= const.physical_constants["Hartree energy"][0] / (
                        const.kilo * const.calorie / const.Avogadro
                    )  # Hartree --> kcal/mol
                if line.startswith("DIPOLE MOMENT:"):
                    tokens = line.split()
                    result.dipole = [
                        float(tokens[2][1:-1]),
                        float(tokens[3][:-1]),
                        float(tokens[4][:-1]),
                        float(tokens[7][:-1]),
                    ]

        mullFile = os.path.join(directory, "scr", "charge_mull.xls")
        if os.path.exists(mullFile):
            result.mulliken = list(np.loadtxt(mullFile, usecols=(2,)))

        if result.energy is None or result.dipole is None:
            result.errored = True

        geomFile = os.path.join(
            directory, "scr", "optim.xyz" if self.optimize else "xyz.xyz"
        )
        if os.path.exists(geomFile):
            result.coords = Molecule(geomFile).coords[:, :, -1][:, :, None]
        else:
            result.errored = True

        if self.esp_points is not None:
            raise NotImplemented("ESP is not available")

        return result


if __name__ == "__main__":

    import sys
    import doctest

    # There is no TeraChem on Travis
    if "TRAVIS" not in os.environ:
        sys.exit(doctest.testmod().failed)
