# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from abc import ABC, abstractmethod
import logging
import os
import re

from moleculekit.periodictable import periodictable
import numpy as np
from protocolinterface import ProtocolInterface, val

from jobqueues.localqueue import LocalCPUQueue
from htmd.queues.playqueue import PlayQueue

logger = logging.getLogger(__name__)


class QMResult:
    """
    Class containg QM calculation results

    Attributes
    ----------
    errored : bool
        If QM failed, it is set to True, overwise False.
    energy: float
        Total QM energy in kcal/mol
    coords : numpy.ndarray
        Atomic coordinates in Angstrom. The array shape is (number_of_atoms, 3, 1).
    dipole : list
        Dipole moment in Debye. The list has 4 elements corresponding to x, y, z conponents, and the total.
    quadrupole: list
        Quadrople moment in Debye. The list has 6 elements.
    mulliken : list
        Mulliken charges in electron charges. The list has an element for each atom.
    esp_points : numpy.ndarray
        Point coordinates (in Angstrom) where ESP values are computed. The array shape is (number_of_points, 3).
    esp_values : numpy.ndarray
        ESP values in elementary_charge/Angstrom. The array shape is (number_of_points,)
    charge : int
        The total charge of the molecule in electron charges.
    """

    def __init__(self):

        self.errored = False
        self.energy = None
        self.coords = None
        self.dipole = None
        self.quadrupole = None
        self.mulliken = None
        self.esp_points = None
        self.esp_values = None
        self.charge = None


class QMBase(ABC, ProtocolInterface):
    """
    Abstract base class to set up and run QM calculations
    """

    THEORIES = (
        "HF",
        "BLYP",
        "PBE",
        "B3LYP",
        "PBE0",
        "B2PLYP",
        "wB97",
        "wB97X",
        "wB97X-D",
    )
    CORRECTIONS = ("none", "D", "D3")
    BASIS_SETS = (
        "3-21G",
        "6-31G",
        "6-31G*",
        "6-31G**",
        "6-31+G",
        "6-31+G*",
        "6-31+G**",
        "6-31++G",
        "6-31++G*",
        "6-31++G**",
        "6-311G",
        "6-311G*",
        "6-311G**",
        "6-311+G",
        "6-311+G*",
        "6-311+G**",
        "6-311++G",
        "6-311++G*",
        "6-311++G**",
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
    )
    SOLVENTS = ("vacuum", "PCM")

    @staticmethod
    def substituteBasisSet(element, basis_set):
        """
        Substitute basis sets, if a given element does not have a requested basis set.

        Pople (6-31G, etc.) and Dunning (cc-pVDZ, etc.) basis sets are replaced with
        Ahlrichs-Karlsruhe (def2-SVP, etc.) basis sets for elements heavier than Ar.

        Arguments
        ---------
        element : str
            Element symbol
        basis_set : str
            Requested basis set name

        Return
        ------
        new_basis_set : str
            Substituted basis set name

        References
        ----------
        https://www.basissetexchange.org
        http://www.psicode.org/psi4manual/master/basissets_byelement.html

        Examples
        --------

        >>> from htmd.qm.base import QMBase
        ffevaluate module is in beta version

        >>> QMBase.substituteBasisSet('F', '3-21G')
        '3-21G'

        >>> QMBase.substituteBasisSet('F', '6-31+G*')
        '6-31+G*'

        >>> QMBase.substituteBasisSet('F', '6-311++G**')
        '6-311++G**'

        >>> QMBase.substituteBasisSet('F', 'cc-pVDZ')
        'cc-pVDZ'

        >>> QMBase.substituteBasisSet('F', 'aug-cc-pVTZ')
        'aug-cc-pVTZ'

        >>> QMBase.substituteBasisSet('Br', '3-21G')
        'def2-SV(P)'

        >>> QMBase.substituteBasisSet('Br', '6-31+G*')
        'def2-SVPD'

        >>> QMBase.substituteBasisSet('Br', '6-311++G**')
        'def2-TZVPD'

        >>> QMBase.substituteBasisSet('Br', 'cc-pVDZ')
        'def2-SVP'

        >>> QMBase.substituteBasisSet('Br', 'aug-cc-pVTZ')
        'def2-TZVPD'
        """

        element = periodictable[element]
        if basis_set not in QMBase.BASIS_SETS:
            raise ValueError(f"Unrecognized basis sets {basis_set}")

        new_basis_set = basis_set

        match_dunning = re.match("^(|aug-)cc-pV(D|T|Q)Z$", basis_set)
        match_pople = re.match("^(3-21|6-31|6-311)([\+]{0,2})G[\*]{0,2}$", basis_set)

        if match_dunning:
            if element.number > 18:
                core = {"D": "S", "T": "TZ", "Q": "QZ"}[match_dunning.group(2)]
                diffuse = "D" if match_dunning.group(1) else ""
                new_basis_set = f"def2-{core}VP{diffuse}"
        elif match_pople:
            if element.number > 18:
                if match_pople.group(1) == "3-21":
                    core = "S"
                    polar = "(P)"
                    diffuse = ""
                elif match_pople.group(1) in ("6-31", "6-311"):
                    core = {"6-31": "S", "6-311": "TZ"}[match_pople.group(1)]
                    polar = "P" if match_pople.group(2) else ""
                    diffuse = (
                        "D" if match_pople.group(2) and match_pople.group(2) else ""
                    )
                else:
                    raise ValueError()
                new_basis_set = f"def2-{core}V{polar}{diffuse}"
        else:
            raise RuntimeError()

        if basis_set != new_basis_set:
            logger.info(
                f"Basis set substitution for {element.symbol}: {basis_set} --> {new_basis_set}"
            )

        return new_basis_set

    def __init__(self):
        from moleculekit.molecule import Molecule

        super().__init__()

        self._arg(
            "molecule",
            ":class: `moleculekit.molecule.Molecule`",
            "Molecule",
            default=None,
            validator=val.Object(Molecule),
            required=True,
        )
        self._arg(
            "charge",
            "int",
            "Charge of the molecule in electron charges",
            default=None,
            validator=val.Number(int, "ANY"),
            required=True,
        )
        self._arg(
            "multiplicity",
            "int",
            "Multiplicity of the molecule",
            default=1,
            validator=val.Number(int, "POS"),
        )
        self._arg(
            "theory",
            "str",
            "Level of theory",
            default="B3LYP",
            validator=val.String(),
            valid_values=self.THEORIES,
        )
        self._arg(
            "correction",
            "str",
            "Empirical dispersion correction",
            default="none",
            validator=val.String(),
            valid_values=self.CORRECTIONS,
        )
        self._arg(
            "basis",
            "str",
            "Basis set",
            default="6-31G*",
            validator=val.String(),
            valid_values=self.BASIS_SETS,
        )
        self._arg(
            "solvent",
            "str",
            "Implicit solvent",
            default="vacuum",
            validator=val.String(),
            valid_values=self.SOLVENTS,
        )
        self._arg(
            "esp_points",
            ":class: `numpy.ndarray`",
            "Point to calculate ESP",
            default=None,
            nargs="*",
        )  # TODO implement validator
        self._arg(
            "optimize",
            "boolean",
            "Optimize geometry",
            default=False,
            validator=val.Boolean(),
        )
        self._arg(
            "restrained_dihedrals",
            ":class: `numpy.ndarray`",
            "List of restrained dihedrals (0-based indices)",
            default=None,
            nargs="*",
        )  # TODO implement validator
        self._arg(
            "queue",
            ":class:`SimQueue <jobqueues.simqueue.SimQueue>` object",
            "Queue object used to run simulations",
            default=LocalCPUQueue(),
        )
        self._arg(
            "directory", "str", "Working directory", default=".", validator=val.String()
        )

    @property
    @abstractmethod
    def _command(self):
        pass

    @abstractmethod
    def _completed(self, directory):
        pass

    @abstractmethod
    def _writeInput(self, directory, iframe):
        pass

    @abstractmethod
    def _readOutput(self, directory):
        pass

    def setup(self):
        """
        Setup QM calculations for all the frames of the molecule
        """

        # Set up the molecule
        # TODO remove molecule coping!
        self._molecule = self.molecule.copy()
        self._nframes = self._molecule.coords.shape[2]
        self._natoms = self._molecule.coords.shape[0]
        if self.charge is None:
            self.charge = int(round(self._molecule.charge.sum()))

        # Set up ESP points
        if self.esp_points is not None:
            # TODO move to a validator
            if self.esp_points.shape[1] != 3:
                raise ValueError("ESP point array must be npoints x 3")
            if self._molecule.coords.shape[2] != 1:
                raise ValueError(
                    "Can only specift ESP point array with a single frame of coords"
                )

        # Set up restrained dihedrals
        self._restrained_dihedrals = None
        if self.restrained_dihedrals is not None:
            self._restrained_dihedrals = (
                self.restrained_dihedrals + 1
            )  # Convert to 1-based indices

        # Create directories and write inputs
        self._directories = []
        for iframe in range(self._nframes):

            # Create a directory
            directory = os.path.join(self.directory, "%05d" % iframe)
            os.makedirs(directory, exist_ok=True)
            self._directories.append(directory)

            if not self._completed(directory):

                # Write input files
                self._writeInput(directory, iframe)

                # Write a point file for ESP
                if self.esp_points is not None:
                    np.savetxt(
                        os.path.join(directory, "grid.dat"), self.esp_points, fmt="%f"
                    )

                # Write a run script
                script = os.path.join(directory, "run.sh")
                with open(script, "w") as f:
                    f.write("#!/bin/sh\n\n%s\n" % self._command)
                os.chmod(script, 0o700)

    def submit(self):
        """
        Submit QM calculations for all the frames of the molecule
        """

        for directory in self._directories:
            if not self._completed(directory):
                self.queue.submit(directory)

    def retrieve(self):
        """
        Retrieve QM calculation results for all the frames of the molecule

        Return
        ------
        results : list
            List of QMResult objects (one for each molecule frames).
        """

        # Wait only if there is something to wait for
        # TODO: queue object should handle this logic
        if self.queue._dirs:
            if isinstance(self.queue, (LocalCPUQueue, PlayQueue)):
                self.queue.wait(sentinel=False)
            else:
                self.queue.wait(sentinel=True)
            self.queue.retrieve()

        # Read output files
        results = [self._readOutput(directory) for directory in self._directories]
        for res in results:
            res.charge = self.charge

        return results

    def run(self):
        """
        Run QM calculations for all the frames of the molecule and return results.

        The method generates input files according to the attributes, submits jobs to the selected queue,
        waits for the calculations to finish, and retrieves the results.

        Example
        -------

        >>> result = qm.run() # doctest: +SKIP

        It is equivalent to:
        >>> qm.setup() # doctest: +SKIP
        >>> qm.submit() # doctest: +SKIP
        >>> result = qm.retrieve() # doctest: +SKIP

        Return
        ------
        results : list
            List of QMResult objects (one for each molecule frames).
        """

        self.setup()  # doctest: +SKIP
        self.submit()  # doctest: +SKIP
        return self.retrieve()  # doctest: +SKIP


if __name__ == "__main__":

    import doctest
    import sys

    import htmd

    sys.exit(doctest.testmod().failed)
