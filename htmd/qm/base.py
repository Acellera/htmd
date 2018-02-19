# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from abc import ABC, abstractmethod

from protocolinterface import ProtocolInterface, val
from htmd.queues.simqueue import SimQueue
from htmd.queues.localqueue import LocalCPUQueue


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


class QMBase(ABC, ProtocolInterface):
    """
    Abstract base class to set up and run QM calculations
    """

    THEORIES = ('HF', 'BLYP', 'PBE', 'B3LYP', 'PBE0', 'B2PLYP')
    CORRECTIONS = ('none', 'D', 'D3')
    BASIS_SETS = ('3-21G',
                  '6-31G',  '6-31G*',  '6-31G**',  '6-31+G',  '6-31+G*',  '6-31+G**',  '6-31++G',  '6-31++G*',  '6-31++G**',
                  '6-311G', '6-311G*', '6-311G**', '6-311+G', '6-311+G*', '6-311+G**', '6-311++G', '6-311++G*', '6-311++G**',
                  'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ')
    SOLVENTS = ('vacuum', 'PCM')

    def __init__(self):
        from htmd.molecule.molecule import Molecule

        super().__init__()

        self._arg('molecule', ':class: `htmd.molecule.molecule.Molecule`', 'Molecule',
                  default=None, validator=val.Object(Molecule), required=True)
        self._arg('netcharge', 'int', 'Net charge', default=None, validator=val.Number(int, 'ANY'), required=True)
        self._arg('multiplicity', 'int', 'Multiplicity of the molecule',
                  default=1, validator=val.Number(int, 'POS'))
        self._arg('theory', 'str', 'Level of theory',
                  default='B3LYP', validator=val.String(), valid_values=self.THEORIES)
        self._arg('correction', 'str', 'Empirical dispersion correction',
                  default='none', validator=val.String(), valid_values=self.CORRECTIONS)
        self._arg('basis', 'str', 'Basis set',
                  default='6-31G*', validator=val.String(), valid_values=self.BASIS_SETS)
        self._arg('solvent', 'str', 'Implicit solvent',
                  default='vacuum', validator=val.String(), valid_values=self.SOLVENTS)
        self._arg('esp_points', ':class: `numpy.ndarray`', 'Point to calculate ESP',
                  default=None, nargs='*')  # TODO implement validator
        self._arg('optimize', 'boolean', 'Optimize geometry',
                  default=False, validator=val.Boolean())
        self._arg('restrained_dihedrals', ':class: `numpy.ndarray`',
                  'List of restrained dihedrals (0-based indices)',
                  default=None, nargs='*')  # TODO implement validator
        self._arg('queue', ':class:`SimQueue <htmd.queues.simqueue.SimQueue>` object',
                  'Queue object used to run simulations',
                  default=LocalCPUQueue())
        self._arg('directory', 'str', 'Working directory',
                  default='.', validator=val.String())

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
        self._charge = self.netcharge

        # Set up ESP points
        if self.esp_points is not None:
            # TODO move to a validator
            if self.esp_points.shape[1] != 3:
                raise ValueError("ESP point array must be npoints x 3")
            if self._molecule.coords.shape[2] != 1:
                raise ValueError("Can only specift ESP point array with a single frame of coords")

        # Set up restrained dihedrals
        self._restrained_dihedrals = None
        if self.restrained_dihedrals is not None:
            self._restrained_dihedrals = self.restrained_dihedrals + 1  # Convert to 1-based indices

        # Create directories and write inputs
        self._directories = []
        for iframe in range(self._nframes):

            # Create a directory
            directory = os.path.join(self.directory, '%05d' % iframe)
            os.makedirs(directory, exist_ok=True)
            self._directories.append(directory)

            if not self._completed(directory):

                # Write input files
                self._writeInput(directory, iframe)

                # Write a point file for ESP
                if self.esp_points is not None:
                    np.savetxt(os.path.join(directory, 'grid.dat'), self.esp_points, fmt='%f')

                # Write a run script
                script = os.path.join(directory, 'run.sh')
                with open(script, 'w') as f:
                    f.write('#!/bin/sh\n\n%s\n' % self._command)
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
            if isinstance(self.queue, LocalCPUQueue):
                self.queue.wait(sentinel=False) # TODO: sentinel is still broken
            else:
                self.queue.wait(sentinel=True)
            self.queue.retrieve()

        # Read output files
        results = [self._readOutput(directory) for directory in self._directories]

        return results

    def run(self):
        """
        Run QM calculations for all the frames of the molecule and return results.

        The method generates input files according to the attributes, submits jobs to the selected queue,
        waits for the calculations to finish, and retrieves the results.

        Example
        -------

        >>> result = qm.run()

        It is equivalent to:
        >>> qm.setup()
        >>> qm.submit()
        >>> result = qm.retrieve()

        Return
        ------
        results : list
            List of QMResult objects (one for each molecule frames).
        """

        self.setup()
        self.submit()
        return self.retrieve()


if __name__ == '__main__':
    pass
