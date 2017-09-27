# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import pickle
import logging
import numpy as np
import nlopt
from scipy.spatial.distance import cdist

from htmd.molecule.util import dihedralAngle
from htmd.qm.base import QMBase, QMResult
from htmd.parameterization.ffevaluate import FFEvaluate

logger = logging.getLogger(__name__)


class FakeQM(QMBase):
    """
    Imitation of QM calculations with MM

    >>> import os
    >>> import numpy as np
    >>> from tempfile import TemporaryDirectory
    >>> from htmd.home import home
    >>> from htmd.molecule.util import dihedralAngle
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> from htmd.qm.fake import FakeQM

    Create a molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2) # doctest: +ELLIPSIS
    Dihedral 0: 2-0-1-3
    ...

    Run a single-point energy and ESP calculation
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = FakeQM()
    ...     qm.molecule = mol
    ...     qm.esp_points = np.array([[1., 1., 1.]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]

    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.fake.FakeQM object at ...>
    >>> result # doctest: +ELLIPSIS
    <htmd.qm.base.QMResult object at ...
    >>> result.errored
    False
    >>> result.energy # doctest: +ELLIPSIS
    8.394800...
    >>> result.esp_points
    array([[ 1.,  1.,  1.]])
    >>> result.esp_values # doctest: +ELLIPSIS
    array([ 0.371352...])
    >>> dihedralAngle(result.coords[[2, 0, 1, 3], :, 0]) # doctest: +ELLIPSIS
    89.999544...

    Run a minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = FakeQM()
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    >>> result.energy # doctest: +ELLIPSIS
    7.729598...
    >>> dihedralAngle(result.coords[[2, 0, 1, 3], :, 0]) # doctest: +ELLIPSIS
    101.44411...

    Run a constrained minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = FakeQM()
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.restrained_dihedrals = np.array([[2, 0, 1, 3]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    >>> result.energy # doctest: +ELLIPSIS
    7.879684...
    >>> dihedralAngle(result.coords[[2, 0, 1, 3], :, 0]) # doctest: +ELLIPSIS
    89.535382...
    """

    # Fake implementations of the abstract methods
    def _command(self): pass
    def _writeInput(self, directory, iframe): pass
    def _readOutput(self, directory): pass

    def _completed(self, directory):
        return os.path.exists(os.path.join(directory, 'data.pkl'))

    def run(self):

        ff = FFEvaluate(self.molecule)

        results = []
        for iframe in range(self.molecule.numFrames):
            self.molecule.frame = iframe

            directory = os.path.join(self.directory, '%05d' % iframe)
            os.makedirs(directory, exist_ok=True)
            pickleFile = os.path.join(directory, 'data.pkl')

            if self._completed(directory):
                with open(pickleFile, 'rb') as fd:
                    result = pickle.load(fd)
                logger.info('Loading QM data from %s' % pickleFile)

            else:
                result = QMResult()
                result.errored = False
                result.coords = self.molecule.coords[:, :, iframe:iframe + 1].copy()

                if self.optimize:
                    optimizer = nlopt.opt(nlopt.LN_COBYLA, result.coords.size)
                    optimizer.set_min_objective(lambda x, _: ff.run(x.reshape((-1, 3)))['total'])
                    if self.restrained_dihedrals is not None:
                        for dihedral in self.restrained_dihedrals:
                            indices = dihedral.copy()
                            ref_angle = np.deg2rad(dihedralAngle(self.molecule.coords[indices, :, iframe]))
                            def constraint(x, _):
                                coords = x.reshape((-1, 3))
                                angle = dihedralAngle(coords[indices])
                                return np.sin(.5*(angle - ref_angle))
                            optimizer.add_equality_constraint(constraint)

                    optimizer.set_xtol_rel(1e-6)
                    optimizer.set_initial_step(1e-3)
                    result.coords = optimizer.optimize(result.coords.ravel()).reshape((-1, 3, 1))

                result.energy = ff.run(result.coords[:, :, 0])['total']
                result.dipole = self.molecule.getDipole()

                if self.optimize:
                    assert optimizer.last_optimum_value() == result.energy # A self-consistency test

                # Compute ESP values
                if self.esp_points is not None:
                    assert self.molecule.numFrames == 1
                    result.esp_points = self.esp_points
                    distances = cdist(result.esp_points, result.coords[:, :, 0]) # Angstrom
                    distances *= 0.52917721067  # Convert Angstrom to a.u. (Bohr) (https://en.wikipedia.org/wiki/Bohr_radius)
                    result.esp_values = np.dot(np.reciprocal(distances), self.molecule.charge) # Hartree/Bohr

                with open(pickleFile, 'wb') as fd:
                    pickle.dump(result, fd)

            results.append(result)

        return results


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
