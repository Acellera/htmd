# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import time
import pickle
import logging

import numpy as np
import nlopt

from htmd.numbautil import dihedralAngle
from htmd.qm.base import QMBase, QMResult

from qmml import QMMLCalculator

logger = logging.getLogger(__name__)


class QMML(QMBase):
    """
    Imitation of QM calculations with MM

    >>> import os
    >>> import numpy as np
    >>> from tempfile import TemporaryDirectory
    >>> from htmd.home import home
    >>> from htmd.numbautil import dihedralAngle
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> from htmd.qm.custom import QMML

    Create a molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2)

    Run a single-point energy and ESP calculation
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = QMML()
    ...     qm.molecule = mol
    ...     qm.esp_points = np.array([[1., 1., 1.]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.custom.QMML object at ...>
    >>> result # doctest: +ELLIPSIS
    <htmd.qm.base.QMResult object at ...
    >>> result.errored
    False
    >>> result.energy # doctest: +ELLIPSIS
    -94970.499...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    89.99...

    Run a minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = QMML()
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> result.energy # doctest: +ELLIPSIS
    -95173.433...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    125.993...

    Run a constrained minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = QMML()
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.restrained_dihedrals = np.array([[2, 0, 1, 3]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> result.energy # doctest: +ELLIPSIS
    -95170.800...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    89.99...
    """

    # Fake implementations of the abstract methods
    def _command(self): pass
    def _writeInput(self, directory, iframe): pass
    def _readOutput(self, directory): pass
    def setup(self): pass
    def submit(self): pass

    def _completed(self, directory):
        return os.path.exists(os.path.join(directory, 'data.pkl'))

    def retrieve(self):

        results = []
        for iframe in range(self.molecule.numFrames):
            self.molecule.frame = iframe

            directory = os.path.join(self.directory, '%05d' % iframe)
            os.makedirs(directory, exist_ok=True)
            pickleFile = os.path.join(directory, 'data.pkl')

            if self._completed(directory):
                with open(pickleFile, 'rb') as fd:
                    result = pickle.load(fd)
                logger.info('Loading data from %s' % pickleFile)

            else:
                start = time.clock()

                result = QMResult()
                result.errored = False
                result.coords = self.molecule.coords[:, :, iframe:iframe + 1].copy()

                calc = QMMLCalculator(maxatoms=self.molecule.numAtoms)

                if self.optimize:
                    opt = nlopt.opt(nlopt.LN_COBYLA, result.coords.size)
                    def objective(x, _):
                        return float(calc.calculate(x.reshape((-1, 3, 1)), self.molecule.element)[0])
                    opt.set_min_objective(objective)
                    if self.restrained_dihedrals is not None:
                        for dihedral in self.restrained_dihedrals:
                            indices = dihedral.copy()
                            ref_angle = dihedralAngle(self.molecule.coords[indices, :, iframe])
                            def constraint(x, _):
                                coords = x.reshape((-1, 3))
                                angle = dihedralAngle(coords[indices])
                                return np.sin(.5*(angle - ref_angle))
                            opt.add_equality_constraint(constraint)
                    opt.set_xtol_abs(1e-3) # Similar to Psi4 default
                    opt.set_maxeval(1000*opt.get_dimension())
                    opt.set_initial_step(1e-3)
                    result.coords = opt.optimize(result.coords.ravel()).reshape((-1, 3, 1))
                    logger.info('Optimization status: %d' % opt.last_optimize_result())

                result.energy = float(calc.calculate(result.coords, self.molecule.element)[0])
                result.dipole = self.molecule.getDipole()

                if self.optimize:
                    assert opt.last_optimum_value() == result.energy # A self-consistency test

                finish = time.clock()
                result.calculator_time = finish - start
                logger.info('QMML calculation time: %f s' % result.calculator_time)

                with open(pickleFile, 'wb') as fd:
                    pickle.dump(result, fd)

            results.append(result)

        return results


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
