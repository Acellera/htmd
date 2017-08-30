# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial.distance import cdist

from htmd.qm.base import QMBase, QMResult


class FakeQM(QMBase):
    """
    Imitation of QM calculations with MM

    >>> import os
    >>> import numpy as np
    >>> from htmd.home import home
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> from htmd.qm.fake import FakeQM

    Create a molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2) # doctest: +ELLIPSIS
    Dihedral 0: 2-0-1-3
    ...

    Create a face QM object
    >>> qm = FakeQM()
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.fake.FakeQM object at ...
    >>> qm.molecule = mol
    >>> qm.esp_points = np.array([[1., 1., 1.]])
    >>> result = qm.run()[0]

    Inspect results
    >>> result # doctest: +ELLIPSIS
    <htmd.qm.base.QMResult object at ...
    >>> result.errored
    False
    >>> result.energy
    0
    >>> result.esp_points
    array([[ 1.,  1.,  1.]])
    >>> result.esp_values
    array([ 0.01507904])
    """

    # Fake implementations of the abstract methods
    def _command(self): pass
    def _completed(self): pass
    def _readOutput(self): pass
    def _writeInput(self): pass

    def run(self):

        assert self.molecule.numFrames == 1
        assert self.optimize == False # TODO

        result = QMResult()
        result.errored = False
        result.coords = self.molecule.coords
        result.energy = 0 # TODO
        result.dipole = self.molecule.getDipole()
        result.esp_points = self.esp_points

        # Compute ESP values
        if self.esp_points is not None:
            reciprocal_distances = 1/cdist(result.esp_points, result.coords[:, :, 0]) # 1/Angstrom
            result.esp_values = np.dot(reciprocal_distances, self.molecule.charge) # Hartree/Angstrom

        return [result]


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
