from htmd.projections.projection import Projection

import logging
import numpy as np


logger = logging.getLogger(__name__)


class MetricNull(Projection):
    """ A dummy metric returning zero for all frames. May be useful as a template for developers,
    or to get the simulation list without doing any processing.

    Parameters
    ----------
    ndim : int
        number of dimensions of the fake projection space

    """

    def __init__(self, ndim):
        logger.info("Initialized")
        self._ndim = ndim


    # Only called if single topology
    def _precalculate(self, mol):
        logger.info("In _precalculate")
        self._precalculation_enabled = True


    def getMapping(self):
        """ Returns optional information.

        Returns
        -------
        map :
            a dummy value
        """

        return "No mapping."


    def project(self, mol):
        """ Compute the actual metric.

        Parameters
        ------------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>`
            A :class:`Molecule <htmd.molecule.molecule.Molecule>` object to project.

        Returns
        -------
        data : np.ndarray
            An array containing the null data.
        """

        trajlen = mol.numFrames
        data = np.zeros((trajlen, self._ndim), dtype=np.float32)

        return data


if __name__ == "__main__":
    # TODO a test
    pass

