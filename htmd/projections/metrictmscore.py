# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
from htmd.molecule.util import molTMscore
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricTMscore(Projection):
    """ Calculates the TMscore of a set of trajectories to a reference structure

    Parameters
    ----------
    refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The reference Molecule to which we want to calculate the TMscore.
    trajtmstr : str
        Atomselection for the trajectories from which to calculate the TMscore
    reftmstr : str, optional
        Atomselection for the reference structure from which to calculate the TMscore. If None, it defaults to `trajrmsdstr`
    centerstr : str, optional
        Atomselection around which to center the wrapping of the trajectories.
    """
    def __init__(self, refmol, trajtmstr, reftmstr=None):
        if reftmstr is None:
            reftmstr = trajtmstr
        self._refmol = refmol
        self._reftmsel = self._refmol.atomselect(reftmstr) & (self._refmol.name == 'CA')
        self._trajtmsel = trajtmstr
        self._pc_trajtmsel = None

    def _precalculate(self, mol):
        self._pc_trajtmsel = mol.atomselect(self._trajtmsel) & (mol.name == 'CA')

    def project(self, mol):
        """ Project molecule.

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>`
            A :class:`Molecule <htmd.molecule.molecule.Molecule>` object to project.

        Returns
        -------
        data : np.ndarray
            An array containing the projected data.
        """
        mol = mol.copy()
        trajtmsel = self._getSelections(mol)

        tm, _ = molTMscore(mol, self._refmol, trajtmsel, self._reftmsel)
        return tm[:, np.newaxis]

    def _getSelections(self, mol):
        if self._pc_trajtmsel is not None:
            trajtmsel = self._pc_trajtmsel
        else:
            trajtmsel = mol.atomselect(self._trajtmsel) & (mol.name == 'CA')
        if np.sum(trajtmsel) == 0:
            raise NameError('RMSD atom selection resulted in 0 atoms.')
        return trajtmsel

    def getMapping(self, mol):
        """ Returns the description of each projected dimension.

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object which will be used to calculate the descriptions of the projected dimensions.

        Returns
        -------
        map : :class:`DataFrame <pandas.core.frame.DataFrame>` object
            A DataFrame containing the descriptions of each dimension
        """
        trajtmsel = self._getSelections(mol)
        from pandas import DataFrame
        types = ['tmscore']
        indexes = [np.where(trajtmsel)[0]]
        description = ['TMscore to reference structure.']
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import numpy as np
    from os import path

    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    ref = mol.copy()
    ref.coords = np.atleast_3d(ref.coords[:, :, 0])
    metr = MetricTMscore(ref, 'protein and name CA')
    data = metr.project(mol)

    lasttm = np.array([0.9633381, 0.96441294, 0.96553609, 0.96088852, 0.96288511, 0.95677591, 0.96544727, 0.96359811,
                       0.95658912, 0.96893117, 0.96623924, 0.96064913, 0.96207041, 0.95947848, 0.96657048, 0.95993426,
                       0.96543296, 0.96806875, 0.96437248, 0.96144066], dtype=np.float32)
    assert np.all(np.abs(data[-20:].flatten() - lasttm) < 0.001), 'Coordinates calculation is broken'


    from htmd.simlist import simlist
    from htmd.projections.metric import Metric
    dd = home(dataDir="adaptive")
    fsims = simlist([path.join(dd, 'data', 'e1s1_1'), path.join(dd, 'data', 'e1s2_1')],
                    path.join(dd, 'generators', '1', 'structure.pdb'))
    ref = Molecule(path.join(dd, 'generators', '1', 'structure.pdb'))

    metr2 = Metric(fsims)
    metr2.set(MetricTMscore(ref, 'protein and name CA'))
    data2 = metr2.project()

    assert data2.trajectories[0].projection.shape == (6, 1)
