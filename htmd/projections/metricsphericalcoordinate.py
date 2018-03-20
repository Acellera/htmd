# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricSphericalCoordinate(Projection):
    """ Creates a MetricSphericalCoordinate object that calculates the spherical coordinates between two centers of
    masses from a set of trajectories.

    Parameters
    ----------
    refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The reference Molecule to which we will align.
    targetcom : str
        Atom selection string from which to calculate the target center of mass.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    refcom : str
        Atom selection string from which to calculate the reference center of mass.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    trajalnstr : str, optional
        Atom selection string for the trajectories from which to align to the reference structure.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    refalnstr : str, optional
        Atom selection string for `refmol` from which to align to the reference structure. If None, it defaults to the
        same as `trajalnstr`. See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    centerstr : str, optional
        Atom selection string around which to wrap the simulation.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    pbc : bool
        Enable or disable coordinate wrapping based on periodic boundary conditions.

    Returns
    -------
    metr : MetricCoordinate object
    """
    def __init__(self, refmol, targetcom, refcom, trajalnstr='protein and name CA', refalnstr=None, centerstr='protein', pbc=True):
        if refalnstr is None:
            refalnstr = trajalnstr
        self._refmol = refmol
        self._refalnsel = self._refmol.atomselect(refalnstr)
        self._trajalnsel = trajalnstr
        self._centersel = centerstr
        self._targetcom = targetcom
        self._refcom = refcom
        self._pc_trajalnsel = None  # pc = Pre-calculated
        self._pc_atomsel = None
        self._pc_centersel = None
        self.pbc = pbc

    def _precalculate(self, mol):
        self._pc_trajalnsel = mol.atomselect(self._trajalnsel)
        self._pc_targetcom = mol.atomselect(self._targetcom)
        self._pc_refcom = mol.atomselect(self._refcom)
        self._pc_centersel = mol.atomselect(self._centersel)

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
        trajalnsel, targetcom, refcom, centersel = self._getSelections(mol)

        mol = mol.copy()
        if self.pbc:
            mol.wrap(self._centersel)
        mol.align(trajalnsel, refmol=self._refmol, refsel=self._refalnsel)

        refcom = np.mean(mol.coords[refcom, :, :], axis=0)
        targetcom = np.mean(mol.coords[targetcom, :, :], axis=0)
        xyz = targetcom - refcom

        r = np.sqrt(np.sum(xyz ** 2, axis=0))
        theta = np.arccos(xyz[2, :] / r)
        phi = np.arctan2(xyz[1, :], xyz[0, :])

        return np.stack((r, theta, phi), axis=1)

    def _getSelections(self, mol):
        if self._pc_trajalnsel is not None and self._pc_refcom is not None and self._pc_targetcom is not None \
                and self._pc_centersel is not None:
            trajalnsel = self._pc_trajalnsel
            targetcom = self._pc_targetcom
            refcom = self._pc_refcom
            centersel = self._pc_centersel
        else:
            trajalnsel = mol.atomselect(self._trajalnsel)
            targetcom = mol.atomselect(self._targetcom)
            refcom = mol.atomselect(self._refcom)
            centersel = mol.atomselect(self._centersel)
        if np.sum(trajalnsel) == 0:
            raise NameError('Alignment selection resulted in 0 atoms.')
        if np.sum(targetcom) == 0:
            raise NameError('Atom selection for `targetcom` resulted in 0 atoms.')
        if np.sum(refcom) == 0:
            raise NameError('Atom selection for `refcom` resulted in 0 atoms.')
        if np.sum(centersel) == 0:
            raise NameError('Centering selection resulted in 0 atoms.')
        return trajalnsel, targetcom, refcom, centersel

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
        trajalnsel, targetcom, refcom, centersel = self._getSelections(mol)
        targetatomidx = np.where(targetcom)[0]
        refatomidx = np.where(refcom)[0]
        from pandas import DataFrame
        types = ['r', 'theta', 'phi']
        indexes = [[targetatomidx, refatomidx]] * 3
        description = ['r', 'theta', 'phi']
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    ref = Molecule(path.join(home(dataDir='metricdistance'), 'filtered.pdb'))
    mol = ref.copy()
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))

    res = MetricSphericalCoordinate(ref, 'resname MOL', 'within 8 of resid 98').project(mol)
    _ = MetricSphericalCoordinate(ref, 'resname MOL', 'within 8 of resid 98').getMapping(mol)

    ref_array = np.load(path.join(home(dataDir='test-metrics'), 'metricsphericalcoordinate', 'res.npy'))
    assert np.allclose(res, ref_array, rtol=0, atol=1e-04)
