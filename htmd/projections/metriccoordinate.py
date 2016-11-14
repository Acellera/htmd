# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricCoordinate(Projection):
    """ Creates a MetricCoordinate object that calculates the atom coordinates from a set of trajectories.

    Parameters
    ----------
    refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The reference Molecule to which we will align.
    atomsel : str
        Atomselection for the atoms whose coordinates we want to retrieve.
    trajalnstr : str, optional
        Atomselection for the trajectories from which to align to the reference structure.
    refalnstr : str, optional
        Atomselection for `refmol` from which to align to the reference structure. If None, it defaults to the same as `trajalnstr`.
    centerstr : str, optional
        Atomselection around which to wrap the simulation.
    pbc : bool
        Enable or disable coordinate wrapping based on periodic boundary conditions.

    Returns
    -------
    metr : MetricCoordinate object
    """
    def __init__(self, refmol, atomsel, trajalnstr='protein and name CA', refalnstr=None, centerstr='protein', pbc=True):
        if refalnstr is None:
            refalnstr = trajalnstr
        self._refmol = refmol
        self._refalnsel = self._refmol.atomselect(refalnstr)
        self._trajalnsel = trajalnstr
        self._centersel = centerstr
        self._atomsel = atomsel
        self._pc_trajalnsel = None  # pc = Pre-calculated
        self._pc_atomsel = None
        self._pc_centersel = None
        self.pbc = pbc

    def _precalculate(self, mol):
        self._pc_trajalnsel = mol.atomselect(self._trajalnsel)
        self._pc_atomsel = mol.atomselect(self._atomsel)
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
        (trajalnsel, atomsel, xxx) = self._getSelections(mol)

        if self.pbc:
            mol.wrap(self._centersel)
        mol.align(trajalnsel, refmol=self._refmol, refsel=self._refalnsel)

        coords = np.concatenate((mol.coords[atomsel, 0, :], mol.coords[atomsel, 1, :], mol.coords[atomsel, 2, :]))
        return np.transpose(coords)  # I need to permute the dimensions here to put frames on the first

    def _getSelections(self, mol):
        if self._pc_trajalnsel is not None and self._pc_atomsel is not None and self._pc_centersel is not None:
            trajalnsel = self._pc_trajalnsel
            atomsel = self._pc_atomsel
            centersel = self._pc_centersel
        else:
            trajalnsel = mol.atomselect(self._trajalnsel)
            atomsel = mol.atomselect(self._atomsel)
            centersel = mol.atomselect(self._centersel)
        if np.sum(trajalnsel) == 0:
            raise NameError('Alignment selection resulted in 0 atoms.')
        if np.sum(atomsel) == 0:
            raise NameError('Atom selection resulted in 0 atoms.')
        if np.sum(centersel) == 0:
            raise NameError('Centering selection resulted in 0 atoms.')
        return trajalnsel, atomsel, centersel

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
        (xxx, atomsel, yyy) = self._getSelections(mol)
        atomidx = np.where(atomsel)[0]
        from pandas import DataFrame
        types = []
        indexes = []
        description = []
        for xyz in ('X', 'Y', 'Z'):
            for i in atomidx:
                types += ['coordinate']
                indexes += [i]
                description += ['{} coordinate of {} {} {}'.format(xyz, mol.resname[i], mol.resid[i], mol.name[i])]
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    ref = mol.copy()
    ref.coords = np.atleast_3d(ref.coords[:, :, 0])
    metr = MetricCoordinate(ref, 'protein and name CA')
    data = metr.project(mol)

    lastcoors = np.array([6.79283285,   5.55226946,   4.49387407,   2.94484425,
                          5.36937141,   3.18590879,   5.75874281,   5.48864174,
                          1.69625032,   1.58790839,   0.57877392,  -2.66498065,
                          -3.70919156,  -3.33702421,  -5.38465405,  -8.43286991,
                          -8.15859032,  -7.85062265, -10.92551327, -13.70733166], dtype=np.float32)
    assert np.all(np.abs(data[-1, -20:] - lastcoors) < 0.001), 'Coordinates calculation is broken'