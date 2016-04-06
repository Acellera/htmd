# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.metric import _OldMetric, _singleMolfile
from htmd.molecule.molecule import Molecule
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

    Returns
    -------
    proj : MetricCoordinate object
    """
    def __init__(self, refmol, atomsel, trajalnstr='protein and name CA', refalnstr=None, centerstr='protein'):
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

    def _precalculate(self, mol):
        self._pc_trajalnsel = mol.atomselect(self._trajalnsel)
        self._pc_atomsel = mol.atomselect(self._atomsel)
        self._pc_centersel = mol.atomselect(self._centersel)

    def project(self, *args, **kwargs):
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
        # -------------- DEPRECATION PATCH --------------
        if isinstance(self, np.ndarray) or isinstance(self, Molecule):
            from warnings import warn
            warn('Static use of the project method will be deprecated in the next version of HTMD. '
                 'Please change your projection code according to the tutorials on www.htmd.org')
            data = _MetricCoordinateOld.project(self, *args, **kwargs)
            logger.warning('Static use of the project method will be deprecated in the next version of HTMD. '
                           'Please change your projection code according to the tutorials on www.htmd.org')
            return data
        # ------------------ CUT HERE -------------------
        mol = args[0]
        (trajalnsel, atomsel, xxx) = self._getSelections(mol)

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
        (xxx, atomsel, yyy) = self._getSelections(mol)
        premap = np.where(atomsel)[0]
        map = np.zeros(len(premap) * 3, dtype=int)
        map[0*len(premap):1*len(premap)] = premap
        map[1*len(premap):2*len(premap)] = premap
        map[2*len(premap):3*len(premap)] = premap


class _MetricCoordinateOld(_OldMetric):
    @staticmethod
    def project(sims, refmol, atomsel, trajalnstr='protein and name CA', refalnstr=None, centerstr='protein', skip=1, update=[]):
        """ Calculates the atom coordinates from a set of trajectories. TODO: Work in progress!

        Parameters
        ----------
        sims : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects or single :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function, or a :class:`Molecule <htmd.molecule.molecule.Molecule>` object
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
        skip : int
            Skip every `skip` frames of the trajectories
        update :
            Not functional yet

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated. First all
            X coordinates of all atoms, then all Y and then all Z coordinates.
        """
        if refalnstr is None:
            refalnstr = trajalnstr
        obj = _MetricCoordinateOld(sims, refmol, trajalnstr, refalnstr, atomsel, centerstr)
        return obj._metrify(sims, skip, update)

    def __init__(self, sims, refmol, trajalnstr, refalnstr, atomsel, centerstr):
        self._refmol = refmol
        self._refalnsel = self._refmol.atomselect(refalnstr)
        self._trajalnsel = trajalnstr
        self._centersel = centerstr
        self._atomsel = atomsel
        self._pc_trajalnsel = None  # pc = Pre-calculated
        self._pc_atomsel = None
        self._pc_centersel = None

        (single, molfile) = _singleMolfile(sims)
        if single:
            mol = Molecule(molfile)
            self._pc_trajalnsel = mol.atomselect(trajalnstr)
            self._pc_atomsel = mol.atomselect(atomsel)
            self._pc_centersel = mol.atomselect(centerstr)

    def _processTraj(self, mol):
        (trajalnsel, atomsel, xxx) = self._getSelections(mol)

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

    def _getMapping(self, mol):
        (xxx, atomsel, yyy) = self._getSelections(mol)
        premap = np.where(atomsel)[0]
        map = np.zeros(len(premap) * 3, dtype=int)
        map[0*len(premap):1*len(premap)] = premap
        map[1*len(premap):2*len(premap)] = premap
        map[2*len(premap):3*len(premap)] = premap


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