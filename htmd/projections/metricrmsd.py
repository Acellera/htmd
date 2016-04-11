# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.metric import _OldMetric, _singleMolfile
from htmd.projections.projection import Projection
from htmd.molecule.molecule import Molecule
from htmd.molecule.util import molRMSD
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricRmsd(Projection):
    """ Calculates the RMSD of a set of trajectories to a reference structure

    Parameters
    ----------
    refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        The reference Molecule to which we want to calculate the RMSD.
    trajrmsdstr : str
        Atomselection for the trajectories from which to calculate the RMSD
    trajalnstr : str, optional
        Atomselection for the trajectories from which to align to the reference structure. If None, it defaults to the same as `trajrmsdstr`.
    refrmsdstr : str, optional
        Atomselection for the reference structure from which to calculate the RMSD. If None, it defaults to `trajrmsdstr`
    refalnstr : str, optional
        Atomselection for the reference structure from which to align to the trajectories. If None, it defaults to `trajalnstr`
    centerstr : str, optional
        Atomselection around which to center the wrapping of the trajectories. CURRENTLY IGNORED!
    """
    def __init__(self, refmol, trajrmsdstr, trajalnstr=None, refrmsdstr=None, refalnstr=None, centerstr='protein', pbc=True):
        if trajalnstr is None:
            trajalnstr = trajrmsdstr
        if refalnstr is None:
            refalnstr = trajalnstr
        if refrmsdstr is None:
            refrmsdstr = trajrmsdstr
        self._refmol = refmol
        self._refalnsel = self._refmol.atomselect(refalnstr)
        self._refrmsdsel = self._refmol.atomselect(refrmsdstr)
        self._trajalnsel = trajalnstr
        self._trajrmsdsel = trajrmsdstr
        self._centersel = centerstr
        self._pc_trajalnsel = None  # pc = Pre-calculated
        self._pc_trajrmsdsel = None
        self._pc_centersel = None
        self._pbc = pbc

    def _precalculate(self, mol):
        self._pc_trajalnsel = mol.atomselect(self._trajalnsel)
        self._pc_trajrmsdsel = mol.atomselect(self._trajrmsdsel)
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
            data = _MetricRmsdOld.project(self, *args, **kwargs)
            logger.warning('Static use of the project method will be deprecated in the next version of HTMD. '
                           'Please change your projection code according to the tutorials on www.htmd.org')
            return data
        # ------------------ CUT HERE -------------------
        mol = args[0]
        mol = mol.copy()
        (trajalnsel, trajrmsdsel, centersel) = self._getSelections(mol)

        if self._pbc:
            mol.wrap()
        #mol.coords = self._wrapPositions(mol.box, mol.coords, centersel)
        mol.align(sel=trajalnsel, refmol=self._refmol, refsel=self._refalnsel)

        return molRMSD(mol, self._refmol, trajrmsdsel, self._refrmsdsel)

    def _getSelections(self, mol):
        if self._pc_trajalnsel is not None and self._pc_trajrmsdsel is not None and self._pc_centersel is not None:
            trajalnsel = self._pc_trajalnsel
            trajrmsdsel = self._pc_trajrmsdsel
            centersel = self._pc_centersel
        else:
            trajalnsel = mol.atomselect(self._trajalnsel)
            trajrmsdsel = mol.atomselect(self._trajrmsdsel)
            centersel = mol.atomselect(self._centersel)
        if np.sum(trajalnsel) == 0:
            raise NameError('Alignment selection resulted in 0 atoms.')
        if np.sum(trajrmsdsel) == 0:
            raise NameError('RMSD atom selection resulted in 0 atoms.')
        if np.sum(centersel) == 0:
            raise NameError('Centering selection resulted in 0 atoms.')
        return trajalnsel, trajrmsdsel, centersel

    def getMapping(self, mol):
        (xxx, trajrmsdsel, yyy) = self._getSelections(mol)
        return np.where(trajrmsdsel)[0]

    def _wrapPositions(self, box, pos, centersel):
        if box is None or np.sum(box) == 0:
            logger.warning('MetricRmsd: The given molecule does not contain box dimensions for wrapping.')
            return pos
        center = np.mean(pos[centersel, :, :], axis=0)
        origin = center - (box / 2)
        pos = pos - origin
        return np.mod(pos, box)


class _MetricRmsdOld(_OldMetric):
    @staticmethod
    def project(sims, refmol, trajrmsdstr, trajalnstr=None, refrmsdstr=None, refalnstr=None, centerstr='protein', skip=1, update=[], pbc=True):
        """ Calculates the RMSD of a set of trajectories to a reference structure

        Parameters
        ----------
        sims : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects or single :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function, or a :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            The reference Molecule to which we want to calculate the RMSD.
        trajrmsdstr : str
            Atomselection for the trajectories from which to calculate the RMSD
        trajalnstr : str, optional
            Atomselection for the trajectories from which to align to the reference structure. If None, it defaults to the same as `trajrmsdstr`.
        refrmsdstr : str, optional
            Atomselection for the reference structure from which to calculate the RMSD. If None, it defaults to `trajrmsdstr`
        refalnstr : str, optional
            Atomselection for the reference structure from which to align to the trajectories. If None, it defaults to `trajalnstr`
        centerstr : str, optional
            Atomselection around which to center the wrapping of the trajectories
        skip : int
            Skip every `skip` frames of the trajectories
        update :
            Not functional yet

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated
        """
        if trajalnstr is None:
            trajalnstr = trajrmsdstr
        if refalnstr is None:
            refalnstr = trajalnstr
        if refrmsdstr is None:
            refrmsdstr = trajrmsdstr
        obj = _MetricRmsdOld(sims, refmol, trajalnstr, refalnstr, trajrmsdstr, refrmsdstr, centerstr, pbc)
        return obj._metrify(sims, skip, update)

    def __init__(self, sims, refmol, trajalnstr, refalnstr, trajrmsdstr, refrmsdstr, centerstr, pbc):
        self._refmol = refmol
        self._refalnsel = self._refmol.atomselect(refalnstr)
        self._refrmsdsel = self._refmol.atomselect(refrmsdstr)
        self._trajalnsel = trajalnstr
        self._trajrmsdsel = trajrmsdstr
        self._centersel = centerstr
        self._pc_trajalnsel = None  # pc = Pre-calculated
        self._pc_trajrmsdsel = None
        self._pc_centersel = None
        self._pbc = pbc

        (single, molfile) = _singleMolfile(sims)
        if single:
            mol = Molecule(molfile)
            self._pc_trajalnsel = mol.atomselect(trajalnstr)
            self._pc_trajrmsdsel = mol.atomselect(trajrmsdstr)
            self._pc_centersel = mol.atomselect(centerstr)

    def _processTraj(self, mol):
        (trajalnsel, trajrmsdsel, centersel) = self._getSelections(mol)

        mol = mol.copy()
        if self._pbc:
            mol.wrap()
        #mol.coords = self._wrapPositions(mol.box, mol.coords, centersel)
        mol.align(sel=trajalnsel, refmol=self._refmol, refsel=self._refalnsel)

        return molRMSD(mol, self._refmol, trajrmsdsel, self._refrmsdsel)

    def _getSelections(self, mol):
        if self._pc_trajalnsel is not None and self._pc_trajrmsdsel is not None and self._pc_centersel is not None:
            trajalnsel = self._pc_trajalnsel
            trajrmsdsel = self._pc_trajrmsdsel
            centersel = self._pc_centersel
        else:
            trajalnsel = mol.atomselect(self._trajalnsel)
            trajrmsdsel = mol.atomselect(self._trajrmsdsel)
            centersel = mol.atomselect(self._centersel)
        if np.sum(trajalnsel) == 0:
            raise NameError('Alignment selection resulted in 0 atoms.')
        if np.sum(trajrmsdsel) == 0:
            raise NameError('RMSD atom selection resulted in 0 atoms.')
        if np.sum(centersel) == 0:
            raise NameError('Centering selection resulted in 0 atoms.')
        return trajalnsel, trajrmsdsel, centersel

    def _getMapping(self, mol):
        (xxx, trajrmsdsel, yyy) = self._getSelections(mol)
        return np.where(trajrmsdsel)[0]

    def _wrapPositions(self, box, pos, centersel):
        if box is None or np.sum(box) == 0:
            logger.warning('MetricRmsd: The given molecule does not contain box dimensions for wrapping.')
            return pos
        center = np.mean(pos[centersel, :, :], axis=0)
        origin = center - (box / 2)
        pos = pos - origin
        return np.mod(pos, box)


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    ref = mol.copy()
    ref.coords = np.atleast_3d(ref.coords[:, :, 0])
    metr = MetricRmsd(ref, 'protein and name CA')
    data = metr.project(mol)

    lastrmsd = np.array([1.30797791,  1.29860222,  1.25042927,  1.31319737,  1.27044261,
                          1.40294552,  1.25354612,  1.30127883,  1.40618336,  1.18303752,
                          1.24414587,  1.34513164,  1.31932807,  1.34282494,  1.2261436 ,
                          1.36359048,  1.26243281,  1.21157813,  1.26476419,  1.29413617], dtype=np.float32)
    assert np.all(np.abs(data[-20:] - lastrmsd) < 0.001), 'Coordinates calculation is broken'