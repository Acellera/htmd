# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.metric import _OldMetric, _singleMolfile
from htmd.projections.projection import Projection
from htmd.molecule.molecule import Molecule
import numpy as np
from htmd.projections.util import pp_calcDistances, pp_calcMinDistances, pp_calcMinDistances_C
import logging
logger = logging.getLogger(__name__)


class MetricDistance(Projection):
    """ Creates a MetricDistance object

    Parameters
    ----------
    sel1 : str
        Atomselection for the first set of atoms
    sel2 : str
        Atomselection for the second set of atoms. If sel1 != sel2 it will calculate inter-set distances. If sel1 == sel2 it will calculate intra-set distances
    groupsel1 : ['all','residue'], optional
        Group all atoms in `sel1` to the single minimum distance. Alternatively can calculate the minimum distance of a residue containing the atoms in sel1.
    groupsel2 : ['all','residue'], optional
        Same as groupsel1 but for `sel2`
    metric : ['distances','contacts'], optional
        Set to 'concacts' to calculate contacts instead of distances
    threshold : float, optional
        The threshold under which a distance is considered in contact
    pbc : bool, optional
        Set to false to disable distance calculations using periodic distances
    truncate : float, optional
        Set all distances larger than `truncate` to `truncate`
    update :
        Not functional yet

    Returns
    -------
    proj : MetricDistance object
    """
    def __init__(self, sel1, sel2, groupsel1=None, groupsel2=None,  metric='distances', threshold=8, pbc=True, truncate=None):
        self.sel1 = sel1
        self.sel2 = sel2
        self.groupsel1 = groupsel1
        self.groupsel2 = groupsel2
        self.metric = metric
        self.threshold = threshold
        self.pbc = pbc
        self.truncate = truncate
        self.precalcsel1 = None
        self.precalcsel2 = None

    def _precalculate(self, mol):
        self.precalcsel1 = self._processSelection(mol, self.sel1, self.groupsel1)
        self.precalcsel2 = self._processSelection(mol, self.sel2, self.groupsel2)

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
            data = _MetricDistanceOld.project(self, *args, **kwargs)
            logger.warning('Static use of the project method will be deprecated in the next version of HTMD. '
                            'Please change your projection code according to the tutorials on www.htmd.org')
            return data
        # ------------------ CUT HERE -------------------
        mol = args[0]
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 1 and np.ndim(sel2) == 1:  # normal distances
            metric = pp_calcDistances(mol, sel1, sel2, self.metric, self.threshold, self.pbc, truncate=self.truncate)
        else:  # minimum distances by groups
            metric = pp_calcMinDistances_C(mol, sel1, sel2)

        return metric

    def getMapping(self, mol):
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 2:
            protatoms = np.ones(len(sel1)) * -1
            for i in range(np.size(sel1, 1)):
                protatoms[i] = np.where(sel1[i] == True)[0][0]
        else:
            protatoms = np.where(sel1)[0]
        if np.ndim(sel2) == 2:
            ligatoms = np.ones(len(sel1)) * -1
            for i in range(np.size(sel2, 1)):
                ligatoms[i] = np.where(sel2[i] == True)[0][0]
        else:
            ligatoms = np.where(sel2)[0]

        numatoms1 = len(protatoms)
        numatoms2 = len(ligatoms)

        if np.array_equal(protatoms, ligatoms):
            map = np.zeros((numatoms1 * (numatoms1-1) / 2, 2), dtype=int)
            start = 0
            for i in range(numatoms1):
                finish = start + numatoms1 - i - 1
                map[start:finish, 0] = protatoms[i]
                map[start:finish, 1] = protatoms[i+1:]
                start = finish
        else:
            map = np.zeros((numatoms1 * numatoms2, 2), dtype=int)
            for i in range(numatoms2):
                start = i * numatoms1
                finish = (i+1) * numatoms1
                map[start:finish, 0] = protatoms
                map[start:finish, 1] = ligatoms[i]
        return map

    def _getSelections(self, mol):
        # If they have been pre-calculated return them.
        if self.precalcsel1 is not None and self.precalcsel2 is not None:
            return self.precalcsel1, self.precalcsel2
        else:  # Otherwise compute them.
            sel1 = self._processSelection(mol, self.sel1, self.groupsel1)
            sel2 = self._processSelection(mol, self.sel2, self.groupsel2)
            return sel1, sel2

    def _processSelection(self, mol, sel, groupsel):
        if isinstance(sel, str):  # If user passes simple string selections
            if groupsel is None:
                sel = mol.atomselect(sel)
            elif groupsel == 'all':
                sel = self._processMultiSelections(mol, [sel])
            elif groupsel == 'residue':
                selCell = self._groupByResidue(mol, sel)
                sel = self._processMultiSelections(mol, selCell)
            else:
                raise NameError('Invalid groupsel argument')
        else:  # If user passes his own sets of groups
            sel = self._processMultiSelections(mol, sel)

        if np.sum(sel) == 0:
            raise NameError('Selection returned 0 atoms')
        return sel

    def _processMultiSelections(self, mol, sel):
        newsel = np.zeros((len(sel), mol.numAtoms), dtype=object)
        for s in range(len(sel)):
            newsel[s, :] = mol.atomselect(sel[s])
        return newsel

    def _groupByResidue(self, mol, sel):
        res = list(set(mol.get('resid', sel=sel)))
        newsel = np.empty(len(res), dtype=object)
        for i in range(len(res)):
            newsel[i] = 'resid ' + str(res[i]) + ' and noh'
        return newsel


class MetricSelfDistance(MetricDistance):
    def __init__(self, sel, groupsel=None,  metric='distances', threshold=8, pbc=True, truncate=None):
        """ Creates a MetricSelfDistance object

        Parameters
        ----------
        sel : str
            Atomselection for the first set of atoms
        groupsel : ['all','residue'], optional
            Group all atoms in `sel1` to the single minimum distance. Alternatively can calculate the minimum distance of a residue containing the atoms in sel1.
        metric : ['distances','contacts'], optional
            Set to 'concacts' to calculate contacts instead of distances
        threshold : float, optional
            The threshold under which a distance is considered in contact
        pbc : bool, optional
            Set to false to disable distance calculations using periodic distances
        truncate : float, optional
            Set all distances larger than `truncate` to `truncate`
        update :
            Not functional yet

        Returns
        -------
        proj : MetricDistance object
        """
        super().__init__(sel, sel, groupsel, groupsel, metric, threshold, pbc, truncate)


class _MetricDistanceOld(_OldMetric):
    @staticmethod
    def project(sims, sel1, sel2=None, groupsel1=None, groupsel2=None, metric='distances', threshold=8,
                pbc=True, gap=None, truncate=None, skip=1, update=[]):
        """ Calculates a set of distances from trajectories

        Parameters
        ----------
        sims : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects or single :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function, or a :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        sel1 : str
            Atomselection for the first set of atoms
        sel2 : str
            Atomselection for the second set of atoms. If sel1 != sel2 it will calculate inter-set distances. If sel1 == sel2 it will calculate intra-set distances
        groupsel1 : ['all','residue'], optional
            Group all atoms in `sel1` to the single minimum distance. Alternatively can calculate the minimum distance of a residue containing the atoms in sel1.
        groupsel2 : ['all','residue'], optional
            Same as groupsel1 but for `sel2`
        metric : ['distances','contacts'], optional
            Set to 'concacts' to calculate contacts instead of distances
        threshold : float, optional
            The threshold under which a distance is considered in contact
        pbc : bool, optional
            Set to false to disable distance calculations using periodic distances
        gap : int, optional
            Not functional yet
        truncate : float, optional
            Set all distances larger than `truncate` to `truncate`
        skip : int
            Skip every `skip` frames of the trajectories
        update :
            Not functional yet

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated
        """
        if sel2 is None:
            sel2 = sel1
        obj = _MetricDistanceOld(sims, sel1, sel2, groupsel1, groupsel2, metric, threshold, pbc, truncate)
        return obj._metrify(sims, skip, update)

    def __init__(self, sims, sel1, sel2, groupsel1, groupsel2,  metric, threshold, pbc, truncate):
        """ Do not call this constructor directly. Use the static `project` method
        """
        self.sel1 = sel1
        self.sel2 = sel2
        self.groupsel1 = groupsel1
        self.groupsel2 = groupsel2
        self.precalcsel1 = None
        self.precalcsel2 = None
        self.metric = metric
        self.threshold = threshold
        self.pbc = pbc
        self.truncate = truncate

        (single, molfile) = _singleMolfile(sims)
        if single:
            mol = Molecule(molfile)
            self.precalcsel1 = self._processSelection(mol, sel1, groupsel1)
            self.precalcsel2 = self._processSelection(mol, sel2, groupsel2)

    def _processTraj(self, mol):
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 1 and np.ndim(sel2) == 1:  # normal distances
            metric = pp_calcDistances(mol, sel1, sel2, self.metric, self.threshold, self.pbc, truncate=self.truncate)
        else:  # minimum distances by groups
            metric = pp_calcMinDistances(mol, sel1, sel2)

        #if not issparse(metric) and not metric.dtype == 'bool' and np.any(np.isreal(metric)) and self.precision != 0:
        #    metric = np.round(metric * (1 / self.precision)) * self.precision
        return metric

    def _getMapping(self, mol):
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 2:
            protatoms = np.ones(len(sel1)) * -1
            for i in range(np.size(sel1, 1)):
                protatoms[i] = np.where(sel1[i] == True)[0][0]
        else:
            protatoms = np.where(sel1)[0]
        if np.ndim(sel2) == 2:
            ligatoms = np.ones(len(sel1)) * -1
            for i in range(np.size(sel2, 1)):
                ligatoms[i] = np.where(sel2[i] == True)[0][0]
        else:
            ligatoms = np.where(sel2)[0]

        numatoms1 = len(protatoms)
        numatoms2 = len(ligatoms)

        if np.array_equal(protatoms, ligatoms):
            map = np.zeros((numatoms1 * (numatoms1-1) / 2, 2), dtype=int)
            start = 0
            for i in range(numatoms1):
                finish = start + numatoms1 - i - 1
                map[start:finish, 0] = protatoms[i]
                map[start:finish, 1] = protatoms[i+1:]
                start = finish
        else:
            map = np.zeros((numatoms1 * numatoms2, 2), dtype=int)
            for i in range(numatoms2):
                start = i * numatoms1
                finish = (i+1) * numatoms1
                map[start:finish, 0] = protatoms
                map[start:finish, 1] = ligatoms[i]
        return map

    def _getSelections(self, mol):
        # If they have been pre-calculated return them.
        if self.precalcsel1 is not None and self.precalcsel2 is not None:
            return self.precalcsel1, self.precalcsel2
        else:  # Otherwise compute them.
            sel1 = self._processSelection(mol, self.sel1, self.groupsel1)
            sel2 = self._processSelection(mol, self.sel2, self.groupsel2)
            return sel1, sel2

    def _processSelection(self, mol, sel, groupsel):
        if isinstance(sel, str):  # If user passes simple string selections
            if groupsel is None:
                sel = mol.atomselect(sel)
            elif groupsel == 'all':
                sel = self._processMultiSelections(mol, sel)
            elif groupsel == 'residue':
                selCell = self._groupByResidue(mol, sel)
                sel = self._processMultiSelections(mol, selCell)
            else:
                raise NameError('Invalid groupSel argument')
        else:  # If user passes his own sets of groups
            sel = self._processMultiSelections(mol, sel)

        if np.sum(sel) == 0:
            raise NameError('Selection returned 0 atoms')
        return sel

    def _processMultiSelections(self, mol, sel):
        newsel = np.zeros((len(sel), mol.numAtoms), dtype=object)
        for s in range(len(sel)):
            newsel[s, :] = mol.atomselect(sel[s])
        return newsel

    def _groupByResidue(self, mol, sel):
        res = list(set(mol.get('resid', sel=sel)))
        newsel = np.empty(len(res), dtype=object)
        for i in range(len(res)):
            newsel[i] = 'resid ' + str(res[i]) + ' and noh'
        return newsel

'''
class MetricDistancePyemma(MetricPyemma):
    @staticmethod
    def project(simList, sel1, sel2, precision=0.1, groupsel1='none', groupsel2='none', metric='distances', threshold=0.8,
                pbc=1, gap=False, truncate=False, skip=1, verbose=1, update=[]):
        obj = MetricDistancePyemma(simList, sel1, sel2, precision, groupsel1, groupsel2, metric, threshold)
        return obj._metrify(simList, skip, verbose, update)

    def __init__(self, sims, sel1, sel2, precision, groupsel1, groupsel2, metric, threshold):
        self.feat = coor.featurizer(sims[0].molfile)
        self.metric = metric
        self.threshold = threshold

        mol = Molecule(sims[0].molfile)
        sel1 = np.where(mol.atomselect(sel1))[0]
        sel2 = np.where(mol.atomselect(sel2))[0]

        #sel1 = self.feat.select(sel1)
        #sel2 = self.feat.select(sel2)
        if np.array_equal(sel1, sel2):
            logging.info('Identical selections, calculating self-distances')
            self.feat.add_distances(sel1)
        else:
            logging.info('Different selections, calculating inter-set distances')
            self.feat.add_distances(sel1, indices2=sel2)

    def _processTraj(self, traj):
        try:
            inp = coor.source(traj, self.feat)
        except:
            logging.warning('Could not read traj: ' + traj[0])
            return np.empty(0), np.empty(0)
        metric = np.squeeze(inp.get_output())

        if self.metric == 'contacts':
            metric = lil_matrix(metric <= self.threshold)

        return metric, inp.trajectory_lengths()

    def _getMapping(self, mol):
        return self.feat.describe()
'''

if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    metr = MetricDistance('protein and name CA', 'resname MOL and noh', metric='contacts')
    data = metr.project(mol)
    contframes = np.array([46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  75,  93,  94,
                           95,  96,  97,  98,  99, 100, 101, 102, 103, 114, 118, 119, 120,
                           121, 122, 123, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134,
                           135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 153, 154,
                           156, 157, 158, 159, 160, 161, 165, 166, 167, 168, 169, 170, 171,
                           172, 173, 174, 175, 176])
    assert np.all(np.unique(np.where(data)[0]) == contframes), 'Contact map calculation is broken'

    metr = MetricDistance('protein and name CA', 'resname MOL and noh', metric='distances')
    data = metr.project(mol)
    lastdists = np.array([ 32.41402817,  35.00286865,  37.82732391,  38.58649445,
                           39.26222229,  36.91499329,  37.93645477,  36.98008728,
                           33.43220901,  33.0732193 ,  30.82616615,  27.90420341,
                           28.2709713 ,  27.94139481,  24.85671616,  23.1801281 ,
                           25.09490013,  24.58997917,  20.71271324], dtype=np.float32)
    assert np.all(np.abs(data[-1, -20:-1] - lastdists) < 0.001), 'Distance calculation is broken'

