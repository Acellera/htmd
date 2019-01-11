# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
from htmd.projections.util import pp_calcDistances, pp_calcMinDistances
from htmd.util import ensurelist
import logging
logger = logging.getLogger(__name__)


class MetricDistance(Projection):
    """ Creates a MetricDistance object

    Parameters
    ----------
    sel1 : str
        Atom selection string for the first set of atoms.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    sel2 : str
        Atom selection string for the second set of atoms. If sel1 != sel2, it will calculate inter-set distances.
        If sel1 == sel2, it will calculate intra-set distances.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    groupsel1 : ['all','residue'], optional
        Group all atoms in `sel1` to the single minimum distance. Alternatively can calculate the minimum distance of a
        residue containing the atoms in sel1.
    groupsel2 : ['all','residue'], optional
        Same as groupsel1 but for `sel2`
    metric : ['distances','contacts'], optional
        Set to 'contacts' to calculate contacts instead of distances
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

    def project(self, mol):
        """ Project molecule.

        Parameters
        ----------
        mol : :class:`Molecule <htmdmol.molecule.Molecule>`
            A :class:`Molecule <htmdmol.molecule.Molecule>` object to project.

        Returns
        -------
        data : np.ndarray
            An array containing the projected data.
        """
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 1 and np.ndim(sel2) == 1:  # normal distances
            metric = pp_calcDistances(mol, sel1, sel2, self.metric, self.threshold, self.pbc, truncate=self.truncate)
        else:  # minimum distances by groups
            metric = pp_calcMinDistances(mol, sel1, sel2, self.metric, self.threshold, pbc=self.pbc)

        return metric

    def getMapping(self, mol):
        """ Returns the description of each projected dimension.

        Parameters
        ----------
        mol : :class:`Molecule <htmdmol.molecule.Molecule>` object
            A Molecule object which will be used to calculate the descriptions of the projected dimensions.

        Returns
        -------
        map : :class:`DataFrame <pandas.core.frame.DataFrame>` object
            A DataFrame containing the descriptions of each dimension
        """
        (sel1, sel2) = self._getSelections(mol)

        if np.ndim(sel1) == 2:
            protatoms = []
            for i in range(sel1.shape[0]):
                protatoms.append(np.where(sel1[i, :] == True)[0])
        else:
            protatoms = np.where(sel1)[0]
        if np.ndim(sel2) == 2:
            ligatoms = []
            for i in range(sel2.shape[0]):
                ligatoms.append(np.where(sel2[i, :] == True)[0])
        else:
            ligatoms = np.where(sel2)[0]

        numatoms1 = len(protatoms)
        numatoms2 = len(ligatoms)

        from pandas import DataFrame
        types = []
        indexes = []
        description = []
        if np.array_equal(sel1, sel2):
            for i in range(numatoms1):
                for j in range(i+1, numatoms1):
                    atm1 = protatoms[i]
                    atm2 = protatoms[j]
                    desc = '{} between {} {} {} and {} {} {}'.format(self.metric[:-1],
                                                                     mol.resname[atm1], mol.resid[atm1], mol.name[atm1],
                                                                     mol.resname[atm2], mol.resid[atm2], mol.name[atm2])
                    types += [self.metric[:-1]]
                    indexes += [[atm1, atm2]]
                    description += [desc]
        else:
            for j in range(numatoms2):
                for i in range(numatoms1):
                    atm1 = protatoms[i]
                    atm2 = ligatoms[j]
                    desc = '{} between {} {} {} and {} {} {}'.format(self.metric[:-1],
                                                                     mol.resname[atm1], mol.resid[atm1], mol.name[atm1],
                                                                     mol.resname[atm2], mol.resid[atm2], mol.name[atm2])
                    types += [self.metric[:-1]]
                    indexes += [[atm1, atm2]]
                    description += [desc]
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})

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
                sel = self._groupByResidue(mol, sel)
            else:
                raise NameError('Invalid groupsel argument')
        else:  # If user passes his own sets of groups
            sel = self._processMultiSelections(mol, sel)

        if np.sum(sel) == 0:
            raise NameError('Selection returned 0 atoms')
        return sel

    def _processMultiSelections(self, mol, sel):
        newsel = np.zeros((len(sel), mol.numAtoms), dtype=bool)
        for s in range(len(sel)):
            newsel[s, :] = mol.atomselect(sel[s])
        return newsel

    def _groupByResidue(self, mol, sel):
        import pandas as pd
        idx = mol.atomselect(sel, indexes=True)
        df = pd.DataFrame({'a': mol.resid[idx]})
        gg = df.groupby(by=df.a).groups  # Grouping by same resids

        newsel = np.zeros((len(gg), mol.numAtoms), dtype=bool)
        for i, res in enumerate(sorted(gg)):
            newsel[i, idx[gg[res]]] = True  # Setting the selected indexes to True which correspond to the same residue
        return newsel


class MetricSelfDistance(MetricDistance):
    def __init__(self, sel, groupsel=None, metric='distances', threshold=8, pbc=True, truncate=None):
        """ Creates a MetricSelfDistance object

        Parameters
        ----------
        sel : str
            Atom selection string for which to calculate the self distance.
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        groupsel : ['all','residue'], optional
            Group all atoms in `sel` to the single minimum distance. Alternatively can calculate the minimum distance
            of a residue containing the atoms in `sel`.
        metric : ['distances','contacts'], optional
            Set to 'contacts' to calculate contacts instead of distances
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


def contactVecToMatrix(vector, atomIndexes):
    from copy import deepcopy
    # Calculating the unique atom groups in the mapping
    uqAtomGroups = []
    atomIndexes = deepcopy(list(atomIndexes))
    for ax in atomIndexes:
        ax[0] = ensurelist(ax[0])
        ax[1] = ensurelist(ax[1])
        if ax[0] not in uqAtomGroups:
            uqAtomGroups.append(ax[0])
        if ax[1] not in uqAtomGroups:
            uqAtomGroups.append(ax[1])
    uqAtomGroups.sort(key=lambda x: x[0])  # Sort by first atom in each atom list
    num = len(uqAtomGroups)

    matrix = np.zeros((num, num), dtype=vector.dtype)
    mapping = np.ones((num, num), dtype=int) * -1
    for i in range(len(vector)):
        row = uqAtomGroups.index(atomIndexes[i][0])
        col = uqAtomGroups.index(atomIndexes[i][1])
        matrix[row, col] = vector[i]
        matrix[col, row] = vector[i]
        mapping[row, col] = i
        mapping[col, row] = i
    return matrix, mapping, uqAtomGroups


def reconstructContactMap(vector, mapping, truecontacts=None, plot=True, figsize=(7, 7), dpi=80, title=None, outfile=None, colors=None):
    """ Plots a given vector as a contact map

    Parameters
    ----------
    vector : np.ndarray or list
        A 1D vector of contacts
    mapping : pd.DataFrame
        A pandas DataFrame which describes the dimensions of the projection
    truecontacts : np.ndarray or list
        A 1D vector of true contacts
    plot : bool
        To plot or not to plot
    figsize : tuple
        The size of the final plot in inches
    dpi : int
        Dots per inch
    outfile : str
        Path of file in which to save the plot

    Returns
    -------
    cm : np.ndarray
        The input vector converted into a 2D numpy array

    Examples
    --------
    >>> reconstructContactMap(contacts, mapping)
    To use it with distances instead of contacts pass ones as the concat vector
    >>> reconstructContactMap(np.ones(dists.shape, dtype=bool), mapping, colors=dists)
    """

    from copy import deepcopy
    from matplotlib import cm as colormaps
    if truecontacts is None:
        truecontacts = np.zeros(len(vector), dtype=bool)
    if len(vector) != len(mapping):
        raise RuntimeError('Vector and map length must match.')

    # Checking if contacts or distances exist in the data
    contactidx = mapping.type == 'contact'
    if not np.any(contactidx):
        contactidx = mapping.type == 'distance'
        if not np.any(contactidx):
            raise RuntimeError(
                'No contacts or distances found in the MetricData object. Check the `.map` property of the object for a description of your projection.')

    # Creating the 2D contact maps
    cm, newmapping, uqAtomGroups = contactVecToMatrix(vector, mapping.atomIndexes)
    cmtrue, _, _ = contactVecToMatrix(truecontacts, mapping.atomIndexes)
    num = len(uqAtomGroups)

    if plot:
        from matplotlib import pylab as plt
        f = plt.figure(figsize=figsize, dpi=dpi)
        plt.imshow(cmtrue / 2, interpolation='none', vmin=0, vmax=1, aspect='equal',
                   cmap='Greys')  # /2 to convert to gray from black

        rows, cols = np.where(cm)
        colorbar = False
        if colors is None:
            truecms = vector & truecontacts
            colors = np.array(['r']*len(vector), dtype=object)
            colors[truecms] = '#ffff00'
        elif isinstance(colors, np.ndarray) and isinstance(colors[0], float):
            mpbl = colormaps.ScalarMappable(cmap=colormaps.jet)
            mpbl.set_array(colors)
            colors = mpbl.to_rgba(colors)
            colorbar = True
        if len(colors) == len(vector):
            colors = colors[newmapping[rows, cols]]

        plt.scatter(rows, cols, s=figsize[0] * 5, marker='o', c=colors, lw=0)
        if colorbar:
            plt.colorbar(mpbl)

        ax = f.axes[0]
        # Major ticks
        ax.set_xticks(np.arange(0, num, 1))
        ax.set_yticks(np.arange(0, num, 1))

        # Labels for major ticks
        ax.set_xticklabels([x[0] for x in uqAtomGroups])
        ax.set_yticklabels([x[0] for x in uqAtomGroups], )

        # Minor ticks
        ax.set_xticks(np.arange(-.5, num, 1), minor=True)
        ax.set_yticks(np.arange(-.5, num, 1), minor=True)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='#969696', linestyle='-', linewidth=1)
        ax.tick_params(axis='both', which='both', length=0)
        plt.xlim([-.5, num - .5])
        plt.ylim([-.5, num - .5])
        plt.xlabel('Atom index')
        plt.ylabel('Atom index')
        if title:
            plt.title(title)
        if outfile is not None:
            plt.savefig(outfile, dpi=dpi, bbox_inches='tight', pad_inches=0.2)
            plt.close()
        else:
            plt.show()
    return cm

import unittest
from htmd.home import home
import os
class TestMetricDistance(unittest.TestCase):
    @classmethod
    def setUpClass(self):
            from htmdmol.molecule import Molecule
            from os import path
            self.mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
            self.mol_skipped = self.mol.copy()

            self.mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
            self.mol_skipped.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'), skip=10)

    def test_contacts(self):
        metr = MetricDistance('protein and name CA', 'resname MOL and noh', metric='contacts')
        data = metr.project(self.mol)
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'contacts.npy'))
        assert np.allclose(data, refdata, atol=1e-3), 'Contact calculation is broken'

    def test_distances(self):
        metr = MetricDistance('protein and name CA', 'resname MOL and noh', metric='distances')
        data = metr.project(self.mol)
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'distances.npy'))
        assert np.allclose(data, refdata, atol=1e-3), 'Distance calculation is broken'

    def test_mindistances(self):
        metr = MetricDistance('protein and noh', 'resname MOL and noh', groupsel1='residue', groupsel2='all')
        data = metr.project(self.mol)
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'mindistances.npy'))
        assert np.allclose(data, refdata, atol=1e-3), 'Minimum distance calculation is broken'

    def test_selfmindistance_manual(self):
        metr = MetricDistance('protein and resid 1 to 50 and noh', 'protein and resid 1 to 50 and noh', groupsel1='residue', groupsel2='residue')
        data = metr.project(self.mol)
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'selfmindistance.npy'))
        assert np.allclose(data, refdata, atol=1e-3), 'Manual self-distance is broken'

    def test_selfmindistance_auto(self):
        metr = MetricSelfDistance('protein and resid 1 to 50 and noh', groupsel='residue')
        data = metr.project(self.mol)
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'selfmindistance.npy'))
        assert np.allclose(data, refdata, atol=1e-3), 'Automatic self-distance is broken'

    def test_mindistances_skip(self):
        metr = MetricSelfDistance('protein and resid 1 to 50 and noh', groupsel='residue')
        data = metr.project(self.mol_skipped)  
        refdata = np.load(os.path.join(home(), 'data', 'metricdistance', 'selfmindistance.npy'))
        assert np.allclose(data, refdata[::10, :], atol=1e-3), 'Minimum distance calculation with skipping is broken'


if __name__ == '__main__':
    unittest.main(verbosity=2)
