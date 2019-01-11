# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.metriccoordinate import MetricCoordinate as _MetricCoordinate
from htmdmol.util import sequenceID
import numpy as np


class MetricFluctuation(_MetricCoordinate):
    """ Creates a MetricFluctuation object that calculates the squared fluctuation of atom positions in trajectories.

    Depending on the `refpos` option the projection either returns the fluctuation from the mean position of the atoms,
    or from coordinates of reference given in `refmol`. This means it calculates for atom coordinates (x,y,z) and
    reference coordinates (xr, yr, yz): (xr-x)**2+(yr-y)**2+(zr-z)**2. If `groupsel` is set to `residue` it will
    calculate the mean of this value over all atoms of the residue.
    To then get the RMSD you need to compute the square root (of the mean) of the results this projection returns over
    the desired atoms.

    Parameters
    ----------
    refmol : :class:`Molecule <htmdmol.molecule.Molecule>` object
        The reference Molecule to which we will align.
    atomsel : str
        Atom selection string for the atoms whose fluctuations we want to calculate.
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
    refpos : str
        Set to 'mean' to calculate the fluctuation from the trajectory mean. Set to 'refmol' to calculate the
        fluctuation from the reference structure.
    groupsel : str or None
        Set to None to get the fluctuation per atom. Set to 'residue' to get the mean fluctuation of the residue by
        grouping all of its atoms given in `atomsel`

    Returns
    -------
    metr : MetricFluctuation object

    Examples
    --------

    """
    def __init__(self, refmol, atomsel, trajalnstr='protein and name CA', refalnstr=None, centerstr='protein', pbc=True, refpos='mean', groupsel=None):
        super().__init__(refmol, atomsel, trajalnstr, refalnstr, centerstr, pbc)
        self._refpos = refpos
        self._groupsel = groupsel

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
        coords = super().project(mol)

        # TODO: Could precalculate some stuff like this to speed it up
        resids = sequenceID(self._refmol.resid)

        if self._refpos == 'mean':
            refcoords = np.mean(coords, axis=0)
        elif self._refpos == 'refmol':
            refcoords = _MetricCoordinate(self._refmol, self._atomsel).project(self._refmol)
        else:
            raise RuntimeError('Wrong refpos option')

        mapping = super().getMapping(self._refmol)
        xyzgroups = mapping.groupby('atomIndexes').groups
        numatoms = len(xyzgroups)

        atomfluct = np.zeros((coords.shape[0], numatoms))
        squarediff = (coords - refcoords) ** 2
        atomresids = np.zeros(numatoms, dtype=int)
        for i, atom in enumerate(sorted(xyzgroups.values(), key=lambda x: x[0])):
            assert len(np.unique(mapping.atomIndexes[atom])) == 1
            atomfluct[:, i] = squarediff[:, atom].sum(axis=1)
            atomresids[i] = resids[int(mapping.atomIndexes[atom[0]])]

        if self._groupsel == 'residue':
            numres = len(np.unique(atomresids))
            meanresfluct = np.zeros((coords.shape[0], numres))
            for i, r in enumerate(np.unique(atomresids)):
                meanresfluct[:, i] = atomfluct[:, atomresids == r].mean(axis=1)
            return meanresfluct
        elif self._groupsel is None:
            return atomfluct

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
        (xxx, atomsel, yyy) = self._getSelections(mol)
        atomidx = np.where(atomsel)[0]
        from pandas import DataFrame
        types = []
        indexes = []
        description = []
        if self._groupsel is None:
            for i in atomidx:
                types += ['fluctuation']
                indexes += [i]
                description += ['Fluctuation of {} {} {}'.format(mol.resname[i], mol.resid[i], mol.name[i])]
        elif self._groupsel == 'residue':
            resids = mol.resid[atomidx]
            for r in np.unique(resids):
                types += ['fluctuation']
                i = atomidx[np.where(resids == r)[0][0]]
                indexes += [i]
                description += ['Mean fluctuation of {} {}'.format(mol.resname[i], mol.resid[i])]

        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})
