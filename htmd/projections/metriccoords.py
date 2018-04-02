# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricCoords(Projection):
    """ Creates a MetricCoords object that seletct the atom coordinates from a set of trajectories.
    Selections are computed once the object is created!

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object, optional
        The  Molecule on which to apply the selections.
    atomsel : str
        Atom selection string for the atoms whose coordinates we want to calculate.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    refmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object, optional
        The reference Molecule to which align the trajectories.
    alignsel : str, optional
        Atom selection string for the trajectories from which to align to the reference structure.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    refsel : str, optional
        Atom selection string for `refmol` from which to align to the reference structure. If None, it defaults to the
        same as `trajalnstr`. See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    pbccentersel : str, optional
        Atom selection string around which to wrap the trajectories for periodic boundary condition, if None is not wrapped
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__

    Returns
    -------
    metr : MetricCoords object
    """
    def __init__(self, mol, atomsel, refmol=None, alignsel=None, refsel=None, pbccentersel=None):
        self._atomsel = atomsel
        self._refmol = refmol
        self._refsel = refsel
        self._alignsel = alignsel
        self._pbccentersel = pbccentersel
        self._align = None
        self._atoms = None
        self._centers = None
        self._refatoms = None 
        #Compute them
        self._computeSelections(mol)

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

        mol = mol.copy() #internal temp copy
        if self._pbccentersel is not None:
            mol.wrap(self._centers)
        if self._alignsel is not None:
            mol.align(self._align, refmol=self._refmol, refsel=self._refsel)

        coords = np.concatenate((mol.coords[self._atoms, 0, :], mol.coords[self._atoms, 1, :], mol.coords[self._atoms, 2, :]))
        return np.transpose(coords)  # I need to permute the dimensions here to put frames on the first

    def _computeSelections(self, mol):
        # _atomsel should always be there
        self._atoms = mol.atomselect(self._atomsel)
        if np.sum(self._atoms) == 0:
            raise NameError('Atom selection resulted in 0 atoms.')

        if self._alignsel is not None:
            self._align = mol.atomselect(self._alignsel)
            if np.sum(self._align) == 0:
                raise NameError('Alignment selection resulted in 0 atoms.')

        if self._pbccentersel is not None: 
            self._centers = mol.atomselect(self._pbccentersel)
            if np.sum(self._centers) == 0:
                raise NameError('Wrapping center selection resulted in 0 atoms.')

        if self._refsel is not None: 
            self._centers = mol.atomselect(self._pbccentersel)
            if np.sum(self._centers) == 0:
                raise NameError('Wrapping center selection resulted in 0 atoms.')

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
        self._computeSelections(mol)
        atoms = self._atoms
        atomidx = np.where(atoms)[0]
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
    metr = MetricCoords(mol, 'protein and name CA',refmol=ref,alignsel='protein and name CA',refsel='protein and name CA',pbccentersel='protein')
    data = metr.project(mol)

    lastcoors = np.array([6.79283285,   5.55226946,   4.49387407,   2.94484425,
                          5.36937141,   3.18590879,   5.75874281,   5.48864174,
                          1.69625032,   1.58790839,   0.57877392,  -2.66498065,
                          -3.70919156,  -3.33702421,  -5.38465405,  -8.43286991,
                          -8.15859032,  -7.85062265, -10.92551327, -13.70733166], dtype=np.float32)
    assert np.all(np.abs(data[-1, -20:] - lastcoors) < 0.001), 'Coords calculation is broken'
