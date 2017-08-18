# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricCmDistance(Projection):

    """ Creates a MetricCmDistance object that calculates the distances of the mass centers of two selections

    Parameters
    ----------
    sel1 : str
        Atomselection for the first set of atoms

    sel2 : str
        Atomselection for the first set of atoms
    
    Returns
    -------
    metr : MetricDistance object
    """


    def __init__(self, sel1, sel2):
        self.sel1 = sel1
        self.sel2 = sel2
        self.precalcsel1 = None
        self.precalcsel2 = None
        
    def _precalculate(self, mol):
        self.precalcsel1 = self._processSelection(mol, self.sel1, None)
        self.precalcsel2 = self._processSelection(mol, self.sel2, None)

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

        (sel1, sel2) = self._getSelections(mol)

        sel1_coords = []
        sel2_coords = []
        for f in range(mol.numFrames):
            mol.frame=f
            sel1_coords.append(mol.get('coords', sel1))
            sel2_coords.append(mol.get('coords', sel2))
        

        cms_sel1 =[ self._getMolCenter(coords) 
                    for coords in sel1_coords ]
        cms_sel2 =[ self._getMolCenter( coords) 
                    for coords in sel2_coords ]


        cms_dist =[ np.linalg.norm(cm_s1 - cm_s2) 
                    for cm_s1, cm_s2 in zip(cms_sel1, cms_sel2) ]
        

        metric = np.atleast_1d( cms_dist)
        

        return metric

    def _getMolCenter(self, coords):
        return np.mean(coords, axis=-2)

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
        (sel1, sel2) = self._getSelections(mol)

        sel1_definition = np.unique(mol.get('resid', sel1))
        sel2_definition = np.unique(mol.get('resid', sel2))


        from pandas import DataFrame
        types = ['cmDist']
        indexes = [ [ np.where(sel1)[0], np.where(sel2)[0] ] ]
        description = ['Mass center distances between resids {} and resids {}'.format(sel1_definition, sel2_definition)]

        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})



    def _getSelections(self, mol):
        # If they have been pre-calculated return them.
        if self.precalcsel1 is not None and self.precalcsel2 is not None:
            return self.precalcsel1, self.precalcsel2
        else:  # Otherwise compute them.
            sel1 = self._processSelection(mol, self.sel1, None)
            sel2 = self._processSelection(mol, self.sel2, None)
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