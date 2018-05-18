# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
from htmd.molecule.molecule import UniqueAtomID
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Dihedral:
    """ Class to store atoms defining a dihedral angle.

    Example
    -------
    >>> # Using the helper functions to construct a dihedral object
    >>> d1 = Dihedral.phi(mol, 5, 6, segid='P0')
    >>> d2 = Dihedral.chi1(mol, 12, segid='P0') # Defining segid
    >>> d3 = Dihedral.chi1(mol, 38, chain='X') # Defining chain
    >>> # Manual construction
    >>> atom1 = {'name': 'N', 'resid': 5, 'segid': 'P'}
    >>> atom2 = {'name': 'CA', 'resid': 3, 'segid': 'P', 'chain': 'A', 'insertion': 'B'}
    >>> atom3 = {'name': 'C', 'resid': 46, 'segid': 'P', 'chain': 'X'}
    >>> atom4 = {'name': 'O', 'resid': 2, 'segid': 'P'}
    >>> d = Dihedral(atom1, atom2, atom3, atom4)
    """
    def __init__(self, atom1, atom2, atom3, atom4, dihedraltype=None):
        self.atoms = [atom1, atom2, atom3, atom4]
        self.dihedraltype = dihedraltype

    def __str__(self):
        descr = ''
        if self.dihedraltype is not None:
            descr = '"{}" dihedral angle including atoms:\n'.format(self.dihedraltype)
        descr += '\n'.join(map(str, self.atoms))
        return descr

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def dihedralsToIndexes(mol, dihedrals):
        """ Converts dihedral objects to atom indexes of a given Molecule
        
        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain atom information
        dihedrals : list
            A single dihedral or a list of Dihedral objects_

        Returns
        -------
        indexes : list of lists
            A list containing a list of atoms that correspond to each dihedral.
            
        Examples
        --------
        >>> dihs = []
        >>> dihs.append(Dihedral.phi(mol, 'resid 1', 'resid 2'))
        >>> dihs.append(Dihedral.psi(mol, 'resid 2', 'resid 3'))
        >>> indexes = Dihedral.dihedralsToIndexes(mol, dihs)
        """
        from htmd.util import ensurelist
        indexes = []
        for dih in ensurelist(dihedrals):
            idx = []
            for a in dih.atoms:
                idx.append(a.selectAtom(mol))
            indexes.append(idx)
        return indexes

    def getIndices(self, mol):
        idx = []
        for a in self.atoms:
            idx.append(a.selectAtom(mol))
        return idx

    @staticmethod
    def _buildDihedral(mol, dihname, dihatoms, sellist):
        from htmd.util import ensurelist
        sellist = ensurelist(sellist)

        atoms = []
        if len(sellist) != 1:
            for sel, dihatomnames in zip(sellist, dihatoms):
                if isinstance(sel, str):
                    sel = mol.atomselect(sel, indexes=True)
                for atomname in dihatomnames:
                    idx = sel[np.where(mol.name[sel] == atomname)[0][0]]
                    atoms.append(UniqueAtomID.fromMolecule(mol, idx=idx))
        else:
            sel = sellist[0]
            if isinstance(sel, str):
                sel = mol.atomselect(sel, indexes=True)
            for atomname in dihatoms:
                idx = sel[np.where(mol.name[sel] == atomname)[0][0]]
                atoms.append(UniqueAtomID.fromMolecule(mol, idx=idx))
        return Dihedral(atoms[0], atoms[1], atoms[2], atoms[3], dihedraltype=dihname)


    @staticmethod
    def proteinDihedrals(mol, sel='protein', dih=('psi', 'phi')):
        """ Returns a list of tuples containing the four resid/atom pairs for each dihedral of the protein

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            Atom selection string to restrict the atoms for which to calculate dihedrals (e.g. only one of many chains).
            See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
        dih : tuple
            A tuple of the dihedral types we want to calculate (phi, psi, omega, chi1, chi2, chi3, chi4, chi5)

        Returns
        -------
        dihedrals : list of :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` objects
            A list of Dihedral objects
        """
        from htmd.builder.builder import sequenceID
        mol = mol.copy()
        mol.filter(sel, _logger=False)

        segments = sequenceID((mol.chain, mol.segid)) # Here I consider segments both chains and segments
        residues = sequenceID((mol.resid, mol.insertion, mol.chain, mol.segid))

        dihedrals = []
        for s in np.unique(segments):
            for r in np.unique(residues[segments == s]):
                res1 = np.where(residues == r)[0]

                if 'psi' in dih and r != residues.max():  # No psi angle for last residue
                    res2 = np.where(residues == (r + 1))[0]
                    dihedrals.append(Dihedral.psi(mol, res1, res2))
                if 'phi' in dih and r != 0:  # No phi angle for first residue
                    res2 = np.where(residues == (r - 1))[0]
                    dihedrals.append(Dihedral.phi(mol, res2, res1))
                if 'omega' in dih and r != residues.max():
                    res2 = np.where(residues == (r + 1))[0]
                    dihedrals.append(Dihedral.omega(mol, res1, res2))
                if 'chi1' in dih:
                    dihedrals.append(Dihedral.chi1(mol, res1))
                if 'chi2' in dih:
                    dihedrals.append(Dihedral.chi2(mol, res1))
                if 'chi3' in dih:
                    dihedrals.append(Dihedral.chi3(mol, res1))
                if 'chi4' in dih:
                    dihedrals.append(Dihedral.chi4(mol, res1))
                if 'chi5' in dih:
                    dihedrals.append(Dihedral.chi5(mol, res1))
        return [d for d in dihedrals if d is not None]

    # Sidechain dihedral atoms taken from
    # http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    @staticmethod
    def phi(mol, sel1, sel2):
        """ Constructs a Dihedral object corresponding to the phi angle of res1 and res2

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        res1 : int
            The resid of the first residue containing the C atom
        res2 : int
            The resid of the second residue containing the N CA C atoms
        segid : str
            The segment id of the residues
        chain : str
            The chain letter of the residues
        insertion1 : str
            The insertion letter of residue 1
        insertion2 : str
            The insertion letter of residue 2

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        phi = (('C'), ('N', 'CA', 'C'))
        return Dihedral._buildDihedral(mol, 'phi', phi, [sel1, sel2])

    @staticmethod
    def psi(mol, sel1, sel2):
        """ Constructs a Dihedral object corresponding to the psi angle of res1 and res2

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        res1 : int
            The resid of the first residue containing the N CA C atoms
        res2 : int
            The resid of the second residue containing the N atom
        segid : str
            The segment id of the residues
        chain : str
            The chain letter of the residues
        insertion1 : str
            The insertion letter of residue 1
        insertion2 : str
            The insertion letter of residue 2

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        psi = (('N', 'CA', 'C'), ('N'))

        return Dihedral._buildDihedral(mol, 'psi', psi, [sel1, sel2])


    @staticmethod
    def omega(mol, sel1, sel2):
        """ Constructs a Dihedral object corresponding to the omega angle of res1 and res2

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        res1 : int
            The resid of the first residue containing the CA C atoms
        res2 : int
            The resid of the second residue containing the N CA atoms
        segid : str
            The segment id of the residues
        chain : str
            The chain letter of the residues
        insertion1 : str
            The insertion letter of residue 1
        insertion2 : str
            The insertion letter of residue 2

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        omega = (('CA', 'C'), ('N', 'CA'))

        return Dihedral._buildDihedral(mol, 'omega', omega, [sel1, sel2])

    @staticmethod
    def chi1(mol, sel):
        """ Constructs a Dihedral object corresponding to the chi1 angle of a residue

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            Atom selection of the residue

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        chi1std = ('N', 'CA', 'CB', 'CG')
        chi1 = {'ARG': chi1std, 'ASN': chi1std, 'ASP': chi1std, 'CYS': ('N', 'CA', 'CB', 'SG'), 'GLN': chi1std,
                'GLU': chi1std, 'HIS': chi1std, 'ILE': ('N', 'CA', 'CB', 'CG1'), 'LEU': chi1std, 'LYS': chi1std,
                'MET': chi1std, 'PHE': chi1std, 'PRO': chi1std, 'SER': ('N', 'CA', 'CB', 'OG'),
                'THR': ('N', 'CA', 'CB', 'OG1'), 'TRP': chi1std, 'TYR': chi1std, 'VAL': ('N', 'CA', 'CB', 'CG1')}

        if isinstance(sel, str):
            sel = mol.atomselect(sel, indexes=True)
        resname = np.unique(mol.resname[sel])[0]
        if resname not in chi1:
            return None
        chi1 = chi1[resname]

        return Dihedral._buildDihedral(mol, 'chi1', chi1, [sel,])

    @staticmethod
    def chi2(mol, sel):
        """ Constructs a Dihedral object corresponding to the chi2 angle of a residue

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            Atom selection of the residue

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        chi2std = ('CA', 'CB', 'CG', 'CD')
        chi2 = {'ARG': chi2std, 'ASN': ('CA', 'CB', 'CG', 'OD1'), 'ASP': ('CA', 'CB', 'CG', 'OD1'), 'GLN': chi2std,
                'GLU': chi2std, 'HIS': ('CA', 'CB', 'CG', 'ND1'), 'ILE': ('CA', 'CB', 'CG1', 'CD1'),
                'LEU': ('CA', 'CB', 'CG', 'CD1'), 'LYS': chi2std,
                'MET': ('CA', 'CB', 'CG', 'SD'), 'PHE': ('CA', 'CB', 'CG', 'CD1'), 'PRO': chi2std,
                'TRP': ('CA', 'CB', 'CG', 'CD1'), 'TYR': ('CA', 'CB', 'CG', 'CD1')}

        if isinstance(sel, str):
            sel = mol.atomselect(sel, indexes=True)
        resname = np.unique(mol.resname[sel])[0]
        if resname not in chi2:
            return None
        chi2 = chi2[resname]

        return Dihedral._buildDihedral(mol, 'chi2', chi2, [sel,])

    @staticmethod
    def chi3(mol, sel):
        """ Constructs a Dihedral object corresponding to the chi3 angle of a residue

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            Atom selection of the residue

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        chi3 = {'ARG': ('CB', 'CG', 'CD', 'NE'), 'GLN': ('CB', 'CG', 'CD', 'OE1'), 'GLU': ('CB', 'CG', 'CD', 'OE1'),
                'LYS': ('CB', 'CG', 'CD', 'CE'), 'MET': ('CB', 'CG', 'SD', 'CE')}

        if isinstance(sel, str):
            sel = mol.atomselect(sel, indexes=True)
        resname = np.unique(mol.resname[sel])[0]
        if resname not in chi3:
            return None
        chi3 = chi3[resname]

        return Dihedral._buildDihedral(mol, 'chi3', chi3, [sel,])

    @staticmethod
    def chi4(mol, sel):
        """ Constructs a Dihedral object corresponding to the chi4 angle of a residue

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            Atom selection of the residue

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object
        """
        chi4 = {'ARG': ('CG', 'CD', 'NE', 'CZ'), 'LYS': ('CG', 'CD', 'CE', 'NZ')}

        if isinstance(sel, str):
            sel = mol.atomselect(sel, indexes=True)
        resname = np.unique(mol.resname[sel])[0]
        if resname not in chi4:
            return None
        chi4 = chi4[resname]

        return Dihedral._buildDihedral(mol, 'chi4', chi4, [sel,])

    @staticmethod
    def chi5(mol, sel):
        """ Constructs a Dihedral object corresponding to the chi5 angle of a residue

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A Molecule object from which to obtain structural information
        sel : str
            The atomselection for the residue

        Returns
        -------
        dihedral : :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
            A Dihedral object

        Examples
        --------
        >>> Dihedral.chi5(mol, 'resid 5 and chain X')
        """
        chi5 = {'ARG': ('CD', 'NE', 'CZ', 'NH1')}

        if isinstance(sel, str):
            sel = mol.atomselect(sel, indexes=True)
        resname = np.unique(mol.resname[sel])[0]
        if resname not in chi5:
            return None
        chi5 = chi5[resname]

        return Dihedral._buildDihedral(mol, 'chi5', chi5, [sel,])


class MetricDihedral(Projection):
    """ Calculates a set of dihedral angles from trajectories

    Parameters
    ----------
    dih : list of :class:`Dihedral <htmd.projections.metricdihedral.Dihedral>` object
        You can provide your own list of Dihedral objects. See example.
    sincos : bool, optional
        Set to True to return the dihedral angles as their sine and cosine components. Makes them periodic.
    protsel : str, optional
        If no `dih` are passed, `MetricDihedral` will automatically calculate the phi/psi dihedral angles of `protsel`.

    Examples
    --------
    >>> mol = Molecule('3PTB')
    >>> mol.filter('not insertion A')
    >>> met = MetricDihedral(protsel='protein and segid 0')
    >>> met.project(mol)
    >>> # More complicated example
    >>> dih = []
    >>> dih.append(Dihedral.chi1(mol, 'resid 45'))
    >>> dih.append(Dihedral.psi(mol, 'resid 29', 'resid 30'))
    >>> met = MetricDihedral(dih)
    >>> met.project(mol)
    >>> met.getMapping(mol)
    """
    def __init__(self, dih=None, sincos=True, protsel='protein'):
        if dih is not None and not isinstance(dih[0], Dihedral):
            raise RuntimeError('Manually passing dihedrals to MetricDihedral requires use of the Dihedral class. Check the example in the documentation')
        self._protsel = protsel
        self._sincos = sincos
        self._dihedrals = dih
        self._pc_dih = None

    def _precalculate(self, mol):
        self._pc_dih = self._dihedralAtomsPrecalc(mol, mol.atomselect(self._protsel))

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
        dih = self._getSelections(mol)
        return self._calcDihedralAngles(mol, dih, sincos=self._sincos)

    def _getSelections(self, mol):
        if self._pc_dih is not None:
            return self._pc_dih
        else:
            return self._dihedralAtomsPrecalc(mol, mol.atomselect(self._protsel))

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
        dihedrals = self._getSelections(mol)

        from pandas import DataFrame
        types = []
        indexes = []
        description = []
        for i, dih in enumerate(dihedrals):
            types += ['dihedral']
            indexes += [dih]
            mapstr = ''
            for d in dih:
                mapstr += '({} {} {} {} {}) '.format(mol.resname[d], mol.resid[d], mol.name[d], mol.segid[d], mol.chain[d])
            if self._sincos:
                description += ['Sine of angle of ' + mapstr]
                description += ['Cosine of angle of ' + mapstr]
                types += ['dihedral']
                indexes += [dih]
            else:
                description += ['Angle of ' + mapstr]
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})

    def _dihedralAtomsPrecalc(self, mol, protsel):
        if self._dihedrals is None:  # Default phi psi dihedrals
            dihedrals = Dihedral.proteinDihedrals(mol, protsel)
        else:
            from htmd.util import ensurelist
            self._dihedrals = ensurelist(self._dihedrals)
            dihedrals = self._dihedrals

        return Dihedral.dihedralsToIndexes(mol, dihedrals)

    def _calcDihedralAngles(self, mol, dihedrals, sincos=True):
        from htmd.numbautil import dihedralAngle
        metric = np.zeros((np.size(mol.coords, 2), len(dihedrals)))

        for i, dih in enumerate(dihedrals):
            metric[:, i] = np.rad2deg(dihedralAngle(mol.coords[dih, :, :]))

        if sincos:
            sc_metric = np.zeros((np.size(metric, 0), np.size(metric, 1) * 2))
            sc_metric[:, 0::2] = np.sin(metric * np.pi / 180.)
            sc_metric[:, 1::2] = np.cos(metric * np.pi / 180.)
            metric = sc_metric
        return metric.astype(np.float32)

    @staticmethod
    def phi(res1, res2):
        raise DeprecationWarning('MetricDihedral.phi has been replaced with Dihedral.phi. Check the latest documentation.')

    @staticmethod
    def psi(res1, res2):
        raise DeprecationWarning('MetricDihedral.psi has been replaced with Dihedral.psi. Check the latest documentation.')

    @staticmethod
    def omega(res1, res2):
        raise DeprecationWarning('MetricDihedral.omega has been replaced with Dihedral.omega. Check the latest documentation.')

    @staticmethod
    def chi1(res, resname):
        raise DeprecationWarning('MetricDihedral.chi1 has been replaced with Dihedral.chi1. Check the latest documentation.')

    @staticmethod
    def chi2(res, resname):
        raise DeprecationWarning('MetricDihedral.chi2 has been replaced with Dihedral.chi2. Check the latest documentation.')

    @staticmethod
    def chi3(res, resname):
        raise DeprecationWarning('MetricDihedral.chi3 has been replaced with Dihedral.chi3. Check the latest documentation.')

    @staticmethod
    def chi4(res, resname):
        raise DeprecationWarning('MetricDihedral.chi4 has been replaced with Dihedral.chi4. Check the latest documentation.')

    @staticmethod
    def chi5(res, resname):
        raise DeprecationWarning('MetricDihedral.chi5 has been replaced with Dihedral.chi5. Check the latest documentation.')


if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from htmd.builder.builder import autoSegment
    from os import path
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    metr = MetricDihedral(protsel='protein')
    data = metr.project(mol)
    dataref = np.load(path.join(home(), 'data', 'metricdihedral', 'ref.npy'))
    assert np.allclose(data, dataref, atol=1e-03), 'Dihedrals calculation gave different results from reference'

    mol = Molecule('5MAT')
    mol.filter('not insertion A and not altloc A B')
    mol = autoSegment(mol)
    data = MetricDihedral().project(mol)
    dataref = np.load(path.join(home(), 'data', 'metricdihedral', '5mat.npy'))
    assert np.allclose(data, dataref, atol=1e-03), 'Dihedrals calculation gave different results from reference'
