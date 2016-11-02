# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricDihedral(Projection):
    """ Calculates a set of dihedral angles from trajectories

    Parameters
    ----------
    dih : list of str
        You can provide your own list of 4 resid/atom pairs for dihedral angles. See example.
    sincos : bool, optional
        Set to True to return the dihedral angles as their sine and cosine components. Makes them periodic.
    protsel : str, optional
        Atomselection for the protein segment for which to calculate dihedral angles. Resids should be unique within
        that segment.

    Examples
    --------
    >>> dih = []
    >>> dih.append(MetricDihedral.chi2(4, 'ARG'))
    >>> dih.append(((15, 'N'), (15, 'CA'), (15, 'C'), (16, 'N')))  # psi angle
    >>> dih.append(MetricDihedral.psi(15, 16))  # psi angle (equivalent to previous command, simple version)
    >>> met = MetricDihedral(dih, 'protein and segid P0')
    """
    def __init__(self, dih=None, sincos=True, protsel='protein'):
        self._protsel = protsel
        self._sincos = sincos
        self._dihquads = dih  # TODO: Calculate the dihedral
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
                mapstr += '({} {} {}) '.format(mol.resname[d], mol.resid[d], mol.name[d])
            if self._sincos:
                description += ['Sine of angle of ' + mapstr]
                description += ['Cosine of angle of ' + mapstr]
                types += ['dihedral']
                indexes += [dih]
            else:
                description += ['Angle of ' + mapstr]
        return DataFrame({'type': types, 'atomIndexes': indexes, 'description': description})

    def _dihedralAtomsPrecalc(self, mol, protsel):
        protatoms = mol.atomselect(protsel)
        resids = np.unique(mol.get('resid', sel=protatoms))

        if self._dihquads is None:  # Default phi psi dihedrals
            dihquads = self._proteinDihedrals(resids)
        else:
            if not isinstance(self._dihquads, list):
                self._dihquads = [self._dihquads]
            dihquads = self._dihquads

        dih = []
        for quad in dihquads:
            dih.append(self._quadAtomIndexes(quad, protatoms, mol))
        return dih

    def _quadAtomIndexes(self, quad, selatoms, mol):
        sels = []
        for p in quad:
            # atomsel = mol.atomselect('resid {} and name {}'.format(p[0], p[1]))
            atomsel = (mol.resid == p[0]) & (mol.name == p[1])
            atomsel = atomsel & selatoms
            if np.sum(atomsel) != 1:
                raise RuntimeError('Expected one atom from atomselection resid {} name {}. Got {} instead.'.format(p[0], p[1], np.sum(atomsel)))
            sels.append(np.where(atomsel)[0][0])
        return sels

    def _calcDihedralAngles(self, mol, dihedrals, sincos=True):
        metric = np.zeros((np.size(mol.coords, 2), len(dihedrals)))

        for i, dih in enumerate(dihedrals):
            metric[:, i] = _angle(mol.coords[dih[0], :, :], mol.coords[dih[1], :, :],
                                  mol.coords[dih[2], :, :], mol.coords[dih[3], :, :])

        if sincos:
            sc_metric = np.zeros((np.size(metric, 0), np.size(metric, 1) * 2))
            sc_metric[:, 0::2] = np.sin(metric * np.pi / 180.)
            sc_metric[:, 1::2] = np.cos(metric * np.pi / 180.)
            metric = sc_metric
        return metric.astype(np.float32)

    @staticmethod
    def _proteinDihedrals(resids, resnames=None, dih=('psi', 'phi')):
        """ Returns a list of tuples containing the four resid/atom pairs for each dihedral of the protein

        Parameters
        ----------
        resids : list
            A list of the protein resids

        Returns
        -------
        quads : list
            A list of quads
        """
        dihquads = []
        for r in range(len(resids)):
            if 'psi' in dih and r != len(resids)-1:  # No psi angle for last residue
                dihquads.append(MetricDihedral.psi(resids[r], resids[r + 1]))
            if 'phi' in dih and r != 0:  # No phi angle for first residue
                dihquads.append(MetricDihedral.phi(resids[r - 1], resids[r]))
            if 'chi1' in dih:
                dihquads.append(MetricDihedral.chi1(resids[r], resnames[r]))
            if 'chi2' in dih:
                dihquads.append(MetricDihedral.chi2(resids[r], resnames[r]))
            if 'chi3' in dih:
                dihquads.append(MetricDihedral.chi3(resids[r], resnames[r]))
            if 'chi4' in dih:
                dihquads.append(MetricDihedral.chi4(resids[r], resnames[r]))
            if 'chi5' in dih:
                dihquads.append(MetricDihedral.chi5(resids[r], resnames[r]))
        return [d for d in dihquads if d is not None]

    # Sidechain dihedral atoms taken from
    # http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    @staticmethod
    def phi(res1, res2):
        """ Get a set of four resid/atom pairs corresponding to the phi angle of res1 and res2

        Parameters
        ----------
        res1 : int
            The first residue containing the C atom
        res2 : int
            The second residue containing the N CA C atoms

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        return (res1, 'C'), (res2, 'N'), (res2, 'CA'), (res2, 'C')

    @staticmethod
    def psi(res1, res2):
        """ Get a set of four resid/atom pairs corresponding to the psi angle of res1 and res2

        Parameters
        ----------
        res1 : int
            The first residue containing the N CA C atoms
        res2 : int
            The second residue containing the N atom

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        return (res1, 'N'), (res1, 'CA'), (res1, 'C'), (res2, 'N')

    @staticmethod
    def omega(res1, res2):
        """ Get a set of four resid/atom pairs corresponding to the omega angle of res1 and res2

        Parameters
        ----------
        res1 : int
            The first residue containing the CA C atoms
        res2 : int
            The second residue containing the N CA atoms

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        return (res1, 'CA'), (res1, 'C'), (res2, 'N'), (res2, 'CA')

    @staticmethod
    def chi1(res, resname):
        """ Get a set of four resid/atom pairs corresponding to the chi1 angle of a residue

        Parameters
        ----------
        res : int
            The resid of the residue

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        chi1std = ('N', 'CA', 'CB', 'CG')
        chi1 = {'ARG': chi1std, 'ASN': chi1std, 'ASP': chi1std, 'CYS': ('N', 'CA', 'CB', 'SG'), 'GLN': chi1std,
                'GLU': chi1std, 'HIS': chi1std, 'ILE': ('N', 'CA', 'CB', 'CG1'), 'LEU': chi1std, 'LYS': chi1std,
                'MET': chi1std, 'PHE': chi1std, 'PRO': chi1std, 'SER': ('N', 'CA', 'CB', 'OG'),
                'THR': ('N', 'CA', 'CB', 'OG1'), 'TRP': chi1std, 'TYR': chi1std, 'VAL': ('N', 'CA', 'CB', 'CG1')}
        if resname not in chi1:
            return None
        return (res, chi1[resname][0]), (res, chi1[resname][1]), (res, chi1[resname][2]), (res, chi1[resname][3])

    @staticmethod
    def chi2(res, resname):
        """ Get a set of four resid/atom pairs corresponding to the chi2 angle of a residue

        Parameters
        ----------
        res : int
            The resid of the residue

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        chi2std = ('CA', 'CB', 'CG', 'CD')
        chi2 = {'ARG': chi2std, 'ASN': ('CA', 'CB', 'CG', 'OD1'), 'ASP': ('CA', 'CB', 'CG', 'OD1'), 'GLN': chi2std,
                'GLU': chi2std, 'HIS': ('CA', 'CB', 'CG', 'ND1'), 'ILE': ('CA', 'CB', 'CG1', 'CD'),
                'LEU': ('CA', 'CB', 'CG', 'CD1'), 'LYS': chi2std,
                'MET': ('CA', 'CB', 'CG', 'SD'), 'PHE': ('CA', 'CB', 'CG', 'CD1'), 'PRO': chi2std,
                'TRP': ('CA', 'CB', 'CG', 'CD1'), 'TYR': ('CA', 'CB', 'CG', 'CD1')}
        if resname not in chi2:
            return None
        return (res, chi2[resname][0]), (res, chi2[resname][1]), (res, chi2[resname][2]), (res, chi2[resname][3])

    @staticmethod
    def chi3(res, resname):
        """ Get a set of four resid/atom pairs corresponding to the chi3 angle of a residue

        Parameters
        ----------
        res : int
            The resid of the residue

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        chi3 = {'ARG': ('CB', 'CG', 'CD', 'NE'), 'GLN': ('CB', 'CG', 'CD', 'OE1'), 'GLU': ('CB', 'CG', 'CD', 'OE1'),
                'LYS': ('CB', 'CG', 'CD', 'CE'), 'MET': ('CB', 'CG', 'SD', 'CE')}
        if resname not in chi3:
            return None
        return (res, chi3[resname][0]), (res, chi3[resname][1]), (res, chi3[resname][2]), (res, chi3[resname][3])

    @staticmethod
    def chi4(res, resname):
        """ Get a set of four resid/atom pairs corresponding to the chi4 angle of a residue

        Parameters
        ----------
        res : int
            The resid of the residue

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        chi4 = {'ARG': ('CG', 'CD', 'NE', 'CZ'), 'LYS': ('CG', 'CD', 'CE', 'NZ')}
        if resname not in chi4:
            return None
        return (res, chi4[resname][0]), (res, chi4[resname][1]), (res, chi4[resname][2]), (res, chi4[resname][3])

    @staticmethod
    def chi5(res, resname):
        """ Get a set of four resid/atom pairs corresponding to the chi5 angle of a residue

        Parameters
        ----------
        res : int
            The resid of the residue

        Returns
        -------
        quad : tuple
            A touple containing four resid/atom pairs
        """
        chi5 = {'ARG': ('CD', 'NE', 'CZ', 'NH1')}
        if resname not in chi5:
            return None
        return (res, chi5[resname][0]), (res, chi5[resname][1]), (res, chi5[resname][2]), (res, chi5[resname][3])


def _angle(A0, A1, A2, A3):
    '''
    http://en.wikipedia.org/wiki/Dihedral_angle#Methods_of_computation
    https://www.cs.unc.edu/cms/publications/dissertations/hoffman_doug.pdf
    http://www.cs.umb.edu/~nurith/cs612/hw4_2012.pdf

    Ua = (A2 - A1) x (A3 - A1)
    Ub = (B2 - B1) x (B3 - B1) = (A3 - A2) x (A4 - A2)
    angle = arccos((Ua * Ub) / (norm(Ua) * norm(Ub)))
    '''
    A10 = np.transpose(A1 - A0)
    A21 = np.transpose(A2 - A1)
    A32 = np.transpose(A3 - A2)
    Ua = np.cross(A10, A21)
    Ub = np.cross(A21, A32)
    #from IPython.core.debugger import Tracer
    #Tracer()()
    if np.any(np.sum(Ua == 0, 1) == 3) or np.any(np.sum(Ub == 0, 1) == 3):
        raise ZeroDivisionError('Two dihedral planes are exactly parallel or antiparallel. There is probably something broken in the simulation.')
    x = np.squeeze(_inner(Ua, Ub)) / (np.squeeze(np.linalg.norm(Ua, axis=1)) * np.squeeze(np.linalg.norm(Ub, axis=1)))

    # Fix machine precision errors (I hope at least that's what I'm doing...)
    if isinstance(x, np.ndarray):
        x[x > 1] = 1
        x[x < -1] = -1
    elif x > 1:
        x = 1
    elif x < -1:
        x = -1

    signV = np.squeeze(np.sign(_inner(Ua, A32)))
    signV[signV == 0] = 1  # Angle sign is 0. Maybe I should handle this better...

    return (np.arccos(x) * 180. / np.pi) * signV


def _inner(A, B):
    return np.sum(A * B, 1)

if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    metr = MetricDihedral(protsel='protein')
    data = metr.project(mol)

    calcdata = np.array([ 0.91631763,  0.40045224,  0.90890707,  0.41699872, -0.99956623,
                          0.02945084,  0.52407037, -0.85167496, -0.67766999, -0.73536616,
                          0.53415969, -0.8453836 , -0.66133656, -0.7500893 ,  0.55669439,
                         -0.83071738, -0.90348715, -0.42861517,  0.5950773 , -0.80366847,
                         -0.5837572 , -0.81192828,  0.71012313, -0.70407751, -0.95668505,
                         -0.29112493,  0.53835619, -0.84271739, -0.9231271 ,  0.38449493,
                         -0.12417973,  0.99225974, -0.93138983,  0.36402332,  0.37667118,
                         -0.92634703, -0.14376672, -0.98961161,  0.37357125, -0.92760149,
                         -0.93655808,  0.35051244,  0.64918191, -0.76063319, -0.93758286,
                         -0.34776195,  0.51787137, -0.8554585 , -0.96970912,  0.2442626 ])

    assert np.all(np.abs(calcdata - data[147, 500:550]) < 0.001), 'Diherdals calculation is broken'