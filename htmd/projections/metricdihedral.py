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


class MetricDihedral(Projection):
    """ Calculates a set of dihedral angles from trajectories

    Parameters
    ----------
    dih : list of str
        You can provide your own list of atomselect strings consisting of 4 atoms each for dihedral angles
    sincos : bool, optional
        To return the dihedral angles as their sine and cosine components. Makes them periodic.
    protsel : str, optional
        Atomselection for the whole protein.
    """
    def __init__(self, dih=None, sincos=True, protsel='protein'):
        self._protsel = protsel
        self._sincos = sincos
        self._dih = dih  # TODO: Calculate the dihedral
        self._pc_dih = None

    def _precalculate(self, mol):
        self._pc_dih = self._dihedralPrecalc(mol, mol.atomselect(self._protsel))

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
            data = _MetricDihedralOld.project(self, *args, **kwargs)
            logger.warning('Static use of the project method will be deprecated in the next version of HTMD. '
                           'Please change your projection code according to the tutorials on www.htmd.org')
            return data
        # ------------------ CUT HERE -------------------
        mol = args[0]
        dih = self._getSelections(mol)
        return self._calcDihedrals(mol, dih, sincos=self._sincos)

    def _getSelections(self, mol):
        if self._pc_dih is not None:
            return self._pc_dih
        else:
            return self._dihedralPrecalc(mol, mol.atomselect(self._protsel))

    def getMapping(self, mol):
        dih = self._getSelections(mol)
        map = []
        for i, d in enumerate(dih):
            resids = mol.resid[d]
            uqresids = np.unique(resids)
            mapstr = ''
            for u in uqresids:
                mapstr += 'resid {} atoms: '.format(u)
                residatoms = mol.name[(mol.resid == u) & d]
                for a in residatoms:
                    mapstr += '{} '.format(a)
            if self._sincos:
                map.append('Sine of angle of ' + mapstr)
                map.append('Cosine of angle of ' + mapstr)
            else:
                map.append('Angle of ' + mapstr)
        return np.array(map)

    def _dihedralPrecalc(self, mol, protsel):
        #logger.info('Precalculating phi and psi angle atom selections')

        # Phi angle: C'-N-Cα-C'
        # Psi angle: N-Cα-C'-N
        protmask = mol.atomselect('protein')
        protmask = protmask & protsel
        residues = np.unique(mol.get('resid', sel=protmask))

        if self._dih is None:
            selstr = []
            # Only psi angle for first residue
            selstr.append('backbone and (resid ' + str(residues[0]) + ' and name N CA C) or (resid ' + str(residues[1]) + ' and name N)')
            for r in range(1, len(residues)-1):
                selstr.append('backbone and (resid ' + str(residues[r]) + ' and name N CA C) or (resid ' + str(residues[r+1]) + ' and name N)')  #psi
                selstr.append('backbone and (resid ' + str(residues[r-1]) + ' and name C) or (resid ' + str(residues[r]) + ' and name N CA C)')  #phi
            # Only phi angle for the last residue
            selstr.append('backbone and (resid ' + str(residues[len(residues)-2]) + ' and name C) or (resid ' + str(residues[len(residues)-1]) + ' and name N CA C)')
        else:
            if isinstance(self._dih, str):
                self._dih = [self._dih]
            selstr = self._dih

        dih = []
        for s in selstr:
            dih.append(mol.atomselect(s))
        #logger.info('Finished precalculating phi and psi.')
        return dih

    def _calcDihedrals(self, mol, dih, sincos=True):
        metric = np.zeros((np.size(mol.coords, 2), len(dih)))

        for i in range(len(dih)):
            metric[:, i] = self._angle(mol.coords[dih[i], :, :])

        if sincos:
            sc_metric = np.zeros((np.size(metric, 0), np.size(metric, 1) * 2))
            sc_metric[:, 0::2] = np.sin(metric * np.pi / 180.)
            sc_metric[:, 1::2] = np.cos(metric * np.pi / 180.)
            metric = sc_metric
        return metric

    def _angle(self, A):
        '''
        http://en.wikipedia.org/wiki/Dihedral_angle#Methods_of_computation
        https://www.cs.unc.edu/cms/publications/dissertations/hoffman_doug.pdf
        http://www.cs.umb.edu/~nurith/cs612/hw4_2012.pdf

        Ua = (A2 - A1) x (A3 - A1)
        Ub = (B2 - B1) x (B3 - B1) = (A3 - A2) x (A4 - A2)
        angle = arccos((Ua * Ub) / (norm(Ua) * norm(Ub)))
        '''
        A10 = np.transpose(A[1, :, :] - A[0, :, :])
        A21 = np.transpose(A[2, :, :] - A[1, :, :])
        A32 = np.transpose(A[3, :, :] - A[2, :, :])
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


class _MetricDihedralOld(_OldMetric):
    @staticmethod
    def project(sims, dih=None, sincos=True, protsel='protein', skip=1, update=[]):
        """ Calculates a set of dihedral angles from trajectories

        Parameters
        ----------
        sims : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects or single :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function, or a :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        dih : list of str
            You can provide your own list of atomselect strings consisting of 4 atoms each for dihedral angles
        sincos : bool, optional
            To return the dihedral angles as their sine and cosine components. Makes them periodic.
        protsel : str, optional
            Atomselection for the whole protein.
        skip : int
            Skip every `skip` frames of the trajectories
        update :
            Not functional yet

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated
        """
        obj = _MetricDihedralOld(sims, protsel, dih, sincos)
        return obj._metrify(sims, skip, update)

    def __init__(self, sims, protsel, dih=None, sincos=True):
        self._protsel = protsel
        self._sincos = sincos
        self._dih = dih  # TODO: Calculate the dihedral
        self._pc_dih = None

        (single, molfile) = _singleMolfile(sims)
        if single:
            mol = Molecule(molfile)
            self._pc_dih = self._dihedralPrecalc(mol, mol.atomselect(protsel))

    def _processTraj(self, mol):
        dih = self._getSelections(mol)
        return self._calcDihedrals(mol, dih, sincos=self._sincos)

    def _getSelections(self, mol):
        if self._pc_dih is not None:
            return self._pc_dih
        else:
            return self._dihedralPrecalc(mol, mol.atomselect(self._protsel))

    def _getMapping(self, mol):
        dih = self._getSelections(mol)
        map = []
        for i, d in enumerate(dih):
            resids = mol.resid[d]
            uqresids = np.unique(resids)
            mapstr = ''
            for u in uqresids:
                mapstr += 'resid {} atoms: '.format(u)
                residatoms = mol.name[(mol.resid == u) & d]
                for a in residatoms:
                    mapstr += '{} '.format(a)
            if self._sincos:
                map.append('Sine of angle of ' + mapstr)
                map.append('Cosine of angle of ' + mapstr)
            else:
                map.append('Angle of ' + mapstr)
        return np.array(map)

    def _dihedralPrecalc(self, mol, protsel):
        logger.info('Precalculating phi and psi angle atom selections')

        # Phi angle: C'-N-Cα-C'
        # Psi angle: N-Cα-C'-N
        protmask = mol.atomselect('protein')
        protmask = protmask & protsel
        residues = np.unique(mol.get('resid', sel=protmask))

        if self._dih is None:
            selstr = []
            # Only psi angle for first residue
            selstr.append('backbone and (resid ' + str(residues[0]) + ' and name N CA C) or (resid ' + str(residues[1]) + ' and name N)')
            for r in range(1, len(residues)-1):
                selstr.append('backbone and (resid ' + str(residues[r]) + ' and name N CA C) or (resid ' + str(residues[r+1]) + ' and name N)')  #psi
                selstr.append('backbone and (resid ' + str(residues[r-1]) + ' and name C) or (resid ' + str(residues[r]) + ' and name N CA C)')  #phi
            # Only phi angle for the last residue
            selstr.append('backbone and (resid ' + str(residues[len(residues)-2]) + ' and name C) or (resid ' + str(residues[len(residues)-1]) + ' and name N CA C)')
        else:
            selstr = self._dih

        dih = []
        for s in selstr:
            dih.append(mol.atomselect(s))
        logger.info('Finished precalculating phi and psi.')
        return dih

    def _calcDihedrals(self, mol, dih, sincos=True):
        metric = np.zeros((np.size(mol.coords, 2), len(dih)))

        for i in range(len(dih)):
            metric[:, i] = self._angle(mol.coords[dih[i], :, :])

        if sincos:
            sc_metric = np.zeros((np.size(metric, 0), np.size(metric, 1) * 2))
            sc_metric[:, 0::2] = np.sin(metric * np.pi / 180.)
            sc_metric[:, 1::2] = np.cos(metric * np.pi / 180.)
            metric = sc_metric
        return metric

    def _angle(self, A):
        '''
        http://en.wikipedia.org/wiki/Dihedral_angle#Methods_of_computation
        https://www.cs.unc.edu/cms/publications/dissertations/hoffman_doug.pdf
        http://www.cs.umb.edu/~nurith/cs612/hw4_2012.pdf

        Ua = (A2 - A1) x (A3 - A1)
        Ub = (B2 - B1) x (B3 - B1) = (A3 - A2) x (A4 - A2)
        angle = arccos((Ua * Ub) / (norm(Ua) * norm(Ub)))
        '''
        A10 = np.transpose(A[1, :, :] - A[0, :, :])
        A21 = np.transpose(A[2, :, :] - A[1, :, :])
        A32 = np.transpose(A[3, :, :] - A[2, :, :])
        Ua = np.cross(A10, A21)
        Ub = np.cross(A21, A32)
        #from IPython.core.debugger import Tracer
        #Tracer()()
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