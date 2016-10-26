from htmd.projections.projection import Projection
from htmd.molecule.util import sequenceID
import numpy as np
import logging
logger = logging.getLogger(__name__)


class MetricSasa(Projection):
    """ Calculate solvent accessible surface area of a molecule.

    Implementation and documentation taken from MDtraj shrake_rupley code.

    Parameters
    ----------
    sel : str
        Atomselection for atoms or residues for which to calculate the SASA
    probeRadius : float
        The radius of the probe, in Angstrom.
    numSpherePoints : int
        The number of points representing the surface of each atom, higher values lead to more accuracy.
    mode : str
        In mode == 'atom', the extracted areas are resolved per-atom. In mode == 'residue', this is consolidated down
        to the per-residue SASA by summing over the atoms in each residue.

    Returns
    -------
    metr : MetricSasa object
    """
    def __init__(self, sel='protein', probeRadius=1.4, numSpherePoints=960, mode='atom'):
        self._probeRadius = probeRadius
        self._numSpherePoints = numSpherePoints
        self._mode = mode
        self._sel = sel
        self._pc_radii = None
        self._pc_atom_mapping = None
        self._pc_sel = None

    def _precalculate(self, mol):
        self._pc_radii, self._pc_atom_mapping, self._pc_sel = self._calcRadiiMapping(mol)

    def _calcRadiiMapping(self, mol):
        sel = mol.atomselect(self._sel)

        _ATOMIC_RADII = {'C': 1.5, 'F': 1.2, 'H': 0.4, 'N': 1.10, 'O': 1.05, 'S': 1.6, 'P': 1.6}
        elements = [n[0] for n in mol.name[sel]]
        atom_radii = np.vectorize(_ATOMIC_RADII.__getitem__)(elements)
        radii = np.array(atom_radii, np.float32) + self._probeRadius

        if self._mode == 'atom':
            atom_mapping = np.arange(np.sum(sel), dtype=np.int32)
        elif self._mode == 'residue':
            atom_mapping = sequenceID((mol.resid[sel], mol.chain[sel], mol.segid[sel])).astype(np.int32)
        else:
            raise ValueError('mode must be one of "residue", "atom". "{}" supplied'.format(self._mode))

        return radii, atom_mapping, sel

    def _getRadiiMapping(self, mol):
        if self._pc_radii is not None and self._pc_atom_mapping is not None and self._pc_sel is not None:
            return self._pc_radii, self._pc_atom_mapping, self._pc_sel
        else:
            return self._calcRadiiMapping(mol)

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
        radii, atom_mapping, sel = self._getRadiiMapping(mol)

        xyz = np.swapaxes(np.swapaxes(np.atleast_3d(mol.coords[sel, :, :]), 1, 2), 0, 1)
        xyz = np.array(xyz.copy(), dtype=np.float32) / 10  # converting to nm

        from mdtraj.geometry._geometry import _sasa as sasa
        out = np.zeros((mol.numFrames, atom_mapping.max() + 1), dtype=np.float32)
        sasa(xyz, radii / 10, int(self._numSpherePoints), atom_mapping, out)  # Divide radii by 10 for nm
        return out

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
        _, atom_mapping, atomsel = self._getRadiiMapping(mol)
        if self._mode == 'atom':
            atomidx = np.where(atomsel)[0]
        elif self._mode == 'residue':
            _, firstidx = np.unique(atom_mapping, return_index=True)
            atomidx = np.where(atomsel)[0][firstidx]
        else:
            raise ValueError('mode must be one of "residue", "atom". "{}" supplied'.format(self._mode))

        from pandas import DataFrame
        types = []
        indexes = []
        description = []
        for i in atomidx:
            types += ['SASA']
            indexes += [i]
            description += ['SASA of {} {} {}'.format(mol.resname[i], mol.resid[i], mol.name[i])]
        return DataFrame({'type': types, 'indexes': indexes, 'description': description})


if __name__ == '__name__':
    from htmd.molecule.molecule import Molecule
    from htmd.projections.metricsasa import MetricSasa
    from htmd.home import home
    from os import path
    import numpy as np
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    mol.dropFrames(keep=[0, 1])  # Keep only two frames cause it's super slow

    metr = MetricSasa(mode='atom')
    sasaA = metr.project(mol)
    metr = MetricSasa(mode='residue')
    sasaR = metr.project(mol)

    sasaA_ref = np.load(path.join(home(), 'data', 'test-metrics', 'metricsasa', 'sasaA.npy'))
    sasaR_ref = np.load(path.join(home(), 'data', 'test-metrics', 'metricsasa', 'sasaR.npy'))

    assert np.array_equal(sasaA, sasaA_ref)
    assert np.array_equal(sasaR, sasaR_ref)
