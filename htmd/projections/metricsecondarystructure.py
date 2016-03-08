# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.metric import _OldMetric, _singleMolfile
from htmd.molecule.util import molRMSD
from htmd.molecule.util import sequenceID
from htmd.projections.projection import Projection
import numpy as np
from mdtraj.geometry._geometry import _dssp as dssp
import logging
logger = logging.getLogger(__name__)


class MetricSecondaryStructure(Projection):
    """ Calculates the secondary structure of the protein. DSSP implementation and documentation taken from MDtraj.

    Parameters
    ----------
    sel : str
        The atomselection for the protein
    simplified: bool
        Uses the simplified 3-letter code
    integer : bool
        Use integers instead of letter codes.

    Notes
    -----
    The DSSP assignment codes are:
       - 'H' : Alpha helix
       - 'B' : Residue in isolated beta-bridge
       - 'E' : Extended strand, participates in beta ladder
       - 'G' : 3-helix (3/10 helix)
       - 'I' : 5 helix (pi helix)
       - 'T' : hydrogen bonded turn
       - 'S' : bend
       - ' ' : Loops and irregular elements
    The simplified DSSP codes are:
       - 'H' : Helix. Either of the 'H', 'G', or 'I' codes. Integer code: 2
       - 'E' : Strand. Either of the 'E', or 'B' codes. Integer code: 1
       - 'C' : Coil. Either of the 'T', 'S' or ' ' codes. Integer code: 0

    A special 'NA' code will be assigned to each 'residue' in the topology which
    isn't actually a protein residue (does not contain atoms with the names
    'CA', 'N', 'C', 'O'), such as water molecules that are listed as 'residue's
    in the topology.

    Returns
    -------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated
    """
    def __init__(self, sel='protein', simplified=True, integer=True):
        self.sel = sel
        self.simplified = simplified
        self.integer = integer
        self._pc_ca_indices = None
        self._pc_nco_indices = None
        self._pc_proline_indices = None
        self._pc_chain_ids = None

    def _precalculate(self, mol):
        self._pc_ca_indices, self._pc_nco_indices, self._pc_proline_indices, self._pc_chain_ids = self._calcarrays(mol)

    def _calcarrays(self, mol):
        mol = mol.copy()
        mol.filter(self.sel, _logger=False)

        def get_or_minus1(idx):
            if len(idx) != 0:
                return idx[0]
            else:
                return -1

        ca_indices = mol.atomselect('name CA', indexes=True)
        chainids = mol.get('chain', sel=ca_indices)
        resnames = mol.get('resname', sel=ca_indices)

        uqchains = np.unique(chainids)
        chain_ids = np.zeros(len(chainids), dtype=np.int32)
        for i, uqc in enumerate(uqchains):
            idx = chainids == uqc
            chain_ids[idx] = i

        nco_indices = []
        uqresid = np.unique(mol.resid)
        for res in uqresid:
            nco_indices.append([get_or_minus1(mol.atomselect('resid {} and backbone and name "N.*"'.format(res), indexes=True)),
                                get_or_minus1(mol.atomselect('resid {} and backbone and name C'.format(res), indexes=True)),
                                get_or_minus1(mol.atomselect('resid {} and backbone and name "O.*"'.format(res), indexes=True))])
        nco_indices = np.array(nco_indices, np.int32)
        proline_indices = np.array(resnames == 'PRO', dtype=np.int32)
        ca_indices = np.array(ca_indices, dtype=np.int32)
        return ca_indices, nco_indices, proline_indices, chain_ids

    def _getSelections(self, mol):
        # If they have been pre-calculated return them.
        if self._pc_ca_indices is not None and self._pc_nco_indices is not None and \
                        self._pc_proline_indices is not None and self._pc_chain_ids is not None:
            return self._pc_ca_indices, self._pc_nco_indices, self._pc_proline_indices, self._pc_chain_ids
        else:  # Otherwise compute them.
            return self._calcarrays(mol)

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
        mol = args[0]

        ca_indices, nco_indices, proline_indices, chain_ids = self._getSelections(mol)
        #print(nco_indices)
        #print(ca_indices)
        #print(proline_indices)
        #print(chain_ids)

        xyz = np.swapaxes(np.swapaxes(np.atleast_3d(mol.coords), 1, 2), 0, 1)
        xyz = np.array(xyz.copy(), dtype=np.float32) / 10  # converting to nm

        data = dssp(xyz, nco_indices, ca_indices, proline_indices, chain_ids)

        if self.simplified:
            trans = str.maketrans('HGIEBTS ', 'HHHEECCC')
            data = data.translate(trans)

        data = np.fromiter(data, dtype=np.dtype('U2'))
        data = data.reshape(mol.numFrames, len(chain_ids))

        if self.integer:
            data[data == 'H'] = 2
            data[data == 'E'] = 1
            data[data == 'C'] = 0
            data = np.array(data, dtype=np.int32)
        #data[:, np.logical_not(protein_indices)] = 'NA'

        return data

    def getMapping(self, mol):
        idx = mol.atomselect('{} and name CA'.format(self.sel), indexes=True)
        return idx


class _MetricSecondaryStructureOld(_OldMetric):
    @staticmethod
    def project(sims, sel='protein', simple=True, skip=1, update=[]):
        """ Calculates the secondary structure of the protein

        Parameters
        ----------
        sims : numpy.ndarray of :class:`Sim <htmd.simlist.Sim>` objects or single :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            A simulation list generated by the :func:`simlist <htmd.simlist.simlist>` function, or a :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        sel : str
            The atomselection for the protein
        simple : bool
            Simple mode. Projects STRIDE assignments onto 3 integers. 0 for no SS, 1 for beta sheet and 2 for helices.

        Returns
        -------
        data : :class:`MetricData <htmd.metricdata.MetricData>` object
            Returns a :class:`MetricData <htmd.metricdata.MetricData>` object containing the metrics calculated
        """
        obj = _MetricSecondaryStructureOld(sims, sel, simple)
        return obj._metrify(sims, skip, update)

    def __init__(self, sims, sel, simple):
        self._pc_sel = None
        self._sel = sel
        self._simple = simple

        (single, molfile) = _singleMolfile(sims)
        if single:
            mol = Molecule(molfile)
            self._pc_sel = mol.atomselect(sel)

    def _processTraj(self, mol):
        sel = self._getSelections(mol)

        mol = mol.copy()
        mol.filter(sel, _logger=False)

        # Ugly conversions to feed it to sensitive C code
        uqresid = np.unique(mol.get('resid'))
        numres = len(uqresid)
        resatmcount = np.bincount(mol.resid)
        resatmcount = resatmcount[uqresid]
        residstr = [str(x) for x in mol.resid]

        atomids = np.ones(mol.numAtoms, dtype=np.int32) * -1
        for i in range(numres):
            resatms = mol.resid == uqresid[i]
            atomids[resatms] = np.arange(np.sum(resatms))

        sslist = stride(mol.coords.ravel(order='F'),
                        np.array([numres], dtype=np.int32),
                        resatmcount.astype(np.int32),
                        np.array(sequenceID(mol.resid) - 1, dtype=np.int32),
                        np.zeros(mol.numAtoms, dtype=np.int32),
                        atomids.astype(np.int32),
                        mol.chain[0],
                        residstr,
                        mol.resname.tolist(),
                        mol.name.tolist(),
                        1,
                        mol.numFrames,
                        mol.numAtoms)

        ss = np.zeros(len(sslist), dtype=object)
        for i in range(len(sslist)):
            ss[i] = np.array(list(sslist[i]))

        if self._simple:
            return _ssmap(ss)
        else:
            return ss

    '''def _processTraj(self, mol):
        from ctypes import cdll, POINTER, c_float, c_int
        from htmd.home import home
        libdir = path.join(home(), 'lib')
        stride = cdll.LoadLibrary(libdir + "/pystride.so")

        sel = self._getSelections(mol)
        mol = mol.copy()
        mol.filter(sel)

        c_coords = mol.coords.ravel(order='F').copy().ctypes.data_as(POINTER(c_float))  # TODO: Make pystride accept C order coordinates

        sslist = stride.stride(c_coords,  # TODO: Fill in the rest here
                               c_int(1),
                               c_int(mol.numFrames),
                               c_int(mol.numAtoms))

        ss = np.zeros(len(sslist), dtype=object)
        for i in range(len(sslist)):
            ss[i] = np.array(list(sslist[i]))

        if self._simple:
            return _ssmap(ss)
        else:
            return ss'''

    def _getSelections(self, mol):
        if self._pc_sel is not None:
            sel = self._pc_sel
        else:
            sel = mol.atomselect(self._sel)
        if np.sum(sel) == 0:
            raise NameError('Protein selection resulted in 0 atoms.')
        return sel

    def _getMapping(self, mol):
        pass


def _ssmap(sschar):
    # Map secondary structures to integers
    #     G = 3-turn helix (310 helix). Min length 3 residues.
    #     H = 4-turn helix (α helix). Min length 4 residues.
    #     I = 5-turn helix (π helix). Min length 5 residues.
    #     T = hydrogen bonded turn (3, 4 or 5 turn)
    #     E = extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
    #     B = residue in isolated β-bridge (single pair β-sheet hydrogen bond formation)
    #     S = bend (the only non-hydrogen-bond based assignment).
    #     C = coil (residues which are not in any of the above conformations).
    ssnum = np.zeros((len(sschar), len(sschar[0])), dtype=np.uint8)

    for i in range(len(sschar)):
        sschar[i][sschar[i] == 'G'] = 'H'
        sschar[i][sschar[i] == 'I'] = 'H'
        sschar[i][sschar[i] == 'B'] = 'E'
        sschar[i][sschar[i] == 'b'] = 'E'
        sschar[i][sschar[i] == 'S'] = 'C'
        sschar[i][sschar[i] == 'T'] = 'C'
        ssnum[i, np.where(sschar[i] == 'C')[0]] = 0
        ssnum[i, np.where(sschar[i] == 'E')[0]] = 1
        ssnum[i, np.where(sschar[i] == 'H')[0]] = 2
    return ssnum

if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    from os import path
    import numpy as np
    mol = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
    #mol.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
    data = MetricSecondaryStructure.project(mol)
    x = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2,
        2, 2, 2, 2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1,
        1, 1, 0, 2, 2, 2, 0, 1, 1, 2, 2, 2, 0, 1, 1, 0, 0, 2, 2, 2, 1, 1,
        1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1,
        1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 2, 2, 2,
        2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0]], dtype=np.uint8)
    assert np.array_equal(data, x), 'MetricSecondaryStructure assertion failed'

