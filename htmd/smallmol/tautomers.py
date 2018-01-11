# The molvs.TautomerCanonicalizer was refactored by:
#           * adding conjugate length score
#           * tautomers filtering based on the scores computed and a given threshold (default: 2)

# Two class were implemented:
#           * MolFragmenter: a class that handles the molecule fragmentations into conjugates fragments and then
#             compute the longest path for electron delocalization
#           * Fragment: a class that is simply an object where to store information of each fragment generated

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import BondType
from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.molvstautomer import TAUTOMER_TRANSFORMS, TAUTOMER_SCORES, MAX_TAUTOMERS, pairwise, memoized_property
import matplotlib.pyplot as plt
from math import ceil
import numpy as np
import logging

logger = logging.getLogger(__name__)

AROMATIC = rdchem.BondType.AROMATIC
SINGLE = rdchem.BondType.SINGLE
DOUBLE = rdchem.BondType.DOUBLE


class MolFragmenter:
    """
    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol - the rdkit molecule object

    Attributes
    ----------

    mol: rdkit.Chem.rdchem.Mol - the rdkit molecule object
    fragments: list - list of Fragment objects
    paths: list - list of list containing fragment objects part of the path


    """
    # conjugate scores. 1pt for each double bond
    SCORE_dict = {'aromatic': 3, 'double': 1}

    def __init__(self, mol):
        if isinstance(mol, SmallMol):
            mol = mol._mol

        self.mol = mol
        self.fragments = []
        self.paths = []

    def fragment(self):
        """
        The method creates and stores Fragment objects for each bond in the molecule and based on the atom type
        """

        # storing all double bonds idences by using rdkit.GetBondType
        all_double = [b.GetIdx() for b in self.mol.GetBonds() if b.GetBondType() == DOUBLE]
        # storing all aromatic bonds idences by using rdkit.GetBondType
        all_aromatics = []
        for ring in self.mol.GetRingInfo().BondRings():
            if sum([1 for bidx in ring if self.mol.GetBondWithIdx(bidx).GetBondType() == AROMATIC]) == len(ring):
                all_aromatics.append(list(ring))
        # storing all others bonds type
        try:  # for aromatics
            all_others = [b.GetIdx() for b in self.mol.GetBonds() if
                          b.GetIdx() not in np.concatenate(all_aromatics) and b.GetIdx() not in all_double]
        except:
            all_others = [b.GetIdx() for b in self.mol.GetBonds() if b.GetIdx() not in all_double]

        # generating Fragment object for each bond in the molecules based on their types
        nameint = 1
        for f in all_double:
            fr = Fragment(str(nameint), 'double')
            bond = self.mol.GetBondWithIdx(f)
            atoms = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            fr.atoms = atoms
            self.fragments.append(fr)
            nameint += 1
        for f in all_aromatics:
            fr = Fragment(str(nameint), 'aromatic')
            atoms = []
            for bidx in f:
                bond = self.mol.GetBondWithIdx(bidx)
                atoms.extend([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            fr.atoms = np.unique(atoms).tolist()
            self.fragments.append(fr)
            nameint += 1
        for f in all_others:
            fr = Fragment(str(nameint), 'single')
            bond = self.mol.GetBondWithIdx(f)
            atoms = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            fr.atoms = atoms
            self.fragments.append(fr)
            nameint += 1

    def get_fragment(self, attr, value):
        """
        The function returns the Fragment objects based on the attribute and value provided

        Parameters
        ----------
        attr: Fragment.attributename -  the Fragment object attribute  (selector)
        value: Fragment.attributevalue - the Fragment object attribute value to match

        Returns
        -------
        fragments: list -  The list of Fragment objects matching the parameters provided

        """
        fragments = [f for f in self.fragments if getattr(f, attr) == value]
        return fragments

    def search_fragment_by_atom(self, atom_idx, fragmenttype=None):
        """
        The function returns the Fragment objects based on the atomIdx  and fragmenttype provided (default: all type)

        Parameters
        ----------
        atom_idx: int - The atom idx
        fragmenttype:  Fragment.ftype value - The fragment type to match

        Returns
        -------
        fragment: list - The list of Fragment objects matching the parameters provided

        """

        if fragmenttype is not None:
            fragments = self.get_fragment('ftype', fragmenttype)
        else:
            fragments = self.fragments

        fragment = [fr for fr in fragments if atom_idx in fr.atoms]

        return fragment

    def get_num_fragments(self):
        """
        Returns: int
                 Total number of fragments.
        """

        return len(self.fragments)

    def get_graphs(self):
        """
        The function inspect for each conjugate Fragment (double aromatic) all the possible pattern
        that will set as connections
        """

        for fr in self.get_fragment('ftype', 'double'):
            graphs = []
            for frsingle in self.get_fragment('ftype', 'single'):
                graphs.extend([frsingle for a in fr.atoms if a in frsingle.atoms])
            fr.connections = [fr.name for fr in graphs]
            fr.graphs = graphs

        for fr in self.get_fragment('ftype', 'aromatic'):
            graphs = []
            for frsingle in self.get_fragment('ftype', 'single'):
                graphs.extend([frsingle for a in fr.atoms if a in frsingle.atoms])
            for fraromatic in self.get_fragment('ftype', 'aromatic'):
                if fr != fraromatic:
                    if len([a for a in fr.atoms if a in fraromatic.atoms]) == 2:
                        graphs.append(fraromatic)

            fr.connections = [fr.name for fr in graphs]
            fr.graphs = graphs

    def compute_paths(self):
        """
        The function inpect and detect all the possible conjugate paths
        """

        fragments = self.get_fragment('ftype', 'aromatic') + self.get_fragment('ftype', 'double')
        for fr in fragments:
            path = [fr]
            for fr2 in fragments:
                if fr == fr2:
                    continue
            path.extend([fr2 for f in path if f.ftype == 'aromatic' and
                         fr2.ftype == 'aromatic' and
                         f.name in fr2.connections])
            path.extend([fr2 for f in path for conn in f.connections if conn in fr2.connections])
        self.paths.append(list(set(path)))


    def get_longest_path(self):
        """
        The function return the longest path

        Returns
        -------
        highest_path: list - The list of fragmaent that represents the longest path

        """
        highest_num_fragments = max([len(p) for p in self.paths])

        longest_path = [p for p in self.paths if len(p) == highest_num_fragments][0]

        return longest_path

    def get_depiction_atoms(self, path):
        """
        The function returns all the atom indeces that are part of a conjugate path

        Paramters
        ---------
        path: list - List of Fragment objects that are part of the same conjugate path

        Returns
        -------
        atoms: list - List of atoms indeces part of a cpnjugate path
        """

        atoms = [a for fr in path for a in fr.atoms]
        return atoms

    def get_conjugate_length(self, path):
        """
        The function returns the length of the conjugate path

        Paramters
        ---------
        path: list - List of Fragment objects that are part of the same conjugate path

        Returns
        -------
        length: int - The lenght of the conjugate path
        """
        length = 0
        for fr in path:
            length += self.SCORE_dict[fr.ftype]
        return length


class Fragment:
    """
    The Fragment class store information about the molecule fragment

    Parameters
    ----------
    name: int - The id of the fragment
    ftype: str - The fragment type (choices: single, aromatic, double)


    Attributes
    ----------
    name: int - The id of the fragment
    ftype: str - The fragment type (choices: single, aromatic, double)
    atoms: list - List of the atom indeces
    graphs: list - List of the Fragment.name of the the Fragment object connected with this Fragment
    connection: list - List of the Fragment.name of the the Fragment object connected with
                       this Fragment (Only double and aromatic.)

    """

    def __init__(self, name, ftype):
        self.name = name
        self.ftype = ftype
        self.atoms = None
        self.graphs = None  # Fragment obj
        self.connections = []

    def get_connection(self):
        """The function returns the connections
        Returns
        -------
        connections: list - List of the Fragment.name of the the Fragment object connected with this Fragment
                            (Only double and aromatic.)
        """
        return self.connections


class TautomerCanonicalizer:
    """ The molvs class refactores"""

    def __init__(self, transforms=TAUTOMER_TRANSFORMS, scores=TAUTOMER_SCORES, max_tautomers=MAX_TAUTOMERS):
        """

        :param transforms: A list of TautomerTransforms to use to enumerate tautomers.
        :param scores: A list of TautomerScores to use to choose the canonical tautomer.
        :param max_tautomers: The maximum number of tautomers to enumerate, a limit to prevent combinatorial explosion.
        """
        self.transforms = transforms
        self.scores = scores
        self.max_tautomers = max_tautomers

    def __call__(self, mol):
        """Calling a TautomerCanonicalizer instance like a function is the same as calling its canonicalize(mol) method."""
        return self.canonicalize(mol)

    def canonicalize(self, mol):
        """Return a canonical tautomer by enumerating and scoring all possible tautomers.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The canonical tautomer.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        # TODO: Overload the mol parameter to pass a list of pre-enumerated tautomers
        tautomers = self._enumerate_tautomers(mol)
        if len(tautomers) == 1:
            return tautomers[0]
        # Calculate score for each tautomer
        highest = None
        for t in tautomers:
            smiles = Chem.MolToSmiles(t, isomericSmiles=True)
            logger.debug('Tautomer: %s', smiles)
            score = 0
            # Add aromatic ring scores
            ssr = Chem.GetSymmSSSR(t)
            for ring in ssr:
                btypes = {t.GetBondBetweenAtoms(*pair).GetBondType() for pair in pairwise(ring)}
                elements = {t.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring}
                if btypes == {BondType.AROMATIC}:
                    logger.debug('Score +100 (aromatic ring)')
                    score += 100
                    if elements == {6}:
                        logger.debug('Score +150 (carbocyclic aromatic ring)')
                        score += 150
            # Add SMARTS scores
            for tscore in self.scores:
                for match in t.GetSubstructMatches(tscore.smarts):
                    logger.debug('Score %+d (%s)', tscore.score, tscore.name)
                    score += tscore.score
            # Add (P,S,Se,Te)-H scores
            for atom in t.GetAtoms():
                if atom.GetAtomicNum() in {15, 16, 34, 52}:
                    hs = atom.GetTotalNumHs()
                    if hs:
                        logger.debug('Score %+d (%s-H bonds)', -hs, atom.GetSymbol())
                        score -= hs
            # Set as highest if score higher or if score equal and smiles comes first alphabetically
            if not highest or highest['score'] < score or (highest['score'] == score and smiles < highest['smiles']):
                logger.debug('New highest tautomer: %s (%s)', smiles, score)
                highest = {'smiles': smiles, 'tautomer': t, 'score': score}
        return highest['tautomer']

    @memoized_property
    def _enumerate_tautomers(self):
        from htmd.smallmol.molvstautomer import TautomerEnumerator
        return TautomerEnumerator(self.transforms, self.max_tautomers)

    def depict_tautomers(self, tautomers, scores, atoms, details=None):
        from rdkit.Chem import Draw
        dep_for_row = len(tautomers) if len(tautomers) < 3 else 3
        legends_tmp = [s for s in scores]
        legends = []
        if details is not None:
            for n, l in enumerate(legends_tmp):
                l = 'TotalScore: {}\n {}'.format(l, "\n".join(['{}: {}'.format(k, v) for k, v in details[n].items()]))
                legends.append(l)

        fig = plt.figure(figsize=(10, 10), dpi=100, tight_layout=True)
        for n, (t, l, a) in enumerate(zip(tautomers, legends, atoms)):
            depcition = Draw.MolToImage(t, subImgSize=(300, 300), highlightAtoms=a)
            ax = fig.add_subplot(ceil(len(tautomers) / dep_for_row), dep_for_row, n + 1)
            ax.set_xlabel(l, fontsize=5)
            ax.imshow(depcition, aspect='equal')
        plt.show()

    def get_conjugate(self, tautomer):
        """The functions called the MolFragmeter and returns the conjugate length  and the atoms member of the path

        Parameters
        ----------
        tatutomer: rdkit.Chem.Molecule - A tautomer object

        Returns
        -------
        n_conj: int - the lenght of the conjugate path
        atoms: list - The list of the atoms part of the path

        """

        mf = MolFragmenter(tautomer)
        mf.fragment()
        mf.get_graphs()
        mf.compute_paths()
        path = mf.get_longest_path()
        n_conj = mf.get_conjugate_length(path)
        atoms = mf.get_depiction_atoms(path)
        return n_conj, atoms

    def compute_score(self, mol, returndetails=False, log=False):
        """
        Return a canonical tautomer by enumerating and scoring all possible tautomers.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The canonical tautomer.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        t_scores = []
        scores_detail = []
        t_depict = []

        if isinstance(mol, SmallMol):
            mol = mol._mol

        # TODO: Overload the mol parameter to pass a list of pre-enumerated tautomers
        tautomers = self._enumerate_tautomers(mol)

        for t in tautomers:
            tmp_score_details = {'ArRing': 0,
                                 'CarbArRing': 0,
                                 'MatchFeature': [0, []],
                                 'Penalty': [0, []],
                                 'Conjugate': 0}
            smiles = Chem.MolToSmiles(t, isomericSmiles=True)
            if log:
                print('Tautomer: %s', smiles)
            score = 0
            # Add aromatic ring scores
            ssr = Chem.GetSymmSSSR(t)
            for ring in ssr:
                btypes = {t.GetBondBetweenAtoms(*pair).GetBondType() for pair in pairwise(ring)}
                elements = {t.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring}
                if btypes == {BondType.AROMATIC}:
                    if log:
                        print('Score +100 (aromatic ring)')
                    score += 100
                    tmp_score_details['ArRing'] += 1
                    if elements == {6}:
                        if log:
                            print('Score +150 (carbocyclic aromatic ring)')
                        score += 150
                        tmp_score_details['CarbArRing'] += 1

            # Add SMARTS scores, Chem.MolToSmiles(t))
            for tscore in self.scores:

                for match in t.GetSubstructMatches(tscore.smarts):
                    if log:
                        print('Score %+d (%s)' % (tscore.score, tscore.name))
                    score += tscore.score
                    tmp_score_details['MatchFeature'][0] += 1
                    tmp_score_details['MatchFeature'][1].append(tscore.name)

            # Add (P,S,Se,Te)-H scores
            for atom in t.GetAtoms():
                if atom.GetAtomicNum() in {15, 16, 34, 52}:
                    hs = atom.GetTotalNumHs()
                    if hs:
                        if log:
                            print('Score %+d (%s-H bonds)' % (-hs, atom.GetSymbol()))
                        score -= hs
                        tmp_score_details['Penalty'][0] += 1
                        tmp_score_details['Penalty'][1].append(atom.GetSymbol())

            # compute the conjuggate system
            n_conjugate, depictionatoms = self.get_conjugate(t)
            t_depict.append(depictionatoms)

            tmp_score_details['Conjugate'] = n_conjugate
            score += n_conjugate * 2

            scores_detail.append(tmp_score_details)
            t_scores.append(score)
        if returndetails:
            return [SmallMol(tautomer) for tautomer in tautomers], t_scores, t_depict, scores_detail
        return tautomers, t_scores

    @staticmethod
    def filter_tautomers(tautomers, scores, threshold=2):
        """ The function returns the tautomers as rdkit molecule objects based on the scores and the threshold

        Parameters
        ----------
        tautomers: list - List of rdkit.Chem.Molecule of the tautomers identified
        scores: list - List of the scores for each tatutomer
        threshold: int - The threshold value to be used as difference from the highest one

        Returns
        -------
        t_filtered: list - List of rdkit.Chem.Molecule of the tautomers filtered
        """

        tautomers = [mol._mol if isinstance(mol, SmallMol) else mol for mol in tautomers]

        tautomers_sorted = [x for _, x in sorted(zip(scores, tautomers), key=lambda pair: pair[0], reverse=True)]
        scores.sort(reverse=True)
        t_filterd = [t for t, s in zip(tautomers_sorted, scores) if s >= max(scores) - threshold]
        return [SmallMol(rdmol) for rdmol in t_filterd]
