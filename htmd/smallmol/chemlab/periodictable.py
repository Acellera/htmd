from rdkit.Chem.rdchem import HybridizationType, BondType, ChiralType
from rdkit.Chem import GetPeriodicTable
from rdkit import rdBase
import numpy as np


_heavy_atoms = ['C', 'N', 'O', 'F',
                'Si', 'P', 'S', 'Cl', 'Br', 'I']

_hetero_atoms = ['N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I']

_halogen_atoms = ['F', 'Cl', 'Br', 'I']


_hybridizations_IdxToType = HybridizationType.values
_hybridizations_StringToType = {'S': HybridizationType.S,
                                'SP': HybridizationType.SP,
                                'SP2': HybridizationType.SP2,
                                'SP3': HybridizationType.SP3}

rdBase.DisableLog('rdApp.error')
atom_mapping = {"Hydrophobe": 0,
                "LumpedHydrophobe": 0,
                "Aromatic": 1,
                "Acceptor": 2,
                "Donor": 3,
                "PosIonizable": 4,
                "NegIonizable": 5}

_chiral_type = {ChiralType.CHI_TETRAHEDRAL_CW: 'clockwise',
                ChiralType.CHI_TETRAHEDRAL_CCW: 'anticlockwise'}
_chiral_type_Dict = ChiralType.values

_bondtypes_IdxToType = BondType.values
_bondtypes_StringToType = {'SINGLE': BondType.SINGLE,
                           'DOUBLE': BondType.DOUBLE,
                           'TRIPLE': BondType.TRIPLE}

_atoms = {'H': ['H', 0, 0, '', 1, 1],
          'C': ['C', 0, 0, '', 4, 1]
          }


class PeriodicTable:
    """
    Class that manages all the chemical information about atoms based on the atom. It exposes the rdkit PeriodicTable
    functionalities.
    """

    def __init__(self):
        self.PeriodicaTable = GetPeriodicTable()

    def _evaluateMissingAtoms(self, valence, btypes):
        """
        Returns the difference bewtween the current valence computed and the expected one.

        Parameters
        ----------
        valence: int
            The valence of the atom
        btypes: list
            A list of the bondtypes connected to the atom

        Returns
        -------
        missatoms: int
            The number of atoms that miss for completing the atom valence

        """
        normbondvalence = sum([1.5 if int(b) == 12 else int(b) for b in btypes])

        missatoms = int(valence - normbondvalence)

        return missatoms

    def getBondRadius(self, element):
        """
        Returns the atom radius when bonded

        Parameters
        ----------
        element: str
            The element of the atom

        Returns
        -------
        bradius: float
            The radius of the bond when bonded

        """

        pT = self.PeriodicaTable
        bradius = pT.GetRb0(element)

        return bradius

    def getMissingValence(self, smallmol, atomidx, formalcharge=0):
        """
        Returns the number of missing atom to complete the atom valence based on the types of bonds and formalcharge

        Parameters
        ----------
        smallmol: htmd.smallmol.smallmol.SmallMol
            The SmallMol object
        atomidx: int
            The index of the atom
        formalcharge: int
            The formalcharge of the atom

        Returns
        -------
        missingatoms: int
            The number of atom that are missing for completing the valence
        """

        pT = self.PeriodicaTable

        element = smallmol.element[atomidx]
        btypes = smallmol.bondtypes[atomidx]
        # steric_number = len(smallmol.neighbors[atomidx])

        valences = list(pT.GetValenceList(element))

        # Some atoms can have multiple valid valences (like S, P). For them I should predict the most reliable one.
        # At the moment not implemented. I will raise an error now.
        if len(valences) == 1:
            valence = valences[0]
        else:
            class Notmplemented(Exception):
                pass
            raise Notmplemented('The missing valence is not implemented for atom with more than one valid valence')

        missingatoms = self._evaluateMissingAtoms(valence,  btypes)

        missingatoms += formalcharge

        return missingatoms

    def getAtom(self, element, attachTo=None, coords=None):
        """
        Returns a dictionary with all the information necessary to create a new atom by using
        htmd.smallmol.chemlab.Builder

        Parameters
        ----------
        element: str
            The element of the atom
        attachTo: int
            The index of the atom you want the new atom will be bonded to
        coords:
            The coordinates of the new atom

        Returns
        -------
        atom_data: dict
            The atom data

        Example
        -------
        >>> PT.getAtom('H')
        {'attachTo': '', 'bondtype': 1, 'charge': 0, 'chiral': '', 'coords': array([[[0],[0],[0]]]), 'element': 'H',
        'formalcharge': 0, 'hybridization': 1}
        """

        if element not in _atoms:
            raise ValueError('The element {} is not available'.format(element))

        keys = ['element', 'charge', 'formalcharge', 'chiral', 'hybridization', 'bondtype']

        coords = coords if coords is not None else np.array([0, 0, 0])
        coords = {'coords': coords.reshape(1, 3, 1)}

        atom_data = {k: v for k, v in zip(keys, _atoms[element])}
        atom_data.update(coords)

        if attachTo is not None:
            atom_data.update({'attachTo': attachTo})
        else:
            atom_data.update({'attachTo': ''})

        return atom_data

    def listHybridizations(self):
        """
        Lists the Hybridization type available

        """

        dict_hyb = HybridizationType.values

        print("Index  --> Hybridization\n")
        for i, hyb in dict_hyb.items():
            print("%5s " % i, '-->', " %-13s" % str(hyb))

    def listBondTypes(self):
        """
        Lists the Bonds type available

        """

        dict_bt = BondType.values

        print("Index  --> BondType\n")
        for i, bt in dict_bt.items():
            print("%5s " % i, '-->', " %-13s" % str(bt))


    def getHybridization(self, index):
        """
        Returns the hybridization type from its index value

        Parameters
        ----------
        index: int
            The index representing the hybridization type

        Returns
        -------
        hybtype: str
            The hybridization type

        """

        dict_hyb = HybridizationType.values

        return str(dict_hyb[index])

    def getBondType(self, index):
        """
        Returns the bond type from its index value

        Parameters
        ----------
        index: int
            The index representing the bond type

        Returns
        -------
        btype: str
            The bond type

        """

        dict_bt = BondType.values

        return str(dict_bt[index])