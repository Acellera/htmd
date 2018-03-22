
from rdkit.Chem.rdchem import HybridizationType, BondType, ChiralType
from rdkit.Chem import GetPeriodicTable
from rdkit import rdBase
import numpy as np


_heavy_atoms = ['C', 'N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']

_hetero_atoms = ['N', 'O', 'F',
                'Si','P', 'S', 'Cl', 'Br', 'I']

_hetero_atoms = ['N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']
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

_chiral_type = {ChiralType.CHI_TETRAHEDRAL_CW:'clockwise',
                ChiralType.CHI_TETRAHEDRAL_CCW:'anitclockwise'}
_chiral_type_Dict = ChiralType.values

_bondtypes_IdxToType = BondType.values
_bondtypes_StringToType = {'SINGLE': BondType.SINGLE,
                           'DOUBLE': BondType.DOUBLE,
                           'TRIPLE': BondType.TRIPLE}

_atoms = {'H': [ 'H', 0, 0, '', 1, 1],
         'C': ['C', 0, 0, '', 4, 1],
        }

class PeriodicTable:

    def __init__(self):
        self.PeriodicaTable = GetPeriodicTable()

    def _evaluateMissingAtoms(self, valence, btypes):
        normbondvalence = sum([ 1.5 if int(b) == 12 else int(b) for b in btypes  ])

        return int(valence - normbondvalence)

    def getRadiusBond(self, element):
        pT = self.PeriodicaTable

        return pT.GetRb0(element)


    def getMissingValence(self, atomidx, smallmol, formalcharge=0):
        pT = self.PeriodicaTable

        element = smallmol.element[atomidx]
        btypes = smallmol.bondtypes[atomidx]
        #steric_number = len(smallmol.neighbors[atomidx])

        valences = list(pT.GetValenceList(element))
        print(valences)
        if len(valences) == 1:
            valence = valences[0]
        else:
            #TODO
            pass

        missingatoms = self._evaluateMissingAtoms(valence,  btypes)

        return missingatoms + formalcharge

    def getAtom(self, element, attachTo=None, coords=None ):

        keys = ['element', 'charge','formalcharge', 'chiral', 'hybridization', 'bondtype']


        coords = coords if coords is not None else np.array([0,0,0])
        coords = {'coords': coords.reshape(1,3,1)}

        atom_data = {k:v for k,v in zip(keys, _atoms[element]) }
        atom_data.update(coords)

        if attachTo is not None:
            atom_data.update({'attachTo':attachTo})
        else:
            atom_data.update({'attachTo': ''})

        return atom_data

