

from rdkit.Chem import GetPeriodicTable
import numpy as np


#{'element':'H', 'charge':0, 'formalcharge':0, 'chiral':'', 'hybridization':1, 'bondtype':0},
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

