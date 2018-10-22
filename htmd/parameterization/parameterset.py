# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import unittest

import numpy as np

_ATOM_TYPE_REG_EX = re.compile('^\S+x\d+$')


def getParameter(type, parameterfield):
    if type in parameterfield:
        return parameterfield[type]
    type = tuple(reversed(type))
    if type in parameterfield:
        return parameterfield[type]
    raise RuntimeError('Could not find parameters for key {}'.format(type))


def getImproperParameter(type, parameters):
    from itertools import permutations
    type = np.array(type)
    perms = np.array([x for x in list(permutations((0, 1, 2, 3))) if x[2] == 2])
    for p in perms:
        if tuple(type[p]) in parameters.improper_types:
            return parameters.improper_types[tuple(type[p])], 'improper_types'
        elif tuple(type[p]) in parameters.improper_periodic_types:
            return parameters.improper_periodic_types[tuple(type[p])], 'improper_periodic_types'
    raise RuntimeError('Could not find improper parameters for key {}'.format(type))


def recreateParameters(mol, originaltypes, parameters):
    from copy import copy
    from parmed.parameters import ParameterSet

    newparams = ParameterSet()
    uqtypes = np.unique(mol.atomtype)

    for type in uqtypes:
        newparams.atom_types[type] = copy(parameters.atom_types[originaltypes[type]])
        if type != originaltypes[type]:
            newparams.atom_types[type].number = newparams.atom_types[type].number + 900 # add a big offset so it doesn't collide with real charm types

    for idx in mol.bonds:
        newkey = tuple(mol.atomtype[idx])
        oldkey = tuple(np.vectorize(originaltypes.get)(newkey))
        newparams.bond_types[newkey] = copy(parameters.bond_types[oldkey])

    for idx in mol.angles:
        newkey = tuple(mol.atomtype[idx])
        oldkey = tuple(np.vectorize(originaltypes.get)(newkey))
        newparams.angle_types[newkey] = copy(parameters.angle_types[oldkey])

    for idx in mol.dihedrals:
        newkey = tuple(mol.atomtype[idx])
        oldkey = tuple(np.vectorize(originaltypes.get)(newkey))
        newparams.dihedral_types[newkey] = copy(parameters.dihedral_types[oldkey])

    for idx in mol.impropers:
        newkey = tuple(mol.atomtype[idx])
        oldkey = np.vectorize(originaltypes.get)(newkey)
        oldval, field = getImproperParameter(oldkey, parameters)
        newparams.__dict__[field][newkey] = copy(oldval)

    return newparams


def createMultitermDihedralTypes(parameters, nterms=6, scee=1.2, scnb=2):
    from parmed.topologyobjects import DihedralTypeList, DihedralType
    from copy import deepcopy

    parameters = deepcopy(parameters)

    for key, val in parameters.dihedral_types.items():
        dihlist = DihedralTypeList()
        for i in range(1, nterms+1):
            found = False
            for d in val: # Check if this term already exists in the parameters.
                if d.per == i:
                    dihlist.append(d)
                    found = True
                    break
            if not found: # Else create an unparametrized term
                dihtype = DihedralType(0, i, 0, scee=scee, scnb=scnb)
                dihlist.append(dihtype)
        parameters.dihedral_types[key] = dihlist

    return parameters


def inventAtomTypes(mol, fit_dihedrals, equivalents):
    """
    Duplicate atom types of the dihedral, so its parameters are unique.
    """
    # TODO check symmetry

    mol = mol.copy()

    alltypes = list(mol.atomtype)
    originaltype = {type: type for type in np.unique(mol.atomtype)}

    # Duplicate the atom types of the dihedral
    for dih in fit_dihedrals:
        for d in dih:
            oldtype = mol.atomtype[d]
            # if the type is already duplicated
            if re.match(_ATOM_TYPE_REG_EX, oldtype):
                continue
            # Create a new atom type name
            i = 0
            while ('{}x{}'.format(oldtype, i)) in alltypes:
                i += 1
            newtype = '{}x{}'.format(oldtype, i)
            alltypes.append(newtype)
            originaltype[newtype] = oldtype

            mol.atomtype[d] = newtype
            # Rename the atom types of the equivalent atoms
            for index in equivalents[1][d]:
                if index != d:
                    assert not re.match(_ATOM_TYPE_REG_EX, mol.atomtype[index])
                    mol.atomtype[index] = newtype

        equivalent_dihedrals = _getEquivalentDihedrals(mol, equivalents, np.array(dih))
        if len(equivalent_dihedrals) > 1:
            raise ValueError("Dihedral term still not unique after duplication. Dihedral {} has {} equivalent "
                             "dihedrals: {}".format(dih, len(equivalent_dihedrals), equivalent_dihedrals))

    return mol, originaltype


def _getEquivalentDihedrals(mol, equivalents, dihedral):
    """
    Find equivalent dihedral angles to the specificied one.
    """
    equivalent_group_by_atom = equivalents[2]
    types = mol.atomtype[dihedral]

    same_type_dihedrals = []
    for dihedral_indices in mol.dihedrals:
        dihedral_types = mol.atomtype[dihedral_indices]
        if np.all(types == dihedral_types) or np.all(types == dihedral_types[::-1]):
            same_type_dihedrals.append(dihedral_indices)

    # Now for each of the uses, remove any which are equivalent
    unique_dihedrals = [dihedral]
    groups = [equivalent_group_by_atom[index] for index in dihedral]
    for dihed in same_type_dihedrals:
        dihedral_groups = [equivalent_group_by_atom[index] for index in dihed]
        if groups != dihedral_groups and groups != dihedral_groups[::-1]:
            unique_dihedrals.append(dihed)

    return unique_dihedrals


class Test(unittest.TestCase):

    def setUp(self):
        from htmd.home import home
        from htmd.parameterization.fftype import fftype
        from htmd.parameterization.util import getEquivalentsAndDihedrals
        from htmd.molecule.molecule import Molecule

        molFile = os.path.join(home('test-param'), 'glycol.mol2')
        mol = Molecule(molFile)
        mol, self.equivalents, all_dihedrals = getEquivalentsAndDihedrals(mol)
        _, mol = fftype(mol, method='GAFF2')
        self.mol = mol

    def test_getEquivalentDihedrals(self):
        self.assertListEqual(_getEquivalentDihedrals(self.mol, self.equivalents, [0, 1, 2, 3]), [[0, 1, 2, 3]])
        self.assertListEqual(_getEquivalentDihedrals(self.mol, self.equivalents, [4, 0, 1, 2]), [[4, 0, 1, 2]])
        self.assertListEqual(_getEquivalentDihedrals(self.mol, self.equivalents, [5, 1, 2, 7]), [[5, 1, 2, 7]])

    def test_inventAtomTypes(self):
        self.assertListEqual(self.mol.atomtype.tolist(), ['oh', 'c3', 'c3', 'oh', 'ho', 'h1', 'h1', 'h1', 'h1', 'ho'])

        mol = self.mol.copy()
        mol, _ = inventAtomTypes(mol, [[0, 1, 2, 3]], self.equivalents)
        self.assertListEqual(mol.atomtype.tolist(), ['ohx0', 'c3x0', 'c3x0', 'ohx0', 'ho', 'h1', 'h1', 'h1', 'h1', 'ho'])

        mol = self.mol.copy()
        mol, _ = inventAtomTypes(mol, [[4, 0, 1, 2]], self.equivalents)
        self.assertListEqual(mol.atomtype.tolist(), ['ohx0', 'c3x0', 'c3x0', 'ohx0', 'hox0', 'h1', 'h1', 'h1', 'h1', 'hox0'])

        mol = self.mol.copy()
        mol, _ = inventAtomTypes(mol, [[5, 1, 2, 7]], self.equivalents)
        self.assertListEqual(mol.atomtype.tolist(), ['oh', 'c3x0', 'c3x0', 'oh', 'ho', 'h1x0', 'h1x0', 'h1x0', 'h1x0', 'ho'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
