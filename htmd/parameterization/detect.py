# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from collections import OrderedDict
import numpy as np
import scipy.sparse.csgraph as sp
import networkx as nx

from htmd.molecule.vdw import elements


def getMolecularGraph(molecule):
    """
    Generate a graph from the topology of molecule, i.e the graph nodes represent atoms and  the graph edges represent
    bonds. Also, the graph nodes store element information.
    """

    graph = nx.Graph()
    for i, element in enumerate(molecule.element):
        graph.add_node(i, element=element, number=elements.index(element))
    graph.add_edges_from(molecule.bonds)

    return graph

def getMolecularTree(graph, source):
    """
    Generate a tree from a molecular graph. The tree starts from source node (atom) and grows along the edges (bonds)
    unrolling all encountered loops. The tree grows untill all the nodes (atoms) of the graph are included.
    """

    assert nx.is_connected(graph)

    tree = nx.DiGraph()
    tree.add_node(0, base=source, element=graph.nodes[source]['element'])
    current_nodes = list(tree.nodes)
    base_nodes = {source}

    while True:
        new_nodes = []
        neighbor_filter = lambda node: node not in base_nodes
        for current_node in current_nodes:
            for neighbor in filter(neighbor_filter, graph.neighbors(tree.nodes[current_node]['base'])):
                new_node = len(tree.nodes)
                tree.add_node(new_node, base=neighbor, element=graph.nodes[neighbor]['element'])
                tree.add_edge(current_node, new_node)
                new_nodes.append(new_node)

        current_nodes = new_nodes
        base_nodes = {base for _, base in tree.nodes.data('base')}
        if base_nodes == set(graph.nodes):
            break

    return tree

def checkIsomorphism(graph1, graph2):
    """
    Check if two graphs are isomorphic based on topology and elements
    """

    return nx.is_isomorphic(graph1, graph2,
                            node_match = lambda node1, node2: node1['element'] == node2['element'])


def detectEquivalentAtoms(molecule):
    """
    Detect topologically equivalent atoms.

    Arguments
    ---------
    molecule : FFMolecule
        Molecule object

    Return
    ------
    equivalent_groups : list of tuples
        List of equivalent atom group. Each element is a tuple contain equivalent atom indices.
    equivalent_atoms : list of tuples
        List of equivalent atom group for each atom. Each element is a tuple contain equivalent atom indices.
    equivalent_group_by_atom : list
        List of equivalent group indices for each atom. The indices corresponds to `equivalent_groups` order.

    Examples
    --------

    >>> import os
    >>> from htmd.home import home
    >>> from htmd.parameterization.ffmolecule import FFMolecule
    >>> from htmd.parameterization.detect import detectEquivalentAtoms

    Get benzamidine
    >>> molFile = os.path.join(home('test-param'), 'benzamidine.mol2')
    >>> mol = FFMolecule(molFile)

    Find the equivalent atoms of bezamidine
    >>> equivalent_groups, equivalent_atoms, equivalent_group_by_atom = detectEquivalentAtoms(mol)
    >>> equivalent_groups
    [(0,), (1, 5), (2, 4), (3,), (6,), (7, 11), (8, 10), (9,), (12, 13), (14, 15, 16, 17)]
    >>> equivalent_atoms
    [(0,), (1, 5), (2, 4), (3,), (2, 4), (1, 5), (6,), (7, 11), (8, 10), (9,), (8, 10), (7, 11), (12, 13), (12, 13), (14, 15, 16, 17), (14, 15, 16, 17), (14, 15, 16, 17), (14, 15, 16, 17)]
    >>> equivalent_group_by_atom
    [0, 1, 2, 3, 2, 1, 4, 5, 6, 7, 6, 5, 8, 8, 9, 9, 9, 9]

    Get dicarbothioic acid
    >>> molFile = os.path.join(home('test-param'), 'dicarbothioic_acid.mol2')
    >>> mol = FFMolecule(molFile)

    Find the equivalent atoms of dicarbothioic acid
    >>> equivalent_groups, equivalent_atoms, equivalent_group_by_atom = detectEquivalentAtoms(mol)
    >>> equivalent_groups
    [(0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,)]
    >>> equivalent_atoms
    [(0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,)]
    >>> equivalent_group_by_atom
    [0, 1, 2, 3, 4, 5, 6, 7]
    """

    graph = getMolecularGraph(molecule)
    trees = [getMolecularTree(graph, node) for node in graph.nodes]

    equivalent_atoms = [tuple([i for i, tree1 in enumerate(trees) if checkIsomorphism(tree1, tree2)]) for tree2 in trees]
    equivalent_groups = sorted(list(set(equivalent_atoms)))
    equivalent_group_by_atom = [equivalent_groups.index(atoms) for atoms in equivalent_atoms]

    return equivalent_groups, equivalent_atoms, equivalent_group_by_atom

def _my_breadth_first(mol, con, start, depth=1):

    con = con.copy()
    N = mol.coords.shape[0]

    ret = []

    search = []
    for j in start:
        for i in range(N):

            if con[j, i] == 1 and j != i:
                # con[j,i] = 0
                con[i, j] = 0
                if i not in search:
                    search.append(i)
    if len(search):
        ret.append(search)
        n = _my_breadth_first(mol, con, search, depth + 1)
        if len(n):
            for i in n:
                ret.append(i)

    # top levelreturn. canonicalise the ordering of the elements at each level
    if depth == 1:
        rr = []
        for r in ret:
            for i in range(len(r)):
                r[i] = mol.element[r[i]]
            r.sort()
            rr.append("!")
            for i in r:
                rr.append(i)
        ret = rr
    return ret


def detectEquivalents(mol):

    bonds = mol.bonds
    natoms = mol.coords.shape[0]

    # Make a connectivity matrix
    conn = np.zeros((natoms, natoms), dtype=np.bool)
    for b in bonds:
        conn[b[0], b[1]] = True
        conn[b[1], b[0]] = True

    paths = []
    uniq_paths = []
    for b in range(natoms):
        a = _my_breadth_first(mol, conn, [b])
        paths.append(a)

    for b in paths:
        found = False
        for c in uniq_paths:
            if c == b:
                found = True
                break
        if not found:
            uniq_paths.append(b)

    equiv_groups = []
    for b in uniq_paths:
        e = []
        i = 0
        for c in paths:
            if b == c:
                e.append(i)
            i += 1
        equiv_groups.append(e)

    equiv_atoms = list(range(natoms))
    equiv_group_by_atom = list(range(natoms))
    i = 0
    for a in equiv_groups:
        if type(a) == int:
            a = [a]
        for b in a:
            equiv_atoms[b] = a
            equiv_group_by_atom[b] = i
        i += 1

    return equiv_groups, equiv_atoms, equiv_group_by_atom


class SoftDihedral:

    def __init__(self, atoms, left=None, right=None, left_1=None, right_1=None, equiv=None):

        self.atoms = atoms
        self.left = [] if left is None else left
        self.right = [] if right is None else right
        self.left_1 = [] if left_1 is None else left_1
        self.right_1 = [] if right_1 is None else right_1
        self.equivalents = [] if equiv is None else equiv


def getMethylGraph():

    methyl = nx.Graph()
    methyl.add_node(0, element='C')
    methyl.add_node(1, element='H')
    methyl.add_node(2, element='H')
    methyl.add_node(3, element='H')
    methyl.add_edge(0 ,1)
    methyl.add_edge(0, 2)
    methyl.add_edge(0, 3)

    return methyl

def filterCores(core_sides):

    _, sides = core_sides

    if len(sides[0]) == 1 or len(sides[1]) == 1:
        return False

    methyl = getMethylGraph()
    if checkIsomorphism(sides[0], methyl) or checkIsomorphism(sides[1], methyl):
        return False

    return True

def detectParameterizableCores(molecule):

    graph = getMolecularGraph(molecule)

    all_core_sides = []
    for core in list(nx.bridges(graph)):
        graph.remove_edge(*core)
        sides = list(nx.connected_component_subgraphs(graph))
        graph.add_edge(*core)

        if core[0] in sides[1]:
            sides = sides[::-1]
        assert core[0] in sides[0] and core[1] in sides[1]

        all_core_sides.append((core, sides))

    all_core_sides = filter(filterCores, all_core_sides)

    return all_core_sides

def chooseTerminal(centre, sideGraph):

    terminals = list(sideGraph.neighbors(centre))

    # Get a subgraph for each terminal
    sideGraph = sideGraph.copy()
    sideGraph.remove_node(centre)
    subGraphs = list(nx.connected_component_subgraphs(sideGraph))
    subGraphs = [[subGraph for subGraph in subGraphs if terminal in subGraph.nodes][0] for terminal in terminals]

    # Compute a score for each terminal
    numberOfAtoms = np.array([len(subGraph.nodes) for subGraph in subGraphs])
    numberOfProtons = np.array([sum(nx.get_node_attributes(subGraph, 'number').values()) for subGraph in subGraphs])
    assert max(numberOfProtons) < 1000000 # No monstrous molecules, please!
    scores = 1000000*numberOfAtoms + numberOfProtons

    terminal = terminals[np.argmax(scores)]

    return terminal

def detectParameterizableDihedrals(molecule):
    """
    Detect parameterizable dihedral angles

    Arguments
    ---------
    molecule : FFMolecule
        Molecule object

    Return
    ------
    dihedrals : list of list of tuples
        List of equivalent dihedral angle groups. Each group is a list of equivalent dihedral angles.
        Each angle is defined as a tuple of four atom indices (0-based).

    Examples
    --------

    >>> import os
    >>> from htmd.home import home
    >>> from htmd.parameterization.ffmolecule import FFMolecule
    >>> from htmd.parameterization.detect import detectParameterizableDihedrals

    Get benzamidine
    >>> molFile = os.path.join(home('test-param'), 'benzamidine.mol2')
    >>> mol = FFMolecule(molFile)

    Find the parameterizable dihedrals of benzamidine
    >>> detectParameterizableDihedrals(mol)
    [[(1, 0, 6, 12)], [(0, 6, 12, 16), (0, 6, 13, 14)]]

    Get glycol
    >>> molFile = os.path.join(home('test-param'), 'glycol.mol2')
    >>> mol = FFMolecule(molFile)

    Find the parameterizable dihedrals of glycol
    >>> detectParameterizableDihedrals(mol)
    [[(2, 1, 0, 4), (1, 2, 3, 9)], [(0, 1, 2, 3)]]
    """

    # Get pameterizable dihedral angles
    dihedrals = []
    for core, sides in detectParameterizableCores(molecule):
        terminals = [chooseTerminal(centre, side) for centre, side in zip(core, sides)]
        dihedrals.append((terminals[0], core[0], core[1], terminals[1]))

    # Get equivantence groups for each atom for each dihedral
    _, _, equivalent_group_by_atom = detectEquivalentAtoms(molecule)
    dihedral_groups = [tuple([equivalent_group_by_atom[atom] for atom in dihedral]) for dihedral in dihedrals]

    # Group equivalent dihedral angles and reverse them that equivalent atoms are matched
    equivalent_dihedrals = OrderedDict()
    for dihedral, groups in zip(dihedrals, dihedral_groups):
        if groups[::-1] < groups:
            groups = groups[::-1]
            dihedral = dihedral[::-1]
        equivalent_dihedrals[groups] = equivalent_dihedrals.get(groups, []) + [dihedral]
    equivalent_dihedrals = list(equivalent_dihedrals.values())

    return equivalent_dihedrals

def detectSoftDihedrals(mol, equivalent_atoms):

    bonds = mol.bonds
    natoms = mol.coords.shape[0]

    # Make a connectivity matrix
    conn = np.zeros((natoms, natoms), dtype=np.bool)
    for b in bonds:
        conn[b[0], b[1]] = True
        conn[b[1], b[0]] = True

    # Iterate over each of the dihedrals, checking to see which partition the graph
    possible_soft = []
    for d in mol.dihedrals:
        a0 = d[0]
        a1 = d[1]
        a2 = d[2]
        a3 = d[3]
        c = conn.copy()
        c[a1, a2] = c[a2, a1] = False
        left = sp.breadth_first_tree(c, a1, directed=False).indices.flatten()
        right = sp.breadth_first_tree(c, a2, directed=False).indices.flatten()
        left = np.unique(left)
        right = np.unique(right)

        c = conn.copy()
        c[a0, a1] = c[a1, a0] = False
        c[a1, a2] = c[a2, a1] = False
        c[a3, a2] = c[a2, a3] = False
        left_1 = sp.breadth_first_tree(c, a0, directed=False)
        left_1 = left_1.indices
        left_1 = left_1.flatten()

        right_1 = sp.breadth_first_tree(c, a3, directed=False)
        right_1 = right_1.indices
        right_1 = right_1.flatten()

        left_1 = np.unique(left_1)
        right_1 = np.unique(right_1)

        if not (a2 in left) and not (a1 in right):
            possible_soft.append(SoftDihedral(d, left, right, left_1, right_1, []))

    final_soft = []
    e = mol.element

    for d in possible_soft:
        a1 = d.atoms[1]
        a2 = d.atoms[2]

        left = d.left
        right = d.right

        # Exclude methyls
        if len(left) == 3:
            if e[a1] == 'C' and e[left[0]] == 'H' and e[left[1]] == 'H' and e[left[2]] == 'H':
                continue
        if len(right) == 3:
            if e[a2] == 'C' and e[right[0]] == 'H' and e[right[1]] == 'H' and e[right[2]] == 'H':
                continue
        found = False

        # check to see if the torsional pair of atoms are already included in the list.
        for g in final_soft:
            f = g.atoms
            if f[1] == a1 and f[2] == a2:
                found = True
                break
            if f[2] == a1 and f[1] == a2:
                found = True
                break
        if not found:
            final_soft.append(d)

    final_soft = remove_equivalents(mol, final_soft, equivalent_atoms)

    return final_soft


def remove_equivalents(mol, soft, equiv):

    equivalent_atom_groups = equiv[0]  # list of groups of equivalent atoms
    equivalent_atoms = equiv[1]  # list of equivalent atoms, indexed by atom
    equivalent_group_by_atom = equiv[2]  # mapping from atom index to equivalent atom group

    # for each of the soft dihedrals, remove any which are equivalent to others through symmetry
    # compare only the middle two atoms since we care about not duplicating X-A-B-X and X-A'-B'-X
    final_soft = []

    for i in range(len(soft)):
        h1 = equivalent_group_by_atom[soft[i].atoms[0]]
        h2 = equivalent_group_by_atom[soft[i].atoms[1]]
        h3 = equivalent_group_by_atom[soft[i].atoms[2]]
        h4 = equivalent_group_by_atom[soft[i].atoms[3]]
        found = None
        for j in range(len(final_soft)):
            g1 = equivalent_group_by_atom[final_soft[j].atoms[0]]
            g2 = equivalent_group_by_atom[final_soft[j].atoms[1]]
            g3 = equivalent_group_by_atom[final_soft[j].atoms[2]]
            g4 = equivalent_group_by_atom[final_soft[j].atoms[3]]

            if h2 == g2 and h3 == g3:
                found = j
            if h3 == g2 and h2 == g3:
                found = j

        if found is not None:
            # check to see which the two -- the one in already in the final list
            # and the one we've just found -- has more "weight"
            # So that we aren't choosing based on the arbitrary graph traversal ordering
            # Weight is the # of atoms in the left_1 + right_1 groups
            already_in_list = final_soft[found]
            candidate = soft[i]
            if (len(already_in_list.left_1) + len(already_in_list.right_1)) < \
                    (len(candidate.left_1) + len(candidate.right_1)):
                final_soft.remove(already_in_list)
                final_soft.append(candidate)

        else:
            final_soft.append(soft[i])

    # now for each of the unique soft dihedrals note the dihedrals that also use the same type
    for s in final_soft:
        a1 = s.atoms[0]
        a2 = s.atoms[1]
        a3 = s.atoms[2]
        a4 = s.atoms[3]

        h1 = equivalent_group_by_atom[a1]
        h2 = equivalent_group_by_atom[a2]
        h3 = equivalent_group_by_atom[a3]
        h4 = equivalent_group_by_atom[a4]

        for d in mol.dihedrals:
            b1 = d[0]
            b2 = d[1]
            b3 = d[2]
            b4 = d[3]

            g1 = equivalent_group_by_atom[b1]
            g2 = equivalent_group_by_atom[b2]
            g3 = equivalent_group_by_atom[b3]
            g4 = equivalent_group_by_atom[b4]

            found = False

            if a1 == b1 and a2 == b2 and a3 == b3 and a4 == b4:
                found = True
            if a4 == b1 and a3 == b2 and a2 == b3 and a1 == b4:
                found = True
            if found is False:  # d is not the soft dihedral itself
                found = False
                if h1 == g1 and h2 == g2 and h3 == g3 and h4 == g4:
                    found = True
                if h4 == g1 and h3 == g2 and h2 == g3 and h1 == g4:
                    found = True
                if found is True:  # this dihedral shares a type with the soft dihedral
                    s.equivalents.append([b1, b2, b3, b4])

    return final_soft


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)