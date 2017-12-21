# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
from collections import OrderedDict
import itertools
import numpy as np
import networkx as nx

from htmd.molecule.vdw import elements

logger = logging.getLogger(__name__)


def getMolecularGraph(molecule):
    """
    Generate a graph from the topology of molecule

    The graph nodes represent atoms and the graph edges represent bonds. Also, the graph nodes store element
    information.
    """

    graph = nx.Graph()
    for i, element in enumerate(molecule.element):
        graph.add_node(i, element=element, number=elements.index(element))
    graph.add_edges_from(molecule.bonds)

    return graph

def getMolecularTree(graph, source):
    """
    Generate a tree from a molecular graph

    The tree starts from source node (atom) and grows along the edges (bonds) unrolling all encountered loops.
    The tree grows until all the nodes (atoms) of the graph are included.
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
    Check if two molecular graphs are isomorphic based on topology (bonds) and elements.
    """

    return nx.is_isomorphic(graph1, graph2, node_match = lambda node1, node2: node1['element'] == node2['element'])


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

    # Generate a molecular tree for each atom
    graph = getMolecularGraph(molecule)
    trees = [getMolecularTree(graph, node) for node in graph.nodes]

    equivalent_atoms = [tuple([i for i, tree1 in enumerate(trees) if checkIsomorphism(tree1, tree2)]) for tree2 in trees]
    equivalent_groups = sorted(list(set(equivalent_atoms)))
    equivalent_group_by_atom = list(map(equivalent_groups.index, equivalent_atoms))

    return equivalent_groups, equivalent_atoms, equivalent_group_by_atom

def getMethylGraph():
    """
    Generate a molecular graph for methyl group
    """

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
    """
    Filter irrelevant dihedral angle cores (central atom pairs)
    """

    _, sideGraphs = core_sides
    methyl = getMethylGraph()

    # Check if a side graph is a methyl group
    if checkIsomorphism(sideGraphs[0], methyl) or checkIsomorphism(sideGraphs[1], methyl):
        return False

    return True

def detectParameterizableCores(molecule):
    """
    Detect parametrizable dihedral angle cores (central atom pairs)

    The cores are detected by looking for bridges (bonds with divide the molecule into two parts) in a molecular graph.
    Terminal cores are skipped.
    """

    graph = getMolecularGraph(molecule)

    all_core_sides = []
    for core in list(nx.bridges(graph)):

        # Get side graphs of the core
        graph.remove_edge(*core)
        sideGraphs = list(nx.connected_component_subgraphs(graph))
        graph.add_edge(*core)

        # Skip terminal bridges, which cannot for dihedral angles
        if len(sideGraphs[0]) == 1 or len(sideGraphs[1]) == 1:
            continue

        # Swap the side graphs to match the order of the core
        if core[0] in sideGraphs[1]:
            sideGraphs = sideGraphs[::-1]
        assert core[0] in sideGraphs[0] and core[1] in sideGraphs[1]

        all_core_sides.append((core, sideGraphs))

    # Filer cores
    all_core_sides = filter(filterCores, all_core_sides)

    return all_core_sides

def chooseTerminals(centre, sideGraph):
    """
    Choose dihedal angle terminals (outer atoms)

    The terminals are chosen by:
    1. Largest number of atoms
    2. Largest molecular mass
    """

    terminals = list(sideGraph.neighbors(centre))

    # Get a subgraph for each terminal
    sideGraph = sideGraph.copy()
    sideGraph.remove_node(centre)
    terminalGraphs = list(nx.connected_component_subgraphs(sideGraph))
    terminalGraphs = [[terminalGraph for terminalGraph in terminalGraphs if terminal in terminalGraph.nodes][0] for terminal in terminals]

    # Compute a score for each terminal
    numberOfAtoms = np.array([len(terminalGraph.nodes) for terminalGraph in terminalGraphs])
    numberOfProtons = np.array([sum(nx.get_node_attributes(terminalGraph, 'number').values()) for terminalGraph in terminalGraphs])
    assert max(numberOfProtons) < 1000000 # No monstrous molecules, please!
    scores = 1000000*numberOfAtoms + numberOfProtons

    # Choose the terminals
    chosen_terminals = []
    refTerminalGraph = None
    for terminal, score, terminalGraph in zip(terminals, scores, terminalGraphs):
        if score < np.max(scores):
            continue

        if not chosen_terminals:
            chosen_terminals.append(terminal)
            refTerminalGraph = terminalGraph
            continue

        if checkIsomorphism(terminalGraph, refTerminalGraph):
            chosen_terminals.append(terminal)
        else:
            logger.warn('Molecular scoring function is not sufficient. '
                        'Dihedal atom selection depends on the atom order!')

    return chosen_terminals

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
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> from htmd.parameterization.detect import detectParameterizableDihedrals

    Find the parameterizable dihedrals of glycol
    >>> molFile = os.path.join(home('test-param'), 'glycol.mol2')
    >>> mol = FFMolecule(molFile)
    >>> detectParameterizableDihedrals(mol)
    [[(2, 1, 0, 4), (1, 2, 3, 9)], [(0, 1, 2, 3)]]

    Find the parameterizable dihedrals of ethanolamine
    >>> molFile = os.path.join(home('test-param'), 'ethanolamine.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2) # Does not work with CGenFF
    >>> detectParameterizableDihedrals(mol)
    [[(2, 1, 0, 4)], [(0, 1, 2, 3)], [(1, 2, 3, 9), (1, 2, 3, 10)]]

    Find the parameterizable dihedrals of benzamidine
    >>> molFile = os.path.join(home('test-param'), 'benzamidine.mol2')
    >>> mol = FFMolecule(molFile)
    >>> detectParameterizableDihedrals(mol)
    [[(1, 0, 6, 12), (1, 0, 6, 13), (5, 0, 6, 12), (5, 0, 6, 13)], [(0, 6, 12, 16), (0, 6, 12, 17), (0, 6, 13, 14), (0, 6, 13, 15)]]

    # Check if the atom swapping does not affect results

    Find the parameterizable dihedrals of chlorethene
    >>> molFile = os.path.join(home('test-param'), 'chlorethene_1.mol2')
    >>> mol = FFMolecule(molFile)
    >>> detectParameterizableDihedrals(mol)
    [[(2, 1, 0, 4), (2, 1, 0, 5)]]

    Find the parameterizable dihedrals of chlorethene (with swapped atoms)
    >>> molFile = os.path.join(home('test-param'), 'chlorethene_2.mol2')
    >>> mol = FFMolecule(molFile)
    >>> detectParameterizableDihedrals(mol)
    [[(3, 1, 0, 4), (3, 1, 0, 5)]]

    # Check the scoring function

    Find the parameterizable dihedrals of dicarbothioic acid
    >>> molFile = os.path.join(home('test-param'), 'dicarbothioic_acid.mol2')
    >>> mol = FFMolecule(molFile)
    >>> detectParameterizableDihedrals(mol)
    [[(3, 1, 0, 6)], [(0, 1, 3, 5)], [(1, 3, 5, 7)]]
    """

    # Get pameterizable dihedral angles
    dihedrals = []
    for core, sides in detectParameterizableCores(molecule):

        # Choose the best terminals for each side
        all_terminals = map(chooseTerminals, core, sides)

        # Generate all terminal combinations
        all_terminals = itertools.product(*all_terminals)

        # Generate new dihedral angles
        new_dihedrals = [(terminals[0], *core, terminals[1]) for terminals in all_terminals]
        dihedrals.extend(new_dihedrals)

    # Get equivalent groups for each atom for each dihedral
    _, _, equivalent_group_by_atom = detectEquivalentAtoms(molecule)
    dihedral_groups = [tuple([equivalent_group_by_atom[atom] for atom in dihedral]) for dihedral in dihedrals]

    # Group equivalent dihedral angles and reverse them that equivalent atoms are matched
    equivalent_dihedrals = OrderedDict()
    for dihedral, groups in zip(dihedrals, dihedral_groups):
        dihedral, groups = (dihedral[::-1], groups[::-1]) if groups[::-1] < groups else (dihedral, groups)
        equivalent_dihedrals[groups] = equivalent_dihedrals.get(groups, []) + [dihedral]
    equivalent_dihedrals = list(equivalent_dihedrals.values())

    return equivalent_dihedrals


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)