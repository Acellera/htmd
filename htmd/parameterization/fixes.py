# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import networkx as nx

logger = logging.getLogger(__name__)


def _getMolecularGraph(molecule):
    """
    Generate a graph from the topology of molecule
    """

    graph = nx.Graph()
    for i in range(molecule.numAtoms):
        graph.add_node(i, element=molecule.element[i])
    for i, bond in enumerate(molecule.bonds):
        graph.add_edge(*bond, index=i)

    return graph

def fixPhosphateTypes(molecule):
    """
    >>> from moleculekit.molecule import Molecule
    >>> from htmd.charge import fitGasteigerCharges

    >>> mol = Molecule('/home/raimis/pm-param.git/param_tests_2/ligForRamis/1a1e_ligand.mol2')

    >>> fitGasteigerCharges(mol) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    AttributeError: 'NoneType' object has no attribute 'GetNumAtoms'

    >>> new_mol = fixPhosphateTypes(mol)

    >>> print(mol.atomtype)
    ['C.2' 'O.2' 'C.3' 'N.am' 'C.3' 'C.2' 'O.2' 'C.3' 'C.ar' 'C.ar' 'C.ar'
     'C.ar' 'C.ar' 'C.ar' 'O.3' 'P.3' 'O.co2' 'O.co2' 'O.co2' 'N.am' 'C.3'
     'C.2' 'O.2' 'C.3' 'C.3' 'C.2' 'O.co2' 'O.co2' 'N.am' 'C.3' 'C.3' 'C.3'
     'C.3' 'C.3' 'C.3' 'C.3' 'C.3' 'C.3' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'
     'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'
     'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H']

    >>> print(new_mol.atomtype)
    ['C.2' 'O.2' 'C.3' 'N.am' 'C.3' 'C.2' 'O.2' 'C.3' 'C.ar' 'C.ar' 'C.ar'
     'C.ar' 'C.ar' 'C.ar' 'O.2' 'P.3' 'O.2' 'O.2' 'O.2' 'N.am' 'C.3' 'C.2'
     'O.2' 'C.3' 'C.3' 'C.2' 'O.co2' 'O.co2' 'N.am' 'C.3' 'C.3' 'C.3' 'C.3'
     'C.3' 'C.3' 'C.3' 'C.3' 'C.3' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'
     'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'
     'H' 'H' 'H' 'H' 'H' 'H']

    >>> print(mol.bondtype)
    ['2' '1' '1' '1' '1' '2' '1' 'ar' 'ar' 'ar' 'ar' '1' 'ar' 'ar' '1' 'ar'
     'ar' 'ar' '1' '1' '1' '1' '1' 'ar' 'ar' '2' '1' '1' '1' '1' '1' '1' '1'
     '1' '1' '1' 'am' 'am' 'am' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1'
     '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1'
     '1' '1' '1' '1' '1' '1']

    >>> print(new_mol.bondtype)
    ['2' '1' '1' '1' '1' '2' '1' 'ar' 'ar' 'ar' 'ar' '1' 'ar' 'ar' '1' '2' '1'
     '1' '1' '1' '1' '1' '1' 'ar' 'ar' '2' '1' '1' '1' '1' '1' '1' '1' '1' '1'
     '1' 'am' 'am' 'am' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1'
     '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1'
     '1' '1' '1' '1']

    >>> print(fitGasteigerCharges(new_mol).charge)
    [ 0.21703024 -0.2754467   0.01484334 -0.34422362  0.11166845  0.24307197
     -0.27269083  0.00167463 -0.04517365 -0.05525651 -0.05525651 -0.01933722
     -0.01933722  0.12386158 -0.4663472   0.11378469 -0.78024143 -0.78024143
     -0.31191045 -0.34249002  0.10829015  0.2446678  -0.27249503 -0.01843089
     -0.01066965  0.04147744 -0.55017024 -0.55017024 -0.34064123  0.02054308
      0.02335262 -0.03511205 -0.04707433 -0.02358754 -0.04876717 -0.05323909
     -0.0561604  -0.0653822   0.03302348  0.03302348  0.03302348  0.16380364
      0.06357127  0.03426235  0.03426235  0.06269906  0.06269906  0.06621166
      0.06621166  0.16388117  0.06323483  0.02967983  0.02967983  0.03294112
      0.03294112  0.04937867  0.04937867  0.04967035  0.04967035  0.02847681
      0.02847681  0.02706841  0.02706841  0.03208783  0.02700245  0.02700245
      0.02665417  0.02665417  0.02637094  0.02637094  0.02303533  0.02303533
      0.02303533]
    """

    molecule = molecule.copy()
    graph = _getMolecularGraph(molecule)

    for node in graph.nodes:

        # Skip not P atoms
        if graph.nodes[node]['element'] != 'P':
            continue

        # Sort neighbor atoms according to descending charge
        # Note: a double bond has to be near the most positive oxygen
        charge = lambda neighbor: molecule.charge[neighbor]
        neighbors = graph.neighbors(node)
        neighbors = sorted(neighbors, key=charge, reverse=True)

        num_double = 0

        for neighbor in neighbors:

            # Skip not O atoms
            if graph.nodes[neighbor]['element'] != 'O':
                continue

            # Change O type
            old = molecule.atomtype[neighbor]
            new = 'O.2'
            if old != new:
                molecule.atomtype[neighbor] = new
                logger.info('Change atom {} type: {} --> {}'.format(neighbor, old, new))

            # Get the current P--O bond type
            bond_index = graph.edges[(node, neighbor)]['index']
            old = molecule.bondtype[bond_index]

            # Get the new P--O bond type
            # Note: phosphate can have just one double bond
            num_bonds = len(list(graph.neighbors(neighbor)))
            if num_bonds == 2:
                new = '1'
            elif num_bonds == 1 and num_double == 0:
                new = '2'
                num_double += 1
            elif num_bonds == 1 and num_double == 1:
                new = '1'
            else:
                raise ValueError()

            # Change the P--O bond type
            if old != new:
                molecule.bondtype[bond_index] = new
                logger.info('Change bond {}--{} type: {} --> {}'.format(*molecule.bonds[bond_index], old, new))

    return molecule


if __name__ == '__main__':

    from moleculekit.molecule import Molecule
    from htmd.charge import fitGasteigerCharges

    import sys
    import doctest

    # Prevent HTMD importing inside doctest to fail if importing gives text output
    from htmd.home import home
    home()

    sys.exit(doctest.testmod().failed)