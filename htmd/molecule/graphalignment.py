# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import networkx as nx
import numpy as np


def createProductGraph(G, H, tolerance, fields):
    # Calculate product graph by creating a node for each feature-matching pair of points
    newnodes = []
    for gn in G.nodes():
        for hn in H.nodes():
            matching = True
            for f in fields:
                if G.node[gn][f] != H.node[hn][f]:
                    matching = False
                    break
            if matching:
                newnodes.append((gn, hn))

    Gprod = nx.Graph()
    Gprod.add_nodes_from(newnodes)

    # Add edges when distances in both graphs agree within tolerance error between the two given nodes
    for np1 in range(len(newnodes)):
        for np2 in range(np1 + 1, len(newnodes)):
            pair1 = newnodes[np1]
            pair2 = newnodes[np2]
            if not G.has_edge(pair1[0], pair2[0]) or not H.has_edge(pair1[1], pair2[1]):
                continue
            dist1 = G.edges[pair1[0], pair2[0]]['distance']
            dist2 = H.edges[pair1[1], pair2[1]]['distance']
            if abs(dist1 - dist2) < tolerance:
                Gprod.add_edge(newnodes[np1], newnodes[np2])
    return Gprod


def compareGraphs(G, H, fields=('element',), tolerance=0.5, returnmatching=False):
    # Comparison algorithm based on:
    # "Chemoisosterism in the Proteome", X. Jalencas, J. Mestres, JCIM 2013
    # http://pubs.acs.org/doi/full/10.1021/ci3002974
    if G == H:
        if returnmatching:
            return 1, len(G), [(x, x) for x in G.nodes()]
        else:
            return 1

    if len(G.edges()) == 0 or len(H.edges()) == 0:
        if returnmatching:
            return 0, 0, []
        else:
            return 0

    Gprod = createProductGraph(G, H, tolerance, fields)

    # Calculate the maximal cliques and return the length of the largest one
    maxcliques = np.array(list(nx.find_cliques(Gprod)))
    cllen = np.array([len(x) for x in maxcliques])
    score = (cllen.max() / max(len(G.nodes()), len(H.nodes())))

    if returnmatching:
        return score, cllen.max(), maxcliques[cllen.argmax()]
    else:
        return score


def makeMolGraph(mol, sel, fields):
    from scipy.spatial.distance import pdist, squareform

    if sel != 'all':
        sel = mol.atomselect(sel, indexes=True)
    else:
        sel = np.arange(mol.numAtoms)

    g = nx.Graph()
    for i in sel:
        props = {f: mol.__dict__[f][i] for f in fields}
        g.add_node(i, **props)

    distances = squareform(pdist(mol.coords[sel, :, mol.frame]))
    nodes = list(g.nodes())
    for i in range(len(g)):
        for j in range(i+1, len(g)):
            g.add_edge(nodes[i], nodes[j], distance=distances[i, j])

    return g


def maximalSubstructureAlignment(mol1, mol2, sel1='all', sel2='all', fields=('element',), tolerance=0.5, visualize=False):
    """ Aligns two molecules on the largest common substructure

    Parameters
    ----------
    mol1 : :class:`Molecule`
        The reference molecule on which to align
    mol2 : :class:`Molecule`
        The second molecule which will be rotated and translated to align on mol1
    sel1 : str
        Atom selection string of the atoms of `mol1` to align.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    sel2 : str
        Atom selection string of the atoms of `mol2` to align.
        See more `here <http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node89.html>`__
    fields : tuple
        A tuple of the fields that are used to match atoms
    tolerance : float
        How different can distances be between to atom pairs for them to match in the product graph
    visualize : bool
        If set to True it will visualize the alignment

    Returns
    -------
    newmol : :class:`Molecule`
        A copy of mol2 aligned on mol1
    """
    mol2 = mol2.copy()
    g1 = makeMolGraph(mol1, sel1, fields)
    g2 = makeMolGraph(mol2, sel2, fields)

    _, _, matching = compareGraphs(g1, g2, fields=fields, tolerance=tolerance, returnmatching=True)

    matchnodes1 = [x[0] for x in matching]
    matchnodes2 = [x[1] for x in matching]

    mol2.align(sel=matchnodes2, refmol=mol1, refsel=matchnodes1)

    if visualize:
        mol1.view(sel='index {}'.format(' '.join(map(str, matchnodes1))), style='CPK', hold=True)
        mol1.view(sel='all', style='Lines')

        mol2.view(sel='index {}'.format(' '.join(map(str, matchnodes2))), style='CPK', hold=True)
        mol2.view(sel='all', style='Lines')

    return mol2