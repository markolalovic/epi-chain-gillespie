# src/contact_networks.py

import itertools
import networkx as nx

def example_contact_network():
    """Example weighted contact network with bipartition A, B and cliques on A and B. 
    - See `figures/contact_networks/g_weighted.svg`
    """
    g = nx.Graph()

    nodes_A = ["a1", "a2", "a3"]
    nodes_B = ["b1", "b2", "b3", "b4", "b5", "b6"]

    g.add_nodes_from(nodes_A, group="A", bipartite=0)
    g.add_nodes_from(nodes_B, group="B", bipartite=1)

    # weighted edges between A and B, weight = number of edges = repeated contacts
    edges_AB = [
        ("a1", "b1", 1), ("a1", "b3", 2), ("a1", "b5", 1),
        ("a2", "b2", 1), ("a2", "b4", 6), ("a2", "b5", 1), ("a2", "b6", 5),
        ("a3", "b1", 1), ("a3", "b3", 2), ("a3", "b5", 2), ("a3", "b6", 2),
    ]
    g.add_weighted_edges_from(edges_AB, weight="weight")

    # cliques induced by A and B, with default weight 1
    g.add_weighted_edges_from((u, v, 1) for u, v in itertools.combinations(nodes_A, 2))
    g.add_weighted_edges_from((u, v, 1) for u, v in itertools.combinations(nodes_B, 2))

    return g

