#!/usr/local/bin/sage -python
# NOTE: you need SageMath installed, 
# you can run it with `sage -python scripts/plot_example_network.py`

from pathlib import Path
import sys
import math

from sage.all import Graph

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

from src.contact_networks import example_contact_network


def nx_to_sage_weighted(nxg, hide_weight_one=False):
    G = Graph()
    G.add_vertices(list(nxg.nodes()))

    for u, v, data in nxg.edges(data=True):
        w = data.get("weight", 1)
        gu = nxg.nodes[u].get("group")
        gv = nxg.nodes[v].get("group")

        if gu != gv:
            if hide_weight_one and w == 1:
                label = ""
            else:
                label = w
        else:
            label = ""
        G.add_edge(u, v, label)

    return G

def circle_positions(nodes, center=(0.0, 0.0), radius=1.0, angle=0.0):
    cx, cy = center
    n = len(nodes)
    pos = {}

    for k, v in enumerate(nodes):
        theta = angle + 2.0 * math.pi * k / n
        x = cx + radius * math.cos(theta)
        y = cy + radius * math.sin(theta)
        pos[v] = (x, y)

    return pos

def build_positions(nxg, angle_A=math.pi/2, angle_B=math.pi/6, clique_distance=7.0,
                    radius_A=1.2, radius_B=2.2):
    nodes_A = sorted([v for v, d in nxg.nodes(data=True) if d.get("group") == "A"])
    nodes_B = sorted([v for v, d in nxg.nodes(data=True) if d.get("group") == "B"])

    center_A = (-clique_distance / 2.0, 0.0)
    center_B = ( clique_distance / 2.0, 0.0)

    pos = {}
    pos.update(circle_positions(nodes_A, center=center_A, radius=radius_A, angle=angle_A))
    pos.update(circle_positions(nodes_B, center=center_B, radius=radius_B, angle=angle_B))

    return pos, nodes_A, nodes_B

if __name__ == "__main__":
    angle_A = math.pi / 3 
    angle_B = math.pi 
    
    clique_distance = 7.5

    radius_A = 1.2
    radius_B = 2.4

    hide_weight_one = True

    nxg = example_contact_network()
    G = nx_to_sage_weighted(nxg, hide_weight_one=hide_weight_one)

    pos, nodes_A, nodes_B = build_positions(
        nxg,
        angle_A=angle_A,
        angle_B=angle_B,
        clique_distance=clique_distance,
        radius_A=radius_A,
        radius_B=radius_B,
    )

    outpath = ROOT / "figures" / "contact_networks" / "g_weighted.svg"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    p = G.plot(
        pos=pos,
        vertex_labels=True,
        vertex_size=700,
        vertex_colors={
            "#8AAED4": nodes_A,
            "#E37E45": nodes_B,
        },
        edge_labels=True,
        edge_thickness=1.3,
        graph_border=False,
        figsize=[11, 6],
    )

    p.save(str(outpath))
    print(f"Saved to: {outpath}")

