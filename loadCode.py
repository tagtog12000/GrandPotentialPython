import matplotlib.pyplot as plt
import networkx as nx
import math

def draw_multidigraph(G_edges, radius=2.0, figsize=(4,4), node_color="orange"):
    # Build graph
    G = nx.MultiDiGraph()
    G.add_edges_from(G_edges)

    # Nodes
    nodes = sorted(set([u for u,v in G_edges] + [v for u,v in G_edges]))
    n = len(nodes)

    # Arrange nodes on a circle
    pos = {}
    for i, node in enumerate(nodes, start=1):
        angle = 2 * math.pi * (i-1) / n
        pos[node] = (radius * math.cos(angle), radius * math.sin(angle))

    # Count multiplicities
    edge_count = Counter(G_edges)

    # Track drawn edges
    drawn = set()

    plt.figure(figsize=figsize)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=250, node_color=node_color, edgecolors=node_color)

    # Draw edges with logic
    for (u, v), count in edge_count.items():
        if (u, v) in drawn:
            continue

        opp = (v, u)
        if opp in edge_count:
            # Opposite directed edges → bend left/right
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                   arrowstyle="->", arrowsize=20,
                                   connectionstyle="arc3,rad=0.2")
            nx.draw_networkx_edges(G, pos, edgelist=[(v, u)],
                                   arrowstyle="->", arrowsize=20,
                                   connectionstyle="arc3,rad=0.2")
            drawn.add((u, v))
            drawn.add((v, u))
        elif count == 1:
            # Single edge → straight
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                   arrowstyle="->", arrowsize=20)
            drawn.add((u, v))
        elif count == 2:
            # Two parallel edges same direction → bend left/right
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                   arrowstyle="->", arrowsize=20,
                                   connectionstyle="arc3,rad=0.2")
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                   arrowstyle="->", arrowsize=20,
                                   connectionstyle="arc3,rad=-0.2")
            drawn.add((u, v))
        else:
            # More than 2 parallels → spread bends
            k = count
            rads = [0] if k % 2 == 1 else []
            for i in range(1, (k//2)+1):
                rads.extend([0.2*i, -0.2*i])
            for r in rads[:k]:
                if r == 0:
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                           arrowstyle="->", arrowsize=20)
                else:
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                           arrowstyle="->", arrowsize=20,
                                           connectionstyle=f"arc3,rad={r}")
            drawn.add((u, v))

    plt.axis("off")
    plt.show()
