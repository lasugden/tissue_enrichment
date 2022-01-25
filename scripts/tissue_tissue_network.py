import matplotlib.pyplot as plt
import networkx as nx

from tissue_enrichment import tscores


def plot_tissues(ts: tscores.TissueScores):
    """Use networkx to plot tissue-tissue distances.
    NOTE: Does not work well due to spring_layout not being a good form of dimensionality reduction

    Args:
        ts (tscores.TissueScores): raw data
    
    """
    tissues = ts.tissues()
    tdists = ts.tissue_distances()

    gx = nx.Graph()
    gx.add_nodes_from(tissues)
    for i, tissue1 in enumerate(tissues):
        for j, tissue2 in enumerate(tissues[i + 1:]):
            gx.add_edge(tissue1, tissue2, weight=tdists[i, j])

    pos = nx.spring_layout(gx)
    nx.draw_networkx_nodes(gx, pos=pos)
    nx.draw_networkx_labels(gx, pos=pos)
    plt.show()


if __name__ == '__main__':
    ts = tscores.load(hugo_aliases=False)
    plot_tissues(ts)
