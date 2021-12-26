import matplotlib.pyplot as plt
import umap

from tissue_enrichment import tscores


def tissue_umap_plot(save_path: str):
    """Plot UMAP results of tissue scores

    Args:
        save_path (str): location to save plot
    
    """
    ts = tscores.load()
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(ts.feature_matrix().transpose())
    plt.scatter(embedding[:, 0], embedding[:, 1])
    for i, tissue in enumerate(ts.tissues()):
        plt.annotate(tissue, (embedding[i, 0], embedding[i, 1]), fontsize=6)
    plt.savefig(save_path)


if __name__ == '__main__':
    # tissue_umap_plot('graphs/tissue_umap_2.pdf')
