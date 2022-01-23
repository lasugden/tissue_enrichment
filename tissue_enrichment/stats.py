import numpy as np
import pandas as pd

from tissue_enrichment import tscores


def tissue_p(ts: tscores.TissueScores,
             gene_weights: list[tuple[str, float]], 
             samples: int = 1_000) -> list[tuple[str, float]]:
    """Search with genes and weights a list of tissues and scores

    Args:
        gene_weights (list[tuple[str, float]]): list of tuples of gene names and weights
        samples (int, optional): number of random samples to compare to. Defaults to 1000.

    Returns:
        list[tuple[str, float]]: a list of tuples of tissue names and scores
    
    """
    

    genes, _ = self._gene_search([v[0] for v in gene_weights])
    gene_weights = list(map({gw[0].lower(): gw for gw in gene_weights}.get, genes))

    assert(len(genes) > 0, 'No genes found')
    assert(max(gene_weights) <= 1, 'Gene weights greater than zero')
    assert(min(gene_weights) >= 0, 'Gene weights less than zero')

    real = self.gene_tissue_scores(gene_weights).to_frame(name='score')

    weights = [gw[1] for gw in gene_weights]
    top_scores = self.random_sample_by_weights(weights, samples)

    real['pvalue'] = real['score'].apply(lambda x: (top_scores > x).sum()/samples)
    return real


if __name__ == '__main__':
    ts = tscores.load(alias=False)
    tw.tissues_by_weight([('hla-a', 0.5), 
                          ('hla-b', 0.6),
                          ('hla-c', 1), 
                          ('hla-dpa1', 1), 
                          ('hla-dpb1', 1), 
                          ('hla-dqa1', 1),
                          ('hla-dqb1', 1), 
                          ('hla-dra', 1), 
                          ('hla-drb1', 1)])
