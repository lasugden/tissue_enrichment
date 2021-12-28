import numpy as np
import pandas as pd


def load(path: str = 'data/gene_median_tpm.csv') -> pd.DataFrame:
    """Load pre-computed edge weights

    Args:
        path (str, optional): Path to precomputed edge weights. Defaults to '../data/gene_median_tpm.csv'.

    Returns:
        pd.DataFrame: Loaded dataframe
    
    """
    df = pd.read_csv(path)
    return df


class TissueWeighter():
    """Compute the tissue specificity of inputs"""
    def __init__(self, weights: pd.DataFrame):
        """Initialize with precomputed weights

        Args:
            weights (pd.DataFrame): precomputed weights
        """
        self.df = weights
        self.df['description'] = self.df['description'].str.lower()
        self.gene_list = self.df['description'].values
        self.numeric_columns = self.df.select_dtypes(exclude='object').columns.drop('tot_tpm')

    def _gene_search(self, genes: list[str], raise_error: bool = False) -> tuple[list[str], list[str]]:
        """Search for a list of gene names

        Args:
            genes (list[str]): list of gene names
            raise_error (bool, optional): raise an error if not found. Defaults to False.

        Returns:
            tuple[list[str], list[str]]: tuple of list of found gene names and unfound gene names
        
        """
        found, unfound = [], []
        for gene in genes:
            if gene.lower() in self.gene_list:
                found.append(gene.lower())
            else:
                if raise_error:
                    raise ValueError(f'Gene name {gene} not found')
                unfound.append(gene)

        return found, unfound

    def gene_tissue_scores(self, gene_weights: list[tuple[str, float]]) -> pd.Series:
        """Compute the tissue scores given genes and gene weights

        Args:
            gene_weights (list[tuple[str, float]]): list of tuples of gene names and weights

        Returns:
            pd.Series: scores per tissue
        
        """
        subdf = self.df.loc[self.df['description'].isin([gw[0] for gw in gene_weights]), :].copy()
        scaled = subdf[self.numeric_columns].multiply(
            subdf['description'].map({gw[0].lower(): gw[1] for gw in gene_weights}), axis=0)
        subdf.loc[:, self.numeric_columns] = scaled
        subdf = subdf[self.numeric_columns].sum(axis=0)
        subdf = subdf.sort_values(ascending=False)
        return subdf

    def random_sample_by_weights(self, weights: list[float], samples: int) -> np.ndarray:
        """Create a random sample of weights

        Args:
            weights (list[float]): list of gene weights as floats
            samples (int): number of random samples to compute

        Returns:
            np.ndarray: array of top scores as floats
        
        """
        rdf = self.df.sample(n=len(weights)*samples, replace=True)[self.numeric_columns].values
        rdf = rdf.transpose()*np.repeat(weights, samples)
        tops = rdf.reshape((rdf.shape[0], len(weights), -1)).sum(axis=1).max(axis=0)
        return tops
        
    def tissues_by_weight(self, 
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
    df = load()
    tw = TissueWeighter(df)
    tw.tissues_by_weight([('hla-a', 0.5), 
                          ('hla-b', 0.6),
                          ('hla-c', 1), 
                          ('hla-dpa1', 1), 
                          ('hla-dpb1', 1), 
                          ('hla-dqa1', 1),
                          ('hla-dqb1', 1), 
                          ('hla-dra', 1), 
                          ('hla-drb1', 1)])
