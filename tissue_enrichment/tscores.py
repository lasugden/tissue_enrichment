from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Iterable


def load(*args, **kwargs) -> TissueScores:
    """Load pre-computed edge weights

    Args:
        path (str, optional): Path to precomputed edge weights. Defaults to '../data/gene_median_tpm.csv'.
        combine_skin (bool, optional): If True, combine two skin tissues as identified in pre-analysis

    Returns:
        TissueScores: Loaded new instance of TissueScores class
    
    """
    ts = TissueScores(*args, **kwargs)
    return ts


class TissueScores():
    def __init__(self, 
                 path: str = 'data/gene_median_tpm.csv', 
                 combine_skin: bool = True,
                 alias: bool = True):
        """Load in tissue scores/pre-computed edge weights

        Args:
            path (str, optional): Path to existing tissue scores. Defaults to 'data/gene_median_tpm.csv'.
            combine_skin (bool, optional): If True, combine two skin tissues as identified in pre-analysis
            alias (bool, optional): If True, load aliases. Defaults to True.

        """
        self.df = pd.read_csv(path)

        if combine_skin:
            self.df['skin'] = self.df[['skin.sun_exposed', 'skin.not_sun_exposed']].mean(axis=1)
            self.df = self.df.drop(columns=['skin.sun_exposed', 'skin.not_sun_exposed'])

        if alias:
            pass

        self.df['description'] = self.df['description'].str.lower()
        self.gene_list = self.df['description'].values
        self.numeric_columns = self.df.select_dtypes(exclude='object').columns.drop('tot_tpm')

    def search(self, genes: list[str], raise_error: bool = False) -> tuple[list[str], list[str]]:
        """Search for a list of gene names

        Args:
            genes (list[str]): list of gene names
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

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
    
    def genes(self, genes: list[str], raise_error: bool = False) -> pd.DataFrame:
        """Return a dataframe containing data only for a list of genes

        Args:
            genes (list[str]): a list of genes to include
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: a dataframe containing only those genes
        
        """
        genes, _ = self.search(genes, raise_error=raise_error)
        assert len(genes) > 0, 'No genes found'
        return self.df.loc[self.df['description'].isin(genes), :].copy()

    def gene_weights(self, gene_weights: list[tuple[str, float]], raise_error: bool = False) -> pd.DataFrame:
        """Return a dataframe containing data only for a list of genes scaled by weights

        Args:
            gene_weights (list[tuple[str, float]]): a list of genes to include
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: a dataframe containing only those genes, scaled by weights
        
        """
        genes, _ = self.search([gw[0].lower() for gw in gene_weights], raise_error=raise_error)
        gene_weights = list(map({gw[0].lower(): gw for gw in gene_weights}.get, genes))

        assert len(genes) > 0, 'No genes found'
        assert max(gene_weights) <= 1, 'Gene weights greater than zero'
        assert min(gene_weights) >= 0, 'Gene weights less than zero'

        subdf = self.df.loc[self.df['description'].isin(genes), :].copy()
        subdf.loc[:, self.numeric_columns] = subdf[self.numeric_columns].multiply(
            subdf['description'].map({gw[0].lower(): gw[1] for gw in gene_weights}), axis=0)
        subdf = subdf[self.numeric_columns].sum(axis=0)
        subdf = subdf.sort_values(ascending=False)

        return subdf

    def feature_matrix(self) -> np.ndarray:
        """Return a feature matrix

        Returns:
            np.ndarray: a 2d feature matrix of genes, tissues

        """
        return self.df[self.numeric_columns].values

    def tissues(self) -> list[str]:
        """Return all tissue names in order

        Returns:
            list[str]: list of tissue names
        
        """
        return list(self.numeric_columns)

    def random_genes(self, 
                     n: int, 
                     samples: int = 1_000, 
                     keep_tissues: bool = False) -> np.ndarray:
        """Create a random sample of weights

        Args:
            n (int): number of genes to include
            samples (int, optional): number of random samples to compute. Defaults to 1000.
            keep_tissues (bool, optional): If True, return results for all tissues. Defaults to False.

        Returns:
            np.ndarray: array of top scores as floats
        
        """
        weights = [1]*n
        return self.random_gene_weights(weights, samples, keep_tissues)

    def random_gene_weights(self, 
                            weights: Iterable[float], 
                            samples: int = 1_000,
                            keep_tissues: bool = False) -> np.ndarray:
        """Create a random sample of weights

        Args:
            weights (Iterable[float]): weights of genes to apply
            samples (int): number of random samples to compute
            keep_tissues (bool, optional): If True, return results for all tissues. Defaults to False.

        Returns:
            np.ndarray: array of top scores as floats
        
        """
        rdf = self.df.sample(n=len(weights)*samples, replace=True)[self.numeric_columns].values
        rdf = rdf.transpose()*np.repeat(weights, samples)
        tops = rdf.reshape((rdf.shape[0], len(weights), -1)).sum(axis=1)
        if keep_tissues:
            return tops
        else:
            return tops.max(axis=0)
