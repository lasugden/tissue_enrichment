from __future__ import annotations
import numpy as np
import pandas as pd


def load(path: str = 'data/gene_median_tpm.csv') -> TissueScores:
    """Load pre-computed edge weights

    Args:
        path (str, optional): Path to precomputed edge weights. Defaults to '../data/gene_median_tpm.csv'.

    Returns:
        pd.DataFrame: Loaded dataframe
    
    """
    ts = TissueScores(path)
    return ts


class TissueScores():
    def __init__(self, path: str = 'data/gene_median_tpm.csv'):
        """Load in tissue scores/pre-computed edge weights

        Args:
            path (str, optional): Path to existing tissue scores. Defaults to 'data/gene_median_tpm.csv'.

        """
        self.df = df = pd.read_csv(path)
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
        return 

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
