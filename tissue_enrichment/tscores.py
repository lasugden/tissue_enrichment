from __future__ import annotations
import enum
import numpy as np
import pandas as pd
from typing import Iterable, Union


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
                 hugo_aliases: bool = True,
                 alias_path: str = 'data/idx.csv',
                 info_path: str = 'data/locus_group.csv'):
        """Load in tissue scores/pre-computed edge weights

        Args:
            path (str, optional): Path to existing tissue scores. Defaults to 'data/gene_median_tpm.csv'.
            combine_skin (bool, optional): If True, combine two skin tissues as identified in pre-analysis
            hugo_aliases (bool, optional): If True, load aliases. Defaults to True.
            alias_path (str, optional): Location of the preprocessed Hugo alias dataframe
            info_path (str, optional): Location of the preprocessed Hugo locus group dataframe

        """
        self.df = pd.read_csv(path).rename(columns={'name': 'ensembl_gene_id', 'description': 'gtex_alias'})
        self.info = None
        self._last_gw = None

        if combine_skin:
            self.df['skin'] = self.df[['skin.sun_exposed', 'skin.not_sun_exposed']].mean(axis=1)
            self.df = self.df.drop(columns=['skin.sun_exposed', 'skin.not_sun_exposed'])

        if hugo_aliases:
            self.aliases = pd.read_csv(alias_path)
            self.df = self.df.merge(pd.read_csv(info_path), how='left', on='ensembl_gene_id')
        else:
            self.aliases = (self.df[['gtex_alias', 'ensembl_gene_id']]
                                .rename(columns={'gtex_alias': 'alias'})
                                .dropna())
            self.aliases['alias'] = self.aliases['alias'].str.lower()

        self.numeric_columns = self.df.select_dtypes(exclude='object').columns.drop('tot_tpm')

    def search(self, genes: list[str], raise_error: bool = False) -> tuple[list[str], list[str]]:
        """Search for a list of gene names

        Args:
            genes (list[str]): list of gene names
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            tuple[list[str], list[str]]: tuple of list of found gene names and unfound gene names
        
        """
        genes = [v.lower() for v in genes]
        found = self.aliases[self.aliases['alias'].isin(genes)]['alias'].to_list()

        if raise_error and len(found) < len(genes):
            missing = list(set(genes).difference(set(found)))
            raise ValueError(f'Gene name(s) {missing} not found')

        unfound = list(set(genes).difference(set(found)))
        return found, unfound

    def search_df(self, genes: pd.DataFrame, raise_error: bool = False) -> pd.DataFrame:
        """Search for a list of gene names and return a dataframe for merging

        Args:
            genes (pd.DataFrame): column of gene names with 'gene' as merge key. Can include other columns.
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: merged dataframe of found ensembl_gene_ids
        
        """
        genes['gene'] = genes['gene'].str.lower()
        found = self.aliases.merge(genes, how='inner', left_on='alias', right_on='gene')

        if raise_error and len(found) < len(genes):
            missing = set(genes['gene'].to_list()).difference(set(found['alias'].to_list()))
            raise ValueError(f'Gene name(s) {missing} not found')

        self._last_gw = found
        return found

    @staticmethod
    def _format_genes(genes: Union[list[str], list[tuple[str, float]]]) -> pd.DataFrame:
        """Convert a list of genes or a list of genes and weights to a dataframe for merging

        Args:
            genes (Union[list[str], list[tuple[str, float]]]): list of genes or list of tuples of genes and weights

        Returns:
            pd.DataFrame: output dataframe with 2 columns of gene and weight
        
        """
        if not isinstance(genes[0], tuple):
            genes = [[v, 1.0] for v in genes]
            
        return pd.DataFrame(genes, columns=['gene', 'weight'])
    
    def genes(self, genes: list[str], raise_error: bool = False) -> pd.DataFrame:
        """Return a dataframe containing data only for a list of genes

        Args:
            genes (list[str]): a list of genes to include
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: a dataframe containing only those genes
        
        """
        gene_df = self.search_df(self._format_genes(genes), raise_error=raise_error)
        assert len(gene_df) > 0, 'No genes found'
        return self.df.merge(gene_df, how='inner', on='ensembl_gene_id')

    def _gene_weights_df(self, gene_df: pd.DataFrame):
        """Run gene weights on gene_df, for use in fast iteration.
        Skips lookup.

        Args:
            gene_weights (pd.DataFrame): dataframe of columns gene, weights
        
        Returns:
            pd.DataFrame: a dataframe containing only those genes, scaled by weights

        """
        subdf = self.df.merge(gene_df, how='inner', on='ensembl_gene_id')
        subdf.loc[:, self.numeric_columns] = subdf[self.numeric_columns].multiply(subdf['weight'], axis=0)
        subdf = subdf[self.numeric_columns].sum(axis=0)
        subdf = subdf.sort_values(ascending=False)

        return subdf

    def gene_weights(self, gene_weights: list[tuple[str, float]], raise_error: bool = False) -> pd.DataFrame:
        """Return a dataframe containing data only for a list of genes scaled by weights

        Args:
            gene_weights (list[tuple[str, float]]): a list of genes to include
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: a dataframe containing only those genes, scaled by weights
        
        """
        gene_df = self.search_df(self._format_genes(gene_weights), raise_error=raise_error)

        assert len(gene_df) > 0, 'No genes found'
        assert gene_df['weight'].max() <= 1, 'Gene weights greater than zero'
        assert gene_df['weight'].min() >= 0, 'Gene weights less than zero'

        return self._gene_weights_df(gene_df)

    def gene_weights_p(self,
                       gene_weights: list[tuple[str, float]], 
                       samples: int = 1_000,
                       within_group: bool = True,
                       raise_error: bool = False) -> pd.DataFrame:
        """Search with genes and weights a list of tissues and scores

        Args:
            gene_weights (list[tuple[str, float]]): list of tuples of gene names and weights
            samples (int, optional): number of random samples to compare to. Defaults to 1000.
            within_group (bool, optional): if True, randomize only within locus_group as defined by hugo
            raise_error (bool, optional): If True, raise an error if at least one gene is not found. Defaults to False.

        Returns:
            pd.DataFrame: return a dataframe of scores and pvalues
        
        """
        real = self.gene_weights(gene_weights, raise_error=raise_error).rename('score').to_frame()
        gene_weights_loci = self._last_gw.merge(self.df, how='left', on='ensembl_gene_id')
        locus_groups = None if not within_group else gene_weights_loci['locus_group'].to_list()
        rand = self.random_gene_weights(gene_weights_loci['weight'].to_list(), 
                                        locus_groups=locus_groups, 
                                        samples=samples)
        real['pvalue'] = real['score'].apply(lambda x: (rand > x).sum()/samples)
        return real

    # def optimal_tissue(self,
    #                    gene_weights: list[tuple[str, float]],
    #                    samples: int = 1_0000,
    #                    within_group: bool = True,
    #                    raise_error: bool = False):

    #     real = self.gene_weights(gene_weights, raise_error=raise_error).rename('score').to_frame()
    #     real['inv_rank'] = 1.0/(np.arange(len(real)) + 1)
    #     tissue_scores = [real['score'].iloc[0]]

    #     gw = self._last_gw
    #     while len(gw) > 1:
    #         subdf = gw.merge(self.df, how='left', on='ensembl_gene_id')
    #         subdf.loc[:, self.numeric_columns] = subdf[self.numeric_columns].multiply(subdf['weight'], axis=0)*real['inv_rank']
    #         subdf['gene_tissue_score'] = subdf[self.numeric_columns].sum(axis=1)
            
    #         gw = gw.drop(subdf['gene_tissue_score'].idxmin()).reset_index(drop=True)
    #         real = self._gene_weights_df(gw).rename('score').to_frame()
    #         real['inv_rank'] = 1.0/(np.arange(len(real)) + 1)
    #         tissue_scores.append(real['score'].iloc[0])

    #     subdf = self.df.merge(gene_df, how='inner', on='ensembl_gene_id')
    #     print(1)

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

    def tissue_distances(self) -> np.ndarray:
        """Return the Euclidean distances between each pair of tissues as a matrix in order of tissues function

        Returns:
            np.ndarray: tissue x tissue matrix of distances
        
        """
        out = np.ones((len(self.numeric_columns), len(self.numeric_columns)))
        for i, tissue1 in enumerate(self.numeric_columns):
            for j, tissue2 in enumerate(self.numeric_columns[i + 1:]):
                out[i, j] = np.linalg.norm(self.df[tissue1].values*self.df['tot_tpm'].values - 
                                           self.df[tissue2].values*self.df['tot_tpm'].values)
                out[j, i] = out[i, j]
        return out

    def ids_and_aliases(self) -> pd.DataFrame:
        """Return ensembl_gene_ids and aliases for integration with other datasets such as Hugo

        Returns:
            pd.DataFrame: A dataframe with two columns, ensembl_gene_id and gtex_alias
        
        """
        return self.df[['ensembl_gene_id', 'gtex_alias']]

    def random_genes_across_tissues(self, n: int = 1000, tissues: Iterable[str] = ('testis', 'pituitary'), all: bool = False) -> pd.DataFrame:
        """Randomly sample genes for specific tissues

        Args:
            n (int, optional): Number of genes to sample. Defaults to 1000.
            tissues (Iterable[str], optional): List of tissues to return. Defaults to ('testis', 'pituitary').
            all (bool, optional): if True, show all genes. Defaults to False.

        Returns:
            pd.DataFrame: Dataframe containing columns of the names of tissues
        
        """
        if all:
            return self.df[tissues]
        else:
            return self.df[tissues].sample(n)

    def random_genes(self, 
                     n: int, 
                     locus_groups: Iterable[str] = None,
                     samples: int = 1_000, 
                     keep_tissues: bool = False) -> np.ndarray:
        """Create a random sample of weights

        Args:
            n (int): number of genes to include
            samples (int, optional): number of random samples to compute. Defaults to 1000.
            locus_groups (str, optional): locus group to randomize within. If None,
                randomize within all locus_groups
            keep_tissues (bool, optional): If True, return results for all tissues. Defaults to False.

        Returns:
            np.ndarray: array of top scores as floats
        
        """
        weights = [1]*n
        return self.random_gene_weights(weights, locus_groups, samples, keep_tissues)

    def random_gene_weights(self, 
                            weights: Iterable[float], 
                            locus_groups: Iterable[str] = None,
                            samples: int = 1_000,
                            keep_tissues: bool = False) -> np.ndarray:
        """Create a random sample of weights

        Args:
            weights (Iterable[float]): weights of genes to apply
            locus_group (str, optional): locus group to randomize within. If None,
                randomize within all locus_groups
            samples (int): number of random samples to compute
            keep_tissues (bool, optional): If True, return results for all tissues. Defaults to False.

        Returns:
            np.ndarray: array of top scores as floats
        
        """
        if locus_groups is None:
            rdf = self.df.sample(n=len(weights)*samples, replace=True)[self.numeric_columns].values
        else:
            lgs = pd.DataFrame({'weights': weights, 'lg': locus_groups})
            lgs = lgs.groupby('lg').agg(list)
            
            rdfs = []
            for lg, weights in lgs['weights'].iteritems():
                rdfs.append(self.df
                                .loc[self.df['locus_group'] == lg, :]
                                .sample(n=len(weights)*samples, replace=True)[self.numeric_columns]
                                .values)
            rdf = np.concatenate(rdfs, axis=0)
        
        rdf = rdf.transpose()*np.repeat(weights, samples)
        tops = rdf.reshape((rdf.shape[0], len(weights), -1)).sum(axis=1)
        if keep_tissues:
            return tops
        else:
            return tops.max(axis=0)


if __name__ == '__main__':
    ts = load()
    ts.optimal_tissue([('hla-a', 0.5), 
                          ('hla-b', 0.6),
                          ('hla-c', 1), 
                          ('hla-dpa1', 1), 
                          ('hla-dpb1', 1), 
                          ('hla-dqa1', 1),
                          ('hla-dqb1', 1), 
                          ('hla-dra', 1), 
                          ('hla-drb1', 1)])
