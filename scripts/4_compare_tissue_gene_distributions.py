import numpy as np
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.stats
from typing import Union

from tissue_enrichment import tscores


def fit_poisson_to_hist(ser: pd.Series, 
                        spacing: float = 0.05, 
                        return_pars: bool = False) -> Union[tuple[np.ndarray, np.ndarray],
                                                            tuple[float]]:
    """Fit a poisson function the histogram of a pandas series

    Args:
        ser (pd.Series): a series of gene expression scores in a tissue
        spacing (float, optional): the spacing of the poisson. Defaults to 0.1
        return_pars (bool, optional): if True, return parameters rather than fits.
            Defaults to False.

    Returns:
        Union[tuple[np.ndarray, np.ndarray],
              tuple[float]]: the x and y coordinates of the best fit
                                or
                             the best-fit parameters, if return_pars is True
    
    """
    def fit_lognorm(x, s, scale):
        return scipy.stats.lognorm.pdf(x, s, 0, scale)
    
    bins = np.arange(0, ser.max()+spacing*0.5, spacing)
    y, bins = np.histogram(ser.values, bins=bins, density=True)
    x = bins[:-1] + spacing/2
    pars = scipy.optimize.curve_fit(fit_lognorm, x, y, p0=[1, 1])

    if return_pars:
        return pars[0]
    else:
        fit_y = fit_lognorm(x, *pars[0])
        return x, fit_y


def compare_tissue_gene_dists(ts: tscores.TissueScores, save_path: str):
    """Draw a histogram for all tissues showing 20 random gene weight distributions

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which histograms should be saved
    
    """
    dist = ts.random_genes(20, keep_tissues=True)
    names = ts.tissues()
    df = pd.DataFrame(dist.transpose(), columns=names)

    for col in df.columns:
        x, y = fit_poisson_to_hist(df[col])
        ax = df[col].hist(bins=20, density=True)
        ax.plot(x, y)
        plt.savefig(os.path.join(save_path, f'gene_hist_20-{col}.pdf'))
        plt.clf()


def plot_best_lognorm_fits(ts: tscores.TissueScores, save_path: str):
    """Draw a histogram for all tissues showing 20 random gene weight distributions

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which histograms should be saved
    
    """
    dist = ts.random_genes(20, keep_tissues=True)
    names = ts.tissues()

    df = pd.DataFrame(dist.transpose(), columns=names)
    mu, sigma, legend = [], [], []
    for col in df.columns:
        pars = fit_poisson_to_hist(df[col], return_pars=True)

        sigma.append(pars[0])
        mu.append(np.log(pars[1]))
        legend.append(col)

    plt.scatter(mu, sigma)
    plt.xlabel('Mu')
    plt.ylabel('Sigma')
    for x, y, tissue in zip(mu, sigma, legend):
        plt.annotate(tissue, (x, y), fontsize=6)
    plt.savefig(os.path.join(save_path, 'gene_hist_20_lognormal_fits.pdf'))
    plt.clf()


def compare_tissue_gene_dists_all(ts: tscores.TissueScores, save_path: str):
    """Draw a histogram for all tissues showing gene weight distributions

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which histograms should be saved
    
    """
    dist = ts.feature_matrix()
    names = ts.tissues()
    df = pd.DataFrame(dist, columns=names)

    for col in df.columns:
        x, y = fit_poisson_to_hist(df[col], spacing=0.005)
        ax = df[col].hist(bins=100, density=True)
        ax.plot(x, y)
        plt.savefig(os.path.join(save_path, f'gene_hist_all-{col}.pdf'))
        plt.clf()


def plot_best_lognorm_fits_all(ts: tscores.TissueScores, save_path: str):
    """Draw a histogram for all tissues showing gene weight distributions

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which histograms should be saved
    
    """
    dist = ts.feature_matrix()
    names = ts.tissues()

    df = pd.DataFrame(dist, columns=names)
    mu, sigma, legend = [], [], []
    for col in df.columns:
        pars = fit_poisson_to_hist(df[col], spacing=0.005, return_pars=True)

        sigma.append(pars[0])
        mu.append(np.log(pars[1]))
        legend.append(col)

    plt.scatter(mu, sigma)
    plt.xlabel('Mu')
    plt.ylabel('Sigma')
    for x, y, tissue in zip(mu, sigma, legend):
        plt.annotate(tissue, (x, y), fontsize=6)
    plt.savefig(os.path.join(save_path, 'gene_hist_all_lognormal_fits.pdf'))


def plot_unique_genes_in_tissues(ts: tscores.TissueScores, save_path: str):
    """Plot the fraction of unique genes in each tissue as a bar plot

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which histograms should be saved
    
    """
    counts = ts.feature_matrix() > 0.99
    counts = pd.DataFrame(counts, columns=ts.tissues())

    counts = counts.sum(axis=0).sort_values(ascending=False)
    counts.plot.bar()
    plt.subplots_adjust(bottom=0.4)
    plt.savefig(os.path.join(save_path, 'tissue_unique_genes.pdf'))


def confirm_y_testis_transcripts(ts: tscores.TissueScores) -> float:
    """Plot the fraction of unique genes in each tissue as a bar plot

    Args:
        ts (tscores.TissueScores): the raw data

    Returns:
        float: the fraction of Y-chromosome transcripts that are uniquely in testis
    
    """
    with open('data/y_chrom_genes.txt') as fp:
        y_chromosome_genes = fp.read().split('\n')

    _, unfound = ts.search(y_chromosome_genes)
    print(f'{len(unfound)} unfound genes')

    subdf = ts.genes(y_chromosome_genes)
    found_genes = len(subdf)
    unique_genes = (subdf['testis'] > 0.99).sum()

    print(subdf)
    print(f'Y chromosome genes are unique in {unique_genes}/{found_genes} cases')
    return subdf.iloc[0]


def compare_testis_versus_reference_expression(ts: tscores.TissueScores, 
                                               save_path: str, 
                                               main_tissue: str = 'testis', 
                                               ref_tissue: str = 'pituitary',
                                               all: bool = False):
    """Plot testis versus reference tissue to look for global inflation

    Args:
        ts (tscores.TissueScores): the raw data
        save_path (str): the path into which scatterplot should be saved
        main_tissue (str, optional): tissue in question. Defaults to 'testis'.
        ref_tissue (str, optional): a particular reference tissue. Defaults to 'pituitary'.
        all (bool, optional): if True, show all genes. Defaults to False.
    
    """
    df = ts.random_genes_across_tissues(tissues=[main_tissue, ref_tissue], all=all)
    df.plot.scatter(x=main_tissue, y=ref_tissue, alpha=(0.05 if all else 0.4))
    plt.xscale('log')
    plt.yscale('log')
    plt.plot([0.00001, 1], [0.00001, 1])
    plt.savefig(os.path.join(save_path, f'tissue_expression_{main_tissue}_{ref_tissue}_{"all" if all else "1000"}.pdf'))


if __name__ == '__main__':
    ts = tscores.load(hugo_aliases=False)
    compare_tissue_gene_dists_all(ts, 'graphs')
    plot_best_lognorm_fits_all(ts, 'graphs')
    plot_unique_genes_in_tissues(ts, 'graphs')
    confirm_y_testis_transcripts(ts)
    compare_testis_versus_reference_expression(ts, 'graphs', all=True)
    compare_testis_versus_reference_expression(ts, 'graphs', all=True, main_tissue='bladder')

