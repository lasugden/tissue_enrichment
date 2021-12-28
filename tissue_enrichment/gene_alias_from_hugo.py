import pandas as pd

from tissue_enrichment import tscores


def read_hugo(path: str, ts: tscores.TissueScores) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read in Hugo file, combine with tissue scores, and return matching and non-matching groups

    Args:
        path (str): path to hugo file
        ts (tscores.TissueScores): instantiated to be used for identifying overlaps

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Dataframe of overlaps with hugo,
                                           Dataframe of non-overlapping entries
    """
    df = (pd.read_csv(path, delimiter='\t')
            .dropna(subset=['ensembl_gene_id']))
    
    # Combine duplicates
    dup_genes = df['ensembl_gene_id'].duplicated(keep=False)
    dups = df.loc[dup_genes, :].dropna(axis=1, how='all').fillna('')  # Drop columns with all NaNs
    dups = (dups.groupby('ensembl_gene_id')
                .agg(lambda x: '|'.join([v for v in x if len(v) > 0]) if (x != x.values[0]).sum() > 0 else x.values[0])
                .reset_index())
    df = pd.concat([df.loc[~dup_genes, :], dups]).reset_index(drop=True)

    # Exctract gene names from Gtex DB for merging
    gene_names = ts.df[['name', 'description']].rename(columns={'name': 'gtex_name', 'description': 'gtex_alias'})
    gene_names['gtex_name'] = gene_names['gtex_name'].str.replace(r'\..*', '', regex=True)
    match = df.merge(gene_names['gtex_name'], how='inner', left_on='ensembl_gene_id', right_on='gtex_name')

    # Identify non-matching entries
    non_match = gene_names.merge(df['ensembl_gene_id'], 
                                 how='left', 
                                 right_on='ensembl_gene_id', 
                                 left_on='gtex_name', 
                                 indicator=True)
    non_match = non_match[(non_match._merge != 'both')]

    return match.reset_index(drop=True), non_match[['gtex_name', 'gtex_alias']].reset_index(drop=True)


def create_idx(df: pd.DataFrame, non_match: pd.DataFrame) -> pd.DataFrame:
    """Create index as Pandas dataframe matching all aliases to ensembl_gene_ids for gtex

    Args:
        df (pd.DataFrame): dataframe of matches between hugo and gtex
        non_match (pd.DataFrame): dataframe of entries only in gtex

    Returns:
        pd.DataFrame: two-column dataframe index of aliases and ensembl_gene_ids

    """
    # Iterate through alias columns ensuring that no duplicates are allowed
    ser = pd.Series([[] for _ in range(len(df))], index=df.index)
    name_ser = df['symbol'].str.lower().str.split('|')
    name_set = name_ser.sum()

    # Iteration ensures column-by-column precedence for symbols
    for col in ['alias_symbol', 'prev_symbol']:
        alias_ser = df[col].str.lower().str.split('|')
        alias_ser[alias_ser.notna()] = alias_ser[alias_ser.notna()].apply(lambda l: list(set(l).difference(name_set)))
        name_ser += alias_ser.fillna(ser)
        name_set = name_ser.sum()

    # Drop any remaining duplicates and create the index
    aliases = sorted([alias for aliases in name_ser for alias in aliases])
    duplicates = set([v for v, vnext in zip(aliases[:-1], aliases[1:]) if v == vnext])
    idx = {key: list(set(val).difference(duplicates)) + [key] 
           for key, val in zip(df['ensembl_gene_id'], name_ser)}

    # Add in indices of those values not found in Hugo
    non_match = non_match.dropna()
    idx.update({key: [val, key] if val not in name_set else [key] 
                for key, val in zip(non_match['gtex_name'], non_match['gtex_alias'])})

    # Invert index and return as dataframe
    return pd.DataFrame([[val, key] for key, vals in idx.items() for val in vals], columns=['alias', 'ensembl_gene_id'])


def save_idx_lgroup(idx: pd.DataFrame, 
                    match: pd.DataFrame, 
                    non_match: pd.DataFrame, 
                    idx_path: str, 
                    info_path: str):
    """Save the name and locus group for each gene

    Args:
        idx (pd.DataFrame): fully created index mapping aliases to ensembl_gene_ids
        match (pd.DataFrame): dataframe containing matches between gtex and hugo
        non_match (pd.DataFrame): dataframe containing gtex-only entries
        idx_path (str): path into which the alias index should be saved
        info_path (str): path into which the additional info should be saved
    
    """
    # Combine non-match info into matches
    for col in ['locus_group', 'locus_type', 'gene_group']:
        non_match[col] = 'uncategorized'

    df = pd.concat([match[['ensembl_gene_id', 'name', 'locus_group', 'locus_type', 'gene_group']],
                    non_match.rename(columns={'gtex_name': 'ensembl_gene_id', 'gtex_alias': 'name'})])
    df.to_csv(info_path, index=False)

    idx.sort_values('alias').to_csv(idx_path, index=False)
    

if __name__ == '__main__':
    import os.path
    ts = tscores.TissueScores(alias=False)

    match, non_match = read_hugo(os.path.join('data', 'non_alt_loci_set.txt'), ts)
    idx = create_idx(match, non_match)
    save_idx_lgroup(idx, match, non_match, os.path.join('data', 'idx.csv'), os.path.join('data', 'locus_group.csv'))
