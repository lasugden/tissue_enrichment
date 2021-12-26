import pandas as pd
import re


def read_gct(path: str) -> pd.DataFrame:
    """
    Read in a GCT file using Pandas.
    cmapPy does not work with Python 3, which is insanity in 2021.

    Args:
        path (str): Path to GCT file

    Returns:
        pd.DataFrame: Output dataframe containing all of the data
    
    """
    df = pd.read_csv(path, delimiter='\t', skiprows=2)
    return df


def tpm_levels_to_edge_weights(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert transcripts-per-million to be normalized across tissues
    Add a total tpm column

    Args:
        df (pd.DataFrame): raw dataframe

    Returns:
        pd.DataFrame: transformed, normalized dataframe
        
    """
    numeric_columns = df.select_dtypes(exclude='object').columns
    tot_tpm = df[numeric_columns].sum(axis=1)  # Unnecessary to include dtypes, but will do so for consistency
    df[numeric_columns] = df[numeric_columns].div(tot_tpm, axis=0)
    df['tot_tpm'] = tot_tpm
    df = df.loc[tot_tpm > 0, :]  # Limit to rows with measured transcripts in any tissue
    return df


def simplify_titles(df: pd.DataFrame) -> pd.DataFrame:
    """Simplify column titles

    Args:
        df (pd.DataFrame): input dataframe from gct

    Returns:
        pd.DataFrame: transformed dataframe with improved column names
    
    """
    transformation = {key: re.sub(r'\([^()]*\)', '', 
                                  key.lower().replace(' - ', '.'))
                                  .strip().replace(' ', '_') 
                                  for key in df.keys()}
    df = df.rename(columns=transformation)
    return df


def combine_similar(df: pd.DataFrame) -> pd.DataFrame:
    """Combine overly similar tissues based on pre-analysis

    Args:
        df (pd.DataFrame): Input dataframe with simplified tissue names

    Returns:
        pd.DataFrame: Output dataframe with overly similar tissues combined

    """
    df['skin'] = df[['skin.sun_exposed', 'skin.not_sun_exposed']].mean(axis=1)
    df = df.drop(columns=['skin.sun_exposed', 'skin.not_sun_exposed'])
    return df


def save_transcripts(df: pd.DataFrame, path: str):
    """Save outputs for efficient reading by backend

    Args:
        df (pd.DataFrame): transformed dataframe
        path (str): location to save file

    """
    df.to_csv(path, index=None)


if __name__ == '__main__':
    import os.path
    df = read_gct(os.path.join(os.path.expanduser('~/Documents/Programs/tissue_enrichment/data'), 
                               'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'))
    df = tpm_levels_to_edge_weights(df)
    df = simplify_titles(df)
    df = combine_similar(df)
    save_transcripts(df, os.path.join(os.path.expanduser('~/Documents/Programs/tissue_enrichment/data'), 
                                      'gene_median_tpm.csv'))
