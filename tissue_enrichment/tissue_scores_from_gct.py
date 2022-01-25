import numpy as np
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
    """Simplify column titles of tissues

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


def fix_ensembl_titles(df: pd.DataFrame) -> pd.DataFrame:
    """Remove version numbers from ensembl titles

    Args:
        df (pd.DataFrame): input dataframe with ensembl titles

    Returns:
        pd.DataFrame: updated dataframe with corrected titles
    
    """
    df['name'] = df['name'].str.replace(r'\..*', '', regex=True)
    df.loc[df['description'].str.lower().duplicated(), 'description'] = pd.NA
    return df


def convert_to_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Convert tissue expression to rank 

    Args:
        df (pd.DataFrame): normalized tissue scores per tissue

    Returns:
        pd.DataFrame: updated dataframe with rank rather than score
    
    """
    numeric_columns = df.select_dtypes(exclude='object').columns.drop('tot_tpm')
    for col in numeric_columns:
        arr = df[col].argsort()
        ranks = np.empty_like(arr)
        ranks[arr] = np.arange(len(arr)).astype(float)/len(arr)
        df[col] = ranks
    return df


def save_transcripts(df: pd.DataFrame, path: str):
    """Save outputs for efficient reading by backend

    Args:
        df (pd.DataFrame): transformed dataframe
        path (str): location to save file

    """
    df.to_csv(path, index=None)


def process_raw(raw_path: str, save_path: str, rank_path: str):
    """Process raw data and save

    Args:
        raw_path (str): raw path to gct file
        save_path (str): path to save data into, usually in data directory in tissue_enrichment repo
        rank_path (str): path to save rank data into, usually in data directory in tissue_enrichment repo
    
    """
    df = read_gct(raw_path)
    df = tpm_levels_to_edge_weights(df)
    df = simplify_titles(df)
    df = fix_ensembl_titles(df)
    save_transcripts(df, save_path)
    convert_to_rank(df)
    save_transcripts(df, rank_path)
