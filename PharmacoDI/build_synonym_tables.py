import os
import re
import glob
import numpy as np
import pandas as pd
from datatable import Frame
from PharmacoDI.combine_pset_tables import write_table


# TODO: add docstrings

def get_metadata(file_name):
    # Check that metadata file exists
    if not os.path.exists(file_name):
        raise ValueError(
            f'No metadata file {file_name} could be found!')
    
    # Read csv file and return df
    return pd.read_csv(file_name)


# --- SYNONYMS TABLES --------------------------------------------------------------------------

def build_cell_synonym_df(cell_file, output_dir):
    # Get metadata file and cell_df
    cell_metadata = get_metadata(cell_file)
    cell_df = pd.read_csv(os.path.join(output_dir, 'cell.csv'))
    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_metadata[[
        col for col in cell_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with cell_df
    cell_synonym_df = melt_and_join(cell_columns, 'unique.cellid', cell_df)
    cell_synonym_df = cell_synonym_df.rename(columns={'id': 'cell_id', 'value': 'cell_name'})

    # Add blank col for dataset_id (TODO)
    cell_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(cell_synonym_df)
    df = write_table(df, 'cell_synonym', output_dir)
    return df
    

def build_tissue_synonym_df(tissue_file, output_dir):
    # Get metadata file and tissue_df (assume taht tissue_df is also in output_dir)
    tissue_metadata = get_metadata(tissue_file)
    tissue_df = pd.read_csv(os.path.join(output_dir, 'tissue.csv'))

    # Find all columns relevant to tissueid
    pattern = re.compile('tissueid')
    tissue_cols = tissue_metadata[[
        col for col in tissue_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with tissue_df
    tissue_synonym_df = melt_and_join(tissue_cols, 'unique.tissueid', tissue_df)
    tissue_synonym_df = tissue_synonym_df.rename(columns={'id': 'tissue_id', 'value': 'tissue_name'})

    # Add blank col for dataset_id (TODO)
    tissue_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(tissue_synonym_df)
    df = write_table(df, 'tissue_synonym', output_dir)
    return df


def build_compound_synonym_df(compound_file, output_dir):
    # Get metadata file and compound_df
    compound_metadata = get_metadata(compound_file)
    compound_df = pd.read_csv(os.path.join(output_dir, 'compound.csv'))

    # Find all columns relevant to compoundid
    # Right now only FDA col is dropped, but may be more metadata in the future
    pattern = re.compile('drugid')
    compound_cols= compound_metadata[[
        col for col in compound_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with compounds_df
    compound_synonym_df = melt_and_join(compound_cols, 'unique.drugid', compound_df)
    compound_synonym_df = compound_synonym_df.rename(columns={'id': 'compound_id', 'value': 'compound_name'})

    # Add blank col for dataset_id (TODO)
    compound_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(compound_synonym_df)
    df = write_table(df, 'compound_synonym', output_dir)
    return df


# Helper function for getting all synonyms related to a certain df
def melt_and_join(meta_df, unique_id, join_df):
    """
    @param meta_df: [`Dask DataFrame`] The DataFrame containing all the synonyms (metadata)
    @param unique_id: [`string`] The name of the column in the metadata containing the unique IDs
    @param join_df: [`Dask DataFrame`] THe DataFrame containing the primary keys that will be used as
        foreign keys in the new synonyms df

    @return [`DataFrame`] The synonys dataframe, with a PK, FK based on join_df, and all unique synonyms
    """
    # Convert wide meta_df to long table
    # Drop 'variable' col (leave only unique ID and synonyms), drop duplicates
    synonyms = pd.melt(meta_df, id_vars=[unique_id])[
        [unique_id, 'value']].drop_duplicates()

    # Drop all rows where value is NA
    synonyms = synonyms[synonyms['value'].notnull()]

    # Join with join_df based on unique_id
    synonyms = pd.merge(synonyms, join_df, left_on=unique_id,
                        right_on='name', how='inner')[['id', 'value']]
    synonyms['id'] = synonyms['id'].astype('int')

    return synonyms

