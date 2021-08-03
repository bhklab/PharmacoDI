# old dependencies
import glob
import pandas as pd
from PharmacoDI.combine_pset_tables import write_table

# new dependenciess
import os
import warnings
import numpy as np
import re
# NOTE: No melt/cast for datatable yet? Using polars instead
from datatable import dt, f, g, join, sort, update, fread
import polars as pl

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_synonym_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)

# --- SYNONYMS TABLES --------------------------------------------------------------------------

@logger.catch
def build_cell_synonym_df(cell_file, output_dir):

    # Get metadata file and cell_df
    cell_metadata = pl.read_csv(cell_file, null_values='NA')
    cell_df = pl.read_csv(os.path.join(output_dir, 'cell.csv'), null_values='NA')

    # Find all columns relevant to cellid
    cell_cols = [col for col in cell_metadata.columns if 
        re.match('.*cellid$', col) and col != 'unique.cellid']

    # Read in which datasets we are working with
    dataset_names = os.listdir('procdata')
    clean_dataset_names = [re.sub('_.*$', '', name) for name in dataset_names]
    dataset_regex = re.compile('|'.join(clean_dataset_names))

    # Filter the cellid columns to only valid datasets
    cell_columns = [name for name in cell_cols if 
        re.match(dataset_regex, name)]

    # Get all unique synonyms and join with cell_df
    cell_meta_long = cell_metadata \
            .melt(id_vars='unique.cellid', value_vars=cell_columns) \
            .drop_nulls() \
            .drop_duplicates() \
            .rename({'value': 'cell_name'})
            #.select(['unique.cellid', 'value']) \
            
    cell_synonym_df = cell_df \
        .join(cell_meta_long, left_on='name', right_on='unique.cellid', how='left') \
        .drop(['tissue_id', 'name']) \
        .drop_duplicates() \
        .rename({'id': 'cell_id'})

    # Add blank col for dataset_id
    cell_synonym_df['dataset_id'] = np.nan * np.empty(cell_synonym_df.shape[0])
    cell_synonym_df['id'] = range(1, cell_synonym_df.shape[0] + 1)
    
    # Convert to datatable.Frame for fast write to disk
    cell_synonym_df.to_csv(os.path.join(output_dir, 'cell_synonym.csv'))
    return cell_synonym_df
    

@logger.catch
def build_tissue_synonym_df(tissue_file, output_dir):
    # Get metadata file and tissue_df (assume taht tissue_df is also in output_dir)
    tissue_metadata = pl.read_csv(tissue_file, null_values='NA')
    tissue_df = pl.read_csv(os.path.join(output_dir, 'tissue.csv'), null_values='NA')

    # Find all columns relevant to tissueid
    tissue_cols = [col for col in tissue_metadata.columns if 
        re.match('.*tissueid$', col) and col != 'unique.tissueid']

    # Read in which datasets we are working with
    dataset_names = os.listdir('procdata')
    clean_dataset_names = [re.sub('_.*$', '', name) for name in dataset_names]
    dataset_regex = re.compile('|'.join(clean_dataset_names))

    # Filter the cellid columns to only valid datasets
    tissue_columns = [name for name in tissue_cols if 
        re.match(dataset_regex, name)]

    # Get all unique synonyms and join with cell_df
    tissue_meta_long = tissue_metadata \
            .melt(id_vars='unique.tissueid', value_vars=tissue_columns) \
            .drop_nulls() \
            .drop_duplicates() \
            .rename({'value': 'tissue_name'})
    
    tissue_synonym_df = tissue_df \
        .join(tissue_meta_long, left_on='name', right_on='unique.tissueid', how='left') \
        .drop_duplicates()

    # Add blank col for dataset_id
    tissue_synonym_df['dataset_id'] = np.nan * np.empty(tissue_synonym_df.shape[0])
    tissue_synonym_df['id'] = range(1, tissue_synonym_df.shape[0] + 1)

    # Convert to datatable.Frame for fast write to disk
    tissue_synonym_df.to_csv(os.path.join(output_dir, 'tissue_synonym.csv'))
    return tissue_synonym_df


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

