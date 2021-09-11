# old dependencies
import glob
import pandas as pd
from PharmacoDI.combine_pset_tables import write_table

# new dependenciess
import os
import warnings
import numpy as np
import re
# NOTE: No melt/pivot for datatable yet? Using polars instead
from datatable import dt, f, g, join, sort, update, fread
import polars as pl
from polars import col

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
    cell_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'cell.jay')).to_arrow()
    )
    dataset_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'dataset.jay')).to_arrow()
    )

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
        .rename({'value': 'cell_name', 'variable': 'dataset_id'})
            
    cell_synonym_df = cell_df \
        .join(cell_meta_long, left_on='name', right_on='unique.cellid', how='left') \
        .drop(['tissue_id', 'name']) \
        .drop_duplicates() \
        .rename({'id': 'cell_id'})

    # Create a map from dataset
    dataset_map = {dct['name']: str(dct['id']) for dct in 
        dataset_df.to_pandas().to_dict(orient='record')}
    # Regex the dataset identifiers to match the dataset map
    cell_synonym_df['dataset_id'] = cell_synonym_df['dataset_id'] \
        .apply(lambda x: re.sub('\.cellid$|[_.].*$', '', x)) \
        .apply(lambda x: re.sub('GDSC$', 'GDSC_v2', x)) \
        .apply(lambda x: re.sub('GDSC1.*$', 'GDSC_v1', x)) \
        .apply(lambda x: dataset_map[x]) \
        .cast(pl.Int64)
    
    cell_synonym_df = cell_synonym_df \
        .drop('cell_uid') \
        .drop_duplicates() \
        .drop_nulls()
    cell_synonym_df['id'] = range(1, cell_synonym_df.shape[0] + 1)
    
    # Convert to datatable.Frame for fast write to disk
    cell_synonym_dt = dt.Frame(cell_synonym_df.to_arrow())
    cell_synonym_dt.to_jay(os.path.join(output_dir, 'cell_synonym.jay'))


@logger.catch
def build_tissue_synonym_df(tissue_file, output_dir):
    # Get metadata file and tissue_df (assume that tissue_df is also in output_dir)
    tissue_metadata = pl.read_csv(tissue_file) # will read NA as string!
    tissue_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'tissue.jay')).to_arrow()
        )
    dataset_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'dataset.jay')).to_arrow()
        )

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
        .rename({'value': 'tissue_name', 'variable': 'dataset_id'})
    
    tissue_synonym_df = tissue_df \
        .join(tissue_meta_long, left_on='name', right_on='unique.tissueid', how='left') \
        .drop_duplicates() \
        .drop_nulls() \
        .rename({'id': 'tissue_id'}) \
        .drop('name')

    # Create a map from dataset
    dataset_map = {dct['name']: str(dct['id']) for dct in 
        dataset_df.to_pandas().to_dict(orient='record')}
    # Regex the dataset identifiers to match the dataset map
    tissue_synonym_df['dataset_id'] = tissue_synonym_df['dataset_id'] \
        .apply(lambda x: re.sub('\.cellid$|[_.].*$', '', x)) \
        .apply(lambda x: re.sub('GDSC$', 'GDSC_v2', x)) \
        .apply(lambda x: re.sub('GDSC1.*$', 'GDSC_v1', x)) \
        .apply(lambda x: dataset_map[x]) \
        .cast(pl.Int64)
    
    tissue_synonym_df = tissue_synonym_df.drop_duplicates()
    tissue_synonym_df['id'] = range(1, tissue_synonym_df.shape[0] + 1)

    # Convert to datatable.Frame for fast write to disk
    tissue_synonym_dt = dt.Frame(tissue_synonym_df.to_arrow())
    tissue_synonym_dt.to_jay(os.path.join(output_dir, 'tissue_synonym.jay'))


def build_compound_synonym_df(compound_file, output_dir):
    # Get metadata file and compound_df
    compound_metadata = pl.read_csv(compound_file, null_values='NA')
    compound_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'compound.jay')).to_arrow()
    )
    dataset_df = pl.from_arrow(
        fread(os.path.join(output_dir, 'dataset.jay')).to_arrow()
    )

    # Find all columns relevant to tissueid
    compound_cols = [col for col in compound_metadata.columns if 
        re.match('.*drugid$', col) and col != 'unique.drugid']
    
    # Read in which datasets we are working with
    dataset_names = os.listdir('procdata')
    clean_dataset_names = [re.sub('_.*$', '', name) for name in dataset_names]
    dataset_regex = re.compile('|'.join(clean_dataset_names))

    # Filter the cellid columns to only valid datasets
    compound_columns = [name for name in compound_cols if 
        re.match(dataset_regex, name)]

    # Get all unique synonyms and join with cell_df
    compound_meta_long = compound_metadata \
            .melt(id_vars='unique.drugid', value_vars=compound_columns) \
            .drop_nulls() \
            .drop_duplicates() \
            .rename({'value': 'compound_name', 'variable': 'dataset_id'}) \
            .filter(col('compound_name') != '')
    
    compound_synonym_df = compound_df \
        .join(compound_meta_long, left_on='name', right_on='unique.drugid', how='left') \
        .rename({'id': 'compound_id'}) \
        .drop(['compound_uid', 'name']) \
        .drop_nulls() \
        .drop_duplicates()

    # Create a map from dataset
    dataset_map = {dct['name']: str(dct['id']) for dct in 
        dataset_df.to_pandas().to_dict(orient='record')}
    # Regex the dataset identifiers to match the dataset map
    compound_synonym_df['dataset_id'] = compound_synonym_df['dataset_id'] \
        .apply(lambda x: re.sub('\.drugid$|[_.].*$', '', x)) \
        .apply(lambda x: re.sub('GDSC2019', 'GDSC_v2', x)) \
        .apply(lambda x: re.sub('GDSC1.*$', 'GDSC_v1', x)) \
        .apply(lambda x: dataset_map[x]) \
        .cast(pl.Int64)
    
    compound_synonym_df = compound_synonym_df.drop_duplicates()
    compound_synonym_df['id'] = range(1, compound_synonym_df.shape[0] + 1)

    # Convert to datatable.Frame for memory mapped output file
    df = dt.Frame(compound_synonym_df.to_arrow())
    df.to_jay(os.path.join(output_dir, 'compound_synonym.jay'))



## ---- DEPRECATED
## TODO:: Write new helper methods to reduce repition

# Helper function for getting all synonyms related to a certain df
def melt_and_join(meta_df, unique_id, join_df):
    """
    @param meta_df: [`polars.DataFrame`] The DataFrame containing all the 
        synonyms (metadata)
    @param unique_id: [`string`] The name of the column in the metadata 
        containing the unique IDs
    @param join_df: [`polars.DataFrame`] The DataFrame containing the primary 
        keys that will be used as foreign keys in the new synonyms df
    @return [`polars.DataFrame`] The synonys dataframe, with a PK, FK based on 
        join_df, and all unique synonyms
    """
    # Convert wide meta_df to long table
    # Drop 'variable' col (leave only unique ID and synonyms), drop duplicates
    synonyms = meta_df \
        .melt(id_vars=unique_id) \
        .select([unique_id, 'value']) \
        .drop_duplicates() \
        .filter(pl.col('value').is_not_null())
    # Join with join_df based on unique_id
    synonyms = synonyms \
        .join(join_df, left_on=unique_id, right_on='name', how='inner') \
        .select(['id', 'value'])
    synonym['id'].cast(pl.Int64)
    return synonyms

