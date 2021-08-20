import os
import glob
import re
import pandas as pd
import numpy as np
from datetime import date
from .build_primary_pset_tables import build_primary_pset_tables, build_cell_df, build_compound_df, build_tissue_df
from .build_experiment_tables import build_experiment_tables, build_experiment_df
from .build_gene_compound_tissue_dataset_tables import build_gene_compound_tissue_dataset_df
from .write_pset_table import write_pset_table
from .build_dataset_join_tables import build_dataset_join_dfs, build_dataset_cell_df

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_all_pset_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def build_all_pset_tables(
    pset_dict: dict, 
    pset_name: str, 
    procdata_dir: str, 
    gene_sig_dir: str
) -> None:
    """
    Build all tables for a dataset and write them to a directory of all processed data.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param procdata_dir: [`string`] The file path to the directory containing processed data
    @param gene_sig_dir: [`string`] The file path to the directory containing gene_compounds data
    @return: [`None`]
    """
    pset_dfs = {}

    # Build primary tables (relating to cells, compounds, tissues, genes)
    logger.info('Building primary tables...')
    pset_dfs = build_primary_pset_tables(pset_dict, pset_name)

    logger.info('Building dataset join tables...')
    pset_dfs = {**pset_dfs, **build_dataset_join_dfs(
        pset_dict, pset_name, pset_dfs)}

    # Build experiment tables
    logger.info('Building experiment tables...')
    # FIX: Modified to use pre-3.9 syntax to ensure backwards compatibility
    pset_dfs = {**pset_dfs, **build_experiment_tables(
        pset_dict, pset_name, pset_dfs['cell'])}

    # Build summary/stats tables
    logger.info('Building mol_cell table...')
    pset_dfs['mol_cell'] = build_mol_cell_df(
        pset_dict, pset_name, pset_dfs['dataset_cell'])
    logger.info('Building dataset_statistics table...')
    pset_dfs['dataset_statistics'] = build_dataset_stats_df(
        pset_dict, pset_name, pset_dfs)

    # Write all tables to CSV files
    for df_name in pset_dfs.keys():
        write_pset_table(pset_dfs[df_name], df_name, pset_name, procdata_dir)

    log_file = open(os.path.join(procdata_dir, pset_name, 
        f'{pset_name}_log.txt'), "w")
    log_file.write('Wrote the following PSet tables: \n')
    log_file.writelines(f'{table}\n' for table in pset_dfs.keys())
    log_file.write(f'on {date.today()}')
    log_file.close()


@logger.catch
def build_mol_cell_df(
    pset_dict: dict, 
    pset_name: str, 
    dataset_cell_df: pd.DataFrame=None, 
    molecularTypes: list=['rna', 'rnaseq', 'cnv', 'mutation']
) -> pd.DataFrame:
    """
    Builds a table that summarizes the number of profiles, per cell line, per molecular data
    type, in this dataset. (Only considers molecular data types for which there are sens stats?)

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param gene_compound_df: [`pd.DataFrame`] The gene_compound table for this PSet
    @param dataset_cell_df: [`pd.DataFrame`] A table containing all the cells in this  
        PSet and the PSet name
    @return: [`pd.DataFrame`] The table with the number of profiles for each cell line, 
        for each molecular data type
    """
    mol_cell_df = pd.DataFrame(
        columns=['cell_id', 'dataset_id', 'mDataType', 'num_prof'])
    if 'molecularProfiles' in pset_dict:
        profiles_dict = pset_dict['molecularProfiles']
        molecularTypes = list(profiles_dict.keys())
    else:
        profiles_dict = None
    if dataset_cell_df is None:
        dataset_cell_df = build_dataset_cell_df(
            pset_dict, pset_name, cell_df=None)
    for mDataType in molecularTypes:
        if isinstance(profiles_dict, dict):
            # Get the number of times each cellid appears in colData for that mDataType
            num_profiles = profiles_dict[mDataType]['colData']['cellid'] \
                .value_counts()
            # Join with datasets cells on cellid
            df = pd.merge(dataset_cell_df, num_profiles,
                left_on='cell_id', right_on=num_profiles.index, how='left')
            # Rename num_profiles column
            df.rename(columns={'cellid': 'num_prof'}, inplace=True)
            # Set mDataType column to the current molecular type
            df['mDataType'] = mDataType
        else:
            # If PSet contains no molecular profiles, set num_prof to 0
            # for all celll lines and all molecular data types
            df = dataset_cell_df.copy()
            df['mDataType'] = mDataType
            df['num_prof'] = 0
        # Append to mol_cell_df
        mol_cell_df = mol_cell_df.append(df)

    # Replace any NaN in the num_profiles column with 0
    mask = mol_cell_df.query('num_prof.isna()').index
    mol_cell_df.loc[mask, 'num_prof'] = 0
    mol_cell_df['num_prof'] = mol_cell_df['num_prof'].astype('int32')

    return mol_cell_df


@logger.catch
def build_dataset_stats_df(pset_dict, pset_name, pset_dfs=None):
    """
    Summarizes how many cell lines, tissues, compounds, and experiments 
    are contained within the dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in 
        the PSet
    @param pset_name: [`string`] The name of the PSet
    @param pset_dfs: [`dict`] A dictionary of tables from the PSet, 
        with table names as the keys
    @return: [`pd.DataFrame`] A one-row table with the summary stats for 
        this PSet
    """
    if pset_dfs is None:
        pset_dfs = {}
    if 'tissue' not in pset_dfs:
        pset_dfs['tissue'] = build_tissue_df(pset_dict)
    if 'cell' not in pset_dfs:
        pset_dfs['cell'] = build_cell_df(pset_dict)
    if 'compound' not in pset_dfs:
        pset_dfs['compound'] = build_compound_df(pset_dict)
    if 'experiment' not in pset_dfs:
        pset_dfs['experiment'] = build_experiment_df(
            pset_dict, pset_name, pset_dfs['cell'])

    return pd.DataFrame({
        'dataset_id': [pset_name],
        'cell_lines': [len(pset_dfs['cell'].index)],
        'tissues': [len(pset_dfs['tissue'].index)],
        'compounds': [len(pset_dfs['compound'].index)],
        'experiments': [len(pset_dfs['experiment'].index)]
    })