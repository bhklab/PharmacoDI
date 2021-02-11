import os
import glob
import re
import pandas as pd
import numpy as np
from PharmacoDI.build_primary_pset_tables import build_primary_pset_tables, build_cell_df, build_drug_df, build_tissue_df
from PharmacoDI.build_experiment_tables import build_experiment_tables, build_experiment_df
from PharmacoDI.build_pset_gene_drugs import build_gene_drug_df
from PharmacoDI.write_pset_table import write_pset_table
from PharmacoDI.build_dataset_join_tables import build_dataset_join_dfs


def build_all_pset_tables(pset_dict, pset_name, procdata_dir, gene_sig_dir):
    """
    Build all tables for a dataset and write them to a directory of all processed data.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param procdata_dir: [`string`] The file path to the directory containing processed data
    @param gene_sig_dir: [`string`] The file path to the directory containing gene_drugs data
    @return: [`None`]
    """
    pset_dfs = {}

    # Build primary tables (relating to cells, drugs, tissues, genes)
    print('Building primary tables...')
    pset_dfs = build_primary_pset_tables(pset_dict, pset_name)

    print('Building dataset join tables...')
    pset_dfs = {**pset_dfs, **build_dataset_join_dfs(
        pset_dict, pset_name, pset_dfs)}

    # Build experiment tables
    print('Building experiment tables...')
    # FIX: Modified to use pre-3.9 syntax to ensure backwards compatibility
    pset_dfs = {**pset_dfs, **build_experiment_tables(
        pset_dict, pset_name, pset_dfs['cell'])}

    # Build gene drugs table
    print('Building gene drug table...')
    pset_dfs['gene_drug'] = build_gene_drug_df(gene_sig_dir, pset_name)
    if not isinstance(pset_dfs['gene_drug'], pd.DataFrame):
        del pset_dfs['gene_drug']

    # Build summary/stats tables
    print('Building mol_cell and dataset_stats tables...')
    if 'gene_drug' in pset_dfs:
        pset_dfs['mol_cell'] = build_mol_cell_df(
            pset_dict, pset_name, pset_dfs['gene_drug'], pset_dfs['dataset_cell'])
    pset_dfs['dataset_statistics'] = build_dataset_stats_df(
        pset_dict, pset_name, pset_dfs)

    # Write all tables to CSV files
    for df_name in pset_dfs.keys():
        write_pset_table(pset_dfs[df_name], df_name, pset_name, procdata_dir)


def build_mol_cell_df(pset_dict, pset_name, gene_drug_df, dataset_cell_df=None):
    """
    Builds a table that summarizes the number of profiles, per cell line, per molecular data
    type, in this dataset. (Only considers molecular data types for which there are sens stats?)

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param gene_drug_df: [`pd.DataFrame`] The gene_drug table for this PSet
    @param dataset_cell_df: [`pd.DataFrame`] A table containing all the cells in this  
        PSet and the PSet name
    @return: [`pd.DataFrame`] The table with the number of profiles for each cell line, 
                                for each molecular data type
    """
    mol_cell_df = pd.DataFrame(
        columns=['cell_id', 'dataset_id', 'mDataType', 'num_prof'])

    molecularTypes = pd.unique(gene_drug_df['mDataType'])

    if 'molecularProfiles' in pset_dict:
        profiles_dict = pset_dict['molecularProfiles']
    else:
        profiles_dict = None

    if dataset_cell_df is None:
        dataset_cell_df = build_dataset_cell_df(
            pset_dict, pset_name, cell_df=None)

    for mDataType in molecularTypes:
        if isinstance(profiles_dict, dict):
            # Get the number of times each cellid appears in colData for that mDataType
            num_profiles = profiles_dict[mDataType]['colData']['cellid'].value_counts(
            )

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


def build_dataset_stats_df(pset_dict, pset_name, pset_dfs=None):
    """
    Summarizes how many cell lines, tissues, drugs, and experiments are contained
    within the dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param pset_dfs: [`dict`] A dictionary of tables from the PSet, with table names 
                                as the keys
    @return: [`pd.DataFrame`] A one-row table with the summary stats for this PSet
    """
    if pset_dfs is None:
        pset_dfs = {}
    if 'tissue' not in pset_dfs:
        pset_dfs['tissue'] = build_tissue_df(pset_dict)
    if 'cell' not in pset_dfs:
        pset_dfs['cell'] = build_cell_df(pset_dict)
    if 'drug' not in pset_dfs:
        pset_dfs['drug'] = build_drug_df(pset_dict)
    if 'experiment' not in pset_dfs:
        pset_dfs['experiment'] = build_experiment_df(
            pset_dict, pset_name, pset_dfs['cell'])

    return pd.DataFrame({
        'dataset_id': [pset_name],
        'cell_lines': [len(pset_dfs['cell'].index)],
        'tissues': [len(pset_dfs['tissue'].index)],
        'drugs': [len(pset_dfs['drug'].index)],
        'experiments': [len(pset_dfs['experiment'].index)]
    })
