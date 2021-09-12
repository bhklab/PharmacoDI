import os
import glob
import pandas as pd
import numpy as np
import polars as pl
import re

from polars import col
from datatable import dt, fread, f, g, by, sort

from .utilities import harmonize_df_columns

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_primary_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def build_primary_pset_tables(pset_dict, pset_name):
    """
    Build the tissue, compound, and gene tables for a PSet and return them
    in a dictionary, with table names as the keys.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`dict`] A dictionary of primary tables, with table names as keys
    """
    pset_dfs = {}
    pset_dfs['dataset'] = pd.Series(pset_name, name='name')
    pset_dfs['tissue'] = build_tissue_df(pset_dict)
    pset_dfs['compound'] = build_compound_df(pset_dict)
    pset_dfs['compound_annotation'] = build_compound_annotation_df(pset_dict)
    pset_dfs['cell'] = build_cell_df(pset_dict)
    # Don't make gene table if there are no molecular profiles
    if 'molecularProfiles' in pset_dict:
        pset_dfs['gene'] = build_gene_df(pset_dict)
        pset_dfs['gene_annotation'] = build_gene_annotation_df(pset_dict)
    return pset_dfs


@logger.catch
def build_gene_df(pset_dict):
    """
    Build a table containing all genes in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The gene table
    """
    gene_df = pd.Series([], name='name', dtype='str')
    for mDataType in pset_dict['molecularProfiles']:
        gene_df = gene_df.append(pd.Series(pd.unique(
            pset_dict['molecularProfiles'][mDataType]['rowData']['.features']),
            name='name', dtype='str'))
    gene_df.replace(r'\.[0-9]*$', '', regex=True, inplace=True)
    gene_df.drop_duplicates(inplace=True)
    return gene_df


@logger.catch
def build_tissue_df(pset_dict):
    """
    Build a table containing all tissues in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The tissue table
    """
    tissue_df = pd.Series(
        pd.unique(pset_dict['cell']['tissueid']),
        name='name'
    )
    tissue_df[tissue_df.name.isna(), 'name'] = 'NA'
    tissue_df.sort_values('name', inplace=True)
    return tissue_df


@logger.catch
def build_compound_df(pset_dict):
    """
    Build a table containing all compounds in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The compound table
    """
    compound_df = pd.Series(pd.unique(pset_dict['drug']['drugid']), name='name')
    return compound_df


@logger.catch
def build_gene_annotation_df(pset_dict):
    """
    Build a table mapping each gene in a dataset to its gene annotations.
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all gene annotations, mapped to genes
    """
    # Extract the all molecular data types for the pSet
    df_list = [pl.from_pandas(pset_dict['molecularProfiles'][mDataType]['rowData'])
        for mDataType in pset_dict['molecularProfiles']]
    # Get columns of interest, add columns needed later
    for i in range(len(df_list)):
        df_list[i] = df_list[i].select(['.features'])
        empty_column = [None for _ in range(len(df_list[i]['.features']))]
        df_list[i]['symbol'] = pl.Series('symbol', empty_column, dtype=pl.Utf8)
        df_list[i]['gene_seq_start'] = pl.Series('gene_seq_start', 
            empty_column, dtype=pl.Int64)
        df_list[i]['gene_seq_end'] = pl.Series('gene_seq_end', 
            empty_column, dtype=pl.Int64)
    # Merge to a single DataFrame
    gene_annotation_df = pl.concat(df_list) \
        .rename({'.features': 'gene_id'})
    # Remove Ensembl gene version
    gene_annotation_df['gene_id'] = gene_annotation_df['gene_id'] \
        .apply(lambda x: re.sub(r'\..*$', '', x))
    gene_annotation_df = gene_annotation_df \
        .drop_duplicates() \
        .to_pandas()
    return gene_annotation_df


@logger.catch
def build_compound_annotation_df(
    pset_dict: dict,
    column_dict: dict={'compound_id': str, 'smiles': str, 'inchikey': str, 
        'pubchem': str, 'FDA': bool},
    rename_dict: dict={'rownames': 'compound_id', 'cid': 'pubchem', 
            'FDA': 'fda_status'}
) -> pd.DataFrame:
    """
    Build a table mapping each compound in a dataset to its compound 
    annotations.

    @param pset_dict: A nested dictionary containing all tables
        in from the PharmacoSet object.
    @return: A table of all compound annotations, 
        mapped to compounds
    """
    compound_annotation_df = pset_dict['drug'].copy()
    compound_annotation_df.rename(columns=rename_dict, inplace=True)
    compound_annotation_df = harmonize_df_columns(
        df=compound_annotation_df,
        column_dict=column_dict
    )
    return compound_annotation_df


# TODO - confirm that you're using the correct cell id
@logger.catch
def build_cell_df(pset_dict):
    """
    Build a table containing all the cells in a dataset, mapped to their tissues.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all cell lines, mapped to tissues
    """
    cell_df = pset_dict['cell'][['cellid', 'tissueid']].copy()
    cell_df.rename(columns={'cellid': 'name',
                            'tissueid': 'tissue_id'}, inplace=True)
    return cell_df
