import os
import glob
import pandas as pd
import numpy as np
import polars as pl
from datatable import dt, fread, f, g, join, by, sort
import re

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_gene_compound_tissue_dataset_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def read_gene_signatures(pset_name, file_path):
    """
    Read all gene signatures for a PSet (to be used in gene_compounds table) 
    from the directory file_path.

    @param pset_name: [`string`] The name of the PSet
    @param file_path: [`string`] The directory that holds all gene signature files
    @return: [`DataFrame`] A dataframe containing all gene signatures for the PSet
    """
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(
        f'{os.path.join(file_path, pset_name)}_gene_sig.parquet')
    if len(pset_file) == 0:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name}_gene_sig.parquet '
            'could be found in {file_path}'
            )
    # Read .parquet file and return df
    return pd.read_parquet(pset_file[0])


@logger.catch
def build_gene_compound_tissue_dataset_df(
        gene_sig_dir: str, 
        pset_name: str, 
        ignore_psets: list=['NCI60', 'PRISM']
) -> pd.DataFrame:
    """
    @param gene_sig_dir: [`string`] The file path to the directory 
        containing the gene signatures for each PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`DataFrame`] The gene_compounds table for this PSet, 
        containing all stats (?)
    """
    # Early return for PSets without gene signatures
    if pset_name in ignore_psets:
        return None
    # Get gene_sig_df from gene_sig_file
    try:
        gene_sig_df = dt.fread(
            os.path.join(gene_sig_dir, 'gene_compound_tissue_dataset.csv'),
            memory_limit=int(1e10) # 10 GBs
        )
    except ValueError:
        return None
    gene_sig_df[f.dataset == pset_name, :]
    # Extract relevant columns
    # gene_compound_tissue_dataset = gctd
    gctd_df = gene_sig_df.loc[:, ['gene', 'compound', 'tissue', 'dataset',
        'estimate', 'lower_analytic', 'upper_analytic', 
        'lower_permutation', 'upper_permutation', 'n', 'pvalue_analytic', 
        'pvalue_permutation', 'df', 'fdr_analytic', 'fdr_permutation',
        'significant_permutation', 'mDataType']]
    # Add missing columns
    gctd_df.loc[:, 'sens_stat'] = 'AAC'
    gctd_df.loc[:, 'permutation_done'] = 0
    gctd_df.loc[~gctd_df['fdr_permutation'].isna(), 'permutation_done'] = 1
    # Rename foreign key columns
    gctd_df.rename(columns={'gene': 'gene_id', 'compound': 'compound_id',
        'tissue': 'tissue_id', 'dataset': 'dataset_id'}, inplace=True)
    gctd_df['gene_id'] = gctd_df['gene_id'] \
        .apply(lambda x: re.sub(r'\..*$', '', x))
    # Reorder columns
    return gctd_df[['gene_id', 'compound_id', 'dataset_id', 'tissue_id',
        'estimate', 'lower_analytic', 'upper_analytic', 
        'lower_permutation', 'upper_permutation', 'n', 'pvalue_analytic', 
        'pvalue_permutation', 'df', 'fdr_analytic', 'fdr_permutation',
        'significant_permutation', 'permutation_done', 'sens_stat', 
        'mDataType']]
