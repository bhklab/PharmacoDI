import os
import glob
import pandas as pd
import numpy as np


def read_gene_signatures(pset_name, file_path):
    """
    Read all gene signatures for a PSet (to be used in gene_compounds table) from the directory file_path.

    @param pset_name: [`string`] The name of the PSet
    @param file_path: [`string`] The directory that holds all gene signature files
    @return: [`DataFrame`] A dataframe containing all gene signatures for the PSet
    """
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(
        f'{os.path.join(file_path, pset_name)}_gene_sig.csv')
    if len(pset_file) == 0:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read .csv file and return df
    return pd.read_csv(pset_file[0])


def build_gene_compound_tissue_dataset_df(gene_sig_dir, pset_name):
    """
    TODO - ask Chris to explain this table again

    @param gene_sig_dir: [`string`] The file path to the directory containing the gene
        signatures for each PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`DataFrame`] The gene_compounds table for this PSet, containing all stats (?)
    """

    # Get gene_sig_df from gene_sig_file
    gene_sig_df = read_gene_signatures(pset_name, gene_sig_dir)

    # Extract relevant columns
    # gene_compound_tissue_dataset = gctd
    gctd_df = gene_sig_df[['gene', 'drug', 'tissue', 'dataset',
        'estimate_analytic', 'lower_analytic', 'upper_analytic', 
        'lower_permutation', 'upper_permutation', 'n', 'pvalue_analytic', 
        'pvalue_permutation', 'df', 'fdr_analytic', 'fdr_permutation',
        'significant_permutation', 'mDataType']].copy()

    # Add missing columns
    gctd_df['sens_stat'] = 'AAC'
    gctd_df['permutation_done'] = 0
    gctd_df[~gctd_df['fdr_analytic'].isna(), 'permutation_done'] = 1

    # Rename foreign key columns
    gctd_df.rename(columns={'gene': 'gene_id', 'drug': 'compound_id',
        'tissue': 'tissue_id', 'dataset': 'dataset_id',
        'estimate_analytic': 'estimate'}, inplace=True)

    # Reorder columns
    return gctd_df[['gene_id', 'compound_id', 'dataset_id', 'tissue_id',
        'estimate', 'lower_analytic', 'upper_analytic', 
        'lower_permutation', 'upper_permutation', 'n', 'pvalue_analytic', 
        'pvalue_permutation', 'df', 'fdr_analytic', 'fdr_permutation',
        'significant_permutation', 'permutation_done', 'sens_stat', 
        'mDataType']]
