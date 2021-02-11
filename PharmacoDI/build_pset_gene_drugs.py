import os
import glob
import pandas as pd
import numpy as np


def read_gene_signatures(pset_name, file_path):
    """
    Read all gene signatures for a PSet (to be used in gene_drugs table) from the directory file_path.

    @param pset_name: [`string`] The name of the PSet
    @param file_path: [`string`] The directory that holds all gene signature files
    @return: [`DataFrame`] A dataframe containing all gene signatures for the PSet
    """
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(
        f'{os.path.join(file_path, pset_name, pset_name)}_gene_sig.csv')
    if len(pset_file) == 0:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read .csv file and return df
    return pd.read_csv(pset_file[0])


def build_gene_drug_df(gene_sig_dir, pset_name):
    """
    TODO - ask Chris to explain this table again

    @param gene_sig_dir: [`string`] The file path to the directory containing the gene
        signatures for each PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`DataFrame`] The gene_drugs table for this PSet, containing all stats (?)
    """
    # If gene signature file doesn't exist, return empty DataFrame
    if not os.path.exists(os.path.join(gene_sig_dir, pset_name)):
        print(
            f'WARNING: gene signature annotations file does not exist for {pset_name} in {gene_sig_dir}')
        return None

    # Get gene_sig_df from gene_sig_file
    gene_sig_df = read_gene_signatures(pset_name, gene_sig_dir)

    # Extract relevant columns
    gene_drug_df = gene_sig_df[[
        'gene', 'drug', 'estimate', 'n', 'pvalue', 'df', 'fdr', 'tissue', 'mDataType']].copy()
    # TODO - cannot find 'se', 'sens_stat' -- is one of these 'significant'???
    # Chris: You will determine significance based on the fdr (false discovery rate) at alpha = 0.05, it will be TRUE or FALSE (or 1 or 0)

    # Chris: 'se' - leave NA/Null/None for now, it will be added as a column to the gene signatures the next time we run them.
    gene_drug_df['se'] = np.nan
    # Chris: 'sens_stat' - I will add this to the function for extracting per PSet gene signatures - for now it is always 'AAC' (Area above dose-response curve)
    gene_drug_df['sens_stat'] = 'AAC'
    # TODO - cannot find 'drug_like_molecule', 'in_clinical_trials'
    # Chris: Have renamed it to tested_in_human_trials, it will indicate a 1 if it has ever been tested in a human clinical trial (even if it failed)
    # Chris: Source for this data will be clinicaltrails.gov
    # TODO - check out API, leave NA for now
    gene_drug_df['tested_in_human_trials'] = np.nan
    gene_drug_df['in_clinical_trials'] = np.nan

    # Rename foreign key columns
    gene_drug_df.rename(
        columns={'gene': 'gene_id', 'drug': 'drug_id', 'tissue': 'tissue_id'}, inplace=True)

    # Add dataset id
    gene_drug_df['dataset_id'] = pset_name

    # Add missing columns (TODO - get this data)
    gene_drug_df['tstat'] = np.nan
    gene_drug_df['fstat'] = np.nan
    gene_drug_df['FWER_genes'] = np.nan
    gene_drug_df['FWER_drugs'] = np.nan
    gene_drug_df['FWER_all'] = np.nan
    gene_drug_df['BF_p_all'] = np.nan
    gene_drug_df['meta_res'] = np.nan

    # Reorder columns
    return gene_drug_df[['gene_id', 'drug_id', 'estimate', 'se', 'n', 'tstat', 'fstat',
                         'pvalue', 'df', 'fdr', 'FWER_genes', 'FWER_drugs', 'FWER_all',
                         'BF_p_all', 'meta_res', 'dataset_id', 'sens_stat', 'tissue_id',
                         'mDataType', 'tested_in_human_trials', 'in_clinical_trials']]
