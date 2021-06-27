import os
import pandas as pd
import numpy as np

file_path = 'rawdata/gene_signatures/metaanalysis/gene_compound_tissue.csv'

def build_gene_compound_tissue_df(file_path):
    """
    Build gene_compound_tissue table (description?)
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'Could not find the {file_path}')

    # gene_compound_tissue df = gct_df
    gct_df = pd.read_csv(file_path)
    gct_df.rename(columns={'Gene': 'gene_id', 'Tissue': 'tissue_id', 'Drug': 'compound_id'}, inplace=True)

    # TODO: we still need these columns (not in CSV)
    gct_df['n'] = np.nan
    gct_df['tstat'] = np.nan
    gct_df['fstat'] = np.nan
    gct_df['df'] = np.nan
    gct_df['fdr'] = np.nan
    gct_df['FWER_drugs'] = np.nan
    gct_df['FWER_all'] = np.nan
    gct_df['BF_p_all'] = np.nan
    gct_df['sens_stat'] = 'AAC'
    gct_df['tested_in_human_trials'] = np.nan
    gct_df['in_clinical_trials'] = np.nan

    # This column is in the CVS but not in the ERD
    gct_df.drop(columns=['hetTestRes'], inplace=True)

    return gct_df


def build_gene_compound_dataset_df():
    pass


def build_gene_compound_df():
    pass
