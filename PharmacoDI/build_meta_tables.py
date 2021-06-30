import os
import pandas as pd
import numpy as np
from datatable import Frame
from PharmacoDI.combine_pset_tables import write_table


def build_gene_compound_tissue_df(gene_compound_tissue_file, output_dir):
    """
    Build gene_compound_tissue table (description?)
    """
    if not os.path.exists(gene_compound_tissue_file):
        raise FileNotFoundError(f'Could not find the {gene_compound_tissue_file}')

    # gene_compound_tissue df = gct_df
    gct_df = pd.read_csv(gene_compound_tissue_file)
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

    # This column is in the CSV but not in the ERD
    gct_df.drop(columns=['hetTestRes'], inplace=True)

    # Convert to datatable.Frame for fast write to disk
    gct_df = Frame(gct_df)
    gct_df = write_table(gct_df, 'gene_compound_tissue', output_dir)
    return gct_df


def build_gene_compound_dataset_df():
    pass


def build_gene_compound_df():
    pass
