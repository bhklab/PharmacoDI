import os
import pandas as pd
import numpy as np
from datatable import Frame
import datatable as dt
from datatable import f, update
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


def build_gene_compound_dataset_df(gene_compound_dataset_file, output_dir):
    """
    Build gene_compound_dataset table. This table contains a pancancer (across all tissues) 
    meta-analysis of compound sensitivity signatures by dataset and molecular data type.

    @param: [`string`] Path to the gene_compound_tissue .csv file.
    @param: [`string`] Path to write the output file to.
    @return: [`None`] Writes a the file 'gene_compound_dataset.csv' to output_dir.
    """
    if not os.path.exists(gene_compound_dataset_file):
        raise FileNotFoundError(f'Could not find the {gene_compound_dataset_file}')
    
    # -- Read in data
    gcd_dt = dt.fread(gene_compound_dataset_file)

    # -- Fix columns to match the ERD


    # -- Join to existing tables to get the proper FK ids
    

    # -- Write to output
    dt.fwrite(gcd_dt, file=os.path.join(output_dir, 'gene_compound_dataset.csv'))


def build_gene_compound_df():
    pass
