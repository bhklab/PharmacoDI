import os
import warnings
import pandas as pd
import numpy as np
from datatable import dt, f, g, join, sort, update, fread


def build_gene_compound_tissue_df(gene_compound_tissue_file, gene_file, compound_file, tissue_file, output_dir):
    """
    Build gene_compound_tissue table (description?)

    @param gene_compound_tissue_file: [`str`] Path to the gene signature .csv file.
    @param gene_file: [`str`] Path to the gene table .csv file.
    @param compound_file: [`str`] Path to the compound table .csv file.
    @param tissue_file: [`str`] Path to the tissue table .csv file.
    @param output_dir: [`str`] Path to write the output file to.

    @return [`datatable.Frame`] Writes the 'gene_compound_tissue.csv' file to 
      output_dir then returns the table.
    """
    # -- Check the input files exist
    for file in [gene_compound_tissue_file, gene_file, compound_file, tissue_file]:
        if not os.path.exists(file):
            raise FileNotFoundError(f'Could not find the {file}')

    # -- Read in mapping tables
    gene_dt = fread(gene_file)
    compound_dt = fread(compound_file)
    tissue_dt = fread(tissue_file)

    # -- Read in gene_compound_tissue table
    gct_dt = fread(gene_compound_tissue_file)

    # -- Fix names and assign missing columns
    gct_dt.names = {'Gene': 'gene_id', 'Tissue': 'tissue_id', 
        'Drug': 'compound_id', 'FWER_genes': 'FWER_gene'}
    # Determine missing columns and assign them, so we don't have to change code 
    #>when new columns are addeds
    gct_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'tissue_id', 
        'estimate', 'lower', 'upper', 'n', 'tstat', 'fstat', 'pvalue', 'df',
        'fdr', 'FWER_gene', 'FWER_compound', 'FWER_all', 'BF_p_all', 'sens_stat', 
        'mDataType', 'tested_in_human_trials', 'in_clinical_trials'))
    gct_missing_columns = np.setdiff1d(gct_table_columns, np.asarray(gct_dt.names))
    for col in gct_missing_columns:
        gct_dt[col] = None
    gct_dt1 = gct_dt[:, list(gct_table_columns)]
    # Sanity check the columns are there
    if not np.all(gct_table_columns == np.asarray(gct_dt1.names)):
        raise ValueError(f'The gene_compound_tissue table',
            ' has missing columns!')

    # -- Map to existing FK ids
    # gene id
    gct_dt1.names = {'gene_id': 'gene_name'}
    gene_dt.names = {'id': 'gene_id', 'name': 'gene_name'}
    gene_dt.key = 'gene_name'
    # NOTE: the g object references the joined tables namespace
    gct_dt1[:, update(gene_id=g.gene_id), join(gene_dt)]

    ## TODO:: rewrite as a helper function
    # regex match failed ids, then assign to the table
    failed_genes = gct_dt1[dt.isna(f.gene_id), 'gene_name'].to_numpy().flatten()
    if len(failed_genes) > 0:
        match_genes = [gene_dt['gene_name'] \
            .to_pandas().gene_name.str.match(f'{gene}.*') for gene in failed_genes]
        failed_match_idx = np.where(~np.asarray([np.any(gene) for gene in match_genes]))[0]
        match_idx = np.where(match_genes)[1]
        # needs to be float64 because Numpy has no NaN for int types... makes no sense!?
        # Pad with NaNs for failed matchess
        gene_ids = gene_dt[match_idx, 'gene_id'].to_pandas().gene_id
        if (len(failed_match_idx) > 1):
            gene_ids = pd.Series(np.insert(gene_ids, failed_match_idx, None), dtype='int32')
        gct_dt1[dt.isna(f.gene_id), update(gene_id=gene_ids)]


    if (np.any(gct_dt1[:, dt.isna(f.gene_id)].to_numpy())):
        warnings.warn('Some gene_ids in gene_compound_tissue are still NA! Dropping'
            'the missing rows...')
        gct_dt1 = gct_dt1[~dt.isna(f.gene_id), :]

    # compound id
    gct_dt1.names = {'compound_id': 'compound_name'}
    compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
    del compound_dt[:, 'compound_uid']
    compound_dt.key = 'compound_name'
    gct_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]

    # tissue id
    gct_dt1.names = {'tissue_id': 'tissue_name'}
    tissue_dt.names = {'id': 'tissue_id', 'name': 'tissue_name'}
    tissue_dt.key = 'tissue_name'
    gct_dt1[:, update(tissue_id=g.tissue_id), join(tissue_dt)]

    ## TODO: Handle failed tissue mappings?

    # -- Sort then assign the primary key column
    ## TODO:: Is there a way to sort by reference? Does setting the key do this
    gct_dt2 = gct_dt1[:, list(gct_table_columns), sort('gene_id', 'compound_id', 'tissue_id')]
    gct_dt2[:, update(id=range(1, gct_dt2.nrows + 1))]

    # Sanity check we didn't lose any rows
    if not gct_dt.nrows == gct_dt2.nrows:
        warnings.warn('The compound_gene_tissue table has lost some rows!')

    dt.fwrite(gct_dt2, file=os.path.join(output_dir, 'compound_gene_tissue.csv'))
    return(gct_dt2)


## FIXME:: This function is almost identical to the previous one, refactor
##>into a helper for gene_compound instead of copy pasting code
def build_gene_compound_dataset_df(gene_compound_dataset_file, gene_file, compound_file, dataset_file, output_dir):
    """
    Build gene_compound_dataset table (description?)

    @param gene_compound_dataset_file: [`str`] Path to the gene signature .csv file.
    @param gene_file: [`str`] Path to the gene table .csv file.
    @param compound_file: [`str`] Path to the compound table .csv file.
    @param dataset_file: [`str`] Path to the tissue table .csv file.
    @param output_dir: [`str`] Path to write the output file to.

    @return [`datatable.Frame`] Writes the 'gene_compound_dataset.csv' file to 
      output_dir the returns the table.
    """
    # -- Check the input files exist
    for file in [gene_compound_dataset_file, gene_file, compound_file, dataset_file]:
        if not os.path.exists(file):
            raise FileNotFoundError(f'Could not find the {file}')

    # -- Read in mapping tables
    gene_dt = fread(gene_file)
    compound_dt = fread(compound_file)
    dataset_dt = fread(dataset_file)

    # -- Read in gene_compound_tissue table
    gct_dt = fread(gene_compound_dataset_file)

    # -- Fix names and assign missing columns
    gct_dt.names = {'Gene': 'gene_id', 'Tissue': 'tissue_id', 
        'Drug': 'compound_id', 'FWER_genes': 'FWER_gene'}
    # Determine missing columns and assign them, so we don't have to change code 
    #>when new columns are addeds
    gct_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'tissue_id', 
        'estimate', 'lower', 'upper', 'n', 'tstat', 'fstat', 'pvalue', 'df',
        'fdr', 'FWER_gene', 'FWER_compound', 'FWER_all', 'BF_p_all', 'sens_stat', 
        'mDataType', 'tested_in_human_trials', 'in_clinical_trials'))
    gct_missing_columns = np.setdiff1d(gct_table_columns, np.asarray(gct_dt.names))
    for col in gct_missing_columns:
        gct_dt[col] = None
    gct_dt1 = gct_dt[:, list(gct_table_columns)]
    # Sanity check the columns are there
    if not np.all(gct_table_columns == np.asarray(gct_dt1.names)):
        raise ValueError(f'The build_gene_compound_tissue table',
            ' has missing columns!')

    # -- Map to existing FK ids
    # gene id
    gct_dt1.names = {'gene_id': 'gene_name'}
    gene_dt.names = {'id': 'gene_id', 'name': 'gene_name'}
    gene_dt.key = 'gene_name'
    # NOTE: the g object references the joined tables namespace
    gct_dt1[:, update(gene_id=g.gene_id), join(gene_dt)]

    ## TODO:: rewrite as a helper function
    # regex match failed ids, then assign to the table
    failed_genes = np.unique(gct_dt1[dt.isna(f.gene_id), 'gene_name'].to_numpy().flatten())
    if len(failed_genes) > 0:
        regex_query = f"{'.*|'.join(failed_genes)}.*"
        match_genes = [gene_dt['gene_name'] \
            .to_pandas().gene_name.str.match(f'{gene}.*') for gene in failed_genes]
        failed_match_idx = np.where(~np.asarray([np.any(gene) for gene in match_genes]))[0]
        match_idx = np.where(match_genes)[1]
        # needs to be float64 because Numpy has no NaN for int types... makes no sense!?
        # Pad with NaNs for failed matchess
        gene_ids = gene_dt[match_idx, 'gene_id'].to_pandas().gene_id
        if (len(failed_match_idx) > 1):
            gene_ids = pd.Series(np.insert(gene_ids, failed_match_idx, None), dtype='int32')
        gct_dt1[dt.isna(f.gene_id), update(gene_id=gene_ids)]


    # if (np.any(gct_dt1[:, dt.isna(f.gene_id)].to_numpy())):
    #     warnings.warn('Some gene_ids in gene_compound_tissue are still NA! Dropping'
    #         'the missing rows...')
    #     gct_dt1 = gct_dt1[~dt.isna(f.gene_id), :]

    # # compound id
    # gct_dt1.names = {'compound_id': 'compound_name'}
    # compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
    # del compound_dt[:, 'compound_uid']
    # compound_dt.key = 'compound_name'
    # gct_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]

    # # tissue id
    # gct_dt1.names = {'dataset_id': 'dataset_name'}
    # dataset_dt.names = {'id': 'dataset_id', 'name': 'dataset_name'}
    # dataset_dt.key = 'dataset_name'
    # gct_dt1[:, update(dataset_id=g.dataset_id), join(dataset_dt)]

    # ## TODO: Handle failed tissue mappings?

    # # -- Sort then assign the primary key column
    # ## TODO:: Is there a way to sort by reference? Does setting the key do this
    # gct_dt2 = gct_dt1[:, list(gct_table_columns), sort('gene_id', 'compound_id', 'dataset_id')]
    # gct_dt2[:, update(id=range(1, gct_dt2.nrows + 1))]

    # # Sanity check we didn't lose any rows
    # if not gct_dt.nrows == gct_dt2.nrows:
    #     warnings.warn('The compound_gene_tissue table has lost some rows!')

    # dt.fwrite(gct_dt2, file=os.path.join(output_dir, 'compound_gene_dataset.csv'))
    # return(gct_dt2)


def build_gene_compound_df():
    pass


  