import os
import warnings
import numpy as np
import re
from datatable import dt, f, g, join, sort, update, fread


def build_gene_compound_tissue_df(gene_compound_tissue_file, gene_file, 
    compound_file, tissue_file, output_dir):
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

    # check for failed genes
    failed_genes = gct_dt1[dt.isna(f.gene_id), 'gene_name'].to_numpy().flatten()
    if len(failed_genes) > 0:
        raise ValueError(f'Could not find the {file}')

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
    ## TODO:: Is there a way to sort by reference?
    gct_dt2 = gct_dt1[:, list(gct_table_columns), 
        sort('gene_id', 'compound_id', 'tissue_id', 'mDataType')]
    gct_dt2[:, update(id=range(1, gct_dt2.nrows + 1))]

    # Sanity check we didn't lose any rows
    if not gct_dt.nrows == gct_dt2.nrows:
        warnings.warn('The compound_gene_tissue table has lost some rows!')

    dt.fwrite(gct_dt2, file=os.path.join(output_dir, 'compound_gene_tissue.csv'))
    return(gct_dt2)


## FIXME:: This function is almost identical to the previous one, refactor
##>into a helper for gene_compound instead of copy pasting code
def build_gene_compound_dataset_df(gene_compound_dataset_file, gene_file, 
    compound_file, dataset_file, output_dir):
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
    gcd_dt = fread(gene_compound_dataset_file)

    # -- Fix names and assign missing columns
    gcd_dt.names = {'gene': 'gene_id', 'drug': 'compound_id', 
        'pSet': 'dataset_id'}
    # Determine missing columns and assign them, so we don't have to change code 
    #>when new columns are addeds
    gcd_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'dataset_id', 
        'estimate', 'lower', 'upper', 'n', 'tstat', 'fstat', 'pvalue', 'df',
        'fdr', 'FWER_gene', 'FWER_compound', 'FWER_all', 'BF_p_all', 'sens_stat', 
        'mDataType', 'tested_in_human_trials', 'in_clinical_trials'))
    gcd_missing_columns = np.setdiff1d(gcd_table_columns, np.asarray(gcd_dt.names))
    for col in gcd_missing_columns:
        gcd_dt[col] = None
    gcd_dt1 = gcd_dt[:, list(gcd_table_columns)]
    # Sanity check the columns are there
    if not np.all(gcd_table_columns == np.asarray(gcd_dt1.names)):
        raise ValueError(f'The build_gene_compound_dataset table',
            ' has missing columns!')

    # -- Map to existing FK ids
    # gene id
    gcd_dt1.names = {'gene_id': 'gene_name'}
    # fix gene ids that still have versions
    gsub = np.vectorize(re.sub)
    gcd_dt1[f.gene_name.re_match('ENS.*[.][0-9]*'), 
        update(gene_name=gsub('[.][0-9]*$', '', gcd_dt1[
            f.gene_name.re_match('ENS.*[.][0-9]*'), :]['gene_name'].to_numpy()))]
    gene_dt.names = {'id': 'gene_id', 'name': 'gene_name'}
    gene_dt.key = 'gene_name'
    # NOTE: the g object references the joined tables namespace
    gcd_dt1[:, update(gene_id=g.gene_id), join(gene_dt)]

    # make sure all genes mapped
    failed_genes = np.unique(gcd_dt1[dt.isna(f.gene_id), 'gene_name'] \
        .to_numpy().flatten())
    if len(failed_genes) > 0:
        warnings.warn(f'The genes: {failed_genes} did not map!')


    if (np.any(gcd_dt1[:, dt.isna(f.gene_id)].to_numpy())):
        warnings.warn('Some gene_ids in gene_compound_dataset are still NA!'
            'Dropping the missing rows...')
        gcd_dt1 = gcd_dt1[~dt.isna(f.gene_id), :]

    # compound id
    gcd_dt1.names = {'compound_id': 'compound_name'}
    compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
    del compound_dt[:, 'compound_uid']
    compound_dt.key = 'compound_name'
    gcd_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]

    # dataset id
    gcd_dt1.names = {'dataset_id': 'dataset_name'}
    dataset_dt.names = {'id': 'dataset_id', 'name': 'dataset_name'}
    dataset_dt.key = 'dataset_name'
    gcd_dt1[:, update(dataset_id=g.dataset_id), join(dataset_dt)]

    # -- Sort then assign the primary key column
    gcd_dt2 = gcd_dt1[:, list(gcd_table_columns), 
        sort('gene_id', 'compound_id', 'dataset_id', 'mDataType')]
    gcd_dt2[:, update(id=range(1, gcd_dt2.nrows + 1))]

    # Sanity check we didn't lose any rows
    if not gcd_dt.nrows == gcd_dt2.nrows:
        warnings.warn('The gene_compound_dataset table has lost some rows!')

    dt.fwrite(gcd_dt2, file=os.path.join(output_dir, 'gene_compound_dataset.csv'))
    return(gcd_dt2)


def build_gene_compound_df(gene_compound_file, gene_file, compound_file, 
    output_dir):
    """
    Build gene_compound table, with pan-cancer and pan-dataset meta-analysis
    results.

    @param gene_compound_file: [`str`] Path to the gene signature .csv file.
    @param gene_file: [`str`] Path to the gene table .csv file.
    @param compound_file: [`str`] Path to the compound table .csv file.
    @param output_dir: [`str`] Path to write the output file to.

    @return [`datatable.Frame`] Writes the 'gene_compound_dataset.csv' file to 
        output_dir the returns the table.
    """
    pass
    # FIXME:: For now this function makes a dummy table since we don't have
    #>real data
    # -- get the all combinations of gene_id x compound_id
    # gene_dt = fread(gene_file)
    # compound_dt = fread(compound_file)
    # gene_id, compound_id = np.meshgrid(gene_dt['id'].to_numpy(), 
    #     compound_dt['id'].to_numpy(), indexing='ij')
    # gene_id = gene_id.flatten()
    # compound_id = compound_id.flatten()

    # # -- build a datatable
    # gc_dt = dt.Frame({
    #     'gene_id': gene_id,
    #     'compound_id': compound_id,
    # })

    # # Determine missing columns and assign them, so we don't have to change code 
    # #>when new columns are addeds
    # gc_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'estimate', 
    #     'lower', 'upper', 'n', 'tstat', 'fstat', 'pvalue', 'df', 'fdr', 
    #     'FWER_gene', 'FWER_compound', 'FWER_all', 'BF_p_all', 'sens_stat', 
    #     'mDataType', 'tested_in_human_trials', 'in_clinical_trials'))
    # gc_missing_columns = np.setdiff1d(gc_table_columns, np.asarray(gc_dt.names))
    # for col in gc_missing_columns:
    #     gc_dt[col] = None
    # gc_dt1 = gc_dt[:, list(gc_table_columns)]
    # # Sanity check the columns are there
    # if not np.all(gc_table_columns == np.asarray(gc_dt1.names)):
    #     raise ValueError(f'The build_gene_compound_dataset table',
    #         ' has missing columns!')

    # gc_dt1['estimate'] = np.random()
    # gc_dt1['n'] =
    # gc_dt1['pvalue'] =


# -- scripts for testing the function
if __name__ == '__main__':
    if not os.getcwd() == 'PharmacoDI_snakemake_pipeline':
        os.chdir('PharmacoDI_snakemake_pipeline')
    sig_data_dir = os.path.join('rawdata/gene_signatures/metaanalysis')
    data_dir = os.path.join('latest')
    gene_compound_tissue_file = os.path.join(sig_data_dir, 'gene_compound_tissue.csv')
    gene_compound_dataset_file = os.path.join(sig_data_dir, 'gene_compound_dataset.csv')
    gene_file = os.path.join(data_dir, 'gene.csv')
    compound_file = os.path.join(data_dir, 'compound.csv')
    tissue_file = os.path.join(data_dir, 'tissue.csv')
    dataset_file = os.path.join(data_dir, 'dataset.csv')
