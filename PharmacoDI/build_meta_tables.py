import os
import warnings
import numpy as np
import re
from datatable import dt, f, g, join, sort, update, fread

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_meta_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
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
    for fl in [gene_compound_tissue_file, gene_file, compound_file, tissue_file]:
        if not os.path.exists(fl):
            raise FileNotFoundError(f'Could not find the {fl}')

    # -- Read in mapping tables
    gene_dt = fread(gene_file)
    compound_dt = fread(compound_file)
    tissue_dt = fread(tissue_file)

    # -- Read in gene_compound_tissue table
    gct_dt = fread(gene_compound_tissue_file)

    # -- Fix names and assign missing columns
    if np.all(np.isin(np.asarray(('Gene', 'Tissue', 'Drug', 'FWER_genes')),
            np.asarray(gct_dt.names))):
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
    gct_dt1 = gct_dt[:, [*gct_table_columns, 'gene_name', 'compound_name', 
        'tissue_name']]
    # Sanity check the columns are there
    if not np.all(np.isin(gct_table_columns, np.asarray(gct_dt1.names))):
        raise ValueError(f'The gene_compound_tissue table',
            ' has missing columns!')

    # -- Map to existing FK ids
    # gene id
    gene_dt.names = {'id': 'gene_id', 'name': 'gene_name'}
    gene_dt.key = 'gene_name'
    # NOTE: the g object references the joined tables namespace
    gct_dt1[:, update(gene_id=g.gene_id), join(gene_dt)]

    # check for failed genes
    failed_genes = gct_dt1[dt.isna(f.gene_id), 'gene_name'].to_numpy().flatten()
    if len(failed_genes) > 0:
        raise ValueError(f'Genes {failed_genes} failed to map!')

    if (np.any(gct_dt1[:, dt.isna(f.gene_id)].to_numpy())):
        warnings.warn('Some gene_ids in gene_compound_tissue are still NA! Dropping'
            'the missing rows...')
        gct_dt1 = gct_dt1[~dt.isna(f.gene_id), :]
    del gct_dt1[:, 'gene_name']

    # compound id
    compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
    del compound_dt[:, 'compound_uid']
    compound_dt.key = 'compound_name'
    gct_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]

    # tissue id
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

    gct_dt2.to_jay(os.path.join(output_dir, 'gene_compound_tissue.jay'))



## FIXME:: This function is almost identical to the previous one, refactor
##>into a helper for gene_compound instead of copy pasting code
@logger.catch
def build_gene_compound_dataset_df(gene_compound_dataset_file, gene_file, 
    compound_file, dataset_file, output_dir, compound_names):
    """
    Build gene_compound_dataset table (description?)

    @param gene_compound_dataset_file: [`str`] Path to the gene signature .csv file.
    @param gene_file: [`str`] Path to the gene table .csv file.
    @param compound_file: [`str`] Path to the compound table .csv file.
    @param dataset_file: [`str`] Path to the tissue table .csv file.
    @param output_dir: [`str`] Path to write the output file to.
    :param compound_name: [`str`] Path to an optional .csv file mapping 
        updated compound names to the dataset. This is to ensure that corrected
        compound annotations still make it into the database without the need
        to rerun all the gene signatures

    @return [`datatable.Frame`] Writes the 'gene_compound_dataset.csv' file to 
        output_dir the returns the table.
    """
    # -- Check the input files exist
    for fl in [gene_compound_dataset_file, gene_file, compound_file, dataset_file]:
        if not os.path.exists(fl):
            raise FileNotFoundError(f'Could not find the {fl}')

    # -- Read in mapping tables
    gene_dt = fread(gene_file)
    compound_dt = fread(compound_file)
    dataset_dt = fread(dataset_file)

    # -- Read in gene_compound_tissue table
    gcd_dt = fread(gene_compound_dataset_file)

    # -- Fix names and assign missing columns
    gcd_dt.names = {'gene': 'gene_id', 'compound': 'compound_id',
        'dataset': 'dataset_id', 'lower': 'lower_analytic',
        'upper': 'upper_analytic', 'pvalue': 'pvalue_analytic',
        'fdr': 'fdr_analytic'}
    del gcd_dt[:, ['significant', 'tissue']]

    # Determine missing columns and assign them, so we don't have to change code 
    #>when new columns are addeds
    gcd_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'dataset_id',
        'estimate', 'lower_analytic', 'upper_analytic', 'lower_permutation',
        'upper_permutation', 'n', 'pvalue_analytic', 'pvalue_permutation', 
        'df', 'fdr_analytic', 'fdr_permutation', 'significant_permutation',
        'permutation_done', 'sens_stat', 'mDataType'))
    gcd_missing_columns = np.setdiff1d(gcd_table_columns, 
        np.asarray(gcd_dt.names))
    for col in gcd_missing_columns:
        gcd_dt[col] = None
    gcd_dt1 = gcd_dt[:, list(gcd_table_columns)]
    # Sanity check the columns are there
    if not np.all(gcd_table_columns == np.asarray(gcd_dt1.names)):
        raise ValueError(f'The build_gene_compound_dataset table',
            ' has missing columns!')

    gcd_dt1[:, update(sens_stat='AAC', permutation_done=0)]

    # -- Map to existing FK ids
    # gene id
    gcd_dt1.names = {'gene_id': 'gene_name'}
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

    # fix compound names 
    ## FIXME:: Remove this when gene signatures are regenerated
    ## START patch
    fix_names_df = dt.fread(compound_names)
    fix_names_df[f.dataset == "GDSC_2020(v1-8.2)", update(dataset="GDSC_v1")]
    fix_names_df[f.dataset == "GDSC_2020(v2-8.2)", update(dataset="GDSC_v2")]
    fix_names_df.names = {"drugid": "compound_name",
        "unique.drugid": "compound_id", "dataset": "dataset_id"}
    fix_names_df.key = ["compound_name", "dataset_id"]
    gcd_dt1.names = {'compound_id': 'compound_name'}
    gcd_dt1[~dt.isna(g.compound_id), update(compound_name=g.compound_id),
        join(fix_names_df)]
    ## END patch

    # compound id
    compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
    del compound_dt[:, 'compound_uid']
    compound_dt.key = 'compound_name'
    gcd_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]

    if np.any(gcd_dt1[:, dt.isna(f.compound_id)].to_numpy()):
        warnings.warn("Some compound_ids in gene_compound_dataset are stll "
            "NA! Dropping the missing rows...")
        gcd_dt1 = gcd_dt1[~dt.isna(f.compound_id)]

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

    gcd_dt2.to_jay(os.path.join(output_dir, 'gene_compound_dataset.jay'))

@logger.catch
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
