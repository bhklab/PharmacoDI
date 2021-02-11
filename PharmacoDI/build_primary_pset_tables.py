import os
import glob
import pandas as pd
import numpy as np


def build_primary_pset_tables(pset_dict, pset_name):
    """
    Build the tissue, drug, and gene tables for a PSet and return them
    in a dictionary, with table names as the keys.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`dict`] A dictionary of primary tables, with table names as keys
    """
    pset_dfs = {}
    pset_dfs['dataset'] = pd.Series(pset_name, name='name')
    pset_dfs['tissue'] = build_tissue_df(pset_dict)
    pset_dfs['drug'] = build_drug_df(pset_dict)
    pset_dfs['drug_annotation'] = build_drug_annotation_df(pset_dict)
    pset_dfs['cell'] = build_cell_df(pset_dict)
    # Don't make gene table if there are no molecular profiles (TODO - check with chris)
    if 'molecularProfiles' in pset_dict:
        pset_dfs['gene'] = build_gene_df(pset_dict)
        pset_dfs['gene_annotation'] = build_gene_annotation_df(pset_dict)
    return pset_dfs


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

    # Many ENSEMBL gene IDs have the version (ex. ENST00000456328.2 instead
    # of ENST00000456328); remove all version numbers
    gene_df.replace('\.[0-9]$', '', regex=True, inplace=True)
    gene_df.drop_duplicates(inplace=True)
    return gene_df


def build_tissue_df(pset_dict):
    """
    Build a table containing all tissues in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The tissue table
    """
    tissue_df = pd.Series(
        pd.unique(pset_dict['cell']['tissueid']), name='name')
    return tissue_df


def build_drug_df(pset_dict):
    """
    Build a table containing all drugs in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The drug table
    """
    drug_df = pd.Series(pd.unique(pset_dict['drug']['drugid']), name='name')
    return drug_df


def build_gene_annotation_df(pset_dict):
    """
    Build a table mapping each gene in a dataset to its gene annotations.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all gene annotations, mapped to genes
    """
    gene_annotation_df = pd.DataFrame(
        columns=['gene_id', 'symbol', 'gene_seq_start', 'gene_seq_end'], dtype='str')

    for mDataType in pset_dict['molecularProfiles']:
        df = pset_dict['molecularProfiles'][mDataType]['rowData'].copy()
        # Get gene annotation columns
        cols = ['.features']
        if 'Symbol' in df.columns:
            cols.append('Symbol')

        df = df[cols]
        df.rename(columns={'.features': 'gene_id',
                           'Symbol': 'symbol'}, inplace=True)
        gene_annotation_df = gene_annotation_df.append(df)

    # Remove all ENSEMBL gene id version numbers (ex. ENST00000456328.2 instead of ENST00000456328)
    gene_annotation_df['gene_id'].replace('\.[0-9]$', '',
                                          regex=True, inplace=True)
    gene_annotation_df.drop_duplicates(subset=['gene_id'], inplace=True)

    return gene_annotation_df


def build_drug_annotation_df(pset_dict):
    """
    Build a table mapping each drug in a dataset to its drug annotations.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all drug annotations, mapped to drugs
    """
    # Make drug_annotations df
    drug_annotation_df = pset_dict['drug'][[
        'rownames', 'smiles', 'inchikey', 'cid', 'FDA']].copy()
    drug_annotation_df.rename(
        columns={'rownames': 'drug_id', 'cid': 'pubchem', 'FDA': 'fda_status'}, inplace=True)

    return drug_annotation_df


# TODO - confirm that you're using the correct cell id
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
