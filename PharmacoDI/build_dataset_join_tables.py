import pandas as pd
from PharmacoDI.build_all_pset_tables import build_cell_df, build_drug_df, build_tissue_df


def build_dataset_join_dfs(pset_dict, pset_name, primary_dfs={}):
    """
    Builds join tables summarizing the cell lines, tissues, and compounds 
    in this PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param primary_dfs: [`pd.DataFrame`] A dictionary containing primary tables (cell,
                                        tissue, compound/drug) from this PSet
    @return: [`dict(string: pd.DataFrame)`] A dictionary of the join tables, with 
                                            table names as keys
    """
    cell_df = primary_dfs['cell'] if 'cell' in primary_dfs else None
    tissue_df = primary_dfs['tissue'] if 'tissue' in primary_dfs else None
    compound_df = primary_dfs['drug'] if 'drug' in primary_dfs else None

    join_dfs = {}
    join_dfs['dataset_cell'] = build_dataset_cell_df(
        pset_dict, pset_name, cell_df)
    join_dfs['dataset_tissue'] = build_dataset_tissue_df(
        pset_dict, pset_name, tissue_df)
    join_dfs['dataset_compound'] = build_dataset_compound_df(
        pset_dict, pset_name, compound_df)
    return join_dfs


def build_dataset_cell_df(pset_dict, pset_name, cell_df=None):
    """
    Builds a join table summarizing the cell lines in this PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param cell_df: [`pd.DataFrame`] The cell table for this PSet
    @return: [`pd.DataFrame`] The join table with all cell lines in this PSet
    """
    if cell_df is None:
        cell_df = build_cell_df(pset_dict)

    dataset_cell_df = pd.DataFrame(
        {'dataset_id': pset_name, 'cell_id': cell_df['name']})

    return dataset_cell_df


def build_dataset_tissue_df(pset_dict, pset_name, tissue_df=None):
    """
    Builds a join table summarizing the tissues in this PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param tissue_df: [`pd.DataFrame`] The tissue table for this PSet
    @return: [`pd.DataFrame`] The join table with all tissues in this PSet
    """
    if tissue_df is None:
        tissue_df = build_tissue_df(pset_dict)

    dataset_tissue_df = pd.DataFrame(
        {'dataset_id': pset_name, 'tissue_id': tissue_df})

    return dataset_tissue_df


def build_dataset_compound_df(pset_dict, pset_name, compound_df=None):
    """
    Builds a join table summarizing the drugs/compounds in this PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param compound_df: [`pd.DataFrame`] The drug/compound table for this PSet
    @return: [`pd.DataFrame`] The join table with all compounds/drugs in this PSet
    """
    if compound_df is None:
        compound_df = build_drug_df(pset_dict)

    dataset_compound_df = pd.DataFrame(
        {'dataset_id': pset_name, 'compound_id': compound_df})

    return dataset_compound_df
