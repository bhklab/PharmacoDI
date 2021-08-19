import pandas as pd
import numpy as np
from .build_primary_pset_tables import build_cell_df
from .utilities import harmonize_df_columns

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_experiment_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def build_experiment_tables(pset_dict, pset_name, cell_df=None):
    """
    Build the experiment, dose response, and profile tables for a PSet 
    and return them in a dictionary, with table names as the keys.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param cell_df: [`pd.DataFrame`] A table of all the cells in the PSet and their tissues
    @return: [`dict`] A dictionary of experiment-related tables
    """
    pset_dfs = {}
    pset_dfs['experiment'] = build_experiment_df(pset_dict, pset_name, cell_df)
    pset_dfs['dose_response'] = build_dose_response_df(pset_dict, pset_name)
    pset_dfs['profile'] = build_profile_df(pset_dict, pset_name)
    return pset_dfs


@logger.catch
def build_experiment_df(pset_dict, pset_name, cell_df=None):
    """
    Build a table with all the experiments of a PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet 
    @param cell_df: [`pd.DataFrame`] A table of all the cells in the PSet and their tissues
    @return: [`pd.DataFrame`] A table containing all experiments in the dataset
    """
    # Build cell_df if not build already
    ## FIX: Truth value of a DataFrame is ambiguous (Chris)
    if cell_df is None:
        cell_df = build_cell_df(pset_dict)

    # Extract relelvant experiment columns
    experiment_df = pset_dict['sensitivity']['info'][[
        '.rownames', 'cellid', 'drugid']].copy()

    # Rename columns
    experiment_df.rename(
        columns={'.rownames': 'experiment_id', 'cellid': 'cell_id',
            'drugid': 'compound_id'}, inplace=True)

    # Add datset_id column
    experiment_df['dataset_id'] = pset_name

    # Add tissue_id column by joining with cells_df
    experiment_df = pd.merge(experiment_df, cell_df[['name', 'tissue_id']],
        left_on='cell_id', right_on='name', how='left')
    experiment_df = experiment_df[[
        'experiment_id', 'cell_id', 'compound_id', 'dataset_id', 'tissue_id']]

    experiment_df.rename(columns={'experiment_id': 'name'}, inplace=True)

    return experiment_df


# TODO:: Do Python functions pass my reference or copy by default?
# - Python neither; it is 'pass-by-object-reference'
# - When you assign a variable in Python, it creates a name in the namespace; that name contains
#   a reference to the value of that variable (You can use the id() function to see the memory address of an object)
# - When you pass a variable to a function in Python, it creates a new variable, which is a copy of the
#   original variable. There are now two DIFFERENT objects, which contain the memory address of the same value
# - Therefore, if you modify that variable BY REFERENCE within the function, then it will modify the actual value;
#   the original variable will then also be modified (because it points to the same memory location)
# - However, if you reassign the variable name within the function scope, it DOES NOT affect the original variable;
#   instead it modifies the copy of the variable with a reference to a different memory location (the new value)
# - Outside of the function, the original variable still holds a reference to the old memory location; therefore
#   the value of that object is not changed
@logger.catch
def build_dose_response_df(pset_dict, pset_name):
    """
    Build a table that, for each experiment in a dataset, lists the 
    compound that was tested, the doses in which that compound was 
    administered, and the viability responses corresponding to all 
    the doses.

    @param pset_dict: [`dict`] A nested dictionary containing all tables 
        in the PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`pd.DataFrame`] A table with all the dose-response 
        mappings for each experiment
    """
    # Get dose and response info from pset
    dose = pset_dict['sensitivity']['raw.Dose']
    response = pset_dict['sensitivity']['raw.Viability']

    # Rename columns so they can be coerced to int later
    dose_pattern = dose.columns[1][:-1]
    rename_dict = {f'{dose_pattern}{n}': str(n) 
        for n in np.arange(1, dose.shape[0])}
    dose.rename(columns=rename_dict, inplace=True)
    response.rename(columns=rename_dict, inplace=True)

    # Reshape the DataFrames using melt to go from 'wide' to 'long'
    dose = dose.melt(id_vars='.exp_id', value_name='dose',
        var_name='dose_id').dropna()
    dose['dose_id'] = dose.dose_id.astype('int')
    response = response.melt(
        id_vars='.exp_id', value_name='response', var_name='dose_id').dropna()
    response['dose_id'] = response.dose_id.astype('int')

    # Set indices for faster joins (~3x)
    dose.set_index(['.exp_id', 'dose_id'], inplace=True)
    response.set_index(['.exp_id', 'dose_id'], inplace=True)

    # Join on index
    dose_response_df = pd.merge(
        dose, response, left_index=True, right_index=True).reset_index()
    dose_response_df.rename(columns={'.exp_id': 'experiment_id'}, inplace=True)
    dose_response_df.drop(columns=['dose_id'], inplace=True)

    # Add dataset_id for joins
    dose_response_df['dataset_id'] = pset_name

    return dose_response_df


@logger.catch
def build_profile_df(
        pset_dict: dict, 
        pset_name: str, 
        column_dict: dict={'experiment_id': str, 'HS': float, 'Einf': float, 
            'EC50': float, 'AAC': float, 'IC50': float, 'DSS1': float,
            'DSS2': float, 'DSS3': float},
        rename_dict: dict={'.rownames': 'experiment_id', 'einf': 'Einf', 
            'E_inf': 'Einf', 'aac_recomputed': 'AAC',
            'ic50_recomputed': 'IC50', 'ec50': 'EC50'
            }
) -> pd.DataFrame:
    """
    Build the molecular profile tables for PharmacoDB for the specified PSet

    @param pset_dict: A nested dictionary containing all tables in the PSet
    @param pset_name: The name of the PSet
    @param column_dict: Dictionary of required columns and their respective
        types. Note that the id column type will not correspond to the
        table schema, since mapping from this column is done later.
    @param rename_dict: Dictionary mapping from PSet specific column names
        to schema defined column names in PharmacoDB.
    @return: A table containing all drug sensitivity summary statistics for 
        in the PSet.
    """
    # Get profiles df and fix column names
    profile_df = pset_dict['sensitivity']['profiles'].copy()
    profile_df.rename(columns=rename_dict, inplace=True)
    if 'HS' not in profile_df.columns:
        profile_df.rename({'slope_recomputed': 'HS'})
    if 'experiment_id' not in profile_df.columns:
        raise ValueError('No experiment_id column in sensitivity profiles!') 
    # Check required columns
    profile_df = harmonize_df_columns(
        df=profile_df,
        column_dict=column_dict
    )
    # Add dataset_id for joins
    profile_df['dataset_id'] = pset_name
    return profile_df
