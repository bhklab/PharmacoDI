import os
import re
import glob
import numpy as np
import pandas as pd
from datatable import Frame
from PharmacoDI.combine_pset_tables import write_table

output_dir = os.path.join('data', 'demo')
metadata_dir = os.path.join('data', 'metadata')

cell_file = "cell_annotation_all.csv"
tissue_file = "cell_annotation_all.csv"
drug_file = "drugs_with_ids.csv"

def get_metadata(file_name, metadata_dir):
    # Find correct metadata annotations CSV file
    annotations_file = glob.glob(
        os.path.join(metadata_dir, file_name))
    if not annotations_file:
        raise ValueError(
            f'No metadata file named {file_name} could be found in {metadata_dir}')
    
    # Read csv file and return df
    return pd.read_csv(annotations_file[0], index_col=[0])


# --- SYNONYMS TABLES --------------------------------------------------------------------------

def build_cell_synonym_df(cell_file, metadata_dir, output_dir):
    # Get metadata file and cell_df
    cell_metadata = get_metadata(cell_file, metadata_dir)
    cell_df = pd.read_csv(os.path.join(output_dir, 'cell.csv'))

    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_metadata[[
        col for col in cell_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with cell_df
    cell_synonym_df = melt_and_join(cell_columns, 'unique.cellid', cell_df)
    cell_synonym_df = cell_synonym_df.rename(columns={'id': 'cell_id', 'value': 'cell_name'})

    # Add blank col for dataset_id (TODO)
    cell_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(cell_synonym_df)
    df = write_table(df, 'cell_synonym', output_dir)
    return df
    

def build_tissue_synonym_df(tissue_file, metadata_dir, output_dir):
    # Get metadata file and tissue_df (assume taht tissue_df is also in output_dir)
    tissue_metadata = get_metadata(tissue_file, metadata_dir)
    tissue_df = pd.read_csv(os.path.join(output_dir, 'tissue.csv'))

    # Find all columns relevant to tissueid
    pattern = re.compile('tissueid')
    tissue_cols = tissue_metadata[[
        col for col in tissue_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with tissue_df
    tissue_synonym_df = melt_and_join(tissue_cols, 'unique.tissueid', tissue_df)
    tissue_synonym_df = tissue_synonym_df.rename(columns={'id': 'tissue_id', 'value': 'tissue_name'})

    # Add blank col for dataset_id (TODO)
    tissue_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(tissue_synonym_df)
    df = write_table(df, 'tissue_synonym', output_dir)
    return df


def build_drug_synonym_df(drug_file, metadata_dir, output_dir):
    # Get metadata file and drug_df
    drug_metadata = get_metadata(drug_file, metadata_dir)
    drug_df = pd.read_csv(os.path.join(output_dir, 'drug.csv'))

    # Find all columns relevant to drugid
    # Right now only FDA col is dropped, but may be more metadata in the future
    pattern = re.compile('drugid')
    drug_cols= drug_metadata[[
        col for col in drug_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with drugs_df
    drug_synonym_df = melt_and_join(drug_cols, 'unique.drugid', drug_df)
    drug_synonym_df = drug_synonym_df.rename(columns={'id': 'drug_id', 'value': 'drug_name'})

    # Add blank col for dataset_id (TODO)
    drug_synonym_df['dataset_id'] = np.nan

    # Convert to datatable.Frame for fast write to disk
    df = Frame(drug_synonym_df)
    df = write_table(df, 'drug_synonym', output_dir)
    return df


# Helper function for getting all synonyms related to a certain df
def melt_and_join(meta_df, unique_id, join_df):
    """
    @param meta_df: [`Dask DataFrame`] The DataFrame containing all the synonyms (metadata)
    @param unique_id: [`string`] The name of the column in the metadata containing the unique IDs
    @param join_df: [`Dask DataFrame`] THe DataFrame containing the primary keys that will be used as
        foreign keys in the new synonyms df

    @return [`DataFrame`] The synonys dataframe, with a PK, FK based on join_df, and all unique synonyms
    """
    # Convert wide meta_df to long table
    # Drop 'variable' col (leave only unique ID and synonyms), drop duplicates
    synonyms = pd.melt(meta_df, id_vars=[unique_id])[
        [unique_id, 'value']].drop_duplicates()

    # Drop all rows where value is NA
    synonyms = synonyms[synonyms['value'].notnull()]

    # Join with join_df based on unique_id
    synonyms = pd.merge(synonyms, join_df, left_on=unique_id,
                        right_on='name', how='inner')[['id', 'value']]
    synonyms['id'] = synonyms['id'].astype('int')

    return synonyms

