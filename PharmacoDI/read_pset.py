import glob
import os
import re
import numpy as np
import pandas as pd
import swifter  # Library to parallelize apply statements automagically

pset_name = 'GDSC_v1'
file_path = os.path.join('data', 'rawdata')
slot_names = ['curation', 'drug', 'molecularProfiles',
              'sensitivity', 'annotation', 'cell']


def read_pset(pset_name, file_path, slot_names=['curation', 'drug', 'molecularProfiles', 'sensitivity', 'annotation', 'cell']):
    """
    Read in all the data associated with a PharmacoSet from the .csv files exported by the writeToCsv method from rPharmacoDI.
    
    @param pset_name: [`string`] Name of the PharmacoSet object as it appears in the directory name for files exported using
        rPharmacoDIs writeToCsv function.
    @param file_path: [`string`] Path to the directory where PharmacoSet data is stored.
    @param slot_names: [`list`] A list of PharmacoSet slots to read in. Defaults to all slots in a PharmacoSet.
    
    @return: [`DataFrame`] A DataFrame with the columns slot, containing the name of the slot data is from, one or more subitem columns, indicating
        the subitem in a slot a data came from, file path, the path data was read from and data, which contains the object associated with the
        selected slot and subitem.
    """
    # Use regex to find the appropriate files
    pset_dir = glob.glob(f'{os.path.join(file_path, pset_name)}_PSet')[0]
    if pset_dir is None:
        raise ValueError(
            f'No PSet directory named {pset_name} could be found in {file_path}')
    
    # List al files for the select PSet, then split on $ to make a DataFrame
    pset_files = pd.Series(os.listdir(pset_dir))
    pset_files_df = pset_files.str.split('$', expand=True)
    
    # Build the file paths to read in data for each row of the DataFrame
    pset_files_df['file_paths'] = [os.path.join(
        pset_dir, file_name) for file_name in pset_files]
    
    # Rename columns
    pset_files_df.columns = [
        'slot', *[f'subitems{i}' for i in range(1, pset_files_df.shape[1] - 1)], 'file_paths']
    
    # Read in PSet data
    pset_files_df['data'] = pset_files_df['file_paths'].swifter.apply(
        read_pset_file)
    
    # Drop file_paths column by reference
    pset_files_df.drop('file_paths', axis='columns', inplace=True)

    # Process id columns to use the proper slot names
    pset_files_df.iloc[:, 0:-1] = pset_files_df.iloc[:, 0:-
                                                     1].apply(lambda col: col.str.replace('.*@|.csv.gz$|.txt', ''))
    
    return pset_files_df


# ---- Helper methods for read_csv
def read_pset_file(file_path):
    """Deal with text files which can't be read in with pd.read_csv"""
    if '.csv.gz' in file_path:
        return pd.read_csv(file_path)
    elif '.txt' in file_path:
        with open(file_path) as f:
            text = [line for line in f]
        return text
    else:
        raise ValueError(
            f'Unsupported file type passed to this function from: {file_path}')


def pset_df_to_nested_dict(df):
    """
    Recrusively turn unique values in the first column of a DataFrame into dictionary keys until only 1 column remains.

    @param df: [`DataFrame`] With one or more ID columns before a data column.
    @return: [`dict`] A nested dict with levels equal to the number of columns of the original DataFrame minus one
    """
    for key_check in pd.unique(df[df.columns[0]]):
        if key_check is not None:
            # Recursively nest the first column as key and remaning columns as values
            return {key: pset_df_to_nested_dict(df.loc[df[df.columns[0]] == key, df.columns[1:]]) if
                    df[1:].shape[1] > 2 else df.loc[df[df.columns[0]] == key, 'data'].values[0] for
                    key in pd.unique(df[df.columns[0]])}
        else:
            # Return the data column if there no no key
            return df.iloc[:, -1].values[0]
