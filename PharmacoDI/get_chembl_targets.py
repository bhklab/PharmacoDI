from chembl_webresource_client.new_client import new_client
import pandas as pd
import os


def get_chembl_targets(target_file):
    """
    Get all ChEMBL targets in humans and write them to a table.

    :target_file: full file path to where the target table should be written
    :return: a DataFrame containing all targets from ChEMBL
    """
    print('Getting all targets from ChEMBL...')
    # Initiate connection to target table
    target = new_client.target

    # Parse all human targets into a DataFrame
    target_result = target.filter(organism__in=['Homo sapiens'])
    results = list(target_result)
    target_df = pd.DataFrame(results)

    # Explode list-like columns into separate rows, duplicating the index
    object_columns = target_df.dtypes[target_df.dtypes ==
                                      'object'].index.values
    for column in object_columns:
        target_df = target_df.explode(column)

    # Drop any targets without cross refs
    target_df = target_df.query("cross_references.notna()").copy()

    # Expand cols with dtype dict into their cols for each key
    for col in ['cross_references', 'target_components']:
        dict_col = pd.json_normalize(target_df[col], max_level=0)
        dict_col.index = target_df.index
        target_df.drop(columns=col, inplace=True)
        target_df = pd.merge(target_df, dict_col,
                             left_index=True, right_index=True)

    # Drop target component cols with dicts (for now; TODO: keep them and also expand them?)
    target_df.drop(columns=['target_component_synonyms',
                            'target_component_xrefs'], inplace=True)
    target_df.drop_duplicates(inplace=True)

    target_df.to_csv(target_file)
    return target_df
