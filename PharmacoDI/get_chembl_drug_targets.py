from chembl_webresource_client.new_client import new_client
import multiprocessing as mp
import numpy as np
import pandas as pd
import os


def get_chembl_drug_target_mappings(drug_annotation_file, target_file, drug_target_file):
    """
    Get drug target mappings for all drugs in the drug_annotation files (using standard
    inchi key to map between the file and ChEMBL) and all targets in the target file.
    Write to drug_target file and return the resulting DataFrame.

    :drug_annotation_file: the full file path to the drug_annotation table
    :target_file: the full file path to the target csv file
    :drug_target_file: the full file path where the drug target file will be written
    :return: the ChEMBL drug target DataFrame
    """
    # Load drug annotations table
    if not os.path.exists(drug_annotation_file):
        raise FileNotFoundError(f"The file {drug_annotation_file} does not exist!")
    drug_df = pd.read_csv(drug_annotation_file)

    # Get inchikeys of all drugs (not smiles bc capitalization in ChEMBL is inconsistent..)
    inchikeys = list(pd.unique(drug_df['inchikey'].dropna()))

    # Get ChEMBL drugs with matching inchikeys
    print('Getting all drugs from ChEMBL...')
    chembl_drug_df = pd.concat(parallelize(
        inchikeys, get_drugs_by_inchikey, 50))
    chembl_drug_df = pd.merge(
        drug_df[['drug_id', 'inchikey']], chembl_drug_df, on='inchikey', how='inner')
    chembl_drug_df.drop(columns='inchikey', inplace=True)
    molecule_ids = list(chembl_drug_df['molecule_chembl_id'])

    # Get targets from target_file
    if not os.path.exists(target_file):
        print(f"ERROR: the ChEMBL target file {target_file} doesn't exist!\n"
              "Call get_chembl_targets to generate this file.")
    target_df = pd.read_csv(target_file, index_col=0)
    target_df = target_df[['target_chembl_id',
                           'pref_name', 'target_type', 'accession']].copy()
    target_df.drop_duplicates(inplace=True)
    target_ids = list(pd.unique(target_df['target_chembl_id']))

    # Get mappings between drugs (molecule_ids) and targets (target_ids)
    # TODO: I'm not sure if parallelization actually helps for this API call
    print('Getting drug-target mappings from ChEMBL...')
    drug_target_mappings = parallelize(
        molecule_ids, get_drug_target_mappings, 50, target_ids)
    drug_target_df = pd.concat(drug_target_mappings).drop_duplicates()
    drug_target_df = pd.merge(
        drug_target_df, chembl_drug_df, on='molecule_chembl_id')
    drug_target_df = pd.merge(drug_target_df, target_df, on='target_chembl_id')

    # Reorder columns and write to .csv
    drug_target_df = drug_target_df[['drug_id', 'molecule_chembl_id', 
        'target_chembl_id', 'pref_name', 'accession', 'target_type']].copy()
    drug_target_df.to_csv(drug_target_file)
    return drug_target_df


# TODO: add progress bar to parallelize
def parallelize(queries, operation, chunksize, *args):
    """
    Splits queries into chunks of chunksize and then uses a pool to 
    parallelize operation on the query chunks.

    :queries: list of arguments passed to function (in tuples)
    :operation: function being parallelized
    :chuksize: integer representing how many queries each process should handle
    :return: list of results for each query chunk
    """
    chunked_queries = [queries[i:i+chunksize]
                       for i in np.arange(0, len(queries), chunksize)]
    pool = mp.Pool(mp.cpu_count())

    # If operation requires extra args, use pool.starmap instead of pool.map
    if args:
        for i in range(len(chunked_queries)):
            chunked_queries[i] = (chunked_queries[i], *args)
        results = pool.starmap(operation, chunked_queries)
    else:
        results = pool.map(operation, chunked_queries)

    pool.close()
    return results


def get_drugs_by_inchikey(inchikeys):
    """
    Get all drugs in the ChEMBL database with matching inchikeys.

    :inchikeys: A list of inchikeys
    :return: A dataframe of drugs (including ChEMBL ID and inchikey)
    """
    chembl_drug_df = pd.DataFrame(columns=['inchikey', 'molecule_chembl_id'])
    # Initiate connection to ChEMBL molecule table
    molecule = new_client.molecule
    molecules = molecule.filter(molecule_structures__standard_inchi_key__in=inchikeys).only(
        ['molecule_chembl_id', 'molecule_structures'])
    for mol in molecules:
        inchikey = mol['molecule_structures']['standard_inchi_key']
        chembl_drug_df = chembl_drug_df.append({'inchikey': inchikey,
                                                'molecule_chembl_id': mol['molecule_chembl_id']},
                                               ignore_index=True)
    return chembl_drug_df


def get_drug_target_mappings(molecule_ids, target_ids):
    """
    Retrieves mapping between drugs specified by ChEMBL molecule_ids and targets
    specified by ChEMBL target_ids.

    :molecule_ids: A list of ChEMBL drug IDs
    :target_ids: A list of ChEMBL target IDs
    :return: A DataFrame of drug target mappings (molecule ChEMBL ID and target ChEMBL ID)
    """
    # Initiate connection to ChEMBL activity table
    activity = new_client.activity
    results = activity.filter(molecule_chembl_id__in=molecule_ids, target_chembl_id__in=target_ids,
                              pchembl_value__isnull=False).only(['molecule_chembl_id', 'target_chembl_id'])
    mappings = list(results)
    df = pd.DataFrame(mappings)
    return df
