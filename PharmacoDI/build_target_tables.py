import os
import requests
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import datatable as dt
from PharmacoDI.combine_pset_tables import write_table, rename_and_key, join_tables
from PharmacoDI.get_chembl_drug_targets import parallelize


drugbank_file = os.path.join(
    "data", "metadata", "drugbank_targets_has_ref_has_uniprot.csv")
chembl_file = os.path.join('data', 'metadata', 'chembl_drug_targets.csv')
output_dir = os.path.join("data", "demo")


def build_target_tables(drugbank_file, chembl_file, output_dir):
    """
    Build the target and drug target tables using data from Drugbank
    and ChEMBL.

    @param drugbank_file: [`string`] The full file path to Drugbank targets
    @param chembl_file: [`string`] The full file path to ChEMBL targets
    @param output_dir: [`string`] The directory of all final PharmacoDB tables
    @return: None
    """
    # Get Drugbank data
    if not os.path.exists(drugbank_file):
        raise FileNotFoundError(f"The file {drugbank_file} doesn't exist!")
    drugbank_df = pd.read_csv(drugbank_file)
    drugbank_df.rename(columns={'polypeptide.external.identifiers.UniProtKB': 'uniprot_id',
                                'drugName': 'drug_name'}, inplace=True)

    # Get ChEMBL data
    if not os.path.exists(chembl_file):
        raise FileNotFoundError(f"The file {chembl_file} doesn't exist!")
    chembl_df = pd.read_csv(chembl_file, index_col=0)
    chembl_df.rename(columns={'pref_name': 'name',
                              'accession': 'uniprot_id'}, inplace=True)

    target_df = build_target_table(chembl_df, drugbank_df, output_dir)
    build_drug_target_table(chembl_df, drugbank_df, target_df, output_dir)
    build_gene_drug_table(chembl_df, drugbank_df, target_df, output_dir)


def build_target_table(chembl_df, drugbank_df, output_dir):
    """
    Using data from the Drugbank and ChEMBL drug target files and
    the UniProt API, build the target table.

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param output_dir: [`string`] The file path to write the final target table
    @return: [`datatable.Frame`] The target table
    """
    # Combine ChEMBL and Drugbank tables to make target table
    target_df = pd.concat([chembl_df[['name']].copy(),
                           drugbank_df[['name']].copy()])
    target_df.drop_duplicates(inplace=True)

    target_df = write_table(dt.Frame(target_df), 'target', output_dir)
    target_df = rename_and_key(target_df, 'target_id')
    return target_df


def build_drug_target_table(chembl_df, drugbank_df, target_df, output_dir):
    """
    Using data from the Drugbank and ChEMBL drug target files and 
    the target table, build the drug target table.

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param target_df: [`datatable.Frame`] The target table, keyed
    @param output_dir: [`string`] The file path with all final PharmacoDB tables
    @return: [`datatable.Frame`] The drug target table
    """
    # Load drug synonym table from output_dir
    drug_synonym_file = os.path.join(output_dir, 'drug_synonym.csv')
    if not os.path.exists(drug_synonym_file):
        raise FileNotFoundError(
            f"There is no drug synonym file in {output_dir}!")
    drug_syn_df = pd.read_csv(drug_synonym_file, dtype={'drug_id': 'int32'})

    # Join drugbank df with drug table (TODO: are we really using drug name to map?)
    drugbank_df = pd.merge(drugbank_df, drug_syn_df, on='drug_name')
    # TODO: from 7521 down to only 122 rows :/

    # Combine ChEMBL and Drugbank tables to make drug target table
    drug_target_df = pd.concat([chembl_df[['name', 'drug_id']].copy(),
                                drugbank_df[['name', 'drug_id']].copy()])
    drug_target_df.rename(columns={'name': 'target_id'}, inplace=True)
    drug_target_df.drop_duplicates(inplace=True)

    # Join with target table
    drug_target_df = dt.Frame(drug_target_df)
    drug_target_df = join_tables(drug_target_df, target_df, 'target_id')
    # Drop rows with no target_id, drop duplicates
    drug_target_df = drug_target_df[dt.f.target_id >= 1, :]
    drug_target_df = drug_target_df[0, :, dt.by(drug_target_df.names)]

    drug_target_df = write_table(
        drug_target_df, 'drug_target', output_dir, add_index=False)
    return drug_target_df


def build_gene_drug_table(chembl_df, drugbank_df, target_df, output_dir):
    """
    Build a join table...

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param target_df: [`datatable.Frame`] The target table, keyed
    @param output_dir: [`string`] The file path with all final PharmacoDB tables
    @return: [`datatable.Frame`] The gene_target table
    """
    # Get target-uniprot mappings from ChEMBL and Drugbank tables
    gene_target_df = pd.concat([chembl_df[['name', 'uniprot_id']].copy(),
                                drugbank_df[['name', 'uniprot_id']].copy()])
    gene_target_df.rename(columns={'name': 'target_id'}, inplace=True)
    gene_target_df.drop_duplicates(inplace=True)

    # Retrieve Uniprot-ENSEMBL gene ID mappings
    uniprot_ids = pd.Series(pd.unique(gene_target_df['uniprot_id']))
    uniprot_ensembl_mappings = pd.concat(
        parallelize(uniprot_ids, map_uniprot_to_ensembl, 1000))
    uniprot_ensembl_mappings.drop_duplicates(inplace=True)

    # Join gene_target table with gene table based on uniprot-ensembl mappings
    gene_target_df = pd.merge(
        gene_target_df, uniprot_ensembl_mappings, on='uniprot_id')
    gene_target_df.drop(columns=['uniprot_id'], inplace=True)

    # Load and key the gene table from output_dir
    gene_file = os.path.join(output_dir, 'gene.csv')
    if not os.path.exists(gene_file):
        raise FileNotFoundError(f"There is no gene file in {output_dir}!")
    gene_df = dt.fread(gene_file, sep=",")
    gene_df = rename_and_key(gene_df, 'gene_id')

    # Join target table with gene table and target table
    gene_target_df = dt.Frame(gene_target_df)
    gene_target_df = join_tables(gene_target_df, gene_df, 'gene_id')
    gene_target_df = join_tables(gene_target_df, target_df, 'target_id')

    # Drop columns that didn't join and drop duplicates
    gene_target_df = gene_target_df[(
        dt.f.target_id >= 1) & (dt.f.gene_id >= 1), :]
    gene_target_df = gene_target_df[0, :, dt.by(gene_target_df.names)]

    gene_target_df = write_table(
        gene_target_df, 'gene_target', output_dir, add_index=False)
    return gene_target_df


def map_uniprot_to_ensembl(uniprot_ids):
    """
    Use the UniProt API to retrieve the ENSEMBL gene IDs
    corresponding to the UniProt IDS.

    @param uniprot_ids: [`list(string)`] A list of UniProt IDs.
    @return: [`pd.DataFrame`] A table mapping UniProt IDs to ENSEMBL gene IDs.
    """
    # Make API call
    params = {
        'from': 'ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': " ".join(uniprot_ids)
    }
    r = requests.get('https://www.uniprot.org/uploadlists/', params=params)

    # Check that request was successful
    r.raise_for_status()

    # Split text into rows (one row per ID) and build df
    # Exclude first row (header)
    gene_id_df = pd.DataFrame(r.text.split("\n")[1:])

    # Split into two columns, UniprotId and EnsemblId
    gene_id_df = gene_id_df[0].str.split(pat="\t", expand=True)

    # Drop empty rows and rename columns
    gene_id_df.dropna(inplace=True)
    gene_id_df.rename(columns={0: 'uniprot_id', 1: 'gene_id'}, inplace=True)

    return gene_id_df
