import os
import requests
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from datatable import dt, fread, join, f, g
import polars as pl
from PharmacoDI.combine_pset_tables import write_table, rename_and_key, join_tables
from PharmacoDI.get_chembl_compound_targets import parallelize

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_target_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def build_target_tables(drugbank_file, chembl_file, output_dir, compound_synonym_file):
    """
    Build the target and drug target tables using data from Drugbank
    and ChEMBL.

    @param drugbank_file: [`string`] The full file path to Drugbank targets
    @param chembl_file: [`string`] The full file path to ChEMBL targets
    @param output_dir: [`string`] The directory of all final PharmacoDB tables
    @param compound_synonym_file: [`string`] The file path to the compound synonym table
    @return: None
    """
    # Get Drugbank data
    if not os.path.exists(drugbank_file):
        raise FileNotFoundError(f"The file {drugbank_file} doesn't exist!")
    drugbank_df = dt.fread(drugbank_file)
    drugbank_df.names = {'polypeptide.external.identifiers.UniProtKB': 'uniprot_id',
        'drugName': 'compound_name'}

    # Get ChEMBL data
    if not os.path.exists(chembl_file):
        raise FileNotFoundError(f"The file {chembl_file} doesn't exist!")
    chembl_df = dt.fread(chembl_file)
    chembl_df.names= {'pref_name': 'name', 'accession': 'uniprot_id'}
    chembl_df.key = chembl_df.names[0]

    target_df = build_target_table(chembl_df, drugbank_df, output_dir)
    build_compound_target_table(chembl_df, drugbank_df, target_df, output_dir, 
        compound_synonym_file)
    build_gene_target_table(chembl_df, drugbank_df, target_df, output_dir)


@logger.catch
def build_target_table(chembl_df, drugbank_df, output_dir):
    """
    Using data from the Drugbank and ChEMBL drug target files and
    the UniProt API, build the target table.

    @param chembl_df: [`dt.Frame`] The ChEMBL drug target table
    @param drugbank_df: [`dt.Frame`] The DrugBank drug target table
    @param output_dir: [`string`] The file path to write the final target table
    @return: [`dt.Frame`] The target table
    """
    # Combine ChEMBL and Drugbank tables to make target table
    target_df = dt.rbind([chembl_df['name'],
        drugbank_df['name']]).to_pandas()
    target_df.drop_duplicates(inplace=True)
    target_df = dt.Frame(target_df)
    target_df = write_table(target_df, 'target', output_dir)
    target_df = rename_and_key(target_df, 'target_id')
    return target_df


@logger.catch
def build_compound_target_table(chembl_df, drugbank_df, target_df, output_dir, compound_synonym_file):
    """
    Using data from the Drugbank and ChEMBL drug target files and 
    the target table, build the drug target table.

    @param chembl_df: [`dt.Frame`] The ChEMBL drug target table
    @param drugbank_df: [`dt.Frame`] The DrugBank drug target table
    @param target_df: [`datatable.Frame`] The target table, keyed
    @param output_dir: [`string`] The file path with all final PharmacoDB tables
    @param compound_synonym_file: [`string`] The file path to the compound synonym table
    @return: [`dt.Frame`] The drug target table
    """
    # Load compound synonym table from output_dir
    if not os.path.exists(compound_synonym_file):
        raise FileNotFoundError(f"The file {compound_synonym_file} doesn't exist!")
    drug_syn_df = dt.fread(compound_synonym_file)
    # Join drugbank df with drug table (TODO: are we really using drug name to map?)
    del drug_syn_df[:, ['dataset_id', 'id']]
    drug_syn_df = dt.Frame(pl.from_arrow(drug_syn_df.to_arrow()).drop_duplicates().to_arrow())
    drug_syn_df.key = 'compound_name'
    drugbank_df = drugbank_df[:, :, join(drug_syn_df)]
    # TODO: from 7521 down to only 122 rows :/
    # Combine ChEMBL and Drugbank tables to make drug target table
    drug_target_df = pd.concat([chembl_df.to_pandas()[['name', 'compound_id']].copy(),
                                drugbank_df.to_pandas()[['name', 'compound_id']].copy()])
    drug_target_df.rename(columns={'name': 'target_id'}, inplace=True)
    drug_target_df.drop_duplicates(inplace=True)
    # Join with target table
    drug_target_df = dt.Frame(drug_target_df)
    drug_target_df = join_tables(drug_target_df, target_df, 'target_id')
    # Drop rows with no target_id, drop duplicates
    drug_target_df = drug_target_df[dt.f.target_id >= 1, :]
    drug_target_df = drug_target_df[0, :, dt.by(drug_target_df.names)]
    drug_target_df = dt.Frame(
        pl.from_arrow(drug_target_df.to_arrow()) \
            .drop_nulls() \
            .to_arrow())
    drug_target_df = write_table(
        drug_target_df, 'compound_target', output_dir, add_index=False)
    return drug_target_df


@logger.catch
def build_gene_target_table(chembl_df, drugbank_df, target_df, output_dir):
    """
    Build a join table...

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param target_df: [`datatable.Frame`] The target table, keyed
    @param output_dir: [`string`] The file path with all final PharmacoDB tables
    @return: [`datatable.Frame`] The gene_target table
    """
    # Get target-uniprot mappings from ChEMBL and Drugbank tables
    gene_target_df = pd.concat([chembl_df.to_pandas()[['name', 'uniprot_id']],
                                drugbank_df.to_pandas()[['name', 'uniprot_id']]])
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
    gene_file = os.path.join(output_dir, 'gene.jay')
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

    gene_target_df.to_jay(os.path.join(output_dir, 'gene_target.jay'))
    return gene_target_df


# TODO: fix this:
"""
2021-08-16T15:08:15.449197+0000 The following gene_ids failed to map:
2021-08-16T15:08:15.450822+0000     | gene_id           
    | str32             
--- + ------------------
  0 | ENSCAFG00000013762
  1 | ENSG00000137332   
  2 | ENSG00000183311   
  3 | ENSG00000196101   
  4 | ENSG00000198457   
  5 | ENSG00000204490   
  6 | ENSG00000206289   
  7 | ENSG00000206298   
  8 | ENSG00000206308   
  9 | ENSG00000206338   
 10 | ENSG00000206376   
 11 | ENSG00000206439   
 12 | ENSG00000206466   
 13 | ENSG00000206511   
 14 | ENSG00000215077   
  … | …                 
124 | ENSG00000288269   
125 | ENSG00000288299   
126 | ENSG00000288359   
127 | ENSG00000288516   
128 | ENSMUSG00000053004
[129 rows x 1 column]
"""

@logger.catch
def map_uniprot_to_ensembl(uniprot_ids):
    """
    Use the UniProt API to retrieve the ENSEMBL gene IDs
    corresponding to the UniProt IDS.

    @param uniprot_ids: [`list(string)`] A list of UniProt IDs.
    @return: [`dt.Frame`] A table mapping UniProt IDs to ENSEMBL gene IDs.
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
