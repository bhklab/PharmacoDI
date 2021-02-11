import os
#import PharmacoDI as di
import pandas as pd
import numpy as np
import urllib
import requests
from io import StringIO
from lxml import etree

def get_target_annotations(pset, annot_dir):
    """
    Annotate a the 'TARGET' in the 'drug' slot of a PSet object using mapping from the UniProt idenitifer
    mapping tool API.

    :param pset:
    :param annot_dir:
    :return:
    """
    # Read in drug target annotations and gene annotations
    drug_targets = pd.read_csv(os.path.join(annot_dir, 'drugbank_drug_targets_all.csv'))
    rnaseq_df = pset.get("molecularProfiles").get("Kallisto_0.46.1.rnaseq").get("elementMetadata")

    # Map genes to drugbank drug ids
    genes_to_drugs = pd.merge(drug_targets.loc[:, ['Name', 'Gene Name', 'Drug IDs']],
                              rnaseq_df.loc[:, ['gene_name', 'gene_id']],
                              left_on='Gene Name', right_on='gene_name')

    # Annotate the genes


    # Expand list columns into rows and annotate drugs
    genes_to_drugs['Drug IDs'] = [str.split(ids, '; ') for ids in genes_to_drugs['Drug IDs'].values]
    genes_to_drugs = genes_to_drugs.explode('Drug IDs')


    # Write to disk if necessary.
    file_path = os.path.join(annot_dir, 'drugbank_drug_to_gene_mappings.csv')
    if not os.isfile(file_path):
        pd.write_csv(genes_to_drugs, file_path)
    pass



def query_uniprot_mapping_api(ids, convert_from="ACC+ID", to="ENSEMBL_ID"):
    """
    Query the UniProt mapping API to convert between different gene/protein identifiers.
    Defaults to converting UniProt AC or ID into ENESBLE gene id. The API, however,
    supports a wide range of identifier conversions such as 'GENENAME', "EMBL", and
    "P_ENTREZGENEID".

    Unmatched ids fail silently and will be excluded from the resulting DataFrame. They
    can be retrieved by redoing your query in manually at https://www.uniprot.org/uploadlists/.

    Documentation for other potential conversions are available at:
        https://www.uniprot.org/help/api_idmapping

    :param ids: [`list`, `tuple` or `ndarray`] An iterable sequence type containing the gene/protein
        identifiers as strings.
    :param convert_from: [`string`] The UniProt abbreviation for a format of gene/protein identifier.
        Defaults to 'ACC+ID'.
    :param to: [`string`] The Uniprot abbreviation for the desired gene/protein identifier. Defaults
        to 'ENSEMBL'.
    :return: [`DataFrame`] With the columns 'From' and 'To', mapping from the current id to the selected
        id type based on annotation in the UniProt database.
    """
    # API URL
    url = 'https://www.uniprot.org/uploadlists/'

    # Build the query
    query = ' '.join([str(id) for id in ids])
    params = {
        'from': convert_from,
        'to': to,
        'format': 'tab',
        'query': query
    }
    query_string = urllib.parse.urlencode(params)

    # Encode query and retrieve get request
    query_string = query_string.encode('utf-8')
    req = urllib.request.Request(url, query_string)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    # Decode query and convert to DataFrame
    table_data = StringIO(response.decode("utf-8"))
    uniprot_mapping_df = pd.read_table(table_data, sep="\t")

    return(uniprot_mapping_df)

