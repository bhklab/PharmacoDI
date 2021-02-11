# import pandas as pd
# import numpy as np
# import glob
# import os
# import re
#
# def pset_to_db_tables(pset, save_dir=os.path.join("..", "data", "procdata"),
#                       annot_dir=os.path.join("..", "data", "metadata"),
#                       api_url="https://www.orcestra.ca/api/psets/canonical"):
#     """
#     Take in a Python dictionary of a PSet and convert it to database tables for PharmacoDB, with a .csv file for
#     each table saved to `save_dir`.
#
#     :param pset: [dict] a nested dictionary containing data from a PSet, as returned by the `convert_pset_to_py`
#         function.
#     :param save_dir: [string] path to the directory where the database table .csv files should be written
#     :param api_url: [string] URL to fetch available PSets from Orecstra, defaults to current URL
#     :return: [None] writes .csv files to `save_dir`
#     """
#     canonical_names = pd.read_json(api_url).name
#     name = re.sub('_', '.*', pset.get('annotation').get('name')[0])
#
#     # ---- Primary tables ----
#
#     ## ---- dataset
#     dataset = pd.DataFrame.from_dict({
#         "id": [np.where(canonical_names.str.match(name))[0][0] + 1], # Wrap single values in list to make 1 row DF
#         "name": [pset.get('annotation').get('name')[0]],
#     })
#
#
#     ## ---- tissue
#     tissue = pd.DataFrame({
#         'id': np.arange(len(pset.get('cell')['tissueid'].unique())) + 1,
#         'name': pset.get('cell')['tissueid'].unique(),
#     })
#
#
#     ## ---- cell
#     # preprocess cell_info for unique cellid:tissueid combinations and map tissue_id to tissue name
#     cell_info = pset.get('cell')
#     tissue_id_map = dict(zip(tissue.name, tissue.id))
#     cell_info['tissue_id'] = [tissue_id_map[tissue_id] for tissue_id in cell_info.tissueid.to_numpy()] # .to_numpy() should speed up iteration
#
#     # build cell db table
#     cell = cell_info[['cellid', 'tissue_id']].copy()
#     cell.insert(0, 'id', np.arange(1, len(cell_info.cellid) + 1))
#     cell.columns = ['id', 'name', 'tissue_id']
#
#
#     ## ---- compound
#     compound = pd.DataFrame({
#         "id": np.arange(1, len(pset.get('drug').drugid) + 1),
#         "name": pset.get('drug').drugid
#     })
#
#
#     ## ---- gene
#     rnaseq_df = pset.get("molecularProfiles").get("Kallisto_0.46.1.rnaseq").get("elementMetadata")
#     rna_df = pset.get("molecularProfiles").get("rna").get("elementMetadata")
#
#     gene_annots = {re.sub(r"^.*/Gencode\.[^\.]*\.|.csv$", "",  file): pd.read_csv(file) for
#                    file in glob.glob(os.path.join(annot_dir, "Gencode.v33.annotation.*.csv"))}
#
#     # Extract individual annotations from dict
#     feature_genes, feature_transcript, tx2gene = gene_annots.values()
#
#     # Coerce pset datatypes to match gene annotations (some weird stuff happens when converting from R object)
#     rnaseq_df = rnaseq_df.astype(feature_genes.dtypes.to_dict())
#
#     gene = pd.DataFrame({
#         "id": np.arange(1, len(rnaseq_df.rownames) + 1),
#         "name": rnaseq_df.rownames
#     })
#
#     ensg = gene.name.apply(lambda name: re.sub(r"\..*$", "", name))
#
#
#     ## ---- target
#
#
#     target = pd.DataFrame({
#         "id": np.arange(1, len(pset.get("drug")['TARGET'])),
#         "name": pset.get("drug")["TARGET"],
#     })
#
#
#     # ---- Annotation tables ----
#
#     ## ---- cellosaurus
#     cellosaurus = pd.DataFrame({
#         'id': np.arange(1, len(cell.id) + 1), ## FIXME: Can't we just use cell_id as PK?
#         'cell_id': cell.id,
#         'identifier': cell_info['COSMIC.identifier'],
#         'accession': cell_info['Cellosaurus.Accession.id'],
#         'as': np.repeat(None, len(cell.id)),
#         'sy': np.repeat(None, len(cell.id)),
#         'dr': np.repeat(None, len(cell.id)),
#         'rx': np.repeat(None, len(cell.id)),
#         'ww': np.repeat(None, len(cell.id)),
#         'cc': np.repeat(None, len(cell.id)),
#         'st': np.repeat(None, len(cell.id)),
#         'di': np.repeat(None, len(cell.id)),
#         'ox': np.repeat(None, len(cell.id)),
#         'hi': np.repeat(None, len(cell.id)),
#         'oi': np.repeat(None, len(cell.id)),
#         'sx': np.repeat(None, len(cell.id)),
#         'ca': np.repeat(None, len(cell.id)),
#     })
#
#     ## TODO: Parse cellosaurus text file into DataFrame
#
#     compound_annotation = pd.DataFrame({
#         "compound_id": compound.id,
#         "name": compound.name,
#         'smiles': pset.get('drug').smiles,
#         'inchikey': pset.get('drug').inchikey,
#         'pubchem': pset.get('drug').cid,
#         'fda': [0 if fda < 0 else 1 for fda in pset.get('drug').FDA]
#     })
#
#     compound_target = pd.DataFrame({
#         "id": 1,
#         "name": []
#     })
#
#     gene_annotations = pd.DataFrame({
#         "id": gene.id,
#         "symbol": rnaseq_df.gene_name,
#         "gene_seq_start": rnaseq_df.start,
#         "gene_seq_end": rnaseq_df.end,
#     })
#
#     dataset_statistics = pd.DataFrame({
#         "id": 1,
#         "name": []
#     })
#
#
#     # ---- Derived tables ----
#
#     compound_drug = pd.read_csv(glob.glob(os.path.join(annot_dir, "gene_signature", "GDSCV1"))
#
#
#     # ---- Join tables ----
#
#     # ---- Write to .csv ----
