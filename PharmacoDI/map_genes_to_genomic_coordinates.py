import datatable as dt
from datatable import update, f
import numpy as np

def map_genes_to_genomic_coordinates(gene_path, gene_annotation_path, gencode_path):
    """
    Reads in the gene and gene_annotation tables along with a gencode annotation file
    and uses the gencode annotations to assign genomic coordinates to genes in gene_annotations
    before writing the updated table back to disk.

    @param gene_path [`string`] Path to the gene table .csv
    @param gene_annotation_path [`string`] Path to the gene_annotation table .csv
    @param gencode_path [`string`] Path tot he genecode annotation table .csv
    
    @return [None] Modifies gene_annotation table and writes .csv to disk
    """

    # -- Load in the required data
    gene = dt.fread(gene_path)
    gene_annot = dt.fread(gene_annotation_path)
    gencode = dt.fread(gencode_path)

    # -- Add gene name back to gene_annotation
    gene.key = 'id'
    # join columns need the same name, stupid...
    gene_annot[:, update(id = f.gene_id)]
    gene_annot.key = 'id'
    gene_a = gene_annot[:, :, dt.join(gene)]

    # -- Prepocess the genomic coordinates
    gencode.names = {'gene_id': 'name'}
    gencode =  gencode[:, ['name', 'start', 'end', 'strand', 'seqnames']]
    gencode.key = 'name'
    gene_a.key = 'name'

    # -- Map coordinates to gene_annotations, check that nothing went wrong
    gene_annotation = gene_a[:, :, dt.join(gencode)].copy()
    # sanity check the mappings didn't get messed up
    if not np.all(gene_annotation['name'].to_numpy() == gene['name'].to_numpy()):
        raise ValueError('The gene_annotation table got mangled while trying to map'
            'genomic coordinates!')

    # -- Clean up the table and write to disk
    gene_annotation[:, update(gene_seq_start=f.start, gene_seq_end=f.end, chr=f.seqnames)]
    del gene_annotation[:, ['name', 'id', 'start', 'end', 'seqnames']]

    dt.fwrite(gene_annotation, file=gene_annotation_path)
