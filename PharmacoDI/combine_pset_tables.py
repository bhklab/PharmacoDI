import glob
import os
import re
import numpy as np
import pandas as pd
from datatable import dt, fread, f, g, join, by, sort, update
import polars as pl
from polars import col

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/combine_all_pset_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)


@logger.catch
def combine_all_pset_tables(
    data_dir: str, 
    output_dir: str, 
    compound_meta_file: str
) -> dict:
    """
    Combine all PSet tables into the final PharmacoDB tables.

    @param data_dir: [`string`] The file path to read the PSet tables
    @param output_dir: [`string`] The file path to write the final tables
    @return: [`dict(string: datatable.Frame)`] A dictionary of all some of the
        final tables, with names as keys, to be used for later joins
    """
    logger.info("Combining all PSet tables...")
    join_dfs = combine_primary_tables(
        data_dir, 
        output_dir, 
        compound_meta_file
    )
    combine_secondary_tables(data_dir, output_dir, join_dfs)
    combine_experiment_tables(data_dir, output_dir, join_dfs)


@logger.catch
def combine_primary_tables(
    data_dir: str, 
    output_dir: str, 
    compound_meta_file: str
) -> dict:
    """
    Build all the primary tables, i.e., tables that require no joins,
    and return them in a dictionary.

    @param data_dir: [`string`] The file path to read the PSet tables
    @param output_dir: [`string`] The file path to write the final tables
    @return: [`dict(string: datatable.Frame)`] A dictionary of all the primary
        tables, with names as keys
    """
    # Load, concatenate, and write primary tables to disk
    tissue_df = load_table("tissue", data_dir)
    tissue_df[dt.isna(f.name), update(name="NA")]
    tissue_df = tissue_df[:, :, sort(f.name)]
    write_table(tissue_df, "tissue", output_dir)
    gene_df = load_join_write("gene", data_dir, output_dir)
    dataset_df = load_join_write("dataset", data_dir, output_dir)

    # Annotate compounds
    compound_df = load_table("compound", data_dir)
    compound_meta_df = dt.fread(compound_meta_file)
    compound_meta_df.names = {"unique.drugid": "name", 
        "PharmacoDB.uid": "compound_uid"}

    # Add compound_uid to compound table and write to disk
    compound_meta_df.key = "name"
    compound_df = compound_df[:, :, dt.join(compound_meta_df)]
    compound_df = write_table(compound_df, "compound", output_dir)

    # Transform tables to be used for joins
    dfs = {}
    dfs["tissue"] = rename_and_key(tissue_df, "tissue_id")
    dfs["compound"] = rename_and_key(compound_df, "compound_id")
    dfs["gene"] = rename_and_key(gene_df, "gene_id")
    dfs["dataset"] = rename_and_key(dataset_df, "dataset_id")
    return dfs


@logger.catch
def combine_secondary_tables(data_dir, output_dir, join_dfs):
    """
    Build all secondary tables, i.e., all tables that have foreign keys corresponding
    to primary keys of primary tables. The function reads PSet tables from 
    data_dir, concatenates and joins them with tables from join_dfs, and 
    writes them to output_dir.

    @param join_dfs: [`dict(string: datatable.Frame)`] A dictionary of all the primary
        tables, with names as keys
    @param data_dir: [`string`] The file path to read the PSet tables
    @param output_dir: [`string`] The file path to write the final tables
    @return: [`dict(string: datatable.Frame)`] The updated dictionary of join tables
    """
    # Build cell table and add to join_dfs dictionary
    cell_df = load_join_write(
        "cell", data_dir, output_dir, ["tissue"], join_dfs)
    join_dfs["cell"] = rename_and_key(cell_df, "cell_id")

    # Build compound annotation table
    load_join_write("compound_annotation", data_dir,
        output_dir, ["compound"], join_dfs, add_index=False)
    # Build gene annotation table
    gene_annot_df = load_table("gene_annotation", data_dir)
    # Join the other way so that genes that got cut out are included back in
    gene_annot_df.key = "gene_id"
    gene_annot_df = join_tables(join_dfs["gene"], gene_annot_df, "gene_id")
    write_table(gene_annot_df, "gene_annotation", output_dir, add_index=False)

    # Build join tables
    load_join_write("dataset_cell", data_dir, output_dir,
                    ["dataset", "cell"], join_dfs, add_index=False)
    load_join_write("dataset_tissue", data_dir, output_dir,
                    ["dataset", "tissue"], join_dfs, add_index=False)
    # TODO: temporary workaround for dataset_compound until we standardize 
    dataset_compound_df = load_table("dataset_compound", data_dir)
    dataset_compound_df = join_tables(
        dataset_compound_df, join_dfs["dataset"], "dataset_id")
    compound_df = join_dfs["compound"].copy()
    compound_df.names = {"compound_id": "compound_id"}
    dataset_compound_df = join_tables(
        dataset_compound_df, compound_df, "compound_id")
    dataset_compound_df = write_table(
        dataset_compound_df, "dataset_compound", output_dir, add_index=False)

    # Build other secondary tables
    load_join_write("mol_cell", data_dir, output_dir,
                    ["cell", "dataset"], join_dfs)
    # mol_cells has Kallisto. not sure why. from CTRPv2 (TODO)
    load_join_write("dataset_statistics", data_dir,
                    output_dir, ["dataset"], join_dfs)
    return join_dfs


@logger.catch
def combine_experiment_tables(data_dir, output_dir, join_dfs):
    """
    Load and process experiment table, then use it to build the dose response
    and profile tables. Drop the "name" column from the experiment table before
    writing to a CSV.

    @param join_dfs: [`dict(string: datatable.Frame)`]
    @param data_dir: [`string`] The file path to the PSet tables
    @param output_dir: [`string`] The file path to the final tables
    @return: [`None`]
    """
    # Load all experiments from PSets
    experiment_df = load_join_write("experiment", data_dir, output_dir, [
                                    "cell", "compound", "dataset", "tissue"], join_dfs)
    # Don"t write the "name" column
    experiment_df[:, ["id", "cell_id", "compound_id", "dataset_id", "tissue_id"]] \
        .to_csv(os.path.join(output_dir, "experiment.jay"))
    # Rename columns and key experiment table based on experiment name and dataset id
    experiment_df.names = {"name": "experiment_id"}
    experiment_df = experiment_df[:, ["id", "experiment_id", "dataset_id"]]
    experiment_df.key = ("dataset_id", "experiment_id")
    join_dfs["experiment"] = experiment_df
    # Nearly the same code as in load_join_write but has special case handling
    for df_name in ["dose_response", "profile"]:
        df = load_table(df_name, data_dir)
        if df_name == "profile":
            df[dt.f.IC50 > 1e54, dt.update(IC50=1e54)]
        for fk in ["dataset", "experiment"]:
            df = join_tables(df, join_dfs[fk], fk+"_id")
        del df[:, "dataset_id"]
        write_table(df, df_name, output_dir,
                    add_index=(df_name == "dose_response"))


@logger.catch
def combine_gene_compound_tissue_dataset_tables(data_dir, output_dir, join_dfs):
    gd_df = load_table("gene_compound_tissue_dataset", data_dir)
    tissue_df = join_dfs["tissue"]

    # Join on gene, compound, dataset, and tissue IDs
    for fk in ["gene", "compound", "dataset", "tissue"]:
        logger.info(f"Joining gene_compound_tissue_dataset table with {fk} table...")
        gd_df = join_tables(gd_df, join_dfs[fk], f"{fk}_id")
    write_table(gd_df, "gene_compound_tissue_dataset", output_dir)


@logger.catch
def load_join_write(name, data_dir, output_dir, foreign_keys=[], join_dfs=None, add_index=True):
    """
    Given the name of a table, load all PSet tables of that name from data_dir,
    join them to any foreign key tables (specified by foreign_keys), and write
    the final combined and joined table to output_dir as a CSV.

    @param name: [`string`] The name of the table
    @param data_dir: [`string`] File path to the directory with all PSet tables
    @param output_dir: [`string`] The file path to the final tables
    @param foreign_keys: [`list(string)`] An optional list of tables that this table
        needs to be joined with
    @param join_dfs: [`dict(string: datatable.Frame)`] An optional dictionary of join
        tables (for building out foreign keys); keys are table names
    @param add_index: [`bool`] Indicates whether or not to add a primary key (1-nrows)
        when writing the final table to a .jay
    @return: [`datatable.Frame`] The final combined and joined table
    """
    df = load_table(name, data_dir)
    for fk in foreign_keys:
        logger.info(f"Joining {name} table with {fk} table...")
        if fk not in join_dfs:
            raise KeyError(f"The {name} table has the foreign key {fk}_id but \
                            there is no {fk} table in the join tables dictionary.")
        df = join_tables(df, join_dfs[fk], f"{fk}_id")
    fk_columns = [f"{fk}_id" for fk in foreign_keys]
    df = df[:, :, sort(fk_columns)]
    df = write_table(df, name, output_dir, add_index)
    return df


@logger.catch
def load_table(name, data_dir):
    """
    Load all PSet tables with name into a datatable, dropping any duplicate 
    rows.
    
    @param name: [`string`] The name of the table
    @param data_dir: [`string`] File path to the directory with all PSet tables
    @return: [`datatable.Frame`] A datatable containing all rows from all PSets
    """
    logger.info(f"Loading PSet-specific {name} tables from {data_dir}...")
    # Get all files
    files = glob.glob(os.path.join(data_dir, "**", f"*{name}.jay"))
    # Filter so that file path are "{data_dir}/{pset}/{pset}_{name}.jay"
    files = [file_name for file_name in files if re.search(
        data_dir + r"/(\w+)/\1_" + name + ".jay$", file_name)]
    # Read and concatenate tables
    df = dt.rbind(*dt.iread(files))
    # Drop duplicates (groups by all columns and
    # selects only the first row from each group)
    df = df[0, :, by(df.names)]
    return df


@logger.catch
def fread_table_for_all_psets(
    table_name: str,
    data_dir: str,
    column_dict: dict,
    rename_dict: dict = None,
    key_columns: list = None
) -> dt.Frame:
    """
    Reads all tables named `table_name` from `data_dir`, using `column_dict`
    to specify the column names, order and types to read in. The resulting
    table iterator is then concatenated using `datatable.rbind` and the
    columns are renamed according to rename dict.
    
    :param table_name:
    :param data_dir:
    :param column_dict:
    :param rename_dict:
    """
    logger.info(f"Loading PSet-specific {table_name} tables from {data_dir}...")
    # Get all files
    files = glob.glob(os.path.join(data_dir, "**", f"*{table_name}.jay"))
    # Filter so that file path are "{data_dir}/{pset}/{pset}_{name}.jay"
    files = [file_name for file_name in files if re.search(
        data_dir + r"/(\w+)/\1_" + table_name + ".jay$", file_name)]
    # Read and concatenate tables
    df = dt.rbind(*dt.iread(files, columns=column_dict), force=True)
    # Drop duplicates (groups by all columns and
    # selects only the first row from each group)
    df = df[0, :, by(df.names)]
    if rename_dict is not None:
        df.names = rename_dict
    if key_columns is not None:
        df = df[0, :, :, by(key_columns)]
    return df


@logger.catch
def rename_and_key(df, join_col, og_col="name"):
    """
    Prepare df to be joined with other tables by renaming the column
    on which it will be joined and by keying it.

    @param df: [`datatable.Frame`] The table to be keyed.
    @param join_col: [`string`] The name of the join column in other tables
                            (ex. "tissue_id", "cell_id", etc.)
    @param og_col: [`string`] The name of the join column in the join table
    @return: [`datatable.Frame`] The keyed and renamed table
    """
    # Rename primary key to match foreign key name (necessary for joins)
    df.names = {og_col: join_col}
    # Only select necessary columns
    df = df[:, ["id", join_col]]
    # Set the key
    df.key = join_col
    return df


@logger.catch
def join_tables(
    df1: dt.Frame, 
    df2: dt.Frame,
    join_col: str, 
    delete_unjoined=True
) -> None:
    """
    Join df2 and df1 based on join_col (left outer join by default).

    @param df1: [`datatable.Frame`] The datatable with the foreign key
    @param df2: [`datatable.Frame`] The join table (ex. tissue datatable)
    @param join_col: [`string`] The name of the columns on which the tables
        will be joined (ex. "tissue_id")
    @param delete_unjoined: [`bool`] An optional parameter (default=True)
        that lets you keep rows in df1 which didn"t join to any rows in df2
    @return [`datatable.Frame`] The new, joined table
    """
    if (join_col not in df1.names) or (join_col not in df2.names):
        logger.info(f"{join_col} is missing from one or both of the datatables " 
            "passed! Make sure you have prepared df2 using rename_and_key().")
        return None
    df = df1[:, :, dt.join(df2)]
    # Check to see if any FKs are null
    if df[dt.isna(df[:, "id"]), :].nrows > 0:
        logger.info(f"The following {join_col}s failed to map:")
        unmatched = df[dt.isna(df[:, "id"]), join_col].copy()
        unmatched = unmatched[0, :, dt.by(join_col)]
        logger.info(unmatched)
        if delete_unjoined:
            logger.info(f"Rows with these {join_col}s will be deleted!")
            del df[dt.isna(df[:, "id"]), :]
    # Rename the join col and drop it
    df.names = {join_col: "drop", "id": join_col}
    del df[:, "drop"]
    return df


@logger.catch
def write_table(df, name, output_dir, add_index=True):
    """
    Add a primary key to df ("id" column) and write it to output_dir
    as a .jay file.

    @param df: [`datatable.Frame`] A PharmacoDB table
    @param name: [`string`] The name of the table
    @param output_dir: [`string`] The directory to write the table to
    @return: [`datatable.Frame`] The indexed PharmacoDB table
    """
    logger.info(f"Writing {name} table to {output_dir}...")
    if add_index:
        # Index datatable
        df = dt.cbind(dt.Frame(id=np.arange(df.nrows) + 1), df)
    df.to_jay(os.path.join(output_dir, f"{name}.jay"))
    return df

