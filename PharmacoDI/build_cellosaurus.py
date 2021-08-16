import pandas as pd
import os
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from datatable import dt, fread, f, g, join
from PharmacoDI.combine_pset_tables import join_tables, write_table, rename_and_key

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_cellosaurus.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)

# Using a default dict because it allows me to append duplicate indexes into a list
# Helper for build_cellosaurus_df
@logger.catch
def build_defaultdict(tuple_list):
    def_dict = defaultdict(list)
    for tup in tuple_list:
        def_dict[tup[0]].append(tup[1])
    return def_dict


@logger.catch
def build_cellosaurus_df(cellosaurus_path, output_dir):
    """
    Build cellosaurus table.

    @param cellosaurus_path: [`string`] Full file path to cellosaurus file
    @param output_dir: [`string`] The directory to write the cellosaurus table
    @param cell_df: [`datatable.Frame`] The cell table; should be renamed, keyed,
                                        and shouldn't have 'tissue_id' column
    @return: [`datatable.Frame`] The cellosaurus table
    """

    with open(cellosaurus_path) as f:
        file = [line for line in f]

    file = file[55:]
    entries = ''.join(file).split('//\n')
    entry_list = [entry.split('\n') for entry in entries]
    entry_split_list = [[item.split('   ')
                         for item in entry] for entry in entry_list]
    entry_tuple_list = [[(item[0], item[1]) for item in entry if len(
        item) > 1] for entry in entry_split_list]

    pool = Pool(cpu_count() - 1)

    dict_list = pool.map(build_defaultdict, entry_tuple_list)
    dict_list = [dict(item) for item in dict_list]
    dict_list = [{key: '|||'.join(value)
                  for key, value in dct.items()} for dct in dict_list]

    cellosaurus_df = pd.DataFrame(dict_list)
    cellosaurus_df.dropna(axis=1, how='all', inplace=True)

    # Always close your pool or you will have a bunch of processes doing nothing
    pool.close()

    # Drop AG and DT columns (age of donor, date)
    cellosaurus_df.drop(columns=['AG', 'DT'], inplace=True)

    # Rename cols and add cell_id column
    rename_dict = {col: col.lower() for col in cellosaurus_df.columns}
    cellosaurus_df.rename(columns=rename_dict, inplace=True)
    cellosaurus_df.rename(
        columns={'id': 'identifier', 'ac': 'accession'}, inplace=True)
    cellosaurus_df['cell_id'] = cellosaurus_df['identifier']

    # Load cell_df
    cell_path = os.path.join(output_dir, 'cell.jay')
    cell_df = rename_and_key(dt.fread(cell_path), 'cell_id')

    # Convert to datatable and join with cell_df
    cellosaurus_df = dt.Frame(cellosaurus_df)
    cellosaurus_df.key = 'cell_id'
    df = cell_df[:, :, join(dt.Frame(cellosaurus_df))]
    df = df[dt.f.id >= 1, :]
    df = df[:, ['cell_id', 'id', 'accession', 'as', 'sy',
                'dr', 'rx', 'ww', 'cc', 'st', 'di', 'ox', 'hi', 'oi', 'sx', 'ca']]
    df.names = {'cell_id': 'identifier', 'id': 'cell_id'}
    df = write_table(df, 'cellosaurus', output_dir)
    return df
