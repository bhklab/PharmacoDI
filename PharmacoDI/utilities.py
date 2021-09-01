import pandas as pd
import numpy as np
from datatable import dt, f, g, join, update


def harmonize_df_columns(
    df: pd.DataFrame,
    column_dict: dict
) -> pd.DataFrame:
    """
    Harmonizes df with the table definition in column dict. 
    
    Specifically, it takes the set difference of the columns in 
    `df` and the `column_dict` keys, then subsets to `df` to 
    columns in `column_dict` or adds them if they are missing, 
    then checks for the correct column types.

    @param df: A pd.DataFrame to pad.
    @param column_dict: A dictionary where required column names are keys
        and respective values are the column numpy type.

    @return df with all the columns specified in column_dict.
    """
    # Check required columns
    columns = np.asarray(list(column_dict.keys()))
    existing_columns = np.asarray(df.columns)
    has_columns = np.isin(columns, existing_columns)
    # Select the columns of interest
    df = df[list(columns[has_columns])].copy()
    # Pad with missing columns
    empty_column = [None for _ in range(df.shape[0])]
    missing_column_dict = {key: value for key, value in column_dict.items() 
        if key in list(columns[~has_columns])}
    for column, dtype in missing_column_dict.items():
        df[column] = pd.Series(empty_column, dtype=dtype)
    # Check column types
    for column, dtype in column_dict.items():
        if not isinstance(df[column][0], dtype):
            try:
                if isinstance(df[column][0], float) and dtype == str:
                    # deal with lack of NA in Pandas int Series
                    # also ensure floats converted to strings don't have decimals
                    df[column] = df[column].astype('Int64').astype(str) \
                        .replace({'<NA>': None})
                else:
                    df[column] = df[column].astype(dtype)
            except TypeError:
                print(f'DataFrame column {column} failed to be coerced to'
                    '{str(dtype)}!')
    return df



def map_foreign_key_to_table(
    primary_df: dt.Frame,
    fk_df: dt.Frame,
    join_column_dict: dict
) -> dt.Frame:
    """
    Performs a left join of `primary_df` to `fk_df` by refence, updating
    the column indicated in `join_column_dict`.

    :primary_df: A `datatable.Frame`. This should be the larger table
        and will ideally be loaded from a .jay file with a `memory_limit`
        specified in `datable.fread`.
    :fk_df: A `datatable.Frame`. This should be a smaller table
        which will be joined to 
    :join_column_dict: A dictionary with keys 'primary_df' and 'fk_df'
        specifying the columns to join the tables on.
    """
    # Check for correct keys in dict
    key_strings = list(join_column_dict.keys())
    if ('primary_df' not in key_strings or 'fk_df' not in key_strings):
        raise ValueError("The join_column_dict item must have keys"
            "'primary_df' and 'fk_df'!")
    # Rename columns
    primary_col = join_column_dict['primary_df']
    fk_col = join_column_dict['fk_df']
    fk_df.names = {fk_col: primary_col}
    fk_df.key = primary_col
    update_expr = {primary_col: g.id}
    # Join, update by reference then coerce to the correct type
    primary_df[:, update(**update_expr), join(fk_df)]
