import os
from datatable import Frame

def write_pset_table(pset_df, df_name, pset_name, df_dir):
    """
    Write a PSet table to a CSV file.

    @param pset_df: [`DataFrame`] A PSet DataFrame
    @param pset_name: [`string`] The name of the PSet
    @param df_dir: [`string`] The name of the directory to hold all the PSet tables
    @return [`None`]
    """
    pset_path = os.path.join(df_dir, pset_name)
    # Make sure directory for this PSet exists
    if not os.path.exists(pset_path):
        os.mkdir(pset_path)

    # Convert to datatable Frame for fast write to disk
    pset_df = Frame(pset_df)

    print(f'Writing {df_name} table to {pset_path}...')
    # Use datatable to convert df to csv
    pset_df.to_csv(os.path.join(pset_path, f'{pset_name}_{df_name}.csv'))
    