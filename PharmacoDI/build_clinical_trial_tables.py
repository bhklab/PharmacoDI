import os
import requests
import pandas as pd
from datatable import dt, fread, f, g, join
from urllib3.exceptions import HTTPError
from .get_chembl_compound_targets import parallelize
from .combine_pset_tables import write_table

# -- Enable logging
from loguru import logger
import sys

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/build_clinical_trails_tables.log", 
            "serialize": True, # Write logs as JSONs
            "enqueue": True}, # Makes logging queue based and thread safe
    ]
}
logger.configure(**logger_config)

@logger.catch
def build_clinical_trial_tables(output_dir: str) -> pd.DataFrame:
    """
    Build the clinical trial and compound trial tables by querying the
    clinicaltrial.gov API. Queries are made by compound names from the compound
    synonyms table.

    @param output_dir: [`string`] The file path to the directory with all
                                    PharmacoDB tables
    @return: None
    """
    # Load compound synonym table
    compound_file = os.path.join(output_dir, 'compound_synonym.jay')
    compound_df = fread(compound_file).to_pandas()[['compound_id', 'compound_name']]

    # Query clinicaltrials.gov API to get clinical trials by compound name
    logger.info('Getting clinical trials from clinicaltrials.gov...')
    all_studies = parallelize(list(compound_df['compound_name']),
        get_clinical_trials_by_compound_names, 50)
    studies_df = pd.concat(all_studies)

    # Explode list-like columns into separate rows, duplicating the index
    # I only use this because all the fields are returned in arrays for some reason
    object_columns = studies_df.dtypes[
        studies_df.dtypes == 'object'
    ].index.values
    for column in object_columns:
        studies_df = studies_df.explode(column)
    # Drop and rename columns
    studies_df.drop(columns='Rank', inplace=True)
    studies_df.rename(
        columns={'NCTId': 'nct', 'SeeAlsoLinkURL': 'link', 
            'OverallStatus': 'status'}, 
        inplace=True
    )

    # Build clinical trials table
    clin_trial_df = studies_df[['nct', 'link', 'status']].copy()
    clin_trial_df.drop_duplicates('nct', inplace=True)
    clin_trial_df.reset_index(inplace=True, drop=True)
    clin_trial_df['clinical_trial_id'] = clin_trial_df.index + 1

    # Build compound trial table
    compound_trial_df = studies_df[['nct', 'compound_name']].copy()
    compound_trial_df.drop_duplicates(inplace=True)    
    compound_trial_df = pd.merge(compound_trial_df, clin_trial_df, on='nct')
    compound_trial_df = pd.merge(compound_trial_df, compound_df, on='compound_name')

    # Write both tables
    write_table(dt.Frame(clin_trial_df), 'clinical_trial', output_dir, add_index=False)
    write_table(dt.Frame(compound_trial_df[['clinical_trial_id', 'compound_id']]), 'compound_trial', output_dir, add_index=False)


@logger.catch
def get_clinical_trials_by_compound_names(compound_names):
    """
    Given a list of compound_names, query the clinicaltrial.gov API iteratively
    to get all trials related to these compounds and return these studies in a table.

    @param compound_names: [`list(string)`] A list of (up to 50) compound names
    @return: [`pd.DataFrame`] A table of all studies, including their rank, study ID,
        NCT id, recruitment status, link, and compound name.
    """
    all_studies = []
    for compound_name in compound_names:
        clean_compound_name = compound_name.replace(' ', '')
        min_rank = 1
        max_rank = 1000
        # Make initial API call for this compound
        failed_compounds = []
        try:
            studies, num_studies_returned, num_studies_found = get_clinical_trials_for_compound(
                clean_compound_name, min_rank, max_rank
            )
        except:
            failed_compounds = [*failed_compounds, compound_name]
            next
        # If not all studies were returned, make additional calls
        while num_studies_found > num_studies_returned:
            min_rank += 1000
            max_rank += 1000
            try:
                more_studies, n_returned, n_found = get_clinical_trials_for_compound(
                    clean_compound_name, min_rank, max_rank
                )
            except:
                next
            studies = pd.concat([studies, more_studies])
            num_studies_returned += n_returned
        studies['compound_name'] = f"{compound_name}"
        all_studies.append(studies)
    if (len(failed_compounds) > 0):
        failed_names = ', '.join(failed_compounds)
        logger.warning(f"Failed compounds: {failed_names}")
    return pd.concat(all_studies)


@logger.catch
def get_clinical_trials_for_compound(
    compound_name: str, 
    min_rank: int, 
    max_rank: int
) -> tuple:
    """
    Given a compound_name, query the clinicaltrial.gov API to get all trials
    for this compound between min_rank and max_rank (inclusive). Return the 
    studies in a table. If the HTTP request fails or if no studies are
    returned, returns an empty DataFrame.

    @param compound_name: [`string`] A compound name
    @param min_rank: [`int`] The minimum rank of a retrieved study
    @param max_rank: [`int`] The maximum rank of a retrieved study
    @return: [`pd.DataFrame`] A table of retrieved studies, including
        their rank, study ID, NCT id, recruitment status, and link.
    """
    # Make API call
    base_url = 'https://clinicaltrials.gov/api/query/study_fields'
    params = {
        'expr': requests.utils.quote(compound_name),
        'fields': 'OrgStudyId,NCTId,OverallStatus,SeeAlsoLinkURL',
        'min_rnk': min_rank,
        'max_rnk': max_rank,
        'fmt': 'json'
    }
    r = requests.get(base_url, params=params)
    studies = pd.DataFrame(columns=['Rank', 'NCTId', 'OverallStatus', 
        'SeeAlsoLinkURL'])
    # Check that request was successful
    if r.status_code != 200 and r.status_code != 404:
        logger.warning(f"API call for clinical trials related to {compound_name} failed!")
        return studies, 0, 0
    elif r.status_code == 404:
        logger.warning(f"Returned HTML 404 page for {compound_name}")
        return studies, 0, 0
    response = r.json()['StudyFieldsResponse']
    if 'StudyFields' in response:
        studies = pd.DataFrame(response['StudyFields'])
    return studies, response['NStudiesFound'], response['NStudiesReturned']
