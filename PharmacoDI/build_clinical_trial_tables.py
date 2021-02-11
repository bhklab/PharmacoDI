import os
import requests
import pandas as pd
import datatable as dt
from urllib3.exceptions import HTTPError
from PharmacoDI.get_chembl_drug_targets import parallelize
from PharmacoDI.combine_pset_tables import join_tables, write_table


# TODO: split into more helpers?
def build_clinical_trial_tables(output_dir):
    """
    Build the clinical trial and drug trial tables by querying the
    clinicaltrial.gov API. Queries are made by drug names from the drug
    synonyms table.

    @param output_dir: [`string`] The file path to the directory with all
                                    PharmacoDB tables
    @return: None
    """
    # Load drug synonym table
    drug_file = os.path.join(output_dir, 'drug_synonym.csv')
    drug_df = pd.read_csv(drug_file)[['drug_id', 'drug_name']]

    # Query clinicaltrials.gov API to get clinical trials by drug name
    print('Getting clinical trials from clinicaltrials.gov...')
    all_studies = parallelize(list(drug_df['drug_name']),
                              get_clinical_trials_by_drug_names, 50)
    studies_df = pd.concat(all_studies)

    # Explode list-like columns into separate rows, duplicating the index
    # I only use this because all the fields are returned in arrays for some reason
    object_columns = studies_df.dtypes[studies_df.dtypes ==
                                       'object'].index.values
    for column in object_columns:
        studies_df = studies_df.explode(column)
    # Drop and rename columns
    studies_df.drop(columns='Rank', inplace=True)
    studies_df.rename(columns={'OrgStudyId': 'clinical_trial_id',
                               'NCTId': 'nct',
                               'SeeAlsoLinkURL': 'link',
                               'OverallStatus': 'status'}, inplace=True)

    # Build clinical trials table
    clin_trial_df = studies_df[['clinical_trial_id',
                                'nct', 'link', 'status']].copy()
    clin_trial_df.drop_duplicates('clinical_trial_id', inplace=True)
    clin_trial_df.reset_index(inplace=True, drop=True)
    write_table(dt.Frame(clin_trial_df), 'clinical_trial',
                output_dir, add_index=False)

    # Build drug trial table
    drug_trial_df = studies_df[['clinical_trial_id', 'drug_name']].copy()
    drug_trial_df.drop_duplicates(inplace=True)
    drug_trial_df = pd.merge(drug_trial_df, drug_df, on='drug_name')
    drug_trial_df.drop(columns='drug_name', inplace=True)
    write_table(dt.Frame(drug_trial_df), 'drug_trial',
                output_dir, add_index=False)


# TODO: shorter names please?
def get_clinical_trials_by_drug_names(drug_names):
    """
    Given a list of drug_names, query the clinicaltrial.gov API iteratively
    to get all trials related to these drugs and return these studies in a table.

    @param drug_names: [`list(string)`] A list of (up to 50) drug names
    @return: [`pd.DataFrame`] A table of all studies, including their rank, study ID,
                                NCT id, recruitment status, link, and drug name.
    """
    all_studies = []
    for drug_name in drug_names:
        min_rank = 1
        max_rank = 1000
        # Make initial API call for this drug
        studies, num_studies_returned, num_studies_found = get_clinical_trials_for_drug(
            drug_name, min_rank, max_rank)

        # If not all studies were returned, make additional calls
        while num_studies_found > num_studies_returned:
            min_rank += 1000
            max_rank += 1000
            more_studies, n_returned, n_found = get_clinical_trials_for_drug(
                drug_name, min_rank, max_rank)
            studies = pd.concat([studies, more_studies])
            num_studies_returned += n_returned
        studies['drug_name'] = drug_name
        all_studies.append(studies)
    return pd.concat(all_studies)


def get_clinical_trials_for_drug(drug_name, min_rank, max_rank):
    """
    Given a drug_name, query the clinicaltrial.gov API to get all trials
    for this drug between min_rank and max_rank (inclusive). Return the 
    studies in a table. If the HTTP request fails or if no studies are
    returned, returns an empty DataFrame.

    @param drug_name: [`string`] A drug name
    @param min_rank: [`int`] The minimum rank of a retrieved study
    @param max_rank: [`int`] The maximum rank of a retrieved study
    @return: [`pd.DataFrame`] A table of retrieved studies, including 
        their rank, study ID, NCT id, recruitment status, and link.
    """
    # Make API call
    base_url = 'https://clinicaltrials.gov/api/query/study_fields'
    params = {
        'expr': drug_name,
        'fields': 'OrgStudyId,NCTId,OverallStatus,SeeAlsoLinkURL',
        'min_rnk': min_rank,
        'max_rnk': max_rank,
        'fmt': 'json'
    }
    r = requests.get(base_url, params=params)

    studies = pd.DataFrame(columns=['Rank', 'OrgStudyId', 'NCTId',
                                    'OverallStatus', 'SeeAlsoLinkURL'])

    # Check that request was successful
    if r.status_code != 200:
        print(f'API call for clinical trials related to {drug_name}',
              ' failed for the following reason:\n', r.json()['Error'])
        return studies, 0, 0

    response = r.json()['StudyFieldsResponse']
    if 'StudyFields' in response:
        studies = pd.DataFrame(response['StudyFields'])

    return studies, response['NStudiesFound'], response['NStudiesReturned']
