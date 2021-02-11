import pandas as pd
from bs4 import BeautifulSoup
import time
from selenium import webdriver
import requests

ensembl_gene_id = 'ENSG00000140465'


def get_gene_targets(ensembl_gene_id):
    """
    Query the APIs of Ensembl and NCBI to fetch information about drugs targeting a given gene.

    :param gene_symbol: [string] The gene symbol for the desired target gene
    :param ensembl_gene_id: [string] The Ensembl gene id for the desired target gene
    :return: [DataFrame] containing annotations for the specified gene
    """
    # Configure requests to automatically raise errors
    http = requests.Session()

    assert_status_hook = lambda response, *args, **kwargs: response.raise_for_status()
    http.hooks["response"] = [assert_status_hook]

    # Query the APIs
    queries = {'genecards_query':
                   f"https://www.genecards.org/cgi-bin/carddisp.pl?id={ensembl_gene_id.upper()}&idtype=ensembl",
               'ncbi_query': f"https://www.ncbi.nlm.nih.gov/gene/?term={ensembl_gene_id.upper()}"},
    f"http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_gene_id.upper()}"

    api_requests = {key: requests.get(query) for key, query in queries.items()}

    # Retry if 403
    is_403 = [req.status_code == "403" for req in api_requests.values()]
    if (any(is_403)):
        to_update = queries.keys()[is_403]
        updates = [requests.get(query) for query in queries.values()[is_403]]
        api_requests.replace(dict(zip(to_update, updates)))  # Map replacement by value to dict key

    parsed_annotation_data = \
        {key: BeautifulSoup(request.text, 'html.parser') if request.status_code == "200" else None
         for key, request in api_requests.items()}


def scrape_genecards(ensembl_gene_id):
    """
    Use headless Firefox browser to make get requests to Genecards despite their anti-scraping software

    WARNING: Requires Firefox and geckodriver be installed to work!

    :param ensembl_gene_id: [string] The ENSEMBL id for the gene you want to query
    :return: [DataFrame] containing the gene card data
    """
    # Configure to run Chrome/Chromium headless
    options = webdriver.FirefoxOptions()
    #options.add_argument("-disable-extensions")
    #options.add_argument("-disable-gpu")
    #options.add_argument("-no-sandbox")
    #options.add_argument("-headless")

    # Initialize Chrome/Chromium
    driver = webdriver.Firefox(options=options)

    # Build the HTTP query and use the browser to get it
    ensembl_query = f"https://www.genecards.org/cgi-bin/carddisp.pl?id={ensembl_gene_id.upper()}&idtype=ensembl"
    driver.get(ensembl_query)
    expand_table_elements = driver.find_elements_by_xpath("//a[@data-role='show-all']")
    for element in expand_table_elements:
        time.sleep(1)
        element.click()
    page_html = driver.page_source
    del driver

    # Parse the page HTML to a DataFrame
    parsed_html = BeautifulSoup(page_html, 'html.parser')
    genecards_dfs = pd.read_html(str(parsed_html.find_all('table')))

