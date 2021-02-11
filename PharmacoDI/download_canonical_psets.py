import pandas as pd
import wget

def download_canonical_psets(save_dir, api_url= "https://www.orcestra.ca/api/psets/canonical"):
    """
    Download the specified canonical PSets from the specified api_url

    :param save_dir: [string] Path to save the PSets in
    :param api_url: [string] URL where available PSets can be retrieved. Defaults to current Orcestra API.
    :return: [None] Downloads PSet .rds files into save_dir using wget.
    """
    pset_df = pd.read_json(api_url)
    url_dict = pset_df.set_index('name').to_dict()['downloadLink']

    for name, url in url_dict.items():
        print("Downloading", name, "from", url, sep=" ")
        wget.download(url, save_dir + "/" + name + '.rds')

    return None
