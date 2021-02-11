import pandas as pd
import wget


def download_psets(names, save_dir, api_url="https://www.orcestra.ca/api/psets/available"):
    """
    Download the specified PSets from the specified api_url

    :param names: [list] Names of PSets to download. Must match the 'names' in the api call JSON.
    :param save_dir: [string] Path to save the PSets in
    :param api_url: [string] URL where available PSets can be retrieved. Defaults to current Orcestra API.
    :return: [None] Downloads PSet .rds files into save_dir using wget.
    """
    pset_df = pd.read_json(api_url)
    names = pd.Series(names)
    if not all(names.isin(pset_df.name)):
        print(names[~names.isin(pset_df.name)] +
                         ' not in names of retrived from api_url')
        raise ValueError(names[~names.isin(pset_df.name)] + 'are not valid pset names')
    pset_df = pset_df[pset_df.name.isin(names)]
    url_dict = pset_df.set_index('name').to_dict()['downloadLink']
    for name, url in url_dict.items():
        print("Downloading", name, "from", url, sep=" ")
        wget.download(url, save_dir + "/" + name + '.rds')

    return None
