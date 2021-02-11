import sys
import os

def download_gene_signatures(user=None, password=None, remote="niagara.computecanada.ca", opts={},
                             remote_path="/scratch/b/bhaibeka/psmirnov/pearson_perm_res/", save_dir="data/rawdata/gene_signatures"):
    """
    Download all precomputed gene signatures from the `remote_path` directory on the `remote` server, excluding those
    already in `save_dir`.

    WARNING: The current implementation relies on rysnc being installed. Thus it will not work on Windows platforms.

    :param user [string] Your username for the remote server. In interactive sessions, if you exclude this argument
        you will be prompted to enter your username.
    :param password [string] Your password for the remote server. To avoid hard coding your password, we recommend
        you read this from an environmental variables using `os.environment`. In interactive sessions, if you exclude
        this argument you will be prompted to enter you password.
    :param remote: [string] Name or IP of remote server
    :param remote_path: [string] Path to the gene signature directory on the remote server
    :param save_dir: [string] Path to save the downloaded gene signature to.
    :return: [None] Syncs `save_dir` with `remote:remote_path`, downloading the most up to date gene signatures for
        PharmacoDB.
    """

    if user == None or password == None:
        if sys.flags.interactive:
            if user == None:
                user = input("Please enter your username for {}: ".format(remote))
            if password == None:
                password = input("Please enter your password for {}: ".format(remote))

        else:
            raise ValueError("You must pass `user` and `password` parameters to this function in non-interative Python sessions")

    os.system(f'rsync -hrtvP {user}@{remote}:{remote_path} {save_dir}')

