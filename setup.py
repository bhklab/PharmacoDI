from setuptools import setup, find_packages

setup(name='PharmacoDI',
      version='0.0.1',
      description="Tools for processing R PharmacoSet objects into .csv files of PharmacoDB database tables.",
      url='https://github.com/bhklab/DataIngestion/tree/master/PharmacoDI',
      install_requires=[
        'dask[dataframe]',
        'swifter',
        'datatable',
        'pandas',
        'chembl_webresource_client',
        'wget',
        'bs4',
        'selenium',
        'lxml'
      ],
      author='Evgeniya Gorobets, Christopher Eeles, Benjamin Haibe-Kains',
      author_email='christopher.eeles@uhnresearch.ca, benjamin.haibe.kains@utoronto.ca',
      license='MIT',
      packages=find_packages(),
      zip_safe=False
      )