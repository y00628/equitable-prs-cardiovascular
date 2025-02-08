import os
import pandas as pd
import numpy as np
import subprocess
import urllib.request


def download_gwas(pops = ['eur', 'eas', 'amr', 'afr']):
    for pop in pops:
        url = f"https://gbmi-sumstats.s3.amazonaws.com/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
        output_file = f"../data/gwas/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"

        try:
            print("Downloading file...")
            urllib.request.urlretrieve(url, output_file)
            print("Download completed successfully.")
        except Exception as e:
            print(f"An error occurred: {e}")