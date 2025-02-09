import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import data
import compute_weights_win
import extract_data
import test
import sys

if __name__ == '__main__':
    args = sys.argv[1:]

    if 'test' in args:
        test.test()

    if 'data' in args:
        if input("Re-download GWAS summary statistics? This may take a while. ([y/n]): ") == 'y':
            print('Starting GBMI GWAS downloads')
            data.download_gwas()
            print('All GWAS files downloaded')
        else:
            print('Skipping GWAS download')
        print('Splitting gene expression files for FUSION.compute_weights.R')
        data.phenotype_split()
        print('Gene expression splitting complete.')
        print('Beginning filtering, and pruning and thresholding for LDREF samples and genotypes')
        data.filter_and_prune_genotype()
        
        
