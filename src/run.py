import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import data
import compute_weights_win
import extract_data
import sys

if __name__ == '__main__':
    args = sys.argv[1:]

    print('\nPress Ctrl-C to cancel at any point.\n')

    if 'data' in args:        
        if input("Download GWAS summary statistics? This may take a while. \n\
Check that the 5 GWAS sumstats files are not already in \
data/gwas/. \n([y/n]): ") == 'y':
            print('Starting GBMI GWAS downloads\n')
            data.download_gwas()
            print('All GWAS files downloaded\n')
        else:
            print('Skipping GWAS download\n')

        
        if input("Split the gene expression files? \nIf this is the \
first time this is being run, enter 'y'. \n([y/n]): ") == 'y':
            print('Splitting gene expression files for FUSION.compute_weights.R')
            data.phenotype_split()
            print('Gene expression splitting complete.')
        else:
            print('Skipping gene expression splitting.')


        if input("Filter and prune the LDREF samples and genotypes? \n If this \
is the first time this is being run, enter 'y'. \n([y/n]): ") == 'y':            
            print('Beginning filtering, and pruning and thresholding for LDREF samples and genotypes')
            data.filter_and_prune_genotype()

        
    if 'compute_weights' in args:
        print('\nWARNING: This code takes incredibly long to run and may take \
\nanywhere from several hours to several days to finish running, depending\
\non your machine and the chromosome chosen. The weights for chromosome 19\
\nhave already been computed if you would like to test further functions.\n')
        ldref_yn = input('Would you like to use the pruned LDREF instead of just the filtered LDREF? ([y/n]): ')
        BFILE_TEMPLATE="../data/LDREF_filtered/1000G.EUR.{chr}"
        if ldref_yn == 'y':
            BFILE_TEMPLATE="../data/LDREF_pruned/1000G.EUR.{chr}"

         
        compute_weights_win.compute_weights_win(
            chromosomes=list(str(input("Which chromosome would you like to compute weights for (1-22): "))), 
            BFILE_TEMPLATE=BFILE_TEMPLATE,
            VERBOSE=0
            )