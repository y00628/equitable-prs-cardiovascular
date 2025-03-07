import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import data
import compute_weights
import association_test
import sys

if __name__ == '__main__':
    args = sys.argv[1:]

    print('\nPress Ctrl-C to terminate the program at any point.\n')

    if not args:
        print('Use the `help` argument to display a list of possible arguments.')

    if 'help' in args:
        print('`help`: List possible arguments')
        print('`data`: Download and process data for model weight computation \
and cross population analyses')
        print('`compute_weights`: Compute model weights for comparison to \
              \nGWAS summary statistics')
        print('`association_test`: Run association tests between calculated SNP weights for European samples \
\n against GWAS sumstats for European, East Asian, African, and American populations.')
        

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
            print('Gene expression splitting complete.\n')
        else:
            print('Skipping gene expression splitting.\n')


        if input("Filter and prune the LDREF samples and genotypes? \n If this \
is the first time this is being run, enter 'y'. \n([y/n]): ") == 'y':            
            print('Beginning filtering, and pruning and thresholding for LDREF samples and genotypes')
            data.filter_and_prune_genotype()
        else:
            print('Skipping filtering and pruning.')

        print('Data processing complete!')

        
    if 'compute_weights' in args:
        print('\nWARNING: The weights have been precomputed for your convenience. \
\nRunning the code to recompute the weights will take incredibly long to run and may take \
\nanywhere from several hours to several days to finish running, depending\
\non your machine.')
        if input("Run weight computation anyway (recommended not to)? \n([y/n]): ") == 'y':            
            print('Beginning filtering, and pruning and thresholding for LDREF samples and genotypes')            
            BFILE_TEMPLATE=os.path.join("..","data","LDREF_pruned","1000G.EUR.{chr}")
            compute_weights.compute_weights(
            chromosomes=list(str(input("Which chromosome would you like to compute weights for (1-22): "))), 
                BFILE_TEMPLATE=BFILE_TEMPLATE,
                VERBOSE=0
            )
            print("Done!")
        else:
            print('Skipping weight computations.')

            
    if 'association_test' in args:
        if input("Gene position files have already been created for FUSION.assoc_test.R for your convenience. \
\nRecreate gene position files anyway? \n This might take a few minutes. \n([y/n]): ") == 'y':
            print("Creating gene position files for FUSION.assoc_test.R...")
            association_test.create_pos_files()
            print('Done!')
        else:
            print('Skipping gene position file creation.')

        print("Generating list of SNPs FUSION.assoc_test.R...")
        association_test.generate_snp_list()
        print('Done!')

        if input("Process GWAS summary statistics for FUSION.assoc_test.R? \n If this \
is the first time this is being run, enter 'y'. \n([y/n]): ") == 'y':
            print('Processing GWAS files for FUSION.assoc_test.R...')
            association_test.process_sumstats()
            print('Done!')

        if input("Run association tests? \n([y/n]): ") == 'y':  
            association_test.association_test()
        else:
            print('Skipping association tests')          
            

