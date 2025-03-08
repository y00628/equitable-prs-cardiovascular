import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import data
import compute_weights
import association_test
import generate_plots
import sys

if __name__ == '__main__':
    args = sys.argv[1:]

    print('\nPress Ctrl-C or Cmd-C to terminate the program at any point.\n')

    commands = ['help', 'data', 'compute_weights', 'association_test', 'generate_plots']
    unrecognized_commands = list(filter(lambda cmd: cmd not in commands, args))

    if 'help' in args or not args:
        print('`help`: List possible arguments.\n')
        print('`data`: Download and process data for model weight computation \
and cross population analyses. \nMust be run before `compute_weights`, `association_test`, and `generate_plots`.\n')
        print('`compute_weights`: Compute model weights for comparison to GWAS summary statistics. \
              \nMust be run before `association_test` and `generate_plots`.\n')
        print('`association_test`: Run association tests between calculated SNP weights for European samples \
\nagainst GWAS sumstats for European, East Asian, African, and American populations. \nMust be run before `generate_plots`.\n')
        print('`generate_plots`: Generate plots from `association_test` results.\n')
        

    if 'data' in args:        
        if input("Download GWAS summary statistics? This may take a while. \
\nEnter 'y' if this is the first time it is being run. \n([y/n]): ") == 'y':
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

        if input("GWAS summary statistics have already by processed for your convenience. \
\nReprocess GWAS summary statistics for FUSION.assoc_test.R? This will take a long time. \n([y/n]): ") == 'y':
            print('Processing GWAS files for FUSION.assoc_test.R...')
            association_test.process_sumstats()
            print('Done!')

        if input("Association test results have been saved for your convenience. \
\nRun association tests again anyway? \n([y/n]): ") == 'y':  
            association_test.association_test()
        else:
            print('Skipping association tests')     

    if 'generate_plots' in args:
        if input("Generate gene loci plot? \n([y/n]): ") == 'y':  
            print('Generating gene loci plot...')
            generate_plots.gene_loci()
            print('Done!')
        if input("Generate gene venn diagram? \n([y/n]): ") == 'y':  
            print('Generating gene venn diagram')
            generate_plots.generate_venn()
            print('Done!')
        if input("Generate Manhattan plots? This might take a couple minutes. \n([y/n]): ") == 'y':  
            print('Creating Manhattan plots...')
            generate_plots.manhattan()
            print('Done!')
        if input("Generate Miami plots? This might take a couple minutes. \n([y/n]): ") == 'y':  
            print('Creating Miami plots...')
            generate_plots.miami()  
            print('Done!')
        if input("Create GO plots? \n([y/n]): ") == 'y':  
            print('Creating GO plots.')
            generate_plots.go_plots()
            print('Done!')
            
    if unrecognized_commands:
        print(f"Unrecognized commands: {unrecognized_commands}\nUse the `help` argument to display a list of commands.")