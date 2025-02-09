import os
import pandas as pd
import numpy as np
import subprocess
import urllib.request

def standardize_phenotype(out='../data/std_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz'):
    '''
    Standardizes and filters gene expressions to have only samples 
    present in the LDREF genotype data.
    
    '''
    print('Loading gene expression data for standardization...')
    all_sample_fid = pd.read_csv('../data/raw/LDREF/1000G.EUR.1.fam', header=None, sep=' ')[0].unique()
    all_genes = pd.read_csv('../data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', sep='\t')
    all_genes = all_genes[['TargetID', 'Chr', 'Coord']+[c for c in all_genes.columns[4:] if c in all_sample_fid]]
    to_std = all_genes[all_genes.columns[4:]].T
    print('Standardizing gene expressions...')
    all_genes[all_genes.columns[4:]] = ((to_std - to_std.mean()) / to_std.std()).T

    all_genes.to_csv(out, sep='\t')


def download_gwas(pops = ['eur', 'eas', 'amr', 'afr', '']):
    '''
    Downloads GBMI GWAS summary statistics.

    Parameters
    pops (list): A list containing the desired population GWAS sumstats to download.
        Possible populations: 
            European: 'eur', 
            East Asian: 'eas', 
            Americas: 'amr', 
            African: 'afr', 
            All: ''
    '''
    count = 1
    for pop in pops:
        if pop:
            url = f"https://gbmi-sumstats.s3.amazonaws.com/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
            output_file = f"../data/gwas/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
        else:
            url = "https://gbmi-sumstats.s3.amazonaws.com/HF_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
            output_file = f"../data/gwas/HF_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"

        try:
            print(f"({count}/{len(pops)}) Downloading {url}")
            urllib.request.urlretrieve(url, output_file)
            print("Download completed successfully.")
        except Exception as e:
            print(f"An error occurred: {e}")


def phenotype_split():
    '''
    Splits the standardized genotypes file into 1 file per gene so that 
    it can be used in FUSION.compute_weights.R

    If the standardized weights do not exist, it will create them.
    '''
    print('Loading gene expression data for splitting...')
    all_sample_fid = pd.read_csv(\
            '../data/LDREF/1000G.EUR.1.fam', header=None, sep=' '
        )[0].unique()
    try:
        all_genes = pd.read_csv(\
            '../data/std_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', \
                sep='\t')
    except Exception as e:
        print(f"An error occurred: {e}")
        if input("Create standardized gene expression file? ([y/n]): ") == 'y':
            standardize_phenotype()
            phenotype_split()
        else:
            return


    all_genes = all_genes[['TargetID', 'Chr', 'Coord'] + \
                          [c for c in all_genes.columns[4:] \
                           if c in all_sample_fid]
                            ]

    for i in all_genes.index:
        temp_gene = all_genes.iloc[i]
        gene_name = temp_gene['TargetID']
        temp_chrom = temp_gene['Chr']
        samples = temp_gene[4:]
        temp_df = pd.DataFrame(samples).reset_index().rename(\
            columns={'index': 'FID', i: 'Expression'})
        temp_df['IID'] = temp_df['FID']
        temp_df = temp_df[['FID', 'IID', 'Expression']]
        os.makedirs(f'../data/gene_expressions/chr_{temp_chrom}/', exist_ok=True)
        temp_df.to_csv(\
            f'../data/gene_expressions/chr_{temp_chrom}/{gene_name}.txt', \
                index=False, header=None, sep='\t')  
    print('Splitting complete')


def filter_and_prune_genotype():
    '''
    Filters the genotype data to keep only samples that exist in both the 
    genotype data and the phenotype (gene expression) data. This data is saved
    to its own directory (data/LDREF_filtered).

    Prunes and thresholds the filtered data to reduce the amount of SNPs in the
    genotype data. This data is saved to its own directory (data/LDREF_pruned).


    '''

    # Path to your input files and desired output folder
    original_file = "../data/raw/LDREF/1000G.EUR."
    keep_file = "../data/samples.txt"
    output_dir_filter = "../data/LDREF_filtered"  # Directory where filtered files will be stored
    output_dir_pruned = "../data/LDREF_pruned"  # Directory where filtered files will be stored

    os.makedirs(output_dir_filter, exist_ok=True)  # Ensure the output directory exists
    os.makedirs(output_dir_pruned, exist_ok=True)  # Ensure the output directory exists

    # LD pruning parameters
    window_size = "50"
    step_size = "5"
    r2_threshold = "0.2"

    # Loop through chromosomes 1 to 22
    for chrom in range(1, 23):
        print(f'Beginning filtering, and pruning and thresholding for chromosome {chrom}')
        filtered_file = f"{output_dir_filter}/1000G.EUR.{chrom}"
        pruned_prefix = f"{output_dir_pruned}/1000G.EUR.{chrom}"

        # Step 1: Filter samples using --keep
        filter_command = [
            "plink",
            "--bfile", f"{original_file}{chrom}",
            "--keep", keep_file,
            "--make-bed",
            "--out", filtered_file
        ]
        
        # Step 2: Perform LD pruning
        prune_command = [
            "plink",
            "--bfile", filtered_file,
            "--indep-pairwise", window_size, step_size, r2_threshold,
            "--out", pruned_prefix
        ]

        # Step 3: Apply pruning (extract only pruned SNPs)
        threshold_command = [
            "plink",
            "--bfile", filtered_file,
            "--extract", f"{pruned_prefix}.prune.in",
            "--make-bed",
            "--out", pruned_prefix
        ]

        # Run the filtering step
        try:
            subprocess.run(filter_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"Successfully filtered chromosome {chrom}")
        except subprocess.CalledProcessError as e:
            print(f"Error filtering chromosome {chrom}: {e.stderr.decode('utf-8')}")
            continue  # Skip to next chromosome if error occurs

        # Run the pruning step
        try:
            subprocess.run(prune_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"Successfully performed LD pruning for chromosome {chrom}")
        except subprocess.CalledProcessError as e:
            print(f"Error in LD pruning for chromosome {chrom}: {e.stderr.decode('utf-8')}")
            continue

        # Run the thresholding step
        try:
            subprocess.run(threshold_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"Successfully applied pruning for chromosome {chrom}")
        except subprocess.CalledProcessError as e:
            print(f"Error applying pruning for chromosome {chrom}: {e.stderr.decode('utf-8')}")

    print("Filtering, LD pruning, and thresholding completed for all chromosomes")
