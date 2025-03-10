import pandas as pd
import numpy as np
import subprocess
import os
from sys import platform

all_platforms = ['linux', 'win32', 'darwin']

plat_os = platform
RSCRIPT_PATH="Rscript"
FUSION_SCRIPT=os.path.join("fusion_twas-master", "FUSION.compute_weights.R")
BFILE_TEMPLATE=os.path.join("..","data","LDREF_pruned","1000G.EUR.{chr}")
TMP_DIR="temp_files"
OUTPUT_DIR=os.path.join("..","data","weights","chr_{chr}")
MODELS="top1,enet,lasso"
PLINK_PATH="plink.exe"
GCTA_PATH=os.path.join("fusion_twas-master","gcta_nr_robust.exe")
if plat_os != 'win32':
    PLINK_PATH = "./plink"
    GCTA_PATH=os.path.join("fusion_twas-master","gcta_nr_robust")
PHENO_DIR=os.path.join("..","data","gene_expressions","chr_{chr}")
VERBOSE=2


def compute_weights(
        chromosomes=[], 
        RSCRIPT_PATH=RSCRIPT_PATH, 
        FUSION_SCRIPT=FUSION_SCRIPT, 
        BFILE_TEMPLATE=BFILE_TEMPLATE, 
        TMP_DIR=TMP_DIR, 
        OUTPUT_DIR=OUTPUT_DIR, 
        MODELS=MODELS, 
        PLINK_PATH=PLINK_PATH, 
        GCTA_PATH=GCTA_PATH, 
        PHENO_DIR=PHENO_DIR,
        VERBOSE=2
        ):
    

    '''
    Runs FUSION.compute_weights.R for LDREF EUR data
    '''
    if not chromosomes:
        print('No chromosomes selected')
        return
    
    os.makedirs(TMP_DIR, exist_ok=True)


    phenotypes = pd.read_csv(os.path.join("..","data","std_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"), compression='gzip', sep='\t')
    if chromosomes:
        filtered_phenotypes = phenotypes[phenotypes['Chr'].isin(chromosomes)]
    else:
        filtered_phenotypes = phenotypes[~phenotypes['Chr'].isin(['X', 'Y', 'M'])]
    to_calc = filtered_phenotypes[['TargetID', 'Chr']]

    def generate_pca(chromosome):
        """Generates PCA components once per chromosome and returns the PCA file path."""
        pca_output_path = os.path.join(TMP_DIR, f"1000G.EUR.{chromosome}.pca10")

        # Check if PCA has already been computed to avoid redundant computation
        if os.path.exists(pca_output_path + ".eigenvec"):
            print(f"PCA already exists for chromosome {chromosome}, skipping recomputation.")
            return pca_output_path + ".eigenvec"

        command_pca = [
            PLINK_PATH, "--bfile", BFILE_TEMPLATE.format(chr=chromosome),
            "--pca", "10",
            "--out", pca_output_path
        ]

        try:
            subprocess.run(command_pca, check=True)
            return pca_output_path + ".eigenvec"  # Return path to PCA file
        except subprocess.CalledProcessError as e:
            print(f"Error generating PCA for chromosome {chromosome}: {e}")
            return None

    # Compute PCA once per chromosome
    chromosome_pca_paths = {} 
    for chromosome in chromosomes:
        cor_path = generate_pca(chromosome) 
        chromosome_pca_paths[chromosome] = cor_path 

    count = 0

    def compute_weights_helper_win(gene, chromosome, cor_path):
        """Runs the FUSION TWAS pipeline for a specific gene and chromosome."""

        # Construct file paths
        bfile = BFILE_TEMPLATE.format(chr=chromosome)
        output_dir = OUTPUT_DIR.format(chr=chromosome)
        pheno_dir = PHENO_DIR.format(chr=chromosome)
        
        os.makedirs(output_dir, exist_ok=True)


        tmp_file_prefix = os.path.join(TMP_DIR, f"test_chr{chromosome}_{gene}")
        output_path = os.path.join(output_dir, f"{gene}_chr{chromosome}")
        pheno_file = os.path.join(pheno_dir, f"{gene}.txt")

        
        # Construct the command
        command = [
            RSCRIPT_PATH, FUSION_SCRIPT,
            "--bfile", bfile,
            "--tmp", tmp_file_prefix,
            "--out", output_path,
            "--models", MODELS,
            "--PATH_plink", PLINK_PATH,
            "--PATH_gcta", GCTA_PATH,
            "--verbose", str(VERBOSE),
            "--pheno", pheno_file, 
            "--covar", cor_path
        ]
        clear_command = ['rm', '-rf', os.path.join(tmp_file_prefix, "*")]

        # Print command for debugging
        print(f"({count}/{to_calc.shape[0]}) Running: {' '.join(command)}")
        
        # Run the command
        try:
            subprocess.run(command, check=True)
            subprocess.run(clear_command, check=True)
            
        except subprocess.CalledProcessError as e:
            print(f"Error running command for gene {gene}, chromosome {chromosome}: {e}")

    for i in to_calc.index:
        count += 1
        gene_data = to_calc.loc[i]
        cor_path = chromosome_pca_paths.get(gene_data['Chr'], None)
        compute_weights_helper_win(gene_data['TargetID'], gene_data['Chr'], cor_path)
