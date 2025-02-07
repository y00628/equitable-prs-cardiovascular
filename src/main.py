import os
import pandas as pd
import numpy as np
import subprocess
import urllib.request

def download_data():
    os.makedirs('data', exist_ok=True)

    # GWAS data for multiple populations

    pops = ['eur', 'eas', 'amr', 'afr']
    for pop in pops:

        url = f"https://gbmi-sumstats.s3.amazonaws.com/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
        os.makedirs('data/gwas', exist_ok=True)

        output_file = f"../data/gwas/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"

        try:
            print("Downloading file...")
            urllib.request.urlretrieve(url, output_file)
            print("Download completed successfully.")
        except Exception as e:
            print(f"An error occurred: {e}")
    return

def compute_weights(chromosomes=[str(i) for i in range(1, 23)], 
                    RSCRIPT_PATH = "Rscript",
                    FUSION_SCRIPT = "fusion_twas-master/FUSION.compute_weights.R",
                    BFILE_TEMPLATE = "../data/LDREF_filtered/1000G.EUR.{chr}_filtered",
                    TMP_DIR = "temp_files",
                    OUTPUT_DIR = "compute_weights_out/",
                    MODELS = "enet",
                    PLINK_PATH = "plink.exe",
                    GCTA_PATH = "fusion_twas-master/gcta_nr_robust.exe",
                    PHENO_DIR = "../data/gene_expressions",
                    VERBOSE = 0
):
    phenotypes = pd.read_csv('../data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', compression='gzip', sep='\t')
    filtered_phenotypes = phenotypes[~phenotypes['Chr'].isin(['X', 'Y', 'M'])]
    to_calc = filtered_phenotypes[['TargetID', 'Chr']]
    to_calc = to_calc[to_calc['Chr'].isin(chromosomes)]
    
    '''
    RSCRIPT_PATH = "Rscript"
    FUSION_SCRIPT = "fusion_twas-master/FUSION.compute_weights.R"
    BFILE_TEMPLATE = "../data/LDREF_filtered/1000G.EUR.{chr}_filtered"  # Placeholder for chromosome
    TMP_DIR = "temp_files"
    OUTPUT_DIR = "compute_weights_out/"
    MODELS = "enet"
    PLINK_PATH = "plink.exe"
    GCTA_PATH = "fusion_twas-master/gcta_nr_robust.exe"
    PHENO_DIR = "../data/gene_expressions"  # Directory with phenotype files
    VERBOSE = 0
    '''

    def compute_weights_helper(gene, chromosome):
        """Runs the FUSION TWAS pipeline for a specific gene and chromosome."""
        # Construct file paths
        bfile = BFILE_TEMPLATE.format(chr=chromosome)
        tmp_file_prefix = os.path.join(TMP_DIR, f"test_chr{chromosome}_{gene}")
        output_path = os.path.join(OUTPUT_DIR, f"{gene}_chr{chromosome}")
        pheno_file = os.path.join(PHENO_DIR, f"{gene}.txt")
        
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
            "--pheno", pheno_file
        ]

        # Print command for debugging
        print(f"Running: {' '.join(command)}")
        
        # Run the command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command for gene {gene}, chromosome {chromosome}: {e}")

    for i in to_calc.index:
        gene_data = to_calc.loc[i]
        compute_weights_helper(gene_data['TargetID'], gene_data['Chr'])

