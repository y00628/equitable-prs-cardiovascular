import pandas as pd
import numpy as np
import subprocess
import os

def compute_weights_win(
        chromosomes=[], 
        RSCRIPT_PATH="Rscript", 
        FUSION_SCRIPT="fusion_twas-master/FUSION.compute_weights.R", 
        BFILE_TEMPLATE="../data/LDREF_filtered/1000G.EUR.{chr}", 
        TMP_DIR="temp_files/", 
        OUTPUT_DIR="../data/weights/chr_{chr}", 
        MODELS="enet", 
        PLINK_PATH="plink.exe", 
        GCTA_PATH="fusion_twas-master/gcta_nr_robust.exe", 
        PHENO_DIR="../data/gene_expressions/chr_{chr}/", 
        VERBOSE=2
        ):
    

    '''
    Runs FUSION.compute_weights.R
    '''
    if not chromosomes:
        print('No chromosomes selected')
        return
    
    os.makedirs(TMP_DIR, exist_ok=True)


    phenotypes = pd.read_csv('../data/std_GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', compression='gzip', sep='\t')
    if chromosomes:
        filtered_phenotypes = phenotypes[phenotypes['Chr'].isin(chromosomes)]
    else:
        filtered_phenotypes = phenotypes[~phenotypes['Chr'].isin(['X', 'Y', 'M'])]
    to_calc = filtered_phenotypes[['TargetID', 'Chr']]
    count = 0
    def compute_weights_helper_win(gene, chromosome):
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
            "--pheno", pheno_file
        ]
        clear_command = ['rm', '-rf', f'{tmp_file_prefix}/*']

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
        compute_weights_helper_win(gene_data['TargetID'], gene_data['Chr'])