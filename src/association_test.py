import os
import pandas as pd
import numpy as np
import subprocess
import urllib.request
from pathlib import Path
import regex as re
from sys import platform
plat_os = platform

PLINK_PATH="plink.exe"
if plat_os != 'win32':
    PLINK_PATH = "./plink"

FUSION_COMPUTE_SCRIPT = os.path.join("fusion_twas-master","FUSION.compute_weights.R")
RSCRIPT_PATH="Rscript"
BFILE_EQTL_TEMPLATE = os.path.join("..","data","LDREF_filtered","1000G.EUR.{chr}_filtered")  # Raw genotype data
BFILE_TEMPLATE = os.path.join("..","data","LDREF_cis_eqtls","1000G.EUR.{chr}")  # Cis-eQTL processed data
TMP_DIR = "temp_files"
GCTA_PATH = os.path.join("fusion_twas-master","gcta64")
PHENO_DIR = os.path.join("..","data","gene_expressions","chr_{chr}")  # Directory with phenotype files


CIS_WINDOW = 500_000  # 500kb cis-window
MAF_THRESHOLD = 0.01  # Minor allele frequency filter
HWE_THRESHOLD = 1e-6  # Hardy-Weinberg filter
LD_PRUNE_WINDOW = 50  # LD pruning window size
LD_PRUNE_STEP = 5  # LD pruning step size
LD_PRUNE_R2 = 0.2  # LD pruning R^2 threshold
VERBOSE = 2
def get_all_Rdat(directory):
  """
  Gets a list of all .dat files in a directory and its subdirectories.

  Args:
    directory: The path to the directory.

  Returns:
    A list of file paths.
  """
  directory = Path(directory)

  # Search recursively for .dat files
  dat_files = list(directory.rglob("*.RDat"))
  return dat_files


def create_pos_files(INPUT_DIR=os.path.join("..", "data", "weights")):
    """
    Creates position file required for FUSION.assoc_test.R
    """
    gene_annot = pd.read_csv(os.path.join("..", "data", "gene_annot.txt.gz"), compression='gzip', sep='\t')
    pos = []
    for chrom in [str(i) for i in range(1, 23)]:
      chrom_weight_dir = os.path.join(INPUT_DIR, f"chr{chrom}")
      weight_files = get_all_Rdat(chrom_weight_dir)  
      for f in weight_files:
          f = str(f)
          sym = re.findall(r'ENSG\d+', f)[0]
          CHR = re.findall(r'chr(\d+)', f)[0]        

          if gene_annot[gene_annot['SYM'].str.contains(sym)].shape[0] > 0:
              ID = gene_annot[gene_annot['SYM'].str.contains(sym)]['ID'].iloc[0]
              P0 = gene_annot[gene_annot['SYM'].str.contains(sym)]['START'].iloc[0]
              P1 = gene_annot[gene_annot['SYM'].str.contains(sym)]['STOP'].iloc[0]
              pos.append([re.findall(r'ENSG.*', f)[0], ID, CHR, P0, P1])
      pos_df = pd.DataFrame(pos, columns=['WGT', 'ID', 'CHR', 'P0', 'P1'])
      pos_df.to_csv(os.path.join(chrom_weight_dir, f"chr{CHR}_weights.pos"), sep='\t', index=False)
      print(f"Created gene position file for chromosome {chrom}")
      

def generate_snp_list_helper(chromosome):
    """
    Creates list of SNPs for a chromosome required for FUSION.assoc_test.R to run
    """
    plink_path = "plink"
    output_fp = os.path.join("..","data","LDREF_filtered",f"snp_list_{chromosome}")
    
    command = [
        plink_path,
        "--bfile", os.path.join("..","data","LDREF_filtered",f"1000G.EUR.{chromosome}"),
        "--write-snplist",
        "--out", output_fp
    ]
    
    subprocess.run(command, check=True)
    print(f"PLINK SNP extraction complete for chromosome {chromosome}.")
    

def generate_snp_list():
    """
    Creates list of SNPs required for FUSION.assoc_test.R to run
    """
    for chrom in [str(i) for i in range(1, 23)]:
        generate_snp_list_helper(chrom)

def process_sumstats_helper(chromosome, input_fp, output_fp):
    """
    Processes and formats a GWAS sumstat file for a single chromosome 
    for FUSION.assoc_test.R
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_fp), exist_ok=True)
    
    # Load data
    df = pd.read_csv(input_fp, compression='gzip', sep='\t', header=0)
    
    # Compute Z-score and CHISQ
    df["Z"] = df["inv_var_meta_beta"] / df["inv_var_meta_sebeta"]
    df["CHISQ"] = df["Z"] ** 2
    
    # Select and rename columns
    df_sumstats = df.rename(columns={
        "rsid": "SNP",
        "ALT": "A1",
        "REF": "A2"
    })[["SNP", "A1", "A2", "N_ctrl", "CHISQ", "Z"]]
    
    # Rename N_ctrl to N
    df_sumstats = df_sumstats.rename(columns={"N_ctrl": "N"})
    
    # Filter for valid SNPs
    df_sumstats = df_sumstats[df_sumstats['SNP'].notna()]
    df_sumstats = df_sumstats[df_sumstats['A1'].apply(len) == 1]
    df_sumstats = df_sumstats[df_sumstats['A2'].apply(len) == 1]
    
    # Load SNP set
    snp_set = set()
    snp_list_fp = os.path.join("..","data","LDREF_filtered",f"snp_list_{chromosome}.snplist")
    with open(snp_list_fp, "r") as f:
        for line in f:
            snp_set.add(line.strip())
    
    print(f"Loaded {len(snp_set):,} SNPs.")
    print("Example SNPs:", list(snp_set)[:5])
    
    # Filter based on SNP set
    df_sumstats = df_sumstats[df_sumstats['SNP'].isin(snp_set)]
    
    # Save processed file
    df_sumstats.to_csv(output_fp, sep="\t", index=False)
    
    print(f"Processed file saved to: {output_fp}")

def process_sumstats():
    pops = ['afr', 'amr', 'eas', 'eur']
    chroms = [str(i) for i in range(1, 23)]
    gwas_base_fp = os.path.join("..", "data", "gwas", "HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz")
    out_base_fp = os.path.join("..", "data", "gwas", "HF_Bothsex_{pop}_{chrom}.sumstats")

    for pop in pops:
        for chrom in chroms:
            process_sumstats_helper(chrom, gwas_base_fp.format(pop=pop), \
                                    out_base_fp.format(pop=pop, chrom=chrom))
            

WEIGHTS_DIR = os.path.join("..", "data", "weights", "chr{chrom}")
FUSION_ASSOC_SCRIPT = os.path.join("fusion_twas-master","FUSION.assoc_test.R")
SUMSTATS_FILE = os.path.join("..","data","gwas","HF_Bothsex_{pop}_{chrom}.sumstats")
WEIGHTS_FILE = os.path.join("..","data","weights","chr{chrom}","chr{chrom}_weights.pos")
OUTPUT_FILE = os.path.join("..","data","gene_associations","chr{chrom}.{pop}.assoc.dat")

def association_test_helper(chrom, pop):
    """Runs the FUSION TWAS association test."""
    command = [
        RSCRIPT_PATH, FUSION_ASSOC_SCRIPT,
        "--sumstats", SUMSTATS_FILE.format(pop=pop, chrom=chrom),
        "--weights", WEIGHTS_FILE.format(chrom=chrom),
        "--weights_dir", WEIGHTS_DIR.format(chrom=chrom),
        "--ref_ld_chr", os.path.join("..","data","LDREF_filtered","1000G.EUR."),
        "--chr", chrom,
        "--out", OUTPUT_FILE.format(chrom=chrom, pop=pop)
    ]
    
    # Print command for debugging
    print(f"Running: {' '.join(command)}")
    
    # Run the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running FUSION TWAS association test: {e}")


def association_test():
    pops = ['afr', 'amr', 'eur', 'eas']
    chroms = [str(i) for i in range(1, 23)]
    for pop in pops:
        for chrom in chroms:
            association_test_helper(chrom, pop)


