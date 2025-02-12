import pandas as pd
import os
import subprocess

def generate_snp_list(chromosome):
    plink_path = "plink"
    output_fp = f"../data/LDREF_pruned/snp_list_{chromosome}"
    
    command = [
        plink_path,
        "--bfile", f"../data/LDREF_pruned/1000G.EUR.{chromosome}",
        "--write-snplist",
        "--out", output_fp
    ]
    
    if not os.path.exists(output_fp + ".snplist"):
        subprocess.run(command, check=True)
        print(f"PLINK SNP extraction complete for chromosome {chromosome}.")
    else:
        print(f"SNP list for chromosome {chromosome} already exists.")

def process_sumstats(chromosome, input_fp, output_fp):
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_fp), exist_ok=True)
    
    # Generate SNP list if it doesn't exist
    generate_snp_list(chromosome)
    
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
    snp_list_fp = f"../data/LDREF_pruned/snp_list_{chromosome}.snplist"
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