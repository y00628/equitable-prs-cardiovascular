import os
from ftplib import FTP

# Base FTP server and directory
ftp_host = "ftp.1000genomes.ebi.ac.uk"
ftp_base_dir = "/vol1/ftp/release/20130502/"

# List of chromosomes to download (1 to 22)
chromosomes = range(1, 23)

# Output directory
output_dir = "1000_genomes_vcf"
os.makedirs(output_dir, exist_ok=True)

# Function to download a file via FTP
def download_file_ftp(ftp, remote_path, local_path):
    try:
        with open(local_path, "wb") as file:
            ftp.retrbinary(f"RETR {remote_path}", file.write)
        print(f"Downloaded: {local_path}")
    except Exception as e:
        print(f"Failed to download {remote_path}: {e}")

# Connect to the FTP server
try:
    ftp = FTP(ftp_host)
    ftp.login()
    ftp.cwd(ftp_base_dir)

    # Loop through all chromosomes and download their VCF and TBI files
    for chrom in chromosomes:
        # Correct file names
        vcf_file = f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        tbi_file = f"{vcf_file}.tbi"

        # Check if the files exist on the server
        if vcf_file in ftp.nlst():
            # Local file paths
            vcf_local_path = os.path.join(output_dir, vcf_file)
            tbi_local_path = os.path.join(output_dir, tbi_file)

            # Download files
            download_file_ftp(ftp, vcf_file, vcf_local_path)
            download_file_ftp(ftp, tbi_file, tbi_local_path)
        else:
            print(f"File not found on FTP server: {vcf_file}")

    ftp.quit()
    print("All downloads complete!")

except Exception as e:
    print(f"FTP connection failed: {e}")