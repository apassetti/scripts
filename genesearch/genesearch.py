#!/usr/bin/env python3

import subprocess
import os
import sys
import shlex

# Define reference and input files
home_dir = os.path.expanduser("~")
ref_file = "gencode.v46.primary_assembly.annotation.gtf.gz"
save_file = os.path.join(home_dir, "Genes")


# Create save directory if it does not exist
if not os.path.exists(save_file):
    os.makedirs(save_file)


def main():
    # Check if reference file exists
    if not os.path.isfile(ref_file):
        command = f"wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/{ref_file}"
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            sys.exit(f"An error occurred while downloading the reference file: {e}")

    # Ask for gene name
    gene_name = input("Enter gene name: ").strip().upper()

    if not gene_name:
        sys.exit("Gene name cannot be empty")

    # Execute bash command
    command = f"gunzip -c {shlex.quote(ref_file)} | grep -v '#' | awk '$3 == \"gene\" && $9 ~ \"gene_name={gene_name}\"' | sort -t ';' -k4,4V | tee {shlex.quote(save_file)}/{shlex.quote(gene_name)}.txt | less -S -x 15"
    try:
        subprocess.run(command, shell=True, check=True)
        save_output = input("Do you want to save the gene information to a file? (y/n): ").strip().lower()
        
        if save_output == "n":
            os.remove(f"{save_file}/{gene_name}.txt")
            print("Gene information not saved")
        else:            
            print(f"Gene information saved to {save_file}/{gene_name}.txt")
    
    except subprocess.CalledProcessError as e:
        sys.exit(f"An error occurred while executing the command: {e}")

    # Ask if user wants to search again
    search_again = input("Do you want to search for another gene? (y/n): ").strip().lower()
    if search_again == "y":
        main()
    else:
        sys.exit("Goodbye!")

if __name__ == "__main__":
    main()
else:
    sys.exit("Goodbye!")