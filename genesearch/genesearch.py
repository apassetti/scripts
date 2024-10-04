#!/usr/bin/env python3

import subprocess
import os
import sys
import shlex

# Define reference and input files
home_dir = os.path.expanduser("~")
ref_file = os.path.join(home_dir, "Reference", "hg38", "gencode.v46.annotation.gff3.gz")
save_fle = os.path.join(home_dir, "Documents")

def main():
    # Check if reference file exists
    if not os.path.isfile(ref_file):
        sys.exit(f"Reference file {ref_file} does not exist")

    # Ask for gene name
    gene_name = input("Enter gene name: ").strip().upper()

    if not gene_name:
        sys.exit("Gene name cannot be empty")

    # Execute bash command
    command = f"gunzip -c {shlex.quote(ref_file)} | grep -v '#' | awk '$3 == \"gene\" && $9 ~ \"gene_name={gene_name}\"' | tee {shlex.quote(save_fle)}/{shlex.quote(gene_name)}.txt | less -S -x 15"
    try:
        subprocess.run(command, shell=True, check=True)
        save_output = input("Do you want to save the gene information to a file? (y/n): ").strip().lower()
        
        if save_output == "n":
            os.remove(f"{save_fle}/{gene_name}.txt")
            sys.exit("Gene information not saved")
        else:            
            print(f"Gene information saved to {save_fle}/{gene_name}.txt")
    
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