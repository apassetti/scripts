#!/usr/bin/env python3

import subprocess
import os
import sys
import shlex

# Define reference and input files
home_dir = os.path.expanduser("~")
ref_file = os.path.join(home_dir, "Reference", "hg38", "gencode.v46.annotation.gff3.gz")

def main():
    # Check if reference file exists
    if not os.path.isfile(ref_file):
        sys.exit(f"Reference file {ref_file} does not exist")

    # Ask for gene name
    gene_name = input("Enter gene name: ").strip()

    if not gene_name:
        sys.exit("Gene name cannot be empty")

    # Execute bash command
    command = f"gunzip -c {shlex.quote(ref_file)} | grep -v '#' | awk '$3 == \"gene\" && $9 ~ \"gene_name={shlex.quote(gene_name)}\"' | less -S -x 15"
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"An error occurred while executing the command: {e}")

if __name__ == "__main__":
    main()