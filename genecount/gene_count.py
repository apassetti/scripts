#!/usr/bin/env python3
import os
import subprocess

# Path to the GFF3 file
gff3_file = "~/src/pseudogenes/data/gencode.v46.annotation.gff3.gz"
output_file = "~/src/pseudogenes/genes_lines.txt"

# Create output file
cmd = f'~/src/pseudogenes/genes_lines.txt'

# Expand the ~ to the full path
output_file = os.path.expanduser(output_file)

# Function to process a chromosome and append the count to the output file
def process_chromosome(chromosome):
    # Build the command
    cmd = f'gunzip -c {gff3_file} | grep -v "#" | awk \'$1 == "chr{chromosome}" && $3 == "gene"\' | wc -l'
    
    # Execute the command
    count = subprocess.check_output(cmd, shell=True).decode('utf-8').strip()
    
    # Append the result to the output file
    with open(output_file, 'a') as f:
        f.write(f"chr{chromosome} {count}\n")

# Loop through chromosomes
for i in range(1, 23):  # For chromosomes 1 to 22
    process_chromosome(i)

# Optional: Add the X and Y chromosomes
process_chromosome("X")
process_chromosome("Y")
