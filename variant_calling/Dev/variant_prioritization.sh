#!/bin/bash

#check directories
home_dir="${HOME}"
folder="${home_dir}/src/CFTR/"
results="${folder}/Results"

mkdir -p ${folder}/Curated

#Extract Funcotator field
cat ${results}/analysis-ready-snps-filteredGT-functotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${folder}/Curated/snp_curated_variants.txt
cat ${results}/analysis-ready-indels-filteredGT-functotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${folder}/Curated/indels_curated_variants.txt

#Extract annotated variant information
cat ${results}/output_snps.table | cut -f 5 | sed 's/|/\t/g' >> ${folder}/Curated/snp_curated_variants.txt
cat ${results}/output_indels.table | cut -f 5 | sed 's/|/\t/g' >> ${folder}/Curated/indels_curated_variants.txt