#!/bin/bash

# Variant calling pipeline 
# GATK best practices

#Folders
home="${HOME}"
folder="${home}/GIAB"
Data="${home}/GIAB/Data"
Ref="/Users/andreapassetti/Reference"
Genome="${Ref}/hg38_chr7.fasta"
Analysis="${folder}/Analysis"

mkdir -p ${Analysis}
mkdir -p ${Analysis}/{Alignment,VariantCalling,VariantFiltering}

if false
then

#Indexing reference genome
bwa index ${Genome}
samtools faidx ${Genome}
samtools dict ${Genome} > ${Ref}/hg38_chr7.dict

fi
### MAPPING AND CLEANUP ###

#Merge reads
cat ${Data}/NIST7035_TAAGGCGA_L001_R1_001_trimmed.fastq ${Data}/NIST7035_TAAGGCGA_L002_R1_001_trimmed.fastq > ${Data}/NIST7035_TAAGGCGA_R1.fastq
cat ${Data}/NIST7035_TAAGGCGA_L001_R2_001_trimmed.fastq ${Data}/NIST7035_TAAGGCGA_L002_R2_001_trimmed.fastq > ${Data}/NIST7035_TAAGGCGA_R2.fastq

#Alignment -> usare hg19
bwa mem -t 8 -R "@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA" ${Genome} ${Data}/NIST7035_TAAGGCGA_R1.fastq ${Data}/NIST7035_TAAGGCGA_R2.fastq > ${Analysis}/Alignment/NA12878_L1_chr7.sam

#fi




