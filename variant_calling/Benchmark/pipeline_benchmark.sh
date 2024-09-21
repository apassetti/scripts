#!/bin/bash

# Variant calling pipeline 
# syntetich benchmark using https://github.com/nh13/DWGSIM

# Folders
folder="/Users/andreapassetti/src/synthetic"
Data="${folder}/Data"
Ref="/Users/andreapassetti/Reference/hg38"
Genome="${Ref}/hg38_chr7.fasta"
Analysis="${folder}/Analysis"
Results="${folder}/Results"

mkdir -p ${Analysis}
mkdir -p ${Analysis}/{Alignment,VariantCalling}

if false
then

# Indexing reference genome
bwa index ${Genome}
samtools faidx ${Genome}
samtools dict ${Genome} > ${Ref}/hg38_chr7.dict

# Map to reference genome using BWA-MEM
bwa mem -t 8 -R "@RG\tID:bench\tSM:sample\tPL:ILLUMINA" ${Genome} ${Data}/output.bwa.read1.fastq ${Data}/output.bwa.read2.fastq > ${Analysis}/Alignment/bench_chr7.sam 
fi

# Mark duplicates (GATK)
# MarkDuplicatesSpark performs both the duplicate marking step and the sort step
gatk MarkDuplicatesSpark -I ${Analysis}/Alignment/bench_chr7.sam -O ${Analysis}/Alignment/bench_chr7_dedup_reads.bam -M ${Analysis}/Alignment/bench_chr7_dedup_metrics.txt

# Base quality recalibration
# 1. Build the model
gatk BaseRecalibrator -I ${Analysis}/Alignment/bench_chr7_dedup_reads.bam \
    -R ${Genome} \
    --known-sites ${Ref}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${Ref}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${Analysis}/Alignment/recal_data.table

# 2. Apply the model
gatk ApplyBQSR -I ${Analysis}/Alignment/bench_chr7_dedup_reads.bam \
    -R ${Genome} \
    --bqsr-recal-file ${Analysis}/Alignment/recal_data.table \
    -O ${Analysis}/Alignment/bench_chr7_sorted_dedup_recal_reads.bam

# 3. Analyze covarities
gatk AnalyzeCovariates -before ${Analysis}/Alignment/recal_data.table \
    -after ${Analysis}/Alignment/bench_chr7_sorted_dedup_recal_reads.bam \
    -plots ${Analysis}/Alignment/recalibration_plots.pdf

# Variant calling

