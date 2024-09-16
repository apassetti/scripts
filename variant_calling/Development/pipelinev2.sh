#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

###################################################### INITIAL STEP ####################################################################

# Get the current home folder dynamically
home_dir=$HOME
echo "Current home directory: $home_dir"

# Make directories
folder="${home_dir}/src/CFTR"
mkdir -p ${folder}/Aligned ${folder}/Data ${folder}/Reads ${folder}/Results ${folder}/Quality 

# Log file path
log_file="${folder}/pipeline.log"
exec > >(tee -i ${log_file}) 2>&1  # Log output and errors to the console and log file

# Error handling function
error_exit() {
    echo "Error on line $1"
    exit 1
}

trap 'error_exit $LINENO' ERR

# Define checkpoint file
CHECKPOINT_FILE="${folder}checkpoint.txt"

# Function to set the checkpoint
set_checkpoint() {
    echo "$1" > "$CHECKPOINT_FILE"
}

# Function to get the current checkpoint
get_checkpoint() {
    if [ -f "$CHECKPOINT_FILE" ]; then
        cat "$CHECKPOINT_FILE"
    else
        echo "0" # Default to 0 if no checkpoint is found
    fi
}

# Get the last checkpoint
checkpoint=$(get_checkpoint)

# Download data
prefetch SRR2136533
fasterq-dump -x SRR2136533 -O ${folder}/Reads

###################################################### VARIANT CALLING STEPS ####################################################################

# directories
ref="${home_dir}/Reference/Homo_sapiens_assembly38.fasta"
known_sites="${home_dir}/Reference/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="${folder}/Aligned"
reads="${folder}/Reads"
results="${folder}/Results"
data="${folder}/Data"
quality="${folder}/Quality"

# -------------------       
# STEP 1: QC - Run fastqc 
# -------------------

if [ "$checkpoint" -le 1 ]; then
    echo "STEP 1: QC - Run fastqc"
    parallel --eta fastqc {} -o ${quality}/ ::: ${reads}/SRR2136533_1.fastq ${reads}/SRR2136533_2.fastq
    
    set_checkpoint 2 
fi


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

if [ "$checkpoint" -le 2 ]; then
    echo "STEP 2: Map to reference using BWA-MEM"
    # BWA index reference 
    bwa index ${ref}
    samtools faidx ${ref}

    # BWA alignment with 8 threads
    bwa mem -t 8 -R "@RG\tID:SRR2136533\tPL:ILLUMINA\tSM:SRR2136533" ${ref} ${reads}/SRR2136533_1.fastq ${reads}/SRR2136533_2.fastq > ${aligned_reads}/SRR2136533.paired.sam
    
    set_checkpoint 3 
fi


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

if [ "$checkpoint" -le 3 ]; then
    echo "STEP 3: Mark Duplicates and Sort - GATK4"
    gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR2136533.paired.sam -O ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam
    
    set_checkpoint4 
fi


# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------

if [ "$checkpoint" -le 4 ]; then
    echo "STEP 4: Base quality recalibration"

    # 1. Build the model
    gatk BaseRecalibrator -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

    # 2. Apply the model to adjust the base quality scores
    gatk ApplyBQSR -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam 
    
    set_checkpoint5
fi


# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

if [ "$checkpoint" -le 5 ]; then
    echo "STEP 5: Collect Alignment & Insert Size Metrics"

    gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam -O ${aligned_reads}/alignment_metrics.txt
    gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf
    set_checkpoint6
fi


# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller (parallelized)
# ----------------------------------------------

if [ "$checkpoint" -le 6 ]; then
    echo "STEP 6: Call Variants - gatk haplotype caller"

    # Parallel variant calling using multiple regions (genome is large, splitting it speeds things up)
    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam \
        -O ${results}/raw_variants.vcf \
        -ERC GVCF \
        -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 \
        --native-pair-hmm-threads 8


    # extract SNPs & INDELS
    gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
    gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
    set_checkpoint7
fi


# -------------------
# STEP 7: Filter Variants - GATK4
# -------------------

if [ "$checkpoint" -le 7 ]; then
    echo "STEP 7: Filtering Variants"
    # Filter SNPs
    gatk VariantFiltration \
        -R ${ref} \
        -V ${results}/raw_snps.vcf \
        -O ${results}/filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"

    # Filter INDELS
    gatk VariantFiltration \
        -R ${ref} \
        -V ${results}/raw_indels.vcf \
        -O ${results}/filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"


    # Select Variants that PASS filters
    gatk SelectVariants --exclude-filtered -V ${results}/filtered_snps.vcf -O ${results}/analysis-ready-snps.vcf
    gatk SelectVariants --exclude-filtered -V ${results}/filtered_indels.vcf -O ${results}/analysis-ready-indels.vcf

    # to exclude variants that failed genotype filters
    grep -v -E "DP_filter|GQ_filter" ${results}/analysis-ready-snps.vcf > ${results}/analysis-ready-snps-filteredGT.vcf
    grep -v -E "DP_filter|GQ_filter" ${results}/analysis-ready-indels.vcf > ${results}/analysis-ready-indels-filteredGT.vcf

    # Annotate VCF with functional information
    snpEff -v GRCh38.99 ${results}/analysis-ready-snps-filteredGT.vcf > ${results}/analysis-ready-snps-annotated.vcf
    snpEff -v GRCh38.99 ${results}/analysis-ready-indels-filteredGT.vcf > ${results}/analysis-ready-indels-annotated.vcf

    # Filter annotated variants
    SnpSift filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" ${results}/analysis-ready-snps-annotated.vcf > ${results}/analysis-ready-snps-highImpact.vcf
    SnpSift filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" ${results}/analysis-ready-indels-annotated.vcf > ${results}/analysis-ready-indels-highImpact.vcf

    set_checkpoint 8 # Final checkpoint (script completed)
fi
echo "Script completed."