#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

# Default parameters (can be overridden via command-line arguments or environment variables)

home_dir="${HOME}"
folder="${home_dir}/src/CFTR"
ref="${home_dir}/Reference/Homo_sapiens_assembly38.fasta"
known_sites="${home_dir}/Reference/Homo_sapiens_assembly38.dbsnp138.vcf"
threads=8  # Default number of threads for parallelization

# Allow user to override the parameters via command-line arguments
while getopts f:r:k:t: flag
do
    case "${flag}" in
        f) folder=${OPTARG};;  # Project folder
        r) ref=${OPTARG};;     # Reference genome
        k) known_sites=${OPTARG};;  # Known sites (dbSNP)
        t) threads=${OPTARG};;  # Number of threads
    esac
done

# Function to log with timestamps
log_with_time() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Input validation function
validate_inputs() {
    if [ ! -f "${ref}" ]; then
        log_with_time "ERROR: Reference genome file ${ref} not found."
        exit 1
    fi

    if [ ! -f "${known_sites}" ]; then
        log_with_time "ERROR: Known sites VCF file ${known_sites} not found."
        exit 1
    fi

    if [ ! -d "${folder}" ]; then
        log_with_time "ERROR: Project folder ${folder} does not exist."
        exit 1
    fi
}

# Call input validation function
validate_inputs

# Print selected configuration
log_with_time "Configuration:"
log_with_time "  Project directory: $folder"
log_with_time "  Reference genome: $ref"
log_with_time "  Known sites VCF: $known_sites"
log_with_time "  Number of threads: $threads"

# Create directories
mkdir -p ${folder}/Aligned ${folder}/Data ${folder}/Reads ${folder}/Results ${folder}/Quality 

# Log file path
log_file="${folder}/pipeline.log"
exec > >(tee -i ${log_file}) 2>&1  # Log output and errors to the console and log file

# Error handling function
error_exit() {
    log_with_time "ERROR on line $1"
    exit 1
}

trap 'error_exit $LINENO' ERR

# Download data
log_with_time "Downloading data with prefetch and fasterq-dump"
prefetch SRR2136533
fasterq-dump -x SRR2136533 -O ${folder}/Reads

###################################################### VARIANT CALLING STEPS ####################################################################

# Directories based on project folder
aligned_reads="${folder}/Aligned"
reads="${folder}/Reads"
results="${folder}/Results"
data="${folder}/Data"
quality="${folder}/Quality"

# -------------------       
# STEP 1: QC - Run fastqc 
# -------------------

log_with_time "STEP 1: QC - Running fastqc"

parallel --eta fastqc {} -o ${quality}/ ::: ${reads}/SRR2136533_1.fastq ${reads}/SRR2136533_2.fastq

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

log_with_time "STEP 2: Mapping to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}
samtools faidx ${ref}

# BWA alignment with parallel threads
bwa mem -t ${threads} -R "@RG\tID:SRR2136533\tPL:ILLUMINA\tSM:SRR2136533" ${ref} ${reads}/SRR2136533_1.fastq ${reads}/SRR2136533_2.fastq > ${aligned_reads}/SRR2136533.paired.sam


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

log_with_time "STEP 3: Marking Duplicates and Sorting"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR2136533.paired.sam -O ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam


# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------

log_with_time "STEP 4: Base quality recalibration"

# 1. Build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam 


# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

log_with_time "STEP 5: Collecting Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam -O ${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller (parallelized)
# ----------------------------------------------

log_with_time "STEP 6: Calling Variants with GATK HaplotypeCaller"

# Parallel variant calling using multithreading
gatk HaplotypeCaller \
    -R ${ref} \
    -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam \
    -O ${results}/raw_variants.vcf \
    --native-pair-hmm-threads ${threads}


# extract SNPs & INDELS
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf


# -------------------
# STEP 7: Filter Variants - GATK4
# -------------------

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