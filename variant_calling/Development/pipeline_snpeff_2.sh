#!/bin/bash

# Script to call germline variants in a human WGS paired-end reads 2 X 100bp
# GATK4 Best Practices Workflow

###################################################### DOWNLOAD STEP ####################################################################

folder="/Users/andreapassetti/src/CFTR"

mkdir -p ${folder}/{Aligned,Data,Reads,Results,Quality}

prefetch SRR2136533
fasterq-dump -x SRR2136533 -O ${folder}/Reads

###################################################### VARIANT CALLING STEPS ####################################################################

# Directories
ref="/Users/andreapassetti/Reference/Homo_sapiens_assembly38.fasta"
known_sites="/Users/andreapassetti/Reference/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="${folder}/Aligned"
reads="${folder}/Reads"
results="${folder}/Results"
data="${folder}/Data"
quality="${folder}/Quality"

# STEP 1: QC
fastqc -t 8 ${reads}/SRR2136533_1.fastq -o ${quality}/
fastqc -t 8 ${reads}/SRR2136533_2.fastq -o ${quality}/

# STEP 2: Map to reference using BWA-MEM
bwa index ${ref}
samtools faidx ${ref}
bwa mem -t 8 -R "@RG\tID:SRR2136533\tPL:ILLUMINA\tSM:SRR2136533" ${ref} ${reads}/SRR2136533_1.fastq ${reads}/SRR2136533_2.fastq > ${aligned_reads}/SRR2136533.paired.sam

# STEP 3: Mark Duplicates and Sort
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR2136533.paired.sam -O ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam

# STEP 4: Base quality recalibration
gatk BaseRecalibrator -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
gatk ApplyBQSR -I ${aligned_reads}/SRR2136533_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam

# STEP 5: Collect Alignment & Insert Size Metrics
gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam -O ${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

# STEP 6: Call Variants
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR2136533_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf

# STEP 7: Filter Variants
gatk VariantFiltration -R ${ref} -V ${results}/raw_snps.vcf -O ${results}/filtered_snps.vcf \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0" \
  -filter-name "SOR_filter" -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
  -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
  -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

gatk VariantFiltration -R ${ref} -V ${results}/raw_indels.vcf -O ${results}/filtered_indels.vcf \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 200.0" \
  -filter-name "SOR_filter" -filter "SOR > 10.0" \
  -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
  -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

gatk SelectVariants --exclude-filtered -V ${results}/filtered_snps.vcf -O ${results}/analysis-ready-snps.vcf
gatk SelectVariants --exclude-filtered -V ${results}/filtered_indels.vcf -O ${results}/analysis-ready-indels.vcf

# Annotate VCF
snpEff -Xmx8g -v GRCh38.99 ${results}/analysis-ready-snps.vcf > ${results}/analysis-ready-snps-annotated.vcf
snpEff -Xmx8g -v GRCh38.99 ${results}/analysis-ready-indels.vcf > ${results}/analysis-ready-indels-annotated.vcf

# Filter Annotated Variants
SnpSift filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" ${results}/analysis-ready-snps-annotated.vcf > ${results}/analysis-ready-snps-highImpact.vcf
SnpSift filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" ${results}/analysis-ready-indels-annotated.vcf > ${results}/analysis-ready-indels-highImpact.vcf