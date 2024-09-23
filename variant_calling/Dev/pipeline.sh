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
Stats="${folder}/Stats"

mkdir -p ${Analysis}
mkdir -p ${Analysis}/{Alignment,VariantCalling,VariantFiltering}


#Indexing reference genome
bwa index ${Genome}
samtools faidx ${Genome}
samtools dict ${Genome} > ${Ref}/hg38_chr7.dict


### MAPPING AND CLEANUP ###

#Merge reads
cat ${Data}/NIST7035_TAAGGCGA_L001_R1_001_trimmed.fastq ${Data}/NIST7035_TAAGGCGA_L002_R1_001_trimmed.fastq > ${Data}/NIST7035_TAAGGCGA_R1.fastq
cat ${Data}/NIST7035_TAAGGCGA_L001_R2_001_trimmed.fastq ${Data}/NIST7035_TAAGGCGA_L002_R2_001_trimmed.fastq > ${Data}/NIST7035_TAAGGCGA_R2.fastq

#Alignment -> usare hg19
bwa mem \
    -t 8 \
    -R "@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA" \
    ${Genome} \
    ${Data}/NIST7035_TAAGGCGA_R1.fastq \
    ${Data}/NIST7035_TAAGGCGA_R2.fastq > ${Analysis}/Alignment/NA12878_L1_chr7.sam

# Mark duplicates (GATK)
# MarkDuplicatesSpark performs both the duplicate marking step and the sort step
gatk MarkDuplicatesSpark \
    -I ${Analysis}/Alignment/NA12878_L1_chr7.sam \
    -O ${Analysis}/Alignment/NA12878_L1_chr7_dedup.bam \
    -M ${Stats}/Alignment/NA12878_L1_chr7_dedup_metrics.txt

# Get some stats
samtools stats ${Analysis}/Alignment/NA12878_L1_chr7_dedup.bam > ${Stats}/Alignment/NA12878_L1_chr7_dedup_stats.txt
mosdepth -x -n --by ${Genome} \
    ${Stats}/Alignment/NA12878_L1_chr7_dedup_stats \
    ${Analysis}/Alignment/NA12878_L1_chr7_dedup.bam

# Base quality recalibration

# 1. Build the pre-recalibration model (before applying recalibration)
gatk BaseRecalibrator -I ${Analysis}/Alignment/NA12878_L1_chr7_dedup.bam \
    -R ${Genome} \
    --known-sites ${Ref}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${Ref}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${Analysis}/Alignment/recal_data_pre.table

# 2. Apply the recalibration to the BAM file (generating recalibrated BAM)
gatk ApplyBQSR -I ${Analysis}/Alignment/NA12878_L1_chr7_dedup.bam \
    -R ${Genome} \
    --bqsr-recal-file ${Analysis}/Alignment/recal_data_pre.table \
    -O ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam

# 3. Build the post-recalibration model (generate recal_data1.table)
gatk BaseRecalibrator -I ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam \
    -R ${Genome} \
    --known-sites ${Ref}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${Ref}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${Analysis}/Alignment/recal_data1.table

# 4. Analyze covariates to compare before and after recalibration
gatk AnalyzeCovariates -before ${Analysis}/Alignment/recal_data_pre.table \
    -after ${Analysis}/Alignment/recal_data1.table \
    -plots ${Analysis}/Alignment/recalibration_plots.pdf \
    -csv ${Analysis}/Alignment/recalibration_plots.csv

# Variant calling
# Variant calling with GATK HaplotypeCaller
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R ${Genome} \
    -I ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam \
    -O ${Analysis}/VariantCalling/NA12878_L1_chr7.vcf.gz \
    --native-pair-hmm-threads 8

# Variant calling with freebayes
freebayes -f ${Genome} ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam > ${Analysis}/VariantCalling/NA12878_L1_chr7_freebayes.vcf

# Variant calling with samtools and bcftools
samtools mpileup -uf ${Genome} \
    ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam \
    | bcftools call -mv > ${Analysis}/VariantCalling/NA12878_L1_chr7_samtools.vcf

# Variant calling with strelka2
strelka --bam ${Analysis}/Alignment/NA12878_L1_chr7_sorted_dedup_recal.bam \
    --referenceFasta ${Genome} \
    --runDir ${Analysis}/VariantCalling/NA12878_L1_chr7_strelka2

# Concatenate all vcf files
vcf-concat ${Analysis}/VariantCalling/NA12878_L1_chr7.vcf.gz \
    ${Analysis}/VariantCalling/NA12878_L1_chr7_freebayes.vcf \
    ${Analysis}/VariantCalling/NA12878_L1_chr7_samtools.vcf \
    ${Analysis}/VariantCalling/NA12878_L1_chr7_strelka2/results/variants/variants.vcf \
    > ${Analysis}/VariantCalling/NA12878_L1_chr7_all.vcf

# Filter variants
vcftools --vcf ${Analysis}/VariantCalling/NA12878_L1_chr7_all.vcf \
    --minQ 30 \
    --minDP 10 \
    --recode \
    --out ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered

# Annotate variants with snpeff
snpeff ann hg38_chr7 ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered.recode.vcf > ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_snpeff.vcf

# Annotate variants with Funcotator
gatk Funcotator \
    -R ${Genome} \
    -V ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered.recode.vcf \
    -O ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_funcotated.vcf \
    --data-sources-path /Users/andreapassetti/Reference/hg38/functotator/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38 \
    --output-file-format VCF

# Annotate variants with VEP
vep -i ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered.recode.vcf \
    -o ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_vep.vcf \
    --cache \
    --dir_cache /Users/andreapassetti/Reference/vep \
    --assembly GRCh38 \
    --fasta ${Genome} \
    --force_overwrite

# Compare vcf files
vcftools --vcf ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_funcotated.vcf \
    --diff ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_vep.vcf \
    --out ${Analysis}/VariantFiltering/differences \
    --diff-site

# Prioritize variants with SnpSift
SnpSift filter "(ANN[0].IMPACT = 'HIGH') & (gnomAD_AF < 0.01)" \
    ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_snpeff.vcf \
    > ${Analysis}/VariantFiltering/NA12878_L1_chr7_filtered_snpeff_prioritized.vcf
