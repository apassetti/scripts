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
bwa mem -t 8 -R "@RG\tID:bench\tSM:sample\tPL:ILLUMINA" ${Genome} ${Data}/chr7_reads.bwa.read1.fastq ${Data}/chr7_reads.bwa.read2.fastq > ${Analysis}/Alignment/bench_chr7.sam 

# Mark duplicates (GATK)
# MarkDuplicatesSpark performs both the duplicate marking step and the sort step
gatk MarkDuplicatesSpark -I ${Analysis}/Alignment/bench_chr7.sam -O ${Analysis}/Alignment/bench_chr7_dedup.bam -M ${Analysis}/Alignment/bench_chr7_dedup_metrics.txt


# Base quality recalibration

# 1. Build the pre-recalibration model (before applying recalibration)
gatk BaseRecalibrator -I ${Analysis}/Alignment/bench_chr7_dedup.bam \
    -R ${Genome} \
    --known-sites ${Ref}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${Ref}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${Analysis}/Alignment/recal_data_pre.table

# 2. Apply the recalibration to the BAM file (generating recalibrated BAM)
gatk ApplyBQSR -I ${Analysis}/Alignment/bench_chr7_dedup.bam \
    -R ${Genome} \
    --bqsr-recal-file ${Analysis}/Alignment/recal_data_pre.table \
    -O ${Analysis}/Alignment/bench_chr7_sorted_dedup_recal.bam

# 3. Build the post-recalibration model (generate recal_data1.table)
gatk BaseRecalibrator -I ${Analysis}/Alignment/bench_chr7_sorted_dedup_recal.bam \
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
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R ${Genome} \
    -I ${Analysis}/Alignment/bench_chr7_sorted_dedup_recal.bam \
    -O ${Analysis}/VariantCalling/bench_chr7.vcf.gz \
    --native-pair-hmm-threads 8
fi

# Annotate variants with functotator
gatk Funcotator \
    -R ${Genome} \
    -V ${Analysis}/VariantCalling/bench_chr7.vcf.gz \
    -O ${Analysis}/VariantCalling/bench_chr7_funcotated.vcf \
    --data-sources-path /Users/andreapassetti/Reference/hg38/functotator/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38 \
    --output-file-format VCF

# Ccompare vcf files
vcftools --vcf ${Analysis}/VariantCalling/bench_chr7_funcotated.vcf \
    --diff ${Data}/chr7_reads.mutations.vcf \
    --out ${Analysis}/VariantCalling/differences \
    --diff-site 


if false
then


### CONTROLLARE ###

# Extract fields from a VCF file to a tab-delimited table
gatk VariantsToTable \
    -V ${Analysis}/VariantCalling/bench_chr7_funcotated.vcf \
    -F CHROM -F POS -F REF -F ALT -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${Analysis}/VariantCalling/output_snps.table

# Extract the values from the VCF file
touch ${Analysis}/VariantCalling/funcotator_values.txt
cat ${Analysis}/VariantCalling/bench_chr7_funcotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${Analysis}/VariantCalling/snp_curated_variants.txt
#cat ${Analysis}/VariantCalling/bench_chr7_funcotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${folder}/Curated/indels_curated_variants.txt

#Extract annotated variant information
cat ${Analysis}/VariantCalling/output_snps.table | cut -f 5 | sed 's/|/\t/g' >> ${Analysis}/VariantCalling/snp_curated_variants.txt
#cat ${results}/output_indels.table | cut -f 5 | sed 's/|/\t/g' >> ${folder}/Curated/indels_curated_variants.txt

fi