#!/bin/bash

#make folders
mkdir -p /Users/andreapassetti/GIAB/Analysis

folder="/Users/andreapassetti/GIAB/Analysis"
mkdir -p ${folder}/{Aligned,Results,Quality}

#directories
Aligned="${folder}/Aligned"
Results="${folder}/Results"
Quality="${folder}/Quality"
Data="/Users/andreapassetti/GIAB/Data"
Reference="/Users/andreapassetti/Reference"


if false
then

### MAPPING AND CLEANUP ###

#create index 
bwa index ${Reference}/hg38_chr7.fasta
samtools faidx ${Reference}/hg38_chr7.fasta

#map to reference genome using BWA-MEM
bwa mem -t 8 -R "@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA" ${Reference}/hg38_chr7.fasta ${Data}/NIST7035_TAAGGCGA_L002_R1_001_trimmed.fastq ${Data}/NIST7035_TAAGGCGA_L002_R2_001_trimmed.fastq > ${Aligned}/NA12878_chr7.sam

#create dictionary
gatk CreateSequenceDictionary -R ${Reference}/hg38_chr7.fasta

#mark duplicates (GATK)
gatk MarkDuplicatesSpark -I ${Aligned}/NA12878_chr7.sam -O ${Aligned}/NA12878_chr7_dedup_reads.bam -M ${Aligned}/NA12878_chr7_dedup_metrics.txt

#sort BAM file (Samtools)
#samtools sort -o ${Aligned}/NA12878_chr7_sorted_dedup_reads.bam ${Aligned}/NA12878_chr7_dedup_reads.bam

#sort BAM file (Picard)
picard SortSam -INPUT ${Aligned}/NA12878_chr7_dedup_reads.bam -OUTPUT ${Aligned}/NA12878_chr7_sorted_dedup_reads.bam -SORT_ORDER coordinate

#indel realignment
gatk BaseRecalibrator \
    -I ${Aligned}/NA12878_chr7_sorted_dedup_reads.bam \
    -R ${Reference}/hg38_chr7.fasta \
    --known-sites ${Reference}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${Reference}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${Data}/recal_data.table

#base recalibration > analysis ready reads
gatk ApplyBQSR -I ${Aligned}/NA12878_chr7_sorted_dedup_reads.bam -R ${Reference}/hg38_chr7.fasta --bqsr-recal-file ${Data}/recal_data.table -O ${Aligned}/NA12878_chr7_sorted_dedup_recal_reads.bam

fi

### VARIANT DISCOVERY ###

#variant calling
gatk --java-options "-Xmx8g" HaplotypeCaller -R ${Reference}/hg38_chr7.fasta -I ${Aligned}/NA12878_chr7_sorted_dedup_recal_reads.bam -O ${Results}/NA12878_chr7.vcf.gz --native-pair-hmm-threads 8

#joint genotyping - combine multiple samples

#variant recalibration separated by variant type > analysis ready variants


### VARIANT EVALUATION ###

#variant annotation and phasing

#variant evaluation