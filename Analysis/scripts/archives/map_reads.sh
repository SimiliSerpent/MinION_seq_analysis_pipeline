#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <path_to_reference_genome.fa> <file_prefixe> <path_to_fastq_files> <path_to_bam_files> <path_to_stats>"
fi

# Input arguments
REFERENCE=$1
PREFIXE=$2
PATH_TO_FASTQ=$3
PATH_TO_BAM=$4
PATH_TO_STATS=$5
PREFIXE_FASTQ=$6

# Output files
SAM_FILE="aligned_${PREFIXE}.sam"
BAM_FILE="aligned_${PREFIXE}.bam"
SORTED_BAM="aligned_${PREFIXE}_sorted.bam"
STATS_FILE="stats_${PREFIXE}.txt"

# Step 1: Align reads to the reference genome
echo "Aligning reads to the reference genome..."
minimap2 -ax map-ont $REFERENCE ${PATH_TO_FASTQ}${PREFIXE_FASTQ}.fastq > ${PATH_TO_BAM}$SAM_FILE

# Step 2: Convert SAM to BAM and sort
echo "Converting SAM to BAM and sorting..."
samtools view -bS ${PATH_TO_BAM}$SAM_FILE > ${PATH_TO_BAM}$BAM_FILE
samtools sort ${PATH_TO_BAM}$BAM_FILE -o ${PATH_TO_BAM}$SORTED_BAM
samtools index ${PATH_TO_BAM}$SORTED_BAM

# Step 3: Generate statistics and extract number of mapped bases
echo "Generating alignment statistics..."
samtools stats ${PATH_TO_BAM}$SORTED_BAM > ${PATH_TO_STATS}$STATS_FILE
MAPPED_BASES=$(grep "^SN" ${PATH_TO_STATS}$STATS_FILE | grep "bases mapped (cigar):" | cut -f 3)

echo "Total bases mapped to the reference genome: $MAPPED_BASES"
