#!/bin/bash
#MSUB -r mapping                 # Requested job name
#MSUB -q normal                                 # Requested partition
#MSUB -N 1                                      # Max and min number of requested nodes
#MSUB -n 20                                     # Requested nb of tasks
#MSUB -T 86400                                  # Time limit (s)
#MSUB -o ./temp/other_logs/mapping.job%I.out   # Job standard output (%I is the job ID)
#MSUB -e ./temp/other_logs/mapping.job%I.err   # Job error output (%I is the job ID)
#MSUB -E "--mem 300G"

module load minimap2
module load samtools

EXP=25070704_MSUP_6
BC=03
# 1_seq_adapt 2_twist_outer 3_twist_primers 4_tso 5_polyA 6_GA
TRIM_STATE=4_tso
REF=SARSCoV2_human

# minimap2 \
# -ax map-ont \
# -t 20 \
# ./Data/references/snakeref/${REF}.fa \
# ./Data/seq_data/${EXP}/trimmed_fastq/barcode${BC}_porechopped_${TRIM_STATE}.fastq \
# > ./Data/seq_data/${EXP}/alignments/barcode${BC}_${TRIM_STATE}_map2_snakeref.sam && \
# samtools sort \
# -t 10 \
# ./Data/seq_data/${EXP}/alignments/barcode${BC}_${TRIM_STATE}_map2_snakeref.sam \
# -o ./Data/seq_data/${EXP}/alignments/barcode${BC}_${TRIM_STATE}_map2_snakeref.bam && \
# samtools index \
# ./Data/seq_data/${EXP}/alignments/barcode${BC}_${TRIM_STATE}_map2_snakeref.bam


### Map barcode raw reads (uncomment to use)

minimap2 \
-ax map-ont \
-t 20 \
./Data/references/snakeref/${REF}.fa \
./Data/seq_data/${EXP}/demultiplexed/barcode${BC}.fastq \
> ./Data/seq_data/${EXP}/alignments/barcode${BC}_raw_map2_snakeref.sam && \
samtools sort \
-t 10 \
./Data/seq_data/${EXP}/alignments/barcode${BC}_raw_map2_snakeref.sam \
-o ./Data/seq_data/${EXP}/alignments/barcode${BC}_raw_map2_snakeref.bam && \
samtools index \
./Data/seq_data/${EXP}/alignments/barcode${BC}_raw_map2_snakeref.bam

### Map unclassified reads (uncomment to use)

# minimap2 \
# -ax map-ont \
# -t 20 \
# ./Data/references/snakeref/${REF}.fa \
# ./Data/seq_data/${EXP}/trimmed_fastq/unclassified_porechopped_${TRIM_STATE}.fastq \
# > ./Data/seq_data/${EXP}/alignments/unclassified_${TRIM_STATE}_map2_snakeref.sam && \
# samtools sort \
# -t 10 \
# ./Data/seq_data/${EXP}/alignments/unclassified_${TRIM_STATE}_map2_snakeref.sam \
# -o ./Data/seq_data/${EXP}/alignments/unclassified_${TRIM_STATE}_map2_snakeref.bam && \
# samtools index \
# ./Data/seq_data/${EXP}/alignments/unclassified_${TRIM_STATE}_map2_snakeref.bam