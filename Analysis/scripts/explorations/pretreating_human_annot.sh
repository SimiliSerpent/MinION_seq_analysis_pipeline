#!/bin/bash
#MSUB -q normal                                 # Requested partition
#MSUB -N 1                                      # Number of requested nodes
#MSUB -n 10                            # Number of requested threads
#MSUB -T 86400                                  # Requested max time (seconds)
#MSUB -r cov_stats                     # Requested job name
#MSUB -o ${CAPASVIR_ANALYSIS}/temp/coverage_comparison/logs/cov.%I.out     # Output log file name (%I is the job ID)
#MSUB -e ${CAPASVIR_ANALYSIS}/temp/coverage_comparison/logs/cov.%I.err     # Error log file name (%I is the job ID)
#MSUB -E'--mem=100G --qos=long'

# # Get list of uniq accession ids
# grep -o '\bNC_[^[:space:]|\:]*' GCF_000001405.gff | uniq | less

# NC_000001.11 -> chr1
# NC_000002.12 -> chr2
# NC_000003.12 -> chr3
# NC_000004.12 -> chr4
# NC_000005.10 -> chr5
# NC_000006.12 -> chr6
# NC_000007.14 -> chr7
# NC_000008.11 -> chr8
# NC_000009.12 -> chr9
# NC_000010.11 -> chr10
# NC_000011.10 -> chr11
# NC_000012.12 -> chr12
# NC_000013.11 -> chr13
# NC_000014.9 -> chr14
# NC_000015.10 -> chr15
# NC_000016.10 -> chr16
# NC_000017.11 -> chr17
# NC_000018.10 -> chr18
# NC_000019.10 -> chr19
# NC_000020.11 -> chr20
# NC_000021.9 -> chr21
# NC_000022.11 -> chr22
# NC_000023.11 -> X
# NC_000024.10 -> Y
# NC_012920.1 -> Mitochondrion

# # Get names of regions in bam
# less ROI/human_regions.txt

# chr1
# chr2
# chr3
# chr4
# chr5
# chr6
# chr7
# chr8
# chr9
# chr10
# chr11
# chr12
# chr13
# chr14
# chr15
# chr16
# chr17
# chr18
# chr19
# chr20
# chr21
# chr22
# chrX
# chrY
# chrM

# # Replace region ids by used names
# sed -i 's/NC_000001.11/chr1/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000002.12/chr2/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000003.12/chr3/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000004.12/chr4/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000005.10/chr5/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000006.12/chr6/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000007.14/chr7/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000008.11/chr8/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000009.12/chr9/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000010.11/chr10/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000011.10/chr11/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000012.12/chr12/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000013.11/chr13/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000014.9/chr14/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000015.10/chr15/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000016.10/chr16/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000017.11/chr17/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000018.10/chr18/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000019.10/chr19/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000020.11/chr20/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000021.9/chr21/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000022.11/chr22/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000023.11/chrX/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_000024.10/chrY/g' Data/references/others/GCF_000001405_mod.gff && \
# sed -i 's/NC_012920.1/chrM/g' Data/references/others/GCF_000001405_mod.gff

# # Limit gff to existing regions
# for i in $(seq 1 22); do grep -P "^chr${i}\t" Data/references/others/GCF_000001405_mod.gff >> Data/references/others/GCF_000001405_mod_restricted.gff; done;
# grep -P "^chrX\t" Data/references/others/GCF_000001405_mod.gff >> Data/references/others/GCF_000001405_mod_restricted.gff
# grep -P "^chrY\t" Data/references/others/GCF_000001405_mod.gff >> Data/references/others/GCF_000001405_mod_restricted.gff
# grep -P "^chrM\t" Data/references/others/GCF_000001405_mod.gff >> Data/references/others/GCF_000001405_mod_restricted.gff

# # Load modules
module load python
module load pysam

# # Produce tsv
# python3 \
# -u Analysis/scripts/get_covered.py \
# -b Data/seq_data/25111906_DSUP1_1/alignments/barcode02_map2_snakeref_sorted.bam \
# -g Data/references/others/GCF_000001405_mod_restricted.gff \
# -o temp/25111906_temp/bc02_cov.tsv \
# --gene_biotype protein_coding

# ..or, for two bam files compared
python3 \
-u Analysis/scripts/get_covered.py \
-b Data/seq_data/25111906_DSUP1_1/alignments/barcode10_map2_snakeref_sorted.bam \
--bam2 Data/seq_data/25111906_DSUP1_1/alignments/barcode02_map2_snakeref_sorted.bam \
-g Data/references/others/GCF_000001405_mod_restricted.gff \
-o temp/25111906_temp/bc02bc10_cov_mRNA.tsv \
--target Data/references/sizes/human_sizes.txt \
--gene_biotype protein_coding \
-v 10

python3 \
-u Analysis/scripts/get_covered.py \
-b Data/seq_data/25111906_DSUP1_1/alignments/barcode10_map2_snakeref_sorted.bam \
--bam2 Data/seq_data/25111906_DSUP1_1/alignments/barcode02_map2_snakeref_sorted.bam \
-g Data/references/others/GCF_000001405_mod_restricted.gff \
-o temp/25111906_temp/bc02bc10_cov_rRNA.tsv \
--target Data/references/sizes/human_sizes.txt \
--gbkey rRNA \
-v 10