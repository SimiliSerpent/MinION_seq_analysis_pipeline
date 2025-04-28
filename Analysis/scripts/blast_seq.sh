#!/bin/bash
#MSUB -q normal                                 # Requested partition
#MSUB -N 1                                      # Number of requested nodes
#MSUB -n 30                            # Number of requested threads
#MSUB -T 86400                                  # Requested max time (seconds)
#MSUB -r blastn                     # Requested job name
#MSUB -o /env/export/cng_n01_scratch/v_scratchE/q_projet_CAPASVIR_880/temp/blasting_unexpected_reads/logs/blastn.%I.out     # Output log file name (%I is the job ID)
#MSUB -e /env/export/cng_n01_scratch/v_scratchE/q_projet_CAPASVIR_880/temp/blasting_unexpected_reads/logs/blastn.%I.err     # Error log file name (%I is the job ID)
#MSUB -E'--mem=300G --qos=long'

module load blast+/2.16.0
test_dir="/env/export/cng_n01_scratch/v_scratchE/q_projet_CAPASVIR_880/temp/blasting_unexpected_reads"
my_seq="test_barcode05_unmapped_selection_no_GA"

# Set path to Blast nt local DB

BLASTDB="/env/export/cng_n01_scratch/v_scratchE/q_projet_CAPASVIR_880/Data/blast_db/nt"
export BLASTDB="/env/export/cng_n01_scratch/v_scratchE/q_projet_CAPASVIR_880/Data/blast_db/nt"

# See options at https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_

blastn \
-task blastn \
-outfmt 5 \
-max_hsps 5 \
-max_target_seqs 300 \
-num_threads 30 \
-db nt \
-query ${test_dir}/${my_seq}.fasta \
-out ${test_dir}/${my_seq}_blast.out
