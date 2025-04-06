#!/bin/bash
#MSUB -r CAPASVIR_assembly                 # Requested job name
#MSUB -q normal                                 # Requested partition
#MSUB -N 1                                      # Max and min number of requested nodes
#MSUB -n 30                                     # Requested nb of tasks
#MSUB -T 86400                                  # Time limit (s)
#MSUB -o ./assembly/spades.job%I.out   # Job standard output (%I is the job ID)
#MSUB -e ./assembly/spades.job%I.err   # Job error output (%I is the job ID)
#MSUB -E "--mem 300G"

module load spades

spades --meta -t 30 -s ./assembly/test_ben_1.fastq -m 298 -k 21,39,59,79,99,119,127 -o assembly/
