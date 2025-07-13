#!/bin/bash
#MSUB -q normal					# Requested partition
#MSUB -N 1					# Number of requested nodes
#MSUB -n ${NTHREADS}				# Number of requested threads
#MSUB -T 86400					# Requested max time (seconds)
#MSUB -r ${EXP_ID}_analysis			# Requested job name
#MSUB -o ../results/${EXP_ID}/logs/snakemake.%I.out	# Output log file name (%I is the job ID)
#MSUB -e ../results/${EXP_ID}/logs/snakemake.%I.err	# Error log file name (%I is the job ID)
#MSUB -E '--mem=300G --qos=default'

module load fastqc
module load gcc/12.2.0
module load minimap2
# module load porechop
module load python
module load samtools
module load snakemake

snakemake \
--cores ${NTHREADS} \
--resources mem_mb=${NTHREADS}000 \
--latency-wait 15 \
--rerun-incomplete \
--keep-going \
--nolock \
--rerun-triggers mtime \
--snakefile Snakefile.py

# Accounting for:

# Available cores
# Available amount of memory
# Wait this amount of second after job execution for file to be produced
# Re-run all jobs the output of which is recognized as incomplete
# Go on with independent jobs if a job fails
# Do not lock the working directory
# Define what triggers the rerunning of a job (here: file modification dates)

# Do not hesitate to use filter_escape_seq.py on outlog like this:
# module load python && python ../scripts/filter_escape_seq.py -f ../results/${EXP_ID}/logs/<outlog_name>
# This allows better visualization of the outlog using the less/more cmds