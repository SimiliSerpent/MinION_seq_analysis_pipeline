#!/bin/bash
#MSUB -r ${EXP_ID}_porechopping			# Requested job name
#MSUB -q normal					# Requested partition
#MSUB -N 1					# Max and min number of requested nodes
#MSUB -n 30					# Requested nb of tasks
#MSUB -T 86400					# Time limit (s)
#MSUB -o ../${EXP_ID}/logs/porechop.job%I.out	# Job standard output (%I is the job ID)
#MSUB -e ../${EXP_ID}/logs/porechop.job%I.err	# Job error output (%I is the job ID)
#MSUB -E "--mem 300G"

. ./trim_long_reads.sh ../../Data/${EXP_ID}/temp/ ../../Data/${EXP_ID}/trimmed_fastq/
