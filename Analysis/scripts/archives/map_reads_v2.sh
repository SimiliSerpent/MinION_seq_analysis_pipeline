#!/bin/bash
#MSUB -r mapping			# Requested job name
#MSUB -q normal					# Requested partition
#MSUB -N 1					# Max and min number of requested nodes
#MSUB -n 20					# Requested nb of tasks
#MSUB -T 86400					# Time limit (s)
#MSUB -o ../${EXP_ID}/logs/mapping.job%I.out	# Job standard output (%I is the job ID)
#MSUB -e ../${EXP_ID}/logs/mapping.job%I.err	# Job error output (%I is the job ID)

module load minimap2
module load samtools

minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode01_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b01_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode02_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b02_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode03_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b03_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode04_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b04_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode05_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b05_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode06_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b06_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/barcode07_porechopped.fastq > ../../Data/${EXP_ID}/alignments/b07_porechopped_map2_multiv3.sam
minimap2 -ax map-ont -t 20 ../../Data/references/multirefv3.fa ../../Data/${EXP_ID}/trimmed_fastq/unclassified_porechopped.fastq > ../../Data/${EXP_ID}/alignments/unc_porechopped_map2_multiv3.sam
