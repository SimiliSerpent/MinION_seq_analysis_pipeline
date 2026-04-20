#!/bin/bash
#MSUB -r ${EXPERIMENT_ID}_${BARCODE}				# Requested job name
#MSUB -o ${EXPERIMENT_ID}/commands_outputs/mapping.%I.out	# Job standard output (%I is the job ID)
#MSUB -e ${EXPERIMENT_ID}/commands_outputs/mapping.%I.err	# Job error output (%I is the job ID)
#MSUB -E "-c 16 -t 120 --mem 160G"				# Requested nb of CPUs, time limite and memory (Gb)

module load minimap2/2.28
module load samtools/1.18

mkdir -p ${EXPERIMENT_ID}/${BARCODE}/bam_files/
mkdir -p ${EXPERIMENT_ID}/${BARCODE}/analysis_results/

. map_reads.sh \
$REFERENCE \
all_reads_${EXPERIMENT_ID}_${BARCODE}_E_coli \
${EXPERIMENT_ID}/${BARCODE}/fastq_files/ \
${EXPERIMENT_ID}/${BARCODE}/bam_files/ \
${EXPERIMENT_ID}/${BARCODE}/analysis_results/ \
${BARCODE}

