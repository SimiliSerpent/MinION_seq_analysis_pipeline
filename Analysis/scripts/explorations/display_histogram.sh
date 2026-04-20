#!/bin/bash

module load samtools

for BARCODE in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07
do
  cat ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n  > ${EXPERIMENT_ID}/${BARCODE}/fastq_files/size_of_all_read_${BARCODE}.txt
  Rscript display_histogram.R ${EXPERIMENT_ID}/${BARCODE}/fastq_files/size_of_all_read_${BARCODE}.txt ${BARCODE} FASTQ
  samtools view ${EXPERIMENT_ID}/${BARCODE}/bam_files/aligned_all_reads_${EXPERIMENT_ID}_${BARCODE}.bam | awk '$5 > 33 || $1 ~ /^@/' | awk '{print length($10)}' > ${EXPERIMENT_ID}/${BARCODE}/bam_files/size_of_all_read_${BARCODE}.txt
  Rscript display_histogram.R ${EXPERIMENT_ID}/${BARCODE}/bam_files/size_of_all_read_${BARCODE}.txt ${BARCODE} BAM_SARS_COV_2
done
