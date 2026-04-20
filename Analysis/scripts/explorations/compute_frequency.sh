#!/bin/bash

for BARCODE in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07
do
  frequence_A=$(grep -v '>' ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | grep -o "[A]" | wc -l)
  frequence_T=$(grep -v '>' ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | grep -o "[T]" | wc -l)
  frequence_C=$(grep -v '>' ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | grep -o "[C]" | wc -l)
  frequence_G=$(grep -v '>' ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | grep -o "[G]" | wc -l)
  frequence_ATCG=$(grep -v '>' ${EXPERIMENT_ID}/${BARCODE}/fastq_files/${BARCODE}.fastq | grep -o "[ATCG]" | wc -l)
  echo ${BARCODE}
  echo "Frequency of A:"
  echo "100*$frequence_A/$frequence_ATCG" | bc -l
  echo "Frequency of T:"
  echo "100*$frequence_T/$frequence_ATCG" | bc -l
  echo "Frequency of C:"
  echo "100*$frequence_C/$frequence_ATCG" | bc -l
  echo "Frequency of G:"
  echo "100*$frequence_G/$frequence_ATCG" | bc -l
  echo "Total number of bases:"
  echo ${frequence_ATCG}
done
