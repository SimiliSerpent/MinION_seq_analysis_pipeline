#!/bin/bash

# Place yourself where you want the files to be imported
# (like in Data/seq_data)
export MINION_PC_EXP_NAME="25_07_07_04_FBA38151"
export INTI_EXP_NAME="25070704_MSUP_1"

mkdir -p ./$INTI_EXP_NAME/raw_fastq
mkdir -p ./$INTI_EXP_NAME/run_reports
cd ./$INTI_EXP_NAME

# Import the basecalled reads
scp -r capasvir@I0017871.illumina.cng.fr:/data/$MINION_PC_EXP_NAME ./minknow_output

# Find reads directory
fastq_pass=$(find . -type d -name "fastq_pass" -print -quit)
fastq_fail=$(find . -type d -name "fastq_fail" -print -quit)

# Group reads
cat $fastq_pass/* > ./raw_fastq/fastq_pass.fastq.gz
cat $fastq_fail/* > ./raw_fastq/fastq_fail.fastq.gz
cat ./raw_fastq/* > ./raw_fastq/all_reads.fastq.gz

# Unzip fastq files
gunzip ./raw_fastq/*

# Extract summary and report
seq_sum=$(find . -type f -name "sequencing_summary_*" -print -quit)
mv $seq_sum ./minknow_sequencing_summary.txt
reports=$(find . -type f -name "report_*" )
for report in $reports; do
    mv $report ./run_reports/.
done
pore_act=$(find . -type f -name "pore_activity_*" -print -quit)
mv $pore_act ./run_reports/.

# Remove directories
rm -rf ./minknow_output
