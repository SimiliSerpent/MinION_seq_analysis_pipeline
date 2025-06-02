#!/bin/bash

# Place yourself where you want the files to be imported
export MINION_PC_EXP_NAME="25_01_24_01_250510_1"
export INTI_EXP_NAME="250124_1_250510_1_FAST"

# Import the basecalled reads
scp -r capasvir@I0017871.illumina.cng.fr:/home/capasvir/basecalling/$MINION_PC_EXP_NAME/basecalling_output ./$INTI_EXP_NAME

# Group reads
cd $INTI_EXP_NAME
for i in $(seq 1 9); do cat barcode0$i/* > barcode0$i.fastq.gz; done && for i in $(seq 10 12); do cat barcode$i/* > barcode$i.fastq.gz; done && cat unclassified/* > unclassified.fastq.gz

# Remove directories
for i in $(seq 1 9); do rm -rf barcode0$i/; done && for i in $(seq 10 12); do rm -rf barcode$i/; done && rm -rf unclassified/

# Unzip fastq files
for i in $(seq 1 9); do gunzip barcode0$i.fastq.gz; done && for i in $(seq 10 12); do gunzip barcode$i.fastq.gz; done && gunzip unclassified.fastq.gz

# Generate raw sequence text file
python3 -u ../../../Analysis/scripts/get_raw_sequence_text.py