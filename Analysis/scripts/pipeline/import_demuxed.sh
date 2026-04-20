#!/bin/bash

# Place yourself in the sequencing data directory
MINION_PC_EXP_LOCATION="26_04_02_08_FBF64370_2/FBF64370/20260402_1508_MN19813_FBF64370_617e1205"
INTI_EXP_NAME="26040208_DSUP2_1"
KIT="TWIST-ALL"

mkdir -p ./$INTI_EXP_NAME/demultiplexed

# Import the basecalled reads
scp -r capasvir@I0017871.illumina.cng.fr:/data/$MINION_PC_EXP_LOCATION/fastq_pass/* \
    ./$INTI_EXP_NAME/demultiplexed/.

# Move read files from their directories to demultiplexed directory
cd ./$INTI_EXP_NAME/demultiplexed
for i in $(seq -w 1 96); do
    if $(test -d barcode$i); then
        mv barcode$i/*.fastq.gz barcode$i.fastq.gz
    fi
done
if $(test -d unclassified); then
    mv unclassified/*.fastq.gz unclassified.fastq.gz
fi

# Remove empty directories
for i in $(seq -w 1 96); do
    if $(test -d barcode$i); then
        rm -rf barcode$i
    fi
done
if $(test -d unclassified); then
    rm -rf unclassified
fi

# Unzip fastq files
gunzip *

# Generate raw sequence text file
# python3 -u ../../../Analysis/scripts/get_raw_sequence_text.py
