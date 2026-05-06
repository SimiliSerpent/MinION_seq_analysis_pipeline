#!/bin/bash

## Import fastq.gz files from distant sequencing device and organize them in a
# way they can be used immediately with the analysis pipeline.
# The sequencing experiment has to comprise barcodes.

# Define the following variables to target the right sequencing data
PROT_GROUP_ID="26_04_20_09_FBF75335" # name given to the exp. in MinKNOW
SAMPLE_ID="FBF75335" # name given to the sample in MinKNOW
START_TIME="20260420_1019" # sequencing starting time
MINION_ID="MN19813" # MinION ID
FLOWCELL_ID="FBF75335" # flowcell ID
SHORT_PROT_GROUP_ID="7044ce5a" # first 8 characters of the protocol_run_id
EXP_NAME="26042009_DSUP1_1" # name given to this analysis

# Check for the required environment variables.
# These can be defined in your shell profile.

: "${CAPASVIR_ANALYSIS:?Error: Environment variable 'CAPASVIR_ANALYSIS' must \
be defined and point towards the MinION analysis repository.}"
echo "CAPASVIR_ANALYSIS is set to: $CAPASVIR_ANALYSIS"

: "${SEQ_USER:?Error: Environment variable 'SEQ_USER' must \
be defined and contain the username to access sequencing data.}"
echo "         SEQ_USER is set to: $SEQ_USER"

: "${SEQ_DEVICE:?Error: Environment variable 'SEQ_DEVICE' must \
be defined and contain the ip adress of the sequencing device.}"
echo "       SEQ_DEVICE is set to: $SEQ_DEVICE"

: "${SEQ_LOCATION:?Error: Environment variable 'SEQ_LOCATION' must \
be defined and contain the path to sequencing data on the sequencing device.}"
echo "     SEQ_LOCATION is set to: $SEQ_LOCATION"

# Build a variable containing the full path to experiments data on seq device
DISTANT_LOC="$SEQ_LOCATION/$PROT_GROUP_ID/$SAMPLE_ID/${START_TIME}_${MINION_ID}"
DISTANT_LOC="${DISTANT_LOC}_${FLOWCELL_ID}_${SHORT_PROT_GROUP_ID}"

# Create the experiment data directory on the cluster
mkdir -p $CAPASVIR_ANALYSIS/Data/seq_data/$EXP_NAME/demultiplexed
cd $CAPASVIR_ANALYSIS/Data/seq_data/$EXP_NAME/demultiplexed

# Import the basecalled reads
scp -r $SEQ_USER@$SEQ_DEVICE:$DISTANT_LOC/fastq_pass/* .

# Move fastq files from their directories to the "demultiplexed" directory
# ...in the case of unbarcoded experiment:
if $(test -f ./$(ls | grep fastq.gz | head -n 1)); then
    mv ./*.fastq.gz ./barcode01.fastq.gz
fi
# ...for barcoded experiment:
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

# To generate raw sequence text file:
# python3 -u ../../../Analysis/scripts/get_raw_sequence_text.py
