#!/bin/bash

# Per-analysis config. Edit before each run.

# Check for the environment variable $CAPASVIR_ANALYSIS
# It should be defined and point towards the directory containg the MinION
# data sequencing analysis repository. It can be defined in your shell profile.
: "${CAPASVIR_ANALYSIS:?Error: Environment variable 'CAPASVIR_ANALYSIS' must \
be defined and point towards the MinION analysis repository.}"
echo "CAPASVIR_ANALYSIS is set to: $CAPASVIR_ANALYSIS"

export EXP_ID="26042009_DSUP1_1"
export NTHREADS=32
export SPECIES="SARSCoV2" # human | mouse | SARSCoV2 (coma-separated)
export BARCODES="TWIST-LRLP-SHv2E" # UNBARCODED | ONT-EXP-PBC001 | TWIST-LRLP-SHv2E

echo "        Using experiment ID: $EXP_ID"
echo "              Nb of threads: $NTHREADS"
echo "             Target species: $SPECIES"
echo "                  Barcoding: $BARCODES"

echo "Creating analysis results directory..."
mkdir -p ../results/${EXP_ID}/logs
