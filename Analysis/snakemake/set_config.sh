#!/bin/bash

# Per-analysis config. Edit before each run.

export EXP_ID="25070704_DSUP2_1"
export NTHREADS=32
export SPECIES="human,SARSCoV2" # human | mouse | SARSCoV2 (coma-separated)
export BARCODES="TWIST-LRLP-SHv2E" # UNBARCODED | ONT-EXP-PBC001 | TWIST-LRLP-SHv2E
export DEDUP_UMIS="True" # true | false
export UMI_PATTERN="5NNNNNNNN"
export UMI_SPLIT_METHOD="density_peaks" # leiden | density_peaks

# Check for the environment variable $CAPASVIR_ANALYSIS
# It should be defined and point towards the directory containg the MinION
# data sequencing analysis repository. It can be defined in your shell profile.
: "${CAPASVIR_ANALYSIS:?Error: Environment variable 'CAPASVIR_ANALYSIS' must \
be defined and point towards the MinION analysis repository.}"
echo "CAPASVIR_ANALYSIS is set to: $CAPASVIR_ANALYSIS"

echo "        Using experiment ID: $EXP_ID"
echo "              Nb of threads: $NTHREADS"
echo "             Target species: $SPECIES"
echo "                  Barcoding: $BARCODES"

if [[ "${DEDUP_UMIS,,}" == "true" ]]; then
    echo "    Dedup UMIs with pattern: $UMI_PATTERN"
    echo "   Split UMIs clusters with: $UMI_SPLIT_METHOD"
fi

echo "Creating analysis results directory..."
mkdir -p ../results/${EXP_ID}/logs
