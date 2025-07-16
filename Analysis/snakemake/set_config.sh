#!/bin/bash

export EXP_ID="25032702_MSUP_1"
export NTHREADS=32
export SPECIES="human,mouse"
export BARCODES="ONT-EXP-PBC001" # UNBARCODED | ONT-EXP-PBC001 | TWIST-LRLP-SHv2E
export SAMPLES="6=Human 100pg,7=Human 100pg sharper cutoff,8=Human 100pg 15' polyA,9=Human 100pg low ATP,10=H2O,11=Mouse 100pg,12=H2O"

mkdir -p ../results/${EXP_ID}/logs
