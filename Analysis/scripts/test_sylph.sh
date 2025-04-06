#!/bin/bash
#MSUB -E "-c 16 -t 120 --mem 160G"                              # Requested nb of CPUs, time limite and memory (Gb)
~/.cargo/bin/sylph profile gtdb-r220-c200-dbv1.syldb 2025_01_24_exp_01/barcode03/fastq_files/barcode03.fastq > profiling.tsv

