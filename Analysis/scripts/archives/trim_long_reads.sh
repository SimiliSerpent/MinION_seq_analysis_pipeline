#!/bin/bash

module load porechop

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_raw_fastq_dir> <path_to_trimmed_fastq_dir>"
fi

# Input argument
PATH_TO_RAW=$1
PATH_TO_TRIMMED=$2

for filename in $PATH_TO_RAW*.fastq; do
	porechop \
	--verbosity 1 \
	--thread 30 \
	--end_size 150 \
	--min_trim_size 4 \
	--extra_end_trim 2 \
	--min_split_read_size 200 \
	-i "$filename" \
	-o ${PATH_TO_TRIMMED}$(basename "$filename" .fastq)_porechopped.fastq
done

