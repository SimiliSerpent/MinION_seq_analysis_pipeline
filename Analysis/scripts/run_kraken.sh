#!/bin/bash

module load kraken
module load python

# Print usage
usage() {
	printf "$(basename "$0") [-h] -f path_to_fastq_dir [-r path_to_report_files] [-k path_to_krona_files]

Classify reads in every fastq file in input directory using
Kraken2, and outputs the output, report and krona visualiza-
tion files to the specified output directory.

options:
    -h  show this help text
    -f  path to directory whith fastq files to classify
    -r  path to directory containing kraken reports (default: working directory)
    -k  path to directory containing krona htmls (default: working directory)
";
	exit;
}

# Set output directory to working directory
PATH_TO_REPORT="./"
PATH_TO_KRONA="./"

# Set abs path to kraken db
KRAKEN_DB="/env/cnrgh/proj/math_stats/scratch/vacus/Projet_CAPASVIR/scratch/Data/kraken2_db"

# Parse options
while getopts ':hf:r:k:-:' option; do
  case "$option" in
    h|help	) usage;;
    f|fastq_dir	) PATH_TO_FASTQ=$OPTARG;;
    r|report_dir) PATH_TO_REPORT=$OPTARG;;
    k|krona_dir	) PATH_TO_KRONA=$OPTARG;;
    :		) printf "missing argument for -%s\n" "$OPTARG" >&2; usage;;
    \?		) printf "illegal option: -%s\n" "$OPTARG" >&2; usage;;
  esac
done
shift $((OPTIND - 1))

if [ -z $PATH_TO_FASTQ ]; then
	printf "missing required option -f\n" >&2; usage;
fi

for filename in $PATH_TO_FASTQ/*.fastq; do
        kraken2 \
        --db $KRAKEN_DB \
        --output ${PATH_TO_REPORT}/$(basename "$filename" .fastq).k2_output \
        --report ${PATH_TO_REPORT}/$(basename "$filename" .fastq).k2_report \
        "$filename"
done

python3 build_krona.py \
-r ${PATH_TO_REPORT}/ \
-o ${PATH_TO_KRONA}/

