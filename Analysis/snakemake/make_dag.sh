module load fastqc
module load gcc
module load minimap2
module load python
module load samtools
module load snakemake

. ./set_config.sh

snakemake -np --filegraph --forceall --snakefile Snakefile.py | dot -Tsvg > dag.svg
