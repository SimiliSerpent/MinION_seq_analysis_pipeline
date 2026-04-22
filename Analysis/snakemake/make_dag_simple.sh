module load fastqc
module load gcc
module load minimap2
module load python
module load samtools
module load snakemake

. ./set_config.sh

snakemake -np --rulegraph --forceall --snakefile Snakefile.py | grep -v 'Building DAG of jobs' | dot -Tsvg > dag_simple.svg
