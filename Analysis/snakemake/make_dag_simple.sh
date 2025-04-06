module load fastqc
module load gcc/12.2.0
module load minimap2
module load python
module load samtools
module load snakemake

. ./set_config.sh

snakemake -np --rulegraph --forceall | dot -Tsvg > dag_simple.svg
