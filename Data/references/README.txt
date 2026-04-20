In order to run the snakemake workflow, provide the species fasta in the public
directory, and the list of regions of interest in the ROI directory, using the
same base name and adding the suffix '_regions.txt'.
The directories work as described in the following sections.

===============================================================================
public - contains public reference sequences

Most references (except SARSCoV2) are named with respect to the following
conventions:
Gspecies.fa   where 'G' is the first letter of the Genus and species is the species name
name.fa       if usage name is not the scientific name

Current names refer to the following public sequences:
human     hs38dh
mouse     grcm39_refseq
SARSCoV2  NC_045512.2
lambda    ENA_J02459
Ecoli     GCF_000005845
Saureus   GCF_000013425

===============================================================================
ROI - contains regions of interest of public sequences

Each <fasta_name>_regions.txt ROI file corresponds to one <fasta_name>.fa public
sequence. It contains the list of regions (=chromosome=contigs) of interest of
that sequence, i.e. often excludes the alternative/unplaced regions and sub-reg-
ions. The format is one region name by line, e.g.:

<start-of-file>
chr1
chr2
MT
<end-of-file>

===============================================================================
sizes - containts the sizes of regions of interest

Each <fasta_name>_sizes.txt sizes file corresponds to one <fasta_name>.fa public
sequence. It contains the list of regions of interest along with their size, us-
ing the following formatting:
<region_name>\t<region_start>\t<region_stop>

For instance:

<start-of-file>
chr1    1     124658
chr2    1     546877
MT      1     46524
<end-of-file>

===============================================================================
snakeref - contains the references used in the snakemake anlysis workflow

===============================================================================
others - contains old stuff and anything that does not fit in the directories d-
escribed previously

===============================================================================
temp - can be deleted at any moment
