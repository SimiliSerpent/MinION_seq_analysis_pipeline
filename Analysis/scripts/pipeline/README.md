# Scripts used for nanopore data analysis pipeline 
This directory contains scripts used in the analysis of Nanopore data generated
in the PhD Project "Pathogen detection in environmental samples using
nucleotide sequencing: Nanopore long-reads sequencing technology for accurate
strain identification".

For most scripts a help has been devised and should be used before running the
script.

---

## deduplicate_umis.py

From a BAM file and a UMI information json, perform UMI deduplication.

Outputs the following files:

  - **reads_stats.json**

Dictionnary containing reads information as follows:

```
{
    "63247054-4dc3-4399-b096-29a12f3d404c": {
        "aligned_len": 958,
        "ref": "NC-045512.2",
        "ref_start": 10174,
        "orientation": "forward"
    },
    "c49a7bb2-abb3-439d-8e16-9d7f846d32ad": {
        "aligned_len": 962,
        "ref": "NC-045512.2",
        "ref_start": 10174,
        "orientation": "forward"
    }, 
    ...
}
```

  - **alignment_heatmap.png**

Heatmap matrix of alignment-based distance between pairs of UMIs.

  - **overlap_heatmap.png**

Heatmap matrix of coverage-overlap-based distance between pairs of UMIs.

  - **umis_stats.tsv**

Tsv file containing information on every identified "raw" UMIs. One raw UMI
corresponds to one contiguous group of reads (in terms of reference coverage)
sharing the exact same UMI sequence at first glance.
Example:

| raw_umi | ref | nb_raw_reads | nb_final_reads | final_umi |
| --- | --- | --- | --- | --- |
| AAAACTCG_chr21 | chr21 | 1 | 2 | ATAAACCG_chr21_1 |
| AAACACGG_chr21 | chr21 | 1 | 434 | TCGATCTG_chr21_2 |
| AAACAGGG_chr2 | chr2 | 1 | 1 | AAACAGGG_chr2 |
| AAACATTC_NC-045512.2 | NC-045512.2 | 1 | 3 | ACAATCTT_NC-045512.2 |
| AAACCGCC_chr21 | chr21 | 1 | 489 | TCGATCTG_chr21_3 |
| ... | ... | ... | ... | ... |

  - **groups.tsv**

Informations on UMI assignment per read. For example:

| read_id | contig | position | gene | umi | umi_count | final_umi | final_umi_count | unique_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 63247054-4dc3-4399-b096-29a12f3d404c | NC-045512.2 | 10174 | NA | GCAGAGGT_NC-045512.2 | 4880 | GCAGAGGT_NC-045512.2 | 5156 | GCAGAGGT_NC-045512.2 |
| c49a7bb2-abb3-439d-8e16-9d7f846d32ad | NC-045512.2 | 10174 | NA | GCAGAGGT_NC-045512.2 | 4880 | GCAGAGGT_NC-045512.2 | 5156 | GCAGAGGT_NC-045512.2 |
| b7c17fd4-daa5-479b-9920-e2c6f5f240c7 | NC-045512.2 | 10174 | NA | GCAGAGGT_NC-045512.2 | 4880 | GCAGAGGT_NC-045512.2 | 5156 | GCAGAGGT_NC-045512.2 |
| d1832d40-7808-405a-9f45-d3b176a413ed | chrM | 8364 | NA | GGGAGACC_chrM | 1 | GGGAGACC_chrM | 1 | GGGAGACC_chrM |
| 46b7c342-085e-47c6-9845-e3b7abdd2032 | NC-045512.2 | 10175 | NA | GAACAACT_NC-045512.2 | 1 | GCAGAGGT_NC-045512.2 | 5156 | GCAGAGGT_NC-045512.2 |
| db60a2dd-0aa5-4bb2-bfde-fdf2d48a0585 | NC-045512.2 | 10175 | NA | GAACAGGG_NC-045512.2 | 1 | GCAGAGGT_NC-045512.2 | 5156 | GCAGAGGT_NC-045512.2 |
| d7cb8354-c1ed-448e-8497-560f713b2feb | chr21 | 8212570 | NA | GAACAGGG_chr21_1 | 1 | TTAATGGG_chr21 | 4 | TTAATGGG_chr21 |
| c0102a23-2f5d-487a-bb4d-d5e9057a8a68_1 | chr21 | 8256779 | NA | GAACAGGG_chr21_2 | 1 | GAACAGGG_chr21_2 | 1 | GAACAGGG_chr21_2 |

  - **dedup.fasta**

Fasta file with deduplicated reads (UMIs trimmed)

Optionally (if errors occur):
  - **deduplication_error_log.txt**
  - **minimap2_log.txt**
  - **racon_log.txt**


