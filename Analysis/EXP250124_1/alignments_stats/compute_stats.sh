#!/bin/bash

module load samtools

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam > b01_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam > b02_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam > b03_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam > b04_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam > b05_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam > b06_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam > b07_porechopped_map2_multiv3.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam > unc_porechopped_map2_multiv3.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b01_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b02_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b03_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b04_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b05_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b06_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > b07_porechopped_map2_multiv3.human.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/hs38dh_human_sizes.txt > unc_porechopped_map2_multiv3.human.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b01_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b02_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b03_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b04_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b05_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b06_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > b07_porechopped_map2_multiv3.mouse.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/grcm39_mouse_sizes.txt > unc_porechopped_map2_multiv3.mouse.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b01_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b02_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b03_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b04_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b05_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b06_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > b07_porechopped_map2_multiv3.lambda.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/lambda_sizes.txt > unc_porechopped_map2_multiv3.lambda.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b01_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b02_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b03_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b04_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b05_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b06_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > b07_porechopped_map2_multiv3.SARSCoV2.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/NC_045512_sarscov2_sizes.txt > unc_porechopped_map2_multiv3.SARSCoV2.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b01_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b02_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b03_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b04_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b05_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b06_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > b07_porechopped_map2_multiv3.Ecoli.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/GCF_000005845_Ecoli_sizes.txt > unc_porechopped_map2_multiv3.Ecoli.samstats

samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b01_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b01_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b02_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b02_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b03_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b03_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b04_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b04_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b05_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b05_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b06_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b06_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/b07_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > b07_porechopped_map2_multiv3.Saureus.samstats
samtools stats --threads 20 ../../../Data/EXP250124_1/alignments/unc_porechopped_map2_multiv3_sorted.sam -t ../../../Data/references/S_aureus_sizes.txt > unc_porechopped_map2_multiv3.Saureus.samstats
