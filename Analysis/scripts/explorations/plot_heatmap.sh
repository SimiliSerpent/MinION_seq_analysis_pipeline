#!/bin/bash
#MSUB -q normal                                 # Requested partition
#MSUB -N 1                                      # Number of requested nodes
#MSUB -n 5                            # Number of requested threads
#MSUB -r temp_job                     # Requested job name
#MSUB -o ./temp.%I.out     # Output log file name (%I is the job ID)
#MSUB -e ./temp.%I.err     # Error log file name (%I is the job ID)
#MSUB -E '-t 700 --mem=50G --qos=default'

module load python

python3 ../scripts/plot_sequence_stats_heatmap.py -p 25102205_DSUP1_1/read_cleaning/pickles -o ../../temp/25102205_temp/25102205_DSUP1_1_size_vs_qual_heatmap.png -s ../../Data/seq_data/25102205_DSUP1_1 -v 5
python3 ../scripts/plot_sequence_stats_heatmap.py -p 25102205_DSUP2_1/read_cleaning/pickles -o ../../temp/25102205_temp/25102205_DSUP2_1_size_vs_qual_heatmap.png -s ../../Data/seq_data/25102205_DSUP2_1 -v 5
python3 ../scripts/plot_sequence_stats_heatmap.py -p 25111906_DSUP1_1/read_cleaning/pickles -o ../../temp/25111906_temp/25111906_DSUP1_1_size_vs_qual_heatmap.png -s ../../Data/seq_data/25111906_DSUP1_1 -v 5
python3 ../scripts/plot_sequence_stats_heatmap.py -p 25111906_DSUP2_1/read_cleaning/pickles -o ../../temp/25111906_temp/25111906_DSUP2_1_size_vs_qual_heatmap.png -s ../../Data/seq_data/25111906_DSUP2_1 -v 5
