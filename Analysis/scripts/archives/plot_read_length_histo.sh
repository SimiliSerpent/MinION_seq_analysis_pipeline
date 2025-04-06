#!/bin/bash

module load python

python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b01_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b02_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b03_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b04_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b05_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b06_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/b07_porechopped.read_lengths
python3 -u plot_read_length_histo.py -d ../EXP250124_1/read_lengths/unc_porechopped.read_lengths
