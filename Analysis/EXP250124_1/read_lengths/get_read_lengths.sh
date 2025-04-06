#!/bin/bash

cat ../alignments_stats/b01_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b01_porechopped.read_lengths
cat ../alignments_stats/b02_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b02_porechopped.read_lengths
cat ../alignments_stats/b03_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b03_porechopped.read_lengths
cat ../alignments_stats/b04_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b04_porechopped.read_lengths
cat ../alignments_stats/b05_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b05_porechopped.read_lengths
cat ../alignments_stats/b06_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b06_porechopped.read_lengths
cat ../alignments_stats/b07_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > b07_porechopped.read_lengths
cat ../alignments_stats/unc_porechopped_map2_multiv3.samstats | grep ^RL | cut -f 2- > unc_porechopped.read_lengths
