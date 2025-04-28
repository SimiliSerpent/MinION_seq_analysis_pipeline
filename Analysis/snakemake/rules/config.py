#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Config file for bioinfo analysis snakemake pipeline.
Imports useful python modules, defines variables, and computes the names of
expected files in order to feed the 'all' rule of the Snakefile.
"""


import importlib.util
import os
import sys

import argparse
import json
import pandas as pd


PROJECT = '/env/cng/scratchE/q_projet_CAPASVIR_880' # path to scratch


# Load utils module using abs path
# (yes, it a pain to import modules from a different directory in python)
spec = importlib.util.spec_from_file_location(
    'utils',
    PROJECT + "/Analysis/scripts/utils.py"
)
utils = importlib.util.module_from_spec(spec)
sys.modules['utils'] = utils
spec.loader.exec_module(utils)


# Get list of target species
species_list = os.environ['SPECIES'].split(',')
species_list.sort()
reference = '_'.join(species_list) + '.fa'

# Get number of threads
nthreads = int(os.environ['NTHREADS'])

# Get files name for output
exp_id = os.environ['EXP_ID']
data_path = PROJECT + '/Data/' + exp_id
analysis_path = PROJECT + '/Analysis/' + exp_id
raw_fastqs = utils.find_file(data_path + '/raw_fastq', '.fastq')
sample_names = [fastq.rstrip('.fastq') for fastq in raw_fastqs]
sample_names.sort()

# List expected output files
out_files = [analysis_path + '/fastqc/' + name + '_porechopped_ONT_TSO_tail_GA_fastqc.html' for name in sample_names] # fastqc
out_files += [analysis_path + '/read_lengths/' + name + '_read_lengths_histo.png' for name in sample_names] # read lengths histogram
out_files += [analysis_path + '/nuc_freq/' + exp_id + '_nuc_freq.png'] # frequency of nucleotides
out_files += [analysis_path + '/' + exp_id + '_samples_composition.png'] # samples composition in targets
out_files += [analysis_path + '/' + exp_id + '_normalized_samples_composition.png'] # samples composition (targets normalized)
out_files += [analysis_path + '/' + exp_id + '_read_size_quality_with_trimming.png'] # evolution of read size vs. quality at each trimming step
