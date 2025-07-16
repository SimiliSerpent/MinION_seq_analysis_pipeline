#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Config file for bioinfo analysis snakemake pipeline.
Imports useful python modules, defines variables, and computes the names of
expected files in order to feed the 'all' rule of the Snakefile.
"""


import datetime
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

# Get type of barcodes
barcode_type = str(os.environ['BARCODES'])
if not barcode_type in ['UNBARCODED', 'ONT-EXP-PBC001', 'TWIST-LRLP-SHv2E']:
    raise ValueError(
        'The selected barcodes type does is not allowed:', barcode_type
    )

# Get list of samples
samples_str = '"' + str(os.environ['SAMPLES']) + '"'

# Get files name for output
exp_id = os.environ['EXP_ID']
data_path = PROJECT + '/Data/seq_data/' + exp_id
analysis_path = PROJECT + '/Analysis/results/' + exp_id
raw_fastqs = utils.find_file(data_path + '/raw_fastq', '.fastq')
barcode_names = [fastq.rstrip('.fastq') for fastq in raw_fastqs]
barcode_names.sort()

# List expected output files
out_files = [
    analysis_path \
    + '/fastqc/' \
    + name \
    + '_porechopped_ONT_TSO_tail_GA_fastqc.html' for name in barcode_names
] # fastqc
out_files += [
    analysis_path \
    + '/read_lengths/' \
    + name \
    + '_read_lengths_histo.png' for name in barcode_names
] # read lengths histogram
out_files += [
    analysis_path + '/nuc_freq/' + exp_id + '_nuc_freq.png'
] # frequency of nucleotides
out_files += [
    analysis_path + '/' + exp_id + '_barcodes_composition.png'
] # samples composition in targets
out_files += [
    analysis_path + '/' + exp_id + '_barcodes_normalized_composition.png'
] # samples composition (targets normalized)
if len(samples_str) > 2:
    out_files += [
        analysis_path + '/' + exp_id + '_samples_composition.png'
    ] # samples composition in targets
    out_files += [
        analysis_path + '/' + exp_id + '_samples_normalized_composition.png'
    ] # samples composition (targets normalized)
out_files += [
    analysis_path + '/' + exp_id + '_read_size_quality_with_trimming.png'
] # evolution of read size vs. quality at each trimming step

# Get date
Y = str(datetime.date.today().year)[2:]
M = str(datetime.date.today().month)
if len(M) == 1:
    M = "0" + M
D = str(datetime.date.today().day)
if len(D) == 1:
    D = "0" + D

# Add config details to the analysis results directory
with open(analysis_path + '/config.txt', 'w') as confile:
    confile.write('Analysis_id\t' + exp_id + '\n')
    confile.write('Nb_threads\t' + str(nthreads) + '\n')
    confile.write('Analysis_date\t' + '/'.join([D, M, Y]) + '\n')
    confile.write('Species\t' + str(species_list) + '\n')
    confile.write('Barcodes\t' + barcode_type + '\n')
    confile.write('Samples\t' + str(
        {i.split('=')[0]: i.split('=')[1] for i in samples_str[1:-1].split(',') if len(i.split('=')) > 1}
    ) + '\n')