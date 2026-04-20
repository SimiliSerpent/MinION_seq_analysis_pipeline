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


BARCODE_TYPES = ['UNBARCODED', 'ONT-EXP-PBC001', 'TWIST-LRLP-SHv2E']


# Load utils module using abs path
# (yes, it a pain to import modules from a different directory in python)
project_dir = os.environ['CAPASVIR_ANALYSIS'] # path to data analysis repo
spec = importlib.util.spec_from_file_location(
    'utils',
    f'{project_dir}/Analysis/scripts/utils.py'
)
utils = importlib.util.module_from_spec(spec)
sys.modules['utils'] = utils
spec.loader.exec_module(utils)


### COMPUTE SOME NEEDED VALUES

# Get list of target species
species_list = os.environ['SPECIES'].split(',')
species_list.sort()
reference = '_'.join(species_list) + '.fa'

# Get number of threads
nthreads = int(os.environ['NTHREADS'])

# Get type of barcodes
barcode_type = str(os.environ['BARCODES'])
if not barcode_type in BARCODE_TYPES:
    raise ValueError(
        f'The selected barcodes type is not allowed: {barcode_type}'
    )

# Get files name for output
exp_id = os.environ['EXP_ID']
data_path = f'{project_dir}/Data/seq_data/{exp_id}'
analysis_path = f'{project_dir}/Analysis/results/{exp_id}'
raw_fastqs = utils.find_file(f'{data_path}/demultiplexed', '.fastq')
barcode_names = [fastq.rstrip('.fastq') for fastq in raw_fastqs]
barcode_names.sort()

# Get list of samples
sample_names_path = utils.find_path(data_path, f'{exp_id}.barcodes')[0]
with open(sample_names_path, 'r') as file:
    line = file.readline()
samples_str = '"' + line.strip() + '"'

### LIST THE OUTPUT FILES TARGETED BY THE WORKFLOW

out_files = [
    f'{analysis_path}/fastqc/{name}_porechopped_6_GA_fastqc.html' \
    for name in barcode_names
] # fastqc
out_files += [
    f'{analysis_path}/read_cleaning/barcodes_cleaning_stats.png'
] # barcodes cleaning stats
out_files += [
    f'{analysis_path}/read_cleaning/samples_cleaning_stats.png'
] # samples cleaning stats
out_files += [
    f'{analysis_path}/read_lengths/{name}_read_lengths_histo.png' \
    for name in barcode_names
] # read lengths histogram
out_files += [
    f'{analysis_path}/nuc_freq/{exp_id}_nuc_freq.png'
] # frequency of nucleotides
out_files += [
    f'{analysis_path}/composition/{exp_id}_barcodes_composition.png'
] # barcodes composition in targets
out_files += [
    f'{analysis_path}/composition/{exp_id}_barcodes_normalized_composition.png'
] # barcodes composition (targets normalized)
out_files += [
    f'{analysis_path}/composition/{exp_id}_samples_composition.png'
] # samples composition in targets
out_files += [
    f'{analysis_path}/composition/{exp_id}_samples_normalized_composition.png'
] # samples composition (targets normalized)
out_files += [
    f'{analysis_path}/size_vs_quality/{exp_id}_barcode_size_vs_quality.png'
] # evolution of read size vs. quality at each trimming step (barcodes)
out_files += [
    f'{analysis_path}/size_vs_quality/{exp_id}_sample_size_vs_quality.png'
] # evolution of read size vs. quality at each trimming step (samples)
out_files += [
    f'{analysis_path}/size_vs_quality/{exp_id}_barcode_heatmap.png'
] # evolution of read size vs. quality for each species reads (barcodes)
out_files += [
    f'{analysis_path}/size_vs_quality/{exp_id}_sample_heatmap.png'
] # evolution of read size vs. quality for each species reads (samples)

### SAVE CONFIGURATION SETTINGS

# Get date
Y = str(datetime.date.today().year)[2:]
M = str(datetime.date.today().month)
if len(M) == 1:
    M = "0" + M
D = str(datetime.date.today().day)
if len(D) == 1:
    D = "0" + D

# Add config details to the analysis results directory
with open(f'{analysis_path}/config.txt', 'w') as config_log_file:
    config_log_file.write(f'Analysis_id\t{exp_id}\n')
    config_log_file.write(f'Nb_threads\t{nthreads}\n')
    config_log_file.write('Analysis_date\t' + '/'.join([D, M, Y]) + '\n')
    config_log_file.write('Species\t' + str(species_list) + '\n')
    config_log_file.write(f'Barcodes\t{barcode_type}\n')