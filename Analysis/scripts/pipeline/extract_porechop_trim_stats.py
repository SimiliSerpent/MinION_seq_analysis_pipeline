#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot porechop trimming ("porechopping") statistics.
Takes as input the path to the snakemake output logs with the escape sequences
filtered. Parse this files and save the number of sequences trimmed by sample.
"""


import os
import sys

import argparse
import json
import matplotlib.patches as pcs
import matplotlib.pyplot as plt
import numpy as np

import utils # import local utility functions


# Define trimming state dictionnary indicating the next trimming step
MY_TRIMMING_STEPS = {
    '1_seq_adapt': 'seq_adapt',
    '2_twist_outer': 'twist_outer',
    '3_twist_primers': 'twist_primers',
    '2_ont_primer_tails': 'ont_primer_tails',
    '3_ont_unknown_seq': 'ont_unknown_seq',
    '4_tso': 'tso',
    '1_twist_full_seq': 'twist_full_seq',
    'ONT': 'ONT_sequences'
}


def parse_trimming_bloc(file, line, porechop_stats):
    """Parse one trimming bloc in a porechop log. Warning: modifies the dict.

    Positionnal arguments:
    file           (file) - file object being parsed
    line           (line) - current line of the parsed file
    porechop_stats (dict) - dictionnary containing the trimmed (by porechop)
        sequences statistics per sample

    Return:
    line            (str) - last read line
    restart        (bool) - boolean indicating whether the bloc was interrupted
        in an unexpected way or not (i.e. it ended normally)
    """

    # Skip lines to adapter trimming log
    while not line.strip().endswith('reads loaded'):
        line = file.readline()
        # Check if no unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)

    input_seq = int(line.strip().split(' ')[0].replace(',', ''))
    adapters_found = False
    start_trim_seq, start_trim_bases = 0, 0
    end_trim_seq, end_trim_bases = 0, 0

    # Skip lines to adapter trimming log
    while not (line.startswith('Trimming adapters from read ends') or \
    line.startswith('No adapters found')):
        line = file.readline()
        # Check if no unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)
    
    if line.startswith('Trimming adapters from read ends'):
        adapters_found = True

    # Skip adapter trimming lines
    while len(line.strip()) != 0:
        line = file.readline()
        # Check if no unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)

    # Skip empty lines
    while len(line.strip()) == 0:
        line = file.readline()
        # Check if no unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)
    
    if adapters_found:

        # Retrieve number of input sequences
        input_seq = int(line.strip().split(' ')[2].replace(',', ''))
        # Retrieve number of sequences trimmed from read start
        start_trim_seq = int(line.strip().split(' ')[0].replace(',', ''))
        # Retrieve number of bases trimmed from read start
        start_trim_bases = int(
            line.strip().split(' ')[-3].replace(',', '').replace('(', '')
        )
        
        # Skip 1 line
        line = file.readline()
        # Check if unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)
        
        # Retrieve number of sequences trimmed from read end
        end_trim_seq = int(line.strip().split(' ')[0].replace(',', ''))
        # Retrieve number of bases trimmed from read end
        end_trim_bases = int(
            line.strip().split(' ')[-3].replace(',', '').replace('(', '')
        )

        # Go 3 lines ahead
        for skip_line in range(3):
            line = file.readline()
            # Check if unexpected trimming interruption
            if not line or line.strip().startswith('Loading'):
                return(line, True)
    
    # Check for read splitting
    if line.strip() == 'Splitting reads containing middle adapters':
        # Go 3 lines ahead
        for skip_line in range(3):
            line = file.readline()
            # Check if unexpected trimming interruption
            if not line or line.strip().startswith('Loading'):
                return(line, True)
        # Retrieve number of split reads
        split = int(line.split(' ')[0].replace(',', '').replace('(', ''))
        # Go 3 lines ahead
        for skip_line in range(3):
            line = file.readline()
            # Check if unexpected trimming interruption
            if not line or line.strip().startswith('Loading'):
                return(line, True)
    else:
        split = None
    
    # Go 2 lines ahead
    for skip_line in range(2):
        line = file.readline()
        # Check if unexpected trimming interruption
        if not line or line.strip().startswith('Loading'):
            return(line, True)
    
    # Retrieve the barcode and the trimming step
    file_name = line.strip().split('/')[-1].split('.fastq')[0]
    barcode = file_name.split('_')[0]
    trim_step_short = file_name.split('_porechopped_')[1]
    trim_step = MY_TRIMMING_STEPS[trim_step_short]
    
    # Add the barcode and the trimming step in the porechopping stats dict
    if not barcode in porechop_stats.keys():
        porechop_stats[barcode] = {}
    porechop_stats[barcode][trim_step] = {}

    # Update the dictionnary with identified values
    porechop_stats[barcode][trim_step]['input_seq'] = input_seq
    porechop_stats[barcode][trim_step]['start_trim_seq'] = start_trim_seq
    porechop_stats[barcode][trim_step]['start_trim_bases'] = start_trim_bases
    porechop_stats[barcode][trim_step]['end_trim_seq'] = end_trim_seq
    porechop_stats[barcode][trim_step]['end_trim_bases'] = end_trim_bases
    if split is not None:
        porechop_stats[barcode][trim_step]['split'] = split

    return(line, False)


def parse_output_logs(path_to_output_log, verbose=0):
    """Extract the porechopping stats from logs file.

    Positionnal arguments:
    path_to_output_log (str) - concatenated snakemake output logs with the 
                               escape sequences filtered.
    verbose            (int) - degree of verbosity (default: 0 = muted).

    Return:
    porechop_stats (dict)    - dictionnary containing the trimmed (by porechop)
                               sequences statistics per sample.
    """

    # Initialize dictionnary with the porechopping stats by sample
    porechop_stats = {}

    # Open logs file and parse it
    with open(path_to_output_log, 'r') as file:
        line  = file.readline()

        # Loop through file lines
        while line :

            # Look for any file porechopping start
            if line.strip() == 'Loading reads':
                utils.send_text('Parsing new porechop output ', verbose, 2, 1)
                line, restart = parse_trimming_bloc(file, line, porechop_stats)

                # If the porechop output is truncated (new output log starts in
                # the middle of the former one, or end of file), skip it
                if restart:
                    continue
                
            line = file.readline()
    
    if verbose > 0:
        utils.print_json(porechop_stats)

    return porechop_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log_files', type=str, required=True,
                        help='Path to the directory containing the snakemake '
                        + 'output logs .txt files with the escape sequences '
                        + 'filtered')
    parser.add_argument('-o', '--output_json', type=str,
                        help='Path to the output file (default: same directory'
                        + ' as input file, same name with json extension)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Find log files
    log_files = utils.find_path(args.log_files, '.txt')
    utils.send_text(f'Found {len(log_files)} porechop log files', v, 1, 0)

    # Retrieve the porechop trimming statistics
    porechop_stats = {}
    for log_file in log_files:
        utils.send_text(f'Parsing {log_file.split("/")[-1]}', v, 2, 0)
        new_stats = parse_output_logs(log_file, verbose=v)
        # Cleanly update the statistics dictionnary
        for bc in new_stats.keys():
            if not bc in porechop_stats.keys():
                porechop_stats[bc] = {}
            for trim_step in new_stats[bc].keys():
                porechop_stats[bc][trim_step] = new_stats[bc][trim_step]

    # Save porechopping statistics to file
    output_path = f'{args.log_files.rstrip("/")}/porechopping_stats.json'
    if args.output_json is not None:
        output_path = args.output_json
    utils.send_text('Saving porechop trim stats to ' + output_path, v, 1, 0)
    utils.save_json(output_path, porechop_stats)

    return 0


if __name__ == "__main__":
    sys.exit(main())
