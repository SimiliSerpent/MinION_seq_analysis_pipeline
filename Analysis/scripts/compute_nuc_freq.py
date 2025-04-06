#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute frequencies of A, T, C and G nucleotides.
Takes as input a fastq file, loop through lines and count the nucleotides.
Outputs the nucleotides frequencies.
"""


import os
import sys

import argparse
import json

import utils # import local utility functions


def increment_counts(count_dict, targets, in_string):
    """Look for substrings and increment a count dictionnary.

    Positionnal arguments:
    count_dict (dict) - dictionnary of targets counts
    targets    (list) - list of target strings
    in_string   (str) - input string to look for targets
    """

    

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', type=str, required=True,
                        help='Path to the input fastq file.')
    parser.add_argument('-o', '--output', type=str,
                        help='Path to the output file (if not specified, only '
                        + 'outputs to stdout)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    nucleotides = ['A', 'T', 'G', 'C']
    count_dict = {charac: 0 for charac in nucleotides}

    # Open the input fastq
    with open(args.fastq, 'r') as fastq:

        # Loop through lines and count nucleotides
        line = fastq.readline()
        read = False
        while line:
            if read:
                for charac in nucleotides:
                    count_dict[charac] += line.count(charac)
                read = False
            if line.startswith('@'):
                read = True
            line = fastq.readline()
    
    # Compute frequencies
    nuc_sum = sum([count_dict[charac] for charac in nucleotides])
    freq_dict = {charac: count_dict[charac]/nuc_sum for charac in nucleotides}

    # Print nucleotides frequencies
    if v > 0:
        for charac in nucleotides:
            print(charac + '\t' + str(100 * freq_dict[charac]) + ' %')

    # Write to output file
    if args.output is not None:
        with open(args.output, 'w') as outfile:
            for charac in nucleotides:
                outfile.write(charac + '\t' + str(100 * freq_dict[charac]) + ' %\n')

    return 0


if __name__ == "__main__":
    sys.exit(main())
