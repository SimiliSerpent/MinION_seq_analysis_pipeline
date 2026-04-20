#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Filter out reads whose lengths have been modified.
From a pair of fastq files, find out what reads have been modified in the
second compared to the first, and output these reads in a separate fastq.
Optionnally outputs untouched reads. Works by comparing reads length.
"""


import os
import sys

import argparse
import contextlib
import numpy as np
import pysam

import utils


def get_lengths(fq_path, v=0):
    """Retrieve reads lengths from fastq
    Uses pysam to parse fastq file and output reads length in a python
    dictionnary.

    Arguments:
    fq_path (str) - Path to input fastq file.
    v       (str) - Level of verbosity (default: 0 = muted)
    
    Return:
    reads_lengths (dict) - Dictionnary containing the reads lengths with reads
                           ids as keys.
    """
    # Initialize dictionnary
    reads_lengths = {}

    # Stream through fastq file
    with pysam.FastxFile(fq_path) as fastq_stream:
        for read in fastq_stream:
            reads_lengths[read.name] = len(read.sequence)
    
    nb_reads = len(reads_lengths)
    utils.send_text(
        f'{nb_reads} reads had their length stored', v, 2, 1
    )
    
    return(reads_lengths)


def filter_reads(fq_path, ref_lengths, mod_path=None, unmod_path=None,
                 unfound_path=None, v=0):
    """Filter reads from fastq with length same/different from reference
    Outputs are optionnal: if none is defined, fastq file is parsed and
    analyzed yet nothing will happen.

    Arguments:
    fq_path      (str)  - Path to input fastq file.
    ref_lengths  (dict) - Dictionnary containing the reads reference lengths
                          with reads ids as keys.
    mod_path     (str)  - Path to the output fastq file containing the reads
                          with modified length relatively to reference length
                          (default: None = reads are discarded).
    unmod_path   (str)  - Path to the output fastq file containing the reads
                          with unmodified length relatively to reference length
                          (default: None = reads are discarded).
    unfound_path (str)  - Path to the output fastq file containing the reads
                          with no reference length (default: None = reads are
                          discarded).
    v            (str)  - Level of verbosity (default: 0 = muted)
    
    Return:
    reads_lengths (dict) - Dictionnary containing the reads lengths with reads
                           ids as keys.
    """
    # Initialize counters to report statistics in std out
    mod_count = 0
    unmod_count = 0
    unfound_count = 0

    # Stream through fastq file if these are defined
    with pysam.FastxFile(fq_path) as fastq_in, (
            open(mod_path, 'w')
            if mod_path is not None
            else contextlib.nullcontext()
        ) as mod_file, (
            open(unmod_path, 'w')
            if unmod_path is not None
            else contextlib.nullcontext()
        ) as unmod_file, (
            open(unfound_path, 'w')
            if unfound_path is not None
            else contextlib.nullcontext()
        ) as unfound_file:
        
        # Loop through reads
        for read in fastq_in:

            # Write read to correct output
            if read.name in ref_lengths.keys():

                if len(read.sequence) == ref_lengths[read.name]:
                    if unmod_path is not None:
                        unmod_file.write(str(read) + '\n')
                        unmod_count += 1
                
                else:
                    if mod_path is not None:
                        mod_file.write(str(read) + '\n')
                        mod_count += 1
            
            else:
                if unfound_path is not None:
                    unfound_file.write(str(read) + '\n')
                    unfound_count += 1
    
    if mod_path is not None:
        utils.send_text(
            f'{mod_count} reads with modified length '
            + f'were written to {mod_path}',
            v, 2, 1
        )
    
    if unmod_path is not None:
        utils.send_text(
            f'{unmod_count} reads with same length '
            + f'were written to {unmod_path}',
            v, 2, 1
        )
    
    if unfound_path is not None:
        utils.send_text(
            f'{unfound_count} reads were not found in former fastq', v, 2, 1
        )
            
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, nargs=2,
                        help='Space-separated pair of paths to fastq files to '
                        + 'be compared. Retained reads are the one from the '
                        + 'second fastq whose length has been altered from '
                        + 'fastq 1. Reads unfound in fastq 1 are discarded.')
    parser.add_argument('-m', '--modified', type=str, required=True,
                        help='Path to output fastq file containing reads with '
                        + 'modified length.')
    parser.add_argument('-u', '--untouched', type=str,
                        help='Path to output fastq file containing reads with '
                        + 'UNmodified length. Optionnal (by default this reads'
                        + ' are not saved).')
    parser.add_argument('-n', '--not_found', type=str,
                        help='Path to output fastq file containing reads from '
                        + 'fastq 2 not found in fastq 1. Optionnal (by default'
                        + ' this reads are not saved).')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Build (huge) dictionnary of fastq 1 reads lengths
    utils.send_text(f'Gathering {args.input[0]} reads lengths', v, 1, 0)
    reads_lengths = get_lengths(args.input[0], v)

    # Filter reads from fastq 2
    utils.send_text(f'Filtering out {args.input[1]}', v, 1, 0)
    filter_reads(args.input[1], reads_lengths, args.modified, args.untouched,
                 args.not_found, v)

    return 0


if __name__ == "__main__":
    sys.exit(main())
