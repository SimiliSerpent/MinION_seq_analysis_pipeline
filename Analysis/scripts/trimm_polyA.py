#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Remove poly-A tails.
Takes as input the fastq file to be trimmed and output the fasta without poly-
A tail. Looks for poly-A tail at the end of reads, and for poly-T tails at the
start of reads. Consider single non-A nucleotides inside a poly-A tail to be a
sequencing/polymerization error. The tail must start with at least 'AAA'.
"""


import os
import sys

import argparse

import utils # import local utility functions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_fastq', type=str, required=True,
                        help='Path to the input fastq file.')
    parser.add_argument('-o', '--out_fastq', type=str,
                        help='Path to the output fastq file (default: <in_name'
                        + '>_no_tail.fastq)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Open input fastq
    with open(args.in_fastq, 'r') as fastq_in:

        # Open output fastq to write to
        if args.out_fastq is None:
            out_path = args.in_fastq.rsplit('.fastq', 1)[0] + '_no_tail.fastq'
        else:
            out_path = args.out_fastq

        with open(out_path, 'w') as fastq_out:

            # Use counter for statistics and booleans for localization
            seq_line = False
            score_line = False
            trimmed_size = 0
            nb_trimmed_polyA = 0
            nb_trimmed_polyT = 0
            nb_trimmed_bases = 0

            # Loop through input lines
            line = fastq_in.readline()
            while line:

                # Copy header and detect coming sequence
                if not score_line and line.startswith('@'):
                    seq_name = line
                    seq_line = True
                
                # Trimm and detect coming seqcuencing score
                elif seq_line:
                    seq = line.strip()
                    seq_size = len(seq)
                    right_trimmed = 0
                    left_trimmed = 0

                    # Trimm poly-A
                    trimmed_polyA = False
                    while seq[-5:].count('A') >= 3:

                        # Remove trailing A-s
                        if seq.endswith('A'):
                            trimmed_polyA = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq.rstrip('A') # Remove A trailing bases
                            right_trimmed += seq_size - len(seq) # Count the number of trimmed bases
                            seq_size = len(seq) # Refresh the size marker
                        
                        # Remove what is thought to be a falsely copied/sequenced A
                        else:
                            trimmed_polyA = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq[:-1] # Remove the last base
                            right_trimmed += 1 # Count the number of trimmed bases
                            seq_size -= 1 # Refresh the size marker

                    # Trimm poly-T
                    trimmed_polyT = False
                    while seq[:5].count('T') >= 3:

                        # Remove trailing T-s
                        if seq.startswith('T'):
                            trimmed_polyT = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq.lstrip('T') # Remove A trailing bases
                            left_trimmed += seq_size - len(seq) # Count the number of trimmed bases
                            seq_size = len(seq) # Refresh the size marker
                        
                        # Remove what is thought to be a falsely copied/sequenced T
                        else:
                            trimmed_polyT = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq[1:] # Remove the first base
                            left_trimmed += 1 # Count the number of trimmed bases
                            seq_size -= 1 # Refresh the size marker
                    
                    # Increment counters
                    if trimmed_polyA: nb_trimmed_polyA += 1
                    if trimmed_polyT: nb_trimmed_polyT += 1
                    nb_trimmed_bases += left_trimmed + right_trimmed

                    # Copy sequence if not empty
                    if seq != '':
                        fastq_out.write(seq_name)
                        fastq_out.write(seq + '\n+\n')
                        seq_line = False
                        score_line = True
                        # Skip '+' line
                        line = fastq_in.readline()
                    
                    else:
                        seq_line = False
                
                # Trim sequencing score accordingly
                elif score_line:
                    score = line.strip()[left_trimmed:]
                    if right_trimmed > 0: score = score[:-right_trimmed]
                    fastq_out.write(score + '\n')
                    score_line = False

                line = fastq_in.readline()
        
        # Output trimming statistics
        if v > 0:
            print('Trimmed', nb_trimmed_polyA, "polyA ending tail.")
            print('Trimmed', nb_trimmed_polyT, "polyT ending tail.")
            print('Total trimmed bases:', nb_trimmed_bases)

    return 0


if __name__ == "__main__":
    sys.exit(main())
