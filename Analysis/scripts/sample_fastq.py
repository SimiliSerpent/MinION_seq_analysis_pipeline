#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute frequencies of A, T, C and G nucleotides.
Takes as input a fastq file, loop through lines and count the nucleotides.
Outputs the nucleotides frequencies.
"""


import os
import sys

import argparse
import random


SANGER_SCORE_OFFSET = 33
INVALID_CHAR_CODE = 200

# Define a mapping allowing ASCII characters translation in Phred scores
q_mapping = bytes(
    (
        letter - SANGER_SCORE_OFFSET
        if SANGER_SCORE_OFFSET <= letter < 94 + SANGER_SCORE_OFFSET
        else INVALID_CHAR_CODE
    )
    for letter in range(256)
)


def get_score(qual_score):

    phred = list(qual_score.encode().translate(q_mapping))

    return(sum(phred)/len(phred))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fastq', type=str, required=True,
                        help='Path to the input fastq file.')
    parser.add_argument('-l', '--min_length', type=int,
                        help='Minimal size of output reads (default: 500)')
    parser.add_argument('-s', '--qual_score', type=int,
                        help='Minimal q score of output reads (default: 500)')
    parser.add_argument('-q', '--output_fastq', type=str,
                        help='Path to the output fastq (if not specified, same'
                        + ' as input with the suffix selection.fastq)')
    parser.add_argument('-a', '--output_fasta', type=str,
                        help='Path to the output fasta (if not specified, same'
                        + ' as input with the suffix selection.fastq)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    fastq = []
    fasta = []
    longer_than_threshold = 0

    keep = ''

    # Open the input fastq
    with open(args.input_fastq, 'r') as fastq_in:

        # Loop through lines and parse reads
        line = fastq_in.readline()
        while line:
            
            read = line # read header
            header = line.strip()[1:].split()[0] # store header
            line = fastq_in.readline() # read sequence
            read += line
            seq = line.strip() # store sequence
            read += fastq_in.readline() # read '+'
            line = fastq_in.readline() # read qual score
            read += line
            qual = line.strip() # store sequence
            qscore = get_score(qual)
            
            # Ask if read should be retained for reads longer than threshold length
            if len(seq) >= args.min_length and qscore >= args.qual_score:

                longer_than_threshold += 1

                print('\nRead ' + header + ' sequence:\n' + seq + '\nThis read has quality: ' + str(qscore))

                tbp_text = 'You have seen ' \
                           + str(longer_than_threshold) \
                           + ' reads longer than ' \
                           + str(args.min_length) \
                           + ' b so far.\nYou have chosen to keep ' \
                           + str(len(fastq)) \
                           + ' of them.\nDo you want to keep this read? ' \
                           + '(y/n/stop) '
                keep = input(tbp_text)

                while not keep in ['y', 'n', 'stop']:
                    print('\nYou must chose one of the following: y|yes n|no' \
                        + ' stop|[stop the process and store the reads in ' \
                        + 'fasta/fastq] qual|[display read quality phred string]')
                    keep = input(tbp_text)

                if keep == 'y':
                    fasta.append('>' + header + '\n' + seq + '\n')
                    fastq.append(read)
            
                if keep == 'stop':
                    break
            
            line = fastq_in.readline()
        
    if keep != '':
        print('You have selected ' + str(len(fastq)) + ' reads.')

    # Write to output file
    if args.output_fastq is not None:
        fq_path = args.output_fastq
    else:
        fq_path = args.input_fastq[:-6] + '_selection.fastq'
    if args.output_fasta is not None:
        fa_path = args.output_fasta
    else:
        fa_path = args.input_fastq[:-6] + '_selection.fasta'
        
    with open(fq_path, 'w') as outfile:
        for read in fastq:
            outfile.write(read)
    with open(fa_path, 'w') as outfile:
        for seq in fasta:
            outfile.write(seq)

    return 0


if __name__ == "__main__":
    sys.exit(main())
