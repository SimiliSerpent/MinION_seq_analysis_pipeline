#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Sub-sample a fastq file.
Takes as input a fastq file, and subsample it keeping only reads matching the
specified criteria (above length & quality threshold). Propose each read to the
user who is free to keep it in output file, or to leave it.
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
    parser.add_argument('-l', '--min_length', type=int, default=400,
                        help='Minimal size of output reads (default: 400)')
    parser.add_argument('-s', '--qual_score', type=int, default=10,
                        help='Minimal q score of output reads (default: 10)')
    parser.add_argument('-n', '--nb_reads', type=int,
                        help='Number of reads to output. If this option is '
                        + 'used, the first n reads matching the specified '
                        + 'condition will be kept in the output file.')
    parser.add_argument('-q', '--output_fastq', type=str,
                        help='Path to the output fastq')
    parser.add_argument('-a', '--output_fasta', type=str,
                        help='Path to the output fasta')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    fastq = []
    fasta = []
    passing_thresholds = 0

    keep = ''

    choose_each_read = True
    if args.nb_reads is not None:
        choose_each_read = False

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
            
            # Check whether read has length and qual above threshold
            if len(seq) >= args.min_length and qscore >= args.qual_score:

                passing_thresholds += 1

                if choose_each_read:
                    # Ask if read should be retained
                    print(
                        '\nRead ' \
                        + header \
                        + ' sequence:\n' \
                        + seq \
                        + '\nThis read has quality: ' \
                        + str(qscore)
                    )

                    tbp_text = 'You have seen '\
                            + str(passing_thresholds)\
                            + ' reads longer than '\
                            + str(args.min_length)\
                            + ' b and with quality over '\
                            + str(args.qual_score)\
                            + ' so far.\nYou have chosen to keep '\
                            + str(len(fastq))\
                            + ' of them.\nDo you want to keep this read? '\
                            + '(y/n/stop) '
                    keep = input(tbp_text)

                    while not keep in ['y', 'n', 'stop']:
                        print('\nYou must chose one of the following: y|yes n|no' \
                            + ' stop|[stop the process and store the reads in ' \
                            + 'fasta/fastq] qual|[display read quality phred string]')
                        keep = input(tbp_text)

                    if keep == 'y':
                        # If read kept by user, add the read to the list
                        fasta.append('>' + header + '\n' + seq + '\n')
                        fastq.append(read)
                
                    if keep == 'stop':
                        # If instructed by user, stop the loop
                        break
                
                else:

                    # Add the read to the list
                    fasta.append('>' + header + '\n' + seq + '\n')
                    fastq.append(read)

                    # If nb_reads found, stop the loop
                    if len(fastq) == args.nb_reads:
                        break
            
            line = fastq_in.readline()
        
    if keep != '':
        print('Number of reads in sample file: ' + str(len(fastq)))

    # Write to output file
    if args.output_fastq is not None:
        with open(args.output_fastq, 'w') as outfile:
            for read in fastq:
                outfile.write(read)
    if args.output_fasta is not None:
        with open(args.output_fasta, 'w') as outfile:
            for seq in fasta:
                outfile.write(seq)

    return 0


if __name__ == "__main__":
    sys.exit(main())
