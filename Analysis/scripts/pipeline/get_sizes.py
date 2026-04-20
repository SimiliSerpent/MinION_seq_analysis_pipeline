#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Write sizes of contigs from multifasta file into txt file.
Takes as input a multifasta file, loop through its lines and count the size of
each contig. Write it in a txt file with following formatting:
contig_name\tstart\tstop
"""


import sys

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='Path to the input fasta file.')
    parser.add_argument('-o', '--output', type=str,
                        help='Path to the output txt file (default: same '
                        + 'directory, same name without fasta extension with '
                        + "suffix '_sizes.txt')")
    parser.add_argument('-r', '--roi', type=str,
                        help='Path to the txt file containing list of regions '
                        + 'of interest: if specified only those regions are '
                        + 'retained.')

    args = parser.parse_args()

    fasta_extension = '.' + args.fasta.split('.')[-1]
    output = args.fasta[:-len(fasta_extension)]
    if args.output is not None:
        output = args.output

    if args.roi is not None:
        with open(args.roi, 'r') as roi_file:
            roi_lines = roi_file.readlines()
        roi_list = [line.strip() for line in roi_lines if len(line.strip()) > 0]

    # Use boolean to evaluate if sequence should be retained
    keepseq = True

    with open(output, 'w') as size_file:
        with open(args.fasta, 'r') as fasta_file:

            # Read first fasta line
            line = fasta_file.readline()

            # Retrieve first sequence name
            seq_name = line.strip().split(' ')[0].lstrip('>')
            # print(line.strip().split(' ')[0].lstrip('>'))
            # size_file.write(
            #     '\n' + line.strip().split(' ')[0].lstrip('>') + '\t1\t'
            # )

            # Initialize first sequence length
            seq_size = 0

            # Loop through fasta lines and retrieve sequences
            line = fasta_file.readline()
            while line:

                # If line contains a new sequence, deal with previous sequence..
                if line[0] == '>':

                    # If ROIs have been specified, check for current seq
                    if args.roi is not None:
                        if not seq_name in roi_list:
                            keepseq = False

                    # Write sequence size in corresponding file
                    if keepseq:
                        size_file.write(
                            seq_name
                            + '\t1\t'
                            + str(seq_size)
                            + '\n'
                        )

                    # Retrieve new sequence name and reinitialize count and bool
                    seq_name = line.strip().split(' ')[0].lstrip('>')
                    seq_size = 0
                    keepseq = True

                # ..else, increment current sequence size
                else:
                    seq_size += len(line.strip())

                line = fasta_file.readline()

            # Eventually, deal with the last sequence
            if args.roi is not None:
                if not seq_name in roi_list:
                    keepseq = False
            if keepseq:
                size_file.write(
                    seq_name
                    + '\t1\t'
                    + str(seq_size)
                    + '\n'
                )

    return 0


if __name__ == "__main__":
    sys.exit(main())
