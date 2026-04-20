#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Modifies a bam to convert all mapped reads to forward orientation
Takes as input a bam file and convert reverse reads to forward orientation.
Returns a bam and (optionnaly) a list of modified reads ids.
"""


import os
import sys

import argparse
import pysam

import utils


def make_all_forward(bam_in, bam_out, ids=None, v=0):
    """Makes the reverse read of a bam forward.
    Parse the sequences in the bam using pysam and updates sequence, quality
    string and flag. Takes a bam file as input and outputs a bam.
    
    bam_in   (str) - Path to the input bam file.
    file_out (str) - Path to the output bam file.
    ids      (str) - Path to output txt file containing the ids of modified
                     reads (optionnal). (default: None = ids not saved)
    v        (str) - Path to the input bam file.
    """
    read_count = 0
    read_rev = 0
    rev_ids = []
    # Open input and output files using `with` to ensure closure
    with pysam.AlignmentFile(bam_in, 'rb') as b_in:
        with pysam.AlignmentFile(bam_out, "wb", header=b_in.header) as b_out:

            # Loop through sequences
            for read in b_in:
                read_count += 1
                if read.is_reverse:
                    read_rev += 1
                    if ids is not None:
                        rev_ids.append(read.query_name)
                    # if read.query_qualities is not None:
                    #     read_qual = read.query_qualities[::-1]
                    # # Reverse-complement sequence
                    # if read.query_sequence is not None:
                    #     read.query_sequence = utils.rev_comp(read.query_sequence)
                    # # Reverse qualities
                    # if read.query_qualities is not None:
                    #     read.query_qualities = read_qual
                    # Clear reverse flag
                    read.flag &= ~16

                b_out.write(read)

    utils.send_text(f'Reversed {read_rev} out of {read_count} reads', v, 2, 2)

    # Index output file
    pysam.index(bam_out)

    # Write out names of modified (reverse) reads
    if ids is not None:
        with open(ids, 'w') as ids_out:
            ids_out.write('\n'.join(rev_ids))
            
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='BAM file containing (some) reverse aligned '
                        + 'reads to convert to forward.')
    parser.add_argument('-o', '--output', type=str,
                        help='Output bam file with only forward reads.')
    parser.add_argument('-n', '--reverse_names', type=str,
                        help='Name of output txt file with ids of reverse '
                        + '(modified) reads. (Optionnal)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    utils.send_text(f'Processing {args.input}', v, 1, 0)
    make_all_forward(args.input, args.output, args.reverse_names, v)

    return 0


if __name__ == "__main__":
    sys.exit(main())
