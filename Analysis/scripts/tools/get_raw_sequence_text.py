#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert sequences in Fastq format to raw text.
Takes as input a fastq file and outputs a txt file with sequences, one by line.
"""


import sys

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, required=True,
                        help='Path to the fastq file to filter.')
    parser.add_argument('-o', '--output_path', type=str,
                        help='Path to the filtered output txt file (default: '
                        + "same as input with '.txt' extension).")

    args = parser.parse_args()

    output_path = args.input_path[:-5] + 'txt'
    if args.output_path is not None:
        output_path = args.output_path

    with open(args.input_path, 'r') as fastq:
        with open(output_path, 'w') as txt:
            line = fastq.readline()
            while line:
                if line.startswith('@'):
                    line = fastq.readline() # go to seq
                    txt.write(line)
                    line = fastq.readline() # go to '+'
                    line = fastq.readline() # go to PHRED
                line = fastq.readline() # go to next sequence

    return 0


if __name__ == "__main__":
    sys.exit(main())
