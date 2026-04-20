#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Place back the UMI at the beginning of the reads after using umi-tools
Takes as input a bam file put back the UMI from the read name to the beginning 
of the sequence. Adds top quality scores for the UMIs and updates CIGAR.
"""


import os
import sys

import argparse
import pysam

import utils


def get_umis_json(json_path):
    """Retrieve reads UMIs from json file.
    
    json_path (str) - Path to the json file with UMIs as keys and list of
                      associated read ids as values.

    Return:
    umis     (dict) - Dictionnary with read ids as keys and associated UMIs as
                      values.
    """
    umis = {}
    umis_dict = utils.load_json(json_path)
    for umi, ids in umis_dict.items():
        for id in ids:
            umis[id] = umi
    return umis


def put_umis_at_start(bam_in, bam_out, umis='end'):
    """Place UMIs back at the beginning of the reads.
    
    bam_in   (str) - Path to the input bam file.
    file_out (str) - Path to the output bam file.
    umis    (dict) - String indicating where to find reads UMIs: can be either
                     'start' (UMI at read id start separated by _); 'end' (UMI
                     at read is end separated by _); or a path to a json files
                     containing UMIs as keys and associated read ids as values.
                     default: 'end'
    """
    # Retrieve UMIs if stored in separated file
    if umis != 'end' and umis != 'start':
        umis_dict = get_umis_json(umis)

    # Open input and output files using `with` to ensure closure
    with pysam.AlignmentFile(bam_in, 'rb') as b_in:
        with pysam.AlignmentFile(bam_out, 'wb', header=b_in.header) as b_out:

            # Loop through sequences
            for read in b_in:

                # Retrieve UMI
                if read.query_name is not None:
                    if umis == 'start':
                        umi = read.query_name.split('_')[0]
                    elif umis == 'end':
                        umi = read.query_name.split('_')[-1]
                    else:
                        umi = umis_dict[read.query_name]

                # Update sequence (store quality beforehand)
                if read.query_qualities is not None:
                    qual = read.query_qualities
                if read.query_sequence is not None:
                    read.query_sequence = umi + read.query_sequence

                # Update quality string
                if read.query_qualities is not None:
                    read.query_qualities = [90] * len(umi) + list(qual)

                # Update CIGAR
                if read.cigartuples is not None:
                    if read.cigartuples[0][0] == 4:
                        read.cigartuples = [
                            (4, read.cigartuples[0][1] + len(umi))
                        ] + read.cigartuples[1:]
                    else:
                        read.cigartuples = [(4, len(umi))] + read.cigartuples
                
                # if read.query_name == 'c27dd2df-e18b-4fb4-b60f-bfc32307235e_TCGATCTG':
                #     L = len(read.query_sequence)
                #     print(f'READ {read.query_name}')
                #     print(f'of len {L}')
                #     print(f'and CIGAR {read.cigartuples}')

                b_out.write(read)

    # Index output file
    pysam.index(bam_out)
            
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input BAM file with extracted UMIs in read '
                        + 'names.')
    parser.add_argument('-o', '--output', type=str,
                        help='Output bam file with UMIs appended at the '
                        + 'beginning of the reads.')
    parser.add_argument('-u', '--umis', type=str, default='end',
                        help="Way to find the UMIs. Can be either 'start' "
                        + "(UMI at read id start separated by _); 'end' (UMI "
                        + "at read is end separated by _); or a path to a json"
                        + " files containing UMIs as keys and associated read "
                        + "ids as values. default: 'end'")
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    utils.send_text(f'Processing {args.input}', v, 1, 0)
    put_umis_at_start(args.input, args.output, args.umis)

    return 0


if __name__ == "__main__":
    sys.exit(main())
