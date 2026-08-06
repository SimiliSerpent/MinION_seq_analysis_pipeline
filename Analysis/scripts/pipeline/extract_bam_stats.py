#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse a BAM file and extract records statistics.
Only takes primary alignements into account. Unmap reads are also absent from
output statistics.
"""


import argparse
import os
import sys
import time

import pysam

import utils


# Set to True if comparing with umi-tools to ensure identical bam input
TEST_COMPARE_UMITOOLS = False

# CIGAR operations used for CIGAR updating
M, I, D, N, S, EQ, X = 0, 1, 2, 3, 4, 7, 8
QUERY_OPS = {M, I, S, EQ, X}
REF_OPS   = {M, D, N, EQ, X}


def retrieve_reads_stats(bam_path, v=0):
    """Retrieve reads alignment stats from a bam file

    Arguments:
    bam_path (str) - Path to input BAM file.
    v        (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    """
    # Initialize dictionnary
    reads_stats = {}

    # Store times for progression status messages
    if v > 0:
        t_zero = time.perf_counter()
        query_count = 0

    # Stream through input fastq file
    with pysam.AlignmentFile(bam_path, 'rb') as b_in:
        
        # Loop through alignments
        for query in b_in:

            if v > 0:
                query_count += 1
                if query_count % 10000 == 0:
                    elapsed = time.perf_counter() - t_zero
                    utils.send_text(
                        f'Processed {query_count} unempty bam records in '
                        + f'{elapsed:.3f} s.',
                        v, 2, 1
                    )

            # One cannot retrieve alignment stats if no sequence is available
            if query.query_sequence is None:
                continue
            # TODO: Deal with case where no reference is accepted
            # (i.e. group reads solely on UMI seq, not mapping proximity) 
            if query.reference_name is None:
                continue
            # Keep only primary alignments
            if query.is_secondary or query.is_supplementary:
                continue

            seq = query.query_sequence
            qual = query.query_qualities
            cigar = query.cigartuples

            # Compute read alignment length (if aligned)
            if cigar is not None:
                aligned_len = sum(
                    length
                    for op, length in query.cigartuples
                    if op in (M, EQ, X)
                )
            else:
                aligned_len = 0

            # TODO: remove this snippet
            if TEST_COMPARE_UMITOOLS:
                # Update name to remove UMI from read id
                umi = query.query_name.split('_')[-1]
                query.query_name = query.query_name[:-(len(umi)+1)]

            # Update reads alignment lengths
            orientation = 'unaligned'
            if cigar is not None:
                orientation = 'reverse' if query.is_reverse else 'forward'
            reads_stats[query.query_name] = {
                'aligned_len': aligned_len,
                'ref': query.reference_name.replace('_', '-') \
                    if cigar is not None \
                    else 'unaligned',
                'ref_start': query.reference_start \
                    if cigar is not None \
                    else -1,
                'orientation': orientation
            }
    
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Processed {query_count} unempty bam records in '
            + f'{elapsed:.3f} s.',
            v, 2, 1
        )
            
    return reads_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input_bam', type=str, required=True,
                        help='Path to input bam file with raw reads aligned to'
                        ' the reference.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output JSON file containing aligned '
                        'reads statistics.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose
    
    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    # Retrieve reads statistics
    utils.send_text(f'Retrieving reads statistics from {args.input_bam}',
                    v, 1, 0)
    reads_stats = retrieve_reads_stats(args.input_bam, v)

    # Save reads data json
    utils.send_text('Saving reads statistics information', v, 1, 0)
    utils.save_json(args.output, reads_stats)

    return 0


if __name__ == "__main__":
    sys.exit(main())
