#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Homemade UMIs extraction tool.
From a BAM file and a UMI sequence, extract UMIs and optionally outputs a BAM
with UMIs extracted. These reads can be re-mapped after UMIs extraction.
"""


import argparse
import contextlib
import itertools
import os
import sys
import time

import pysam

import utils


# CIGAR operations used for CIGAR updating
M, I, D, N, S, EQ, X = 0, 1, 2, 3, 4, 7, 8
QUERY_OPS = {M, I, S, EQ, X}
REF_OPS   = {M, D, N, EQ, X}


def update_cigar(cigar, to_be_removed, v=0):
    """Update CIGAR by removing bases at given indices

    Arguments:
    cigar      (list of tuples) - CIGAR as output by pysam while parsing a bam.
    to_be_removed (list of int) - List of indices of deleted nucleotides.
    v                     (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    new_cigar  (list of tuples) - Updated CIGAR in same format as input.
    """
    target_indices = set(to_be_removed)
    query_pos = 0
    new_cigar = []

    # Stream through CIGAR operation x number of corresponding bases
    for op, length in cigar:

        # Keep CIGAR parts that do not correspond to the updated queried read
        if op not in QUERY_OPS:
            new_cigar.append((op, length))
            continue

        # Otherwise, remove letters that corresponds to removed indices...
        kept = []
        for i in range(length):
            # ...except if those positions are matching the reference !
            # In that case, replace the letter by a deletion instead.
            if query_pos in to_be_removed:
                if op in REF_OPS:
                    kept.append(D)
            else:
                kept.append(op)
            query_pos += 1
            
        for k, grp in itertools.groupby(kept):
            new_cigar.append((k, len(list(grp))))
    
    # Reconstruct the tuple-like CIGAR structure as returned by pysam
    # query.cigartuples, remove 'empty' tuples and aggregate adjacent identical
    # operations
    temp = []
    for op, length in new_cigar:
        if length > 0:
            if len(temp) > 0 and op == temp[-1][0]:
                temp[-1] = (temp[-1][0], temp[-1][1] + length)
            else:
                temp.append((op, length))
    new_cigar = temp

    return new_cigar


def retrieve_umis(bam_path, pattern, out_bam=None, out_UMIs=None, v=0):
    """Retrieve UMIs from a bam file
    UMIs are retrieved following the specified pattern. A BAM file with UMIs
    cut is optionally produced.

    Arguments:
    bam_path (str) - Path to input BAM file.
    pattern  (str) - Pattern telling how to look for UMI at read start. First
                     pattern character must be S|E|3|5 to indicate where to
                     look for the UMI. S = always at read start; E = always at
                     read end; 3 = always at read 3' end; 5 = always at read 5'
                     end. Immediately after should be a sequence composed with
                     N|X: N bases are UMI bases, X bases will be ignored
                     (eg SXXNNNXN). Pattern sequence should always be given in
                     the 5'->3' direction.
    out_bam  (str) - Path to output BAM file. (Optional)
    out_UMIs (str) - Path to output UMIs statistics file. (Optional)
    v        (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    umis        (dict) - Dictionnary containing the UMIs as keys and and the
                         ids of the reads associated with each UMI.
    """
    # Initialize dictionnary
    umis = {}

    # Retrieve read end where to look for UMI
    side, pattern = pattern[0], pattern[1:]
    start = (side == 'S')
    end = (side == 'E')
    five_p = (side == '5')
    three_p = (side == '3')

    # Retrieve indices of UMI bases
    base_umi_inds = [ind for ind, char in enumerate(pattern) if char == 'N']

    # Store times for progression status messages
    if v > 0:
        t_zero = time.perf_counter()
        query_count = 0

    # Stream through input fastq file
    with pysam.AlignmentFile(bam_path, 'rb') as b_in, (
            pysam.AlignmentFile(out_bam, 'wb', header=b_in.header)
            if out_bam is not None
            else contextlib.nullcontext()
        ) as b_out:
        
        # Loop through alignments
        for query in b_in:

            if v > 0:
                query_count += 1
                if query_count % 10000 == 0:
                    elapsed = time.perf_counter() - t_zero
                    utils.send_text(
                        f'Processed {query_count} unempty bam records in '
                        + f'{elapsed:.1f} s.',
                        v, 1, 1
                    )

            # One cannot update the alignment if no sequence is available
            if query.query_sequence is None:
                continue
            # If no CIGAR information is available (unaligned read), then the
            # bam record can only be searched for UMIs if alignment orientation
            # does not matter, i.e. if UMIs is searched for at read start/end
            if query.reference_name is None and (five_p or three_p):
                continue
            # Keep only primary alignments
            if query.is_secondary or query.is_supplementary:
                continue

            seq = query.query_sequence
            qual = query.query_qualities
            cigar = query.cigartuples

            # Update UMI indices depending on UMI position and read
            # orientation and retrieve UMI
            umi_inds = base_umi_inds
            if query.is_reverse and (start or end):
                if start:
                    umi_inds = [len(seq)-i-1 for i in umi_inds[::-1]]
                if end:
                    umi_inds = [len(pattern)-i-1 for i in umi_inds[::-1]]
                umi = utils.rev_comp(''.join([seq[i] for i in umi_inds]))
            if three_p or (end and not query.is_reverse):
                umi_inds = [[len(seq)-len(pattern)+i for i in umi_inds]]
                umi = ''.join([seq[i] for i in umi_inds])
            if five_p or (start and not query.is_reverse):
                umi = ''.join([seq[i] for i in umi_inds])
            
            # Update sequence, quality string and CIGAR (if read aligned)
            new_seq = ''
            new_qual = []
            to_be_removed = set(umi_inds)
            for ind, base in enumerate(seq):
                if ind not in to_be_removed:
                    new_seq += base
                    new_qual.append(qual[ind])
            query.query_sequence = new_seq
            query.query_qualities = new_qual
            if cigar is not None:
                query.cigartuples = update_cigar(cigar, umi_inds)
            
            # Classify read
            if umi in umis.keys():
                umis[umi].append(query.query_name)
            else :
                umis[umi] = [query.query_name]

            # Write out query
            if out_bam is not None and len(query.query_sequence) > 0:
                b_out.write(query)
    
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Processed {query_count} bam records in '
            + f'{elapsed:.1f} s.',
            v, 1, 1
        )
    
    # Write out UMIs
    if out_UMIs is not None:
        header = ['UMI_sequence', 'Nb_of_reads', 'List_of_reads_ids']
        with open(out_UMIs, 'w') as file_out:
            file_out.write('\t'.join(header) + '\n')
            umis_list = list(umis.keys())
            umis_list.sort()
            for umi in umis_list:
                fields =  [
                    umi,
                    str(len(umis[umi])),
                    ','.join(umis[umi])
                ]
                file_out.write('\t'.join(fields) + '\n')
            
    return umis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input_bam', type=str, required=True,
                        help='Path to input bam file with raw reads aligned to'
                        + ' the reference.')
    parser.add_argument('-p', '--pattern', type=str, required=True,
                        help='Pattern telling how to look for UMI at read '
                        + 'start. First pattern character must be S|E|3|5 to '
                        + 'indicate where to look for the UMI. S = always at '
                        + 'read start; E = always at read end; 3 = always at '
                        + "read 3' end; 5 = always at read 5' end. Immediately"
                        + ' after should be a sequence composed with N|X: N '
                        + 'bases are UMI bases, X bases will be ignored (eg '
                        + 'SXXNNNXN). Pattern sequence should always be given '
                        + "in the 5'->3' direction.")
    parser.add_argument('-o', '--output_prefix', type=str, default=None,
                        help='Prefix for output files. Default is the path to '
                        + 'the input bam without the .bam extension. Note '
                        + 'that this can be used to redirect outputs to a '
                        + 'given directory. Output files are: _UMI_trimmed.bam'
                        + ' (bam file with the UMI withdrawn from records, '
                        + 'primary alignments only, unaligned reads discarded '
                        + "if looking for UMI pattern at 5' or 3' end), "
                        + '_raw_umis.tsv (tsv file with raw UMIs and their '
                        + 'associated reads ids), _raw_umis.json')
    parser.add_argument('-ob', '--output_bam', type=str, default=None,
                        help='Prefix for output bam files with UMIs trimmed. '
                        'Replaces {output_prefix}_UMI_trimmed.bam if used.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Define output path
    if args.output_prefix is not None:
        output = args.output_prefix
    else:
        output = args.input_bam.replace('.bam', '')
    # Ensure directory exists
    output_dir = os.path.split(output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    out_bam_path = f'{output}_UMI_trimmed.bam'
    if args.output_bam is not None:
        out_bam_path = args.output_bam
        bam_out_dir = os.path.split(out_bam_path)[0]
        if len(bam_out_dir) > 0 and not os.path.isdir(bam_out_dir):
            os.mkdir(bam_out_dir)
        
    # Retrieve UMIs from reads, along with reads statistics
    utils.send_text(f'Retrieving raw UMIs from {args.input_bam}', v, 1, 0)
    umis = retrieve_umis(
        args.input_bam,
        args.pattern,
        out_bam_path, 
        f'{output}_raw_umis.tsv',
        v
    )
    # Save UMI and reads data json for future time save
    utils.send_text('Saving raw UMIs json', v, 1, 0)
    utils.save_json(f'{output}_raw_umis.json', umis)

    return 0


if __name__ == "__main__":
    sys.exit(main())
