#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Split reads sharing the same raw UMI in coverage-connected components.
From reads alignment information and associations between raw UMIs sequences
and reads ids, devise new UMIs groups by splitting previous groups in new,
smaller groups where all reads (i) share the exact same UMI sequence and (ii)
have a connex coverage profile when aligned to the reference genome.
"""


import argparse
import os
import sys
import time

import utils


def update_umis(old_umis, reads_stats):
    """Update UMIs with only reads in the provided statistics

    Arguments:
    old_umis    (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    
    Return:
    umis        (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    """
    umis = {}
    for umi, old_read_ids in old_umis.items():
        read_ids = [
            r_id for r_id in old_read_ids if r_id in reads_stats.keys()
        ]
        if len(read_ids) > 0:
            umis[umi] = read_ids
    return umis


def split_connex_coverage(umis, reads_stats, v=0):
    """Split UMIs according to associated reads coverage connected components
    Reads associated with a given UMIs are placed in a graph. Edges are drawn
    between reads if and only if they overlap (once aligned on the reference).
    If the graph has multiple connected components, the UMI is split in
    several sub read-lists. In practice, use a 1-D line rather than a graph for
    faster execution

    Arguments:
    umis        (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    v            (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    cc_umis     (dict) - Updated dictionnary containing the UMIs.
    """
    n = len(umis)
    # Initiate new UMIs dictionnary
    by_ref_umis = {}
    cc_umis = {}

    # Store times for progression status messages
    if v > 0:
        t_zero = time.perf_counter()
        query_count = 0
    
    # Split UMIs on references of their associated reads
    for umi, read_ids in umis.items():
        references = {}
        for read_id in read_ids:
            read_ref = reads_stats[read_id]['ref']
            if read_ref in references.keys():
                references[read_ref].append(read_id)
            else:
                references[read_ref] = [read_id]
        for ref_name, new_read_ids in references.items():
            by_ref_umis[f'{umi}_{ref_name}'] = new_read_ids
    by_ref_n = len(by_ref_umis)
    
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Split {n} raw UMIs in {by_ref_n} reference-aware UMIs in '
            + f'{elapsed:.3f} s.',
            v, 2, 1
        )
        t_zero = time.perf_counter()

    # Split UMIs on connected components
    for umi, read_ids in by_ref_umis.items():
        if umi.endswith('unaligned'):
            cc_umis[umi] = read_ids
            continue
        # Build reads interval list
        intervals = [
            (
                reads_stats[r]['ref_start'],
                reads_stats[r]['ref_start'] + reads_stats[r]['aligned_len'] -1,
                r,
            ) for r in read_ids if reads_stats[r]['ref'] != 'unaligned'
        ]
        intervals.sort()
        # Retrieve components by scanning intervals
        components = []
        current_comp = [intervals[0][2]]
        current_start, current_end = intervals[0][0], intervals[0][1]
        for start, end, read_id in intervals[1:]:
            if start <= current_end:  # overlap
                current_comp.append(read_id)
                current_end = max(current_end, end)
            else:
                components.append(current_comp)
                current_comp = [read_id]
                current_start, current_end = start, end
        components.append(current_comp)
        # If only one component keep umi unchanged
        if len(components) == 1:
            cc_umis[umi] = read_ids
        else:
            for i, component in enumerate(components):
                cc_umis[f'{umi}_{i+1}'] = component
        # TODO: Remove
        # if '9173c2f4-e258-4693-96b9-9c2caba906fe_1' in read_ids \
        #     or 'ccef40b9-ff6b-41ae-8351-d09b6efe0a16_1' in read_ids:
        #     print(components)
        #     print(intervals)
        #     print(umi)
    by_cc_n = len(cc_umis)
    
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Split {by_ref_n} ref-aware UMIs in {by_cc_n} coverage-connected '
            + f'UMIs in {elapsed:.3f} s.',
            v, 2, 1
        )
    
    return cc_umis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to output JSON file with UMIs/reads connex '
                        + 'clusters.')
    parser.add_argument('-r', '--raw_umis', type=str, required=True,
                        help='Path to input json file containing raw UMIs '
                        + 'information.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        + 'alignment statistics.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load raw UMIs
    utils.send_text(f'Loading raw UMIs statistics', v, 1, 0)
    raw_umis = utils.load_json(args.raw_umis)

    # Load reads statistics
    utils.send_text(f'Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)

    # Update UMIs to filter out reads that aren't in the statistics
    utils.send_text(f'Filtering out unmapped reads from UMIs', v, 1, 0)
    raw_umis = update_umis(raw_umis, reads_stats)
    
    # Divide UMI-associated read groups that have distinct connex components
    # in their coverage profile
    utils.send_text(f'Splitting umis on coverage connex components', v, 1, 0)
    cc_umis = split_connex_coverage(raw_umis, reads_stats, v)

    # Save updated UMIs
    utils.send_text('Saving split UMIs', v, 1, 0)
    utils.save_json(args.output, cc_umis)

    return 0


if __name__ == "__main__":
    sys.exit(main())
