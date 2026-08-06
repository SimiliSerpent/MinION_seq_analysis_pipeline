#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Add a color-tag to BAM alignment records to visualize UMIs clustering.
"""


import argparse
import os
import sys
import time

import pysam

import utils


# Set to True if comparing with umi-tools to ensure identical bam input
TEST_COMPARE_UMITOOLS = False


def tag_bam_with_clusters(in_bam_path, out_bam_path, clustered_umis, v=0):
    """Colors reads with respect to their clustering.
    Write a copy of the input BAM where each primary read carries:
      - YC: hex color shared by all reads of the same cluster
      - XC: cluster representative UMI (string id)
    Open the output in IGV and choose 'Color alignments by → tag → YC'
    (or group by XC).

    Arguments:
    in_bam_path     (str) - Path to input BAM file.
    out_bam_path    (str) - Path to output BAM file.
    clustered_umis (dict) - Dictionnary containing the UMIs as keys and the
                            ids of the reads associated with each UMI.
    v               (int) - Level of verbosity (default: 0 = muted)
    """
    # read_id -> representative UMI
    read_to_cluster = {}
    for rep, read_ids in clustered_umis.items():
        for r in read_ids:
            read_to_cluster[r] = rep

    # Stable, well-spread palette by hashing the representative name
    def color_for(rep):
        h = hash(rep) & 0xFFFFFF
        r = (h >> 16) & 0xFF
        g = (h >> 8) & 0xFF
        b = h & 0xFF
        return f'{r},{g},{b}'
    cluster_color = {rep: color_for(rep) for rep in clustered_umis}

    if v > 0:
        t_zero = time.perf_counter()
    with pysam.AlignmentFile(in_bam_path, 'rb') as b_in, \
         pysam.AlignmentFile(out_bam_path, 'wb', header=b_in.header) as b_out:
        for query in b_in:
            if query.is_secondary or query.is_supplementary:
                continue
            # The bam was built before umi-trimming so query_name still has
            # the UMI suffix when TEST_COMPARE_UMITOOLS is True. Strip it the
            # same way retrieve_reads_stats does.
            name = query.query_name
            if TEST_COMPARE_UMITOOLS and name not in read_to_cluster:
                stripped = name.rsplit('_', 1)[0]
                if stripped in read_to_cluster:
                    name = stripped
            rep = read_to_cluster.get(name)
            if rep is not None:
                query.set_tag('XC', rep, value_type='Z')
                query.set_tag('YC', cluster_color[rep], value_type='Z')
            b_out.write(query)
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Wrote cluster-tagged BAM in {elapsed:.3f} s.', v, 2, 1
        )
    
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input_bam', type=str, required=True,
                        help='Path to input bam file with raw reads aligned to'
                        + ' the reference.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output bam file where records have a'
                        ' UMI cluster tag allowing colored cluster '
                        'visualization in IGV.')
    parser.add_argument('-u', '--clustered_umis', type=str, required=True,
                        help='Path to input json file containing clustered '
                        'UMIs information.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    # Load clustered UMIs
    utils.send_text(f'Loading clustered UMIs', v, 1, 0)
    clustered_umis =  utils.load_json(args.clustered_umis)

    # Write BAM file with cluster-tag for clusters alignment visualization
    utils.send_text('Writing cluster-tagged BAM', v, 1, 0)
    res = tag_bam_with_clusters(
        args.input_bam,
        args.output,
        clustered_umis,
        v=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
