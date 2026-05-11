#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Cluster reads associated to different UMIs using UMI sequence and alignment.
Use alignment- and coverage-based distances between groups of reads / UMIs and
potentially merge groups of reads considering they come from the same original
molecule and ended up having slightly different UMIs because of polymerisation
or sequencing errors.
"""


import argparse
import os
import pickle
import sys
import time

import networkx as nx
import numpy as np

import utils



def cluster_umis(umi_names, umis, alignment_mat, overlap_mat,
                 alignment_threshold=3, overlap_threshold=0.75,
                 umi_stats_path=None, reads_stats=None, groups_path=None, v=0):
    """Cluster UMIs based on sequence and coverage proximity

    Arguments:
    umi_names  (list of str) - List of UMI names (used to preserve UMIs order).
    umis              (dict) - Dictionnary containing the UMIs as keys and
                               the ids of the reads associated with each UMI.
    alignment_mat (np.array) - Sequence-alignment based distance matrix.
    overlap_mat   (np.array) - Coverage-overlap based distance matrix.
    alignment_threshold (float) - Threshold for alignment-based distance
                                  (these distances are typically positive
                                  integers) (default: 3).
    overlap_threshold   (float) - Threshold for overlap-based distance
                                  (these distances are typically float
                                  comprised between 0 and 1) (default: 0.7).
    umi_stats_path        (str) - Path to tsv file with umis statistics
                                  (optional).
    reads_stats          (dict) - Dictionnary containing the reads ids as keys
                                  and alignment statistics of the reads.
    groups_path           (str) - Path to tsv file with per-read umis
                                  statistics (optional).
    v                     (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    clustered_umis   (np.array) - Updated dictionnary containing the UMIs.
    umis_graph (networkx graph) - Graph of UMIs with edges connecting clustered
                                  UMIs. 
    """
    # Initialize graph
    n = len(umis)
    umis_graph = nx.Graph()
    # Add nodes with their raw read count as an attribute
    for umi in umi_names:
        umis_graph.add_node(umi, n_reads=len(umis[umi]))

    if v > 0:
        t_zero = time.perf_counter()

    # Connect UMIs based on pairwise distances and thresholds
    # if they are aligned to the same reference contig (or both unaligned).
    # Store the distances as edge attributes for further statistics computing.
    umi_refs = {umi: umi.split('_')[1] for umi in umi_names}
    for i in range(n):
        for j in range(i + 1, n):
            if umi_refs[umi_names[i]] != umi_refs[umi_names[j]]:
                continue
            a_dist = alignment_mat[i, j]
            o_dist = overlap_mat[i, j]
            if a_dist <= alignment_threshold:
                if o_dist <= (1 - overlap_threshold):
                    umis_graph.add_edge(
                        umi_names[i], umi_names[j],
                        alignment_distance=int(a_dist),
                        overlap_distance=float(o_dist),
                    )
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Connected UMIs using distance in {elapsed:.3f} s.',
            v, 1, 1
        )
    
    # Cluster UMIs
    clustered_umis = {}
    old_to_new = {}
    for umis_cluster in nx.connected_components(umis_graph):
        # Retain name of the UMI with the highest number of associated reads
        nb_reads = [(len(umis[umi]), umi) for umi in umis_cluster]
        representative = max(nb_reads)[1]
        clustered_umis[representative] = []
        for umi in umis_cluster:
            clustered_umis[representative] += umis[umi]
            old_to_new[umi] = representative
        
    # Write out UMIs statistics
    if umi_stats_path is not None:
        with open(umi_stats_path, 'w') as tsv_out:
            tsv_out.write(
                'raw_umi\tref\tnb_raw_reads\tnb_final_reads\tfinal_umi\n'
            )
            for umi in umi_names:
                ref = umi.split('_')[1]
                rpz = old_to_new[umi]
                n_ini = len(umis[umi])
                n_final = len(clustered_umis[rpz])
                tsv_out.write(f'{umi}\t{ref}\t{n_ini}\t{n_final}\t{rpz}\n')
    # Write out groups umi-tools-like
    if groups_path is not None:
        with open(groups_path, 'w') as tsv_out:
            tsv_out.write(
                'read_id\tcontig\tposition\tgene\tumi\tumi_count\tfinal_umi'
                + '\tfinal_umi_count\tunique_id\n'
            ) # write header
            for umi, read_ids in umis.items():
                for read_id in read_ids:
                    contig = reads_stats[read_id]['ref']
                    position = str(reads_stats[read_id]['ref_start'])
                    gene = 'NA'
                    umi_count = str(len(read_ids))
                    final_umi = old_to_new[umi]
                    final_umi_count = str(len(clustered_umis[final_umi]))
                    tsv_out.write(
                        '\t'.join([read_id, contig, position, gene, umi,
                                   umi_count, final_umi, final_umi_count,
                                   final_umi]) + '\n'
                    )

    return clustered_umis, umis_graph


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_prefix', type=str, required=True,
                        help='Prefix for output files.')
    parser.add_argument('-u', '--cc_umis', type=str, required=True,
                        help='Path to input json file containing coverage '
                        'connex UMIs information.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        'alignment statistics.')
    parser.add_argument('-am', '--dist_mat_alignment', type=str, required=True,
                        help='Path to alignment-based distance matrix pickle '
                        'file.')
    parser.add_argument('-om', '--dist_mat_overlap', type=str, required=True,
                        help='Path to overlap-based distance matrix pickle '
                        'file.')
    parser.add_argument('-at', '--alignment_threshold', type=int, default=3,
                        help='Two clusters of reads associated with two '
                        'distinct UMI sequences must have their UMIs distant '
                        'from at most this threshold to be merged. Using '
                        'Levenshtein distance, a substitution costs 2 and a '
                        'translation costs 2. Default value is 3, meaning that'
                        ' two clusters must have their UMIs differ from at '
                        'most 1 translation or 1 substitution to be able to '
                        'merge together. Keep in mind that overlap between the'
                        " two cluster's reads coverage profiles must be "
                        'sufficient for two clusters to be merged, and that '
                        'because clusters get merged from cluster to cluster, '
                        'two clusters having their UMIs differing from more '
                        'than this threshold can still end up merged '
                        'together.')
    parser.add_argument('-ot', '--overlap_threshold', type=float, default=0.75,
                        help='Two clusters of reads associated with two '
                        'distinct UMI sequences must see the overlap between '
                        'their two coverage profiles below this threshold to '
                        'be merged. Default value is 0.75, meaning that the '
                        'intersection between the two coverage profiles must '
                        'be greater than or equal to 3 fourths of the smallest'
                        ' of the two coverage silhouettes. Keep in mind that '
                        "the sequence similarity between the two cluster's "
                        'UMIs must be sufficient for two clusters to be '
                        'merged, and that because clusters get merged from '
                        'cluster to cluster, two clusters not overlapping '
                        'enough can still end up merged together.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output files do not start with an underscore
    output = args.output_prefix
    if output.endswith('/'):
        output += 'clust_out'
    # Ensure output directory exists
    output_dir = os.path.split(output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load reads statistics
    utils.send_text(f'Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)

    # Load coverage connex UMIs
    utils.send_text(f'Loading coverage connex UMIs', v, 1, 0)
    cc_umis = utils.load_json(args.cc_umis)
    cc_umis_list = list(cc_umis.keys())
    cc_umis_list.sort()

    # Load pairwise distances between UMIs (coverage connected read sets)
    utils.send_text(f'Loading distance matrices', v, 1, 0)
    alignment_dist_matrix = np.load(args.dist_mat_alignment)['align_dist_mat']
    overlap_dist_matrix = np.load(args.dist_mat_overlap)['overlap_dist_mat']

    # Cluster UMIs
    utils.send_text(f'Cluster UMIs using alignment/overlap distances', v, 1, 0)
    clustered_umis, umis_graph = cluster_umis(
        umi_names=cc_umis_list,
        umis=cc_umis,
        alignment_mat=alignment_dist_matrix,
        overlap_mat=overlap_dist_matrix,
        alignment_threshold=args.alignment_threshold,
        overlap_threshold=args.overlap_threshold,
        umi_stats_path=f'{output}_umis_stats.tsv',
        reads_stats=reads_stats,
        groups_path=f'{output}_groups.tsv',
        v=v
    )

    # Save clustered UMIs
    utils.send_text('Saving clustered UMIs', v, 1, 0)
    utils.save_json(f'{output}_clustered_umis.json', clustered_umis)

    # Save UMIs graph
    utils.send_text('Saving UMIs graph', v, 1, 0)
    with open(f'{output}_UMIs_graph.pkl', 'wb') as graph_pkl_file:
        pickle.dump(
            umis_graph,
            graph_pkl_file,
            pickle.HIGHEST_PROTOCOL
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
