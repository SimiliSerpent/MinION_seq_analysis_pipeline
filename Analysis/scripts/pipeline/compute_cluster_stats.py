#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute statistics on connected subgraph from a graph of UMIs.
"""


import argparse
import os
import pickle
import sys
import time

import networkx as nx

import utils


def weighted_diameter(graph, weight):
    """Compute the longest weighted shortest path in a graph.
    
    Arguments:
    graph (networkx graph) - Networkx graph object.
    weight           (str) - Type of weight to use to compute diameter.
    
    Return:
    diameter   (int|float) - Longest weighted shortest path, or +inf if
                             disconnected, or 0 if fewer than 2 nodes.
    """
    if graph.number_of_nodes() < 2:
        return 0
    if not nx.is_connected(graph):
        return float('inf')
    diameter = 0
    for source in graph:
        lengths = nx.single_source_dijkstra_path_length(
            graph, source, weight=weight
        )
        if lengths:
            diameter = max(diameter, max(lengths.values()))
    return diameter


def compute_cluster_stats(graph, umis, reads_stats, out_path, v=0):
    """Compute per-cluster (connected component) statistics from the UMI
    graph and write them to a TSV.

    Cluster stats:
      - representative, reference
      - n_nodes, n_edges
      - n_reads_total
      - n_unique_umis              (nb of unque split UMIs in the cluster)
      - rep_reads_proportion       (reads in representative / total)
      - avg_aligned_len            (avg alignments length for cluster's reads)
      - coverage_breadth           (bp on reference covered by any read)
      - coverage_avg_depth         (sum of read lengths / breadth)
      - n_centroids                (nodes whose n_reads >= all neighbours)
      - alignment_diameter         (weighted by alignment_distance)
      - overlap_diameter           (weighted by overlap_distance)
      - avg_degree
      - density                    (2|E| / n(n-1) if n>1 else 0)

    Arguments:
    graph (networkx graph) - Networkx graph object.
    umis            (dict) - Dictionnary containing the UMIs as keys and the
                             ids of the reads associated with each UMI.
    reads_stats     (dict) - Dictionnary containing the reads ids as keys and
                             alignment statistics of the reads.
    out_path         (str) - Path to tsv file with clusters statistics.
    v                (int) - Level of verbosity (default: 0 = muted)
    """
    if v > 0:
        t_zero = time.perf_counter()

    headers = [
        'representative', 'reference', 'n_nodes', 'n_edges', 'n_reads_total',
        'n_unique_umis', 'rep_reads_proportion', 'avg_aligned_len',
        'coverage_breadth', 'coverage_avg_depth', 'n_centroids',
        'alignment_diameter', 'overlap_diameter', 'avg_degree', 'density'
    ]

    components = sorted(
        nx.connected_components(graph),
        key=lambda comp: sum(graph.nodes[n]['n_reads'] for n in comp),
        reverse=True,
    )

    with open(out_path, 'w') as tsv_out:
        tsv_out.write('\t'.join(headers) + '\n')

        for component in components:
            sub = graph.subgraph(component)
            n_nodes = sub.number_of_nodes()
            n_edges = sub.number_of_edges()

            # Read-count stats
            node_reads = {n: sub.nodes[n]['n_reads'] for n in sub}
            total_reads = sum(node_reads.values())
            representative = max(node_reads, key=node_reads.get)
            rep_prop = node_reads[representative] / total_reads
            ref = representative.split('_')[1]

            # Unique raw UMI sequences in the cluster. Each node name is
            # '<umi_seq>_<ref>[_<cc_idx>]', so the UMI sequence is the first
            # token.
            n_unique_umis = len({node.split('_')[0] for node in sub})

            # Centroids: n_reads > all direct neighbours 
            n_centroids = sum(
                1 for node in sub
                if all(node_reads[node] > node_reads[nb]
                       for nb in sub.neighbors(node))
            )
            n_centroids = max(n_centroids, 1)

            # Diameters (weighted shortest paths)
            align_diam = weighted_diameter(sub, weight='alignment_distance')
            ovlp_diam = weighted_diameter(sub, weight='overlap_distance')

            # Average degree and density
            if n_nodes > 0:
                avg_degree = sum(d for _, d in sub.degree()) / n_nodes
            else:
                avg_degree = 0.0
            if n_nodes > 1:
                density = 2 * n_edges / (n_nodes * (n_nodes - 1))
            else:
                density = 0.0

            # Coverage breadth / avg depth from all cluster reads
            intervals = []
            read_lengths = []
            for node in sub:
                for read_id in umis[node]:
                    read_stat = reads_stats[read_id]
                    if read_stat['ref'] == 'unaligned':
                        continue
                    intervals.append(
                        (read_stat['ref_start'],
                         read_stat['ref_start'] + read_stat['aligned_len'])
                    )
                    read_lengths.append(read_stat['aligned_len'])
            if intervals:
                intervals.sort()
                merged = [list(intervals[0])]
                for start, end in intervals[1:]:
                    if start <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], end)
                    else:
                        merged.append([start, end])
                breadth = sum(e - s for s, e in merged)
                total_bases = sum(e - s for s, e in intervals)
                avg_depth = total_bases / breadth if breadth > 0 else 0.0
            else:
                breadth = 0
                avg_depth = 0.0
            avg_aligned_len = (sum(read_lengths) / len(read_lengths)
                               if read_lengths else 0.0)

            row = [
                representative, ref, n_nodes, n_edges, total_reads,
                n_unique_umis, f'{rep_prop:.4f}', f'{avg_aligned_len:.1f}',
                breadth, f'{avg_depth:.3f}', n_centroids, align_diam,
                f'{ovlp_diam:.4f}', f'{avg_degree:.3f}', f'{density:.4f}'
            ]
            tsv_out.write('\t'.join(str(v) for v in row) + '\n')

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Computed cluster stats in {elapsed:.3f} s.', v, 2, 1
        )
    
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output TSV file containing UMIs '
                        'clusters statistics.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        'alignment statistics.')
    parser.add_argument('-u', '--cc_umis', type=str, required=True,
                        help='Path to input json file containing coverage '
                        'connex UMIs information.')
    parser.add_argument('-g', '--graph_pickle', type=str, required=True,
                        help='Path to input pickle file containing UMIs '
                        'clustering graph.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load reads statistics
    utils.send_text(f'Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)

    # Load coverage connex UMIs
    utils.send_text(f'Loading coverage connex UMIs', v, 1, 0)
    cc_umis = utils.load_json(args.cc_umis)

    # Load UMIs graph
    utils.send_text(f'Loading UMIs graphs', v, 1, 0)
    with open(args.graph_pickle, 'rb') as graph_pkl_file:
        umis_graph = pickle.load(graph_pkl_file)

    # Compute clusters statistics
    utils.send_text('Computing cluster statistics', v, 1, 0)
    res = compute_cluster_stats(
        umis_graph,
        cc_umis,
        reads_stats,
        args.output,
        v=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
