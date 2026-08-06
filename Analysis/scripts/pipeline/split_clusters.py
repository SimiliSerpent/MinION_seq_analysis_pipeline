#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Split over-merged UMI clusters into sub-clusters using Leiden community
detection. Iterates over every connected component of the UMIs graph produced
by cluster_umis.py, re-partitions it with Leiden (standard resolution,
geometric-mean node-weight boost by default), shrinks each detected community
to its heaviest node (with tied-for-max nodes pooled), and discards
sub-clusters whose pooled read count falls below a user-defined threshold.
"""


import argparse
import math
import os
import pickle
import sys
import time

import networkx as nx

import utils
import detect_cluster_modules


def per_cluster_threshold(n_reads_total, floor, log_scale):
    """Per-component min-reads threshold.
    Returns ceil(floor + log_scale * log10(max(1, n_reads_total))), so the
    threshold rises (slowly) with the component's total read count. With
    log_scale=0 the result reduces to the constant `floor` for all
    components.

    Arguments:
    n_reads_total (int) - Total read count of the connected component.
    floor         (int) - Threshold floor (= threshold applied to tiny
                          components and the intercept of the log formula).
    log_scale   (float) - Multiplier on log10(n_reads_total). 0 = flat
                          threshold = floor. 1.0 = +1 to threshold per decade
                          of cluster size.
    
    Return:
    threshold     (int) - Per-component min-reads threshold.
    """
    if log_scale == 0:
        return floor
    return int(math.ceil(
        floor + log_scale * math.log10(max(1, n_reads_total))
    ))


def split_clusters(graph, cc_umis, method='leiden', resolution=1.0,
                   node_weight_mode='geomean', log_prominence=1.0,
                   alignment_threshold=3, overlap_threshold=0.75,
                   min_reads_floor=5, min_reads_log_scale=0.0, v=0):
    """Split every connected component of the UMIs graph into sub-clusters.
    Each component is re-partitioned with the chosen community-detection
    method (Leiden or density-peaks), every community shrinked to its
    heaviest node(s), and sub-clusters whose pooled read count falls below
    the per-component threshold (cf. per_cluster_threshold) are dropped.

    Arguments:
    graph    (networkx graph) - Full UMIs graph as output by cluster_umis.
    cc_umis            (dict) - Dictionnary containing the coverage-connex UMIs
                                as keys and the ids of the reads associated
                                with each UMI (= cluster_umis input, NOT
                                output: the per-node read lists are needed to
                                be re-pooled after community detection).
    method              (str) - Community detection method. 'leiden' uses
                                modularity-based detection with edge weights;
                                'density_peaks' uses a topographic / persistent
                                -homology approach with log-prominence
                                filtering on node read counts (default:
                                'leiden').
    resolution        (float) - Leiden RB resolution. Higher values yield more
                                (smaller) communities. Unused when
                                method='density_peaks' (default: 1.0).
    node_weight_mode    (str) - Multiplicative boost applied to edge weights
                                using endpoint node weights. See
                                detect_cluster_modules.node_weight_boost.
                                Unused when method='density_peaks'
                                (default: 'geomean').
    log_prominence    (float) - Log10-prominence threshold for density_peaks.
                                A peak survives if log10(peak_count /
                                saddle_count) >= threshold. Unused for leiden
                                (default: 1.0 = peak must be >=10x its saddle).
    alignment_threshold (int) - See detect_cluster_modules.edge_similarity
                                (default: 3). Should match the value used at
                                the clustering step. Unused for density_peaks.
    overlap_threshold (float) - See detect_cluster_modules.edge_similarity
                                (default: 0.75). Should match the value used
                                at the clustering step. Unused for
                                density_peaks.
    min_reads_floor     (int) - Threshold floor for the adaptive min-reads
                                filter (default: 5). See per_cluster_threshold.
    min_reads_log_scale (float) - Log-scale slope for the adaptive min-reads
                                  filter (default: 0.0 = flat threshold).
                                  See per_cluster_threshold.
    v                   (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    new_umis           (dict) - Dictionnary containing the sub-cluster
                                representatives as keys and the pooled ids of
                                the reads associated with each sub-cluster.
    stats_rows         (list) - One tuple per connected component:
                                (orig_rep, n_nodes, n_reads, threshold_used,
                                n_communities, n_kept, n_reads_kept,
                                n_reads_dropped).
    """
    if v > 0:
        t_zero = time.perf_counter()

    new_umis = {}
    # Counters for the global summary
    n_cc = 0
    n_cc_singletons = 0
    n_communities_total = 0
    n_kept_total = 0
    n_reads_kept = 0
    n_reads_dropped = 0

    # Accumulator for stats
    stats_rows = []

    for component in nx.connected_components(graph):
        n_cc += 1
        sub = graph.subgraph(component).copy()
        n_reads_orig = sum(sub.nodes[n]['n_reads'] for n in sub)
        threshold = per_cluster_threshold(
            n_reads_orig, min_reads_floor, min_reads_log_scale
        )

        # Singleton component: no graph structure to leverage; either keep
        # the node as-is or drop it depending on the read count. Cheap fast
        # path - UMI data typically has many singleton components.
        if sub.number_of_nodes() == 1:
            n_cc_singletons += 1
            (node,) = sub.nodes
            n_reads = sub.nodes[node]['n_reads']
            n_communities_total += 1
            if n_reads >= threshold:
                new_umis[node] = list(cc_umis[node])
                n_kept_total += 1
                n_reads_kept += n_reads
                kept_here, reads_kept_here, reads_dropped_here = 1, n_reads, 0
            else:
                n_reads_dropped += n_reads
                kept_here, reads_kept_here, reads_dropped_here = 0, 0, n_reads
            stats_rows.append((
                node, 1, n_reads, threshold, 1,
                kept_here, reads_kept_here, reads_dropped_here,
            ))
            continue

        # Multi-node component: run Leiden community detection
        if method == 'leiden':
            communities = detect_cluster_modules.run_leiden(
                sub,
                alignment_threshold=alignment_threshold,
                overlap_threshold=overlap_threshold,
                resolution=resolution,
                node_weight_mode=node_weight_mode,
                v=0,  # mute per-component chatter; we summarize globally
            )
        elif method == 'density_peaks':
            communities, _persist = detect_cluster_modules.run_density_peaks(
                sub,
                log_prominence=log_prominence,
                v=0,  # mute per-component chatter; we summarize globally
            )
        else:
            raise ValueError(
                f"Unknown method: {method!r}. "
                "Must be 'leiden' or 'density_peaks'."
            )
        # Shrink each community to its representative node(s), pooling ties
        # (keep_tied=True: ambiguous max -> co-representatives, reads merged)
        reps_per_comm = detect_cluster_modules.restrict_communities_to_reps(
            communities, sub, keep_tied=True
        )

        n_communities_total += len(communities)
        kept_here = 0
        reads_kept_here = 0
        reads_dropped_here = 0
        for reps in reps_per_comm:
            pooled = []
            for r in reps:
                pooled.extend(cc_umis[r])
            n_reads = len(pooled)
            if n_reads < threshold:
                reads_dropped_here += n_reads
                continue
            # Pick the rep name deterministically: highest n_reads, ties
            # broken by node name (max() over the names tuple). Reads from
            # every tied co-rep are pooled under that single name.
            rep = max(reps, key=lambda n: (sub.nodes[n]['n_reads'], n))
            new_umis[rep] = pooled
            kept_here += 1
            reads_kept_here += n_reads

        n_kept_total += kept_here
        n_reads_kept += reads_kept_here
        n_reads_dropped += reads_dropped_here

        # Original "representative" of the component, à la cluster_umis
        orig_rep = max(
            sub.nodes,
            key=lambda n: (sub.nodes[n]['n_reads'], n),
        )
        stats_rows.append((
            orig_rep,
            sub.number_of_nodes(),
            n_reads_orig,
            threshold,
            len(communities),
            kept_here,
            reads_kept_here,
            reads_dropped_here,
        ))

    # Display processing time and global summary
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Split {n_cc} connected components ({n_cc_singletons} '
            + f'singletons) using {method} into {n_communities_total} '
            + f'communities, of which {n_kept_total} kept after min-reads '
            + f'filter (floor={min_reads_floor}, '
            + f'log_scale={min_reads_log_scale}: {n_reads_kept} reads kept, '
            + f'{n_reads_dropped} dropped) in {elapsed:.3f} s.',
            v, 1, 1
        )

    return new_umis, stats_rows


def plot_top_clusters(graph, cc_umis, reads_stats, fastq_path, top_rows,
                      plots_dir, method='leiden', resolution=1.0,
                      node_weight_mode='geomean', log_prominence=1.0,
                      alignment_threshold=3, overlap_threshold=0.75, v=0):
    """Generate detect_cluster_modules-style diagnostic plots for the top-N
    biggest components. For each row of `top_rows`, calls
    detect_cluster_modules.investigate_cluster with restrict_to_reps=True and
    keep_tied_reps=True, so the dissimilarity matrix and coverage profile
    reflect the rep-pooled reads that split_clusters actually keeps. When
    method='density_peaks', the plot additionally shows a persistence diagram
    (peak birth/death in log10 read-count space). Note that
    investigate_cluster does NOT apply the min_reads filter itself, so
    sub-threshold communities still appear on the plot - they are visible by
    their low read count in the legend, and cross-referenceable with the
    stats TSV.

    Arguments:
    graph    (networkx graph) - Full UMIs graph.
    cc_umis            (dict) - Coverage-connex UMIs (cluster_umis input).
    reads_stats        (dict) - Reads alignment statistics.
    fastq_path          (str) - Path to fastq file with raw reads (used for
                                SPOA consensus calling).
    top_rows           (list) - Rows from split_clusters' stats_rows, in the
                                order to be plotted. Only the first element
                                of each row (= original representative) and
                                the third (= n_reads, used for log progress)
                                are read.
    plots_dir           (str) - Path to output plots directory. Created if
                                it does not exist.
    method              (str) - Community detection method, forwarded to
                                investigate_cluster.
    resolution        (float) - Leiden RB resolution. Should match the value
                                used by split_clusters (default: 1.0).
    node_weight_mode    (str) - Edge-weight node boost mode. Unused for
                                density_peaks (default: 'geomean').
    log_prominence    (float) - Log10-prominence threshold for density_peaks.
                                Should match split_clusters. Unused for
                                leiden (default: 1.0).
    alignment_threshold (int) - Alignment-distance threshold for edge
                                similarity (default: 3). Unused for
                                density_peaks.
    overlap_threshold (float) - Overlap pass threshold for edge similarity
                                (default: 0.75). Unused for density_peaks.
    v                   (int) - Level of verbosity (default: 0 = muted)
    """
    os.makedirs(plots_dir, exist_ok=True)
    if v > 0:
        t_zero = time.perf_counter()

    for i, row in enumerate(top_rows):
        orig_rep = row[0]
        n_reads = row[2]
        # Zero-padded rank prefix so file listing sorts by size desc
        out_path = f'{plots_dir}/{i+1:03d}_{orig_rep}.png'
        utils.send_text(
            f'Plotting cluster {i+1}/{len(top_rows)} around {orig_rep} '
            + f'({n_reads} reads)',
            v, 2, 1
        )
        detect_cluster_modules.investigate_cluster(
            graph,
            cc_umis,
            reads_stats,
            fastq_path,
            orig_rep,
            out_path,
            method=method,
            alignment_threshold=alignment_threshold,
            overlap_threshold=overlap_threshold,
            resolution=resolution,
            log_prominence=log_prominence,
            node_weight_mode=node_weight_mode,
            restrict_to_reps=True,
            keep_tied_reps=True,
            v=0,  # mute investigate_cluster's per-step chatter
        )

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Plotted {len(top_rows)} top clusters in {elapsed:.3f} s.',
            v, 1, 1
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_prefix', type=str, required=True,
                        help='Prefix for output files.')
    parser.add_argument('-g', '--graph_pickle', type=str, required=True,
                        help='Path to input pickle file containing the UMIs '
                        'clustering graph.')
    parser.add_argument('-u', '--cc_umis', type=str, required=True,
                        help='Path to input json file containing coverage '
                        'connex UMIs information.')
    parser.add_argument('-s', '--reads_stats', type=str, default=None,
                        help='Path to input json file containing reads '
                        'alignment statistics. Required only when --plots_dir '
                        'is set.')
    parser.add_argument('-f', '--input_fastq', type=str, default=None,
                        help='Path to input fastq file with raw reads, used '
                        'to build SPOA community consensuses for the top-N '
                        'plots. Required only when --plot_top_n is set.')
    parser.add_argument('-m', '--method', type=str, default='leiden',
                        choices=['leiden', 'density_peaks'],
                        help='Community detection method used to split each '
                        'connected component. leiden uses modularity-based '
                        'detection with edge weights; density_peaks uses a '
                        'topographic / persistent-homology approach on node '
                        'read counts filtered by log-prominence (see '
                        '--log_prominence). Default: leiden.')
    parser.add_argument('-lp', '--log_prominence', type=float, default=1.0,
                        help='Log10-prominence threshold for density_peaks. '
                        'A peak survives if log10(peak_count / saddle_count) '
                        '>= threshold. 1.0 means "peak must be >=10x its '
                        'saddle"; 0.5 means ">=~3x"; 2.0 means ">=100x". '
                        'Unused when --method is leiden. Default: 1.0.')
    parser.add_argument('-at', '--alignment_threshold', type=int, default=3,
                        help='Alignment-distance threshold used to normalize '
                        'edge similarity weights for Leiden. Should match the '
                        'value passed to cluster_umis to reflect the edge '
                        'geometry actually used during clustering. '
                        'Default: 3.')
    parser.add_argument('-ot', '--overlap_threshold', type=float, default=0.75,
                        help='Overlap pass threshold (similarity, not '
                        'distance), used to normalize edge similarity weights '
                        'for Leiden. Should match the value passed to '
                        'cluster_umis. Default: 0.75.')
    parser.add_argument('--resolution', type=float, default=1.0,
                        help='Leiden RB resolution parameter. Higher values '
                        'yield more (smaller) communities. Default: 1.0.')
    parser.add_argument('-nwm', '--node_weight_mode', type=str,
                        default='geomean',
                        choices=['none', 'log_max', 'geomean', 'asymmetric'],
                        help='Multiplicative boost applied to edge weights '
                        'using endpoint node weights. none = read counts '
                        'information unused; log_max = log10(max(c_u, c_v) + '
                        '1), boosts edges touching any heavy node; geomean = '
                        'sqrt(c_u * c_v), aggressive; asymmetric = 1 + log10('
                        "max/min), continuous analog of umi-tools' 2N-1 rule. "
                        'Default: geomean.')
    parser.add_argument('-mr', '--min_reads', type=int, default=5,
                        help='Floor of the adaptive min-reads threshold (= '
                        'threshold applied to tiny components, and the '
                        'intercept of the log formula): sub-clusters whose '
                        'pooled read count is strictly below this threshold '
                        'are discarded. With --min_reads_log_scale 0 (the '
                        'default), the threshold is flat == this value. '
                        'Default: 5.')
    parser.add_argument('-mls', '--min_reads_log_scale', type=float,
                        default=0.0,
                        help='Log-scale slope for the adaptive min-reads '
                        'threshold: per-component threshold = ceil('
                        'min_reads + log_scale * log10(component_total_reads))'
                        '. 0.0 (default) gives a flat threshold == min_reads. '
                        '1.0 adds +1 to the threshold per decade of cluster '
                        'size (e.g. floor=2, scale=1 -> threshold of 3 at '
                        '10 reads, 5 at 1,000, 7 at 100,000). Default: 0.0.')
    parser.add_argument('-n', '--plot_top_n', type=int, default=0,
                        help='Write a per-component splitting summary TSV for '
                        'the N biggest (by read count) original connected '
                        'components. 0 = no stats written. Default: 0.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Validate plot-related args
    if args.plot_top_n > 0:
        missing = [name for name, val in [
            ('--reads_stats', args.reads_stats),
            ('--input_fastq', args.input_fastq),
        ] if val is None]
        if missing:
            parser.error(
                f'--plot_top_n > 0 requires the following arg(s) to be set: '
                + ', '.join(missing)
            )

    # Ensure output files do not start with an underscore
    output = args.output_prefix
    if output.endswith('/'):
        output += 'split_out'
    # Ensure output directory exists
    output_dir = os.path.split(output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load graph
    utils.send_text(f'Loading UMIs graph', v, 1, 0)
    with open(args.graph_pickle, 'rb') as graph_pkl_file:
        graph = pickle.load(graph_pkl_file)

    # Load coverage connex UMIs
    utils.send_text(f'Loading coverage connex UMIs', v, 1, 0)
    cc_umis = utils.load_json(args.cc_umis)

    # Split over-merged clusters using community detection
    utils.send_text(
        f'Splitting clusters using {args.method} community detection',
        v, 1, 0
    )
    new_umis, stats_rows = split_clusters(
        graph=graph,
        cc_umis=cc_umis,
        method=args.method,
        resolution=args.resolution,
        node_weight_mode=args.node_weight_mode,
        log_prominence=args.log_prominence,
        alignment_threshold=args.alignment_threshold,
        overlap_threshold=args.overlap_threshold,
        min_reads_floor=args.min_reads,
        min_reads_log_scale=args.min_reads_log_scale,
        v=v,
    )

    # Save updated clustered UMIs (drop-in for deduplicate_umis.py)
    utils.send_text('Saving split clustered UMIs', v, 1, 0)
    utils.save_json(f'{output}_clustered_umis.json', new_umis)

    # Sort & write stats, optionally plot top-N
    # Sort by original cluster read count desc, keep top N in top_rows
    stats_rows.sort(key=lambda row: row[2], reverse=True)
    top_rows = stats_rows[:args.plot_top_n]

    stats_path = f'{output}_splitting_stats.tsv'
    with open(stats_path, 'w') as tsv_out:
        tsv_out.write(
            'original_representative\tn_nodes\tn_reads\tthreshold_used\t'
            + 'n_communities\tn_communities_kept\tn_reads_kept\t'
            + 'n_reads_dropped\n'
        )
        for row in stats_rows:
            tsv_out.write('\t'.join(str(x) for x in row) + '\n')
    utils.send_text(
        f'Wrote splitting stats to {stats_path}',
        v, 1, 0
    )

    if args.plot_top_n > 0:
        utils.send_text(
            f'Loading reads statistics (for top-N plots)', v, 1, 0
        )
        reads_stats = utils.load_json(args.reads_stats)
        utils.send_text(
            f'Plotting top-{args.plot_top_n} components in '
            + f'{output}_top_{args.plot_top_n}_clusters',
            v, 1, 0
        )
        plot_top_clusters(
            graph=graph,
            cc_umis=cc_umis,
            reads_stats=reads_stats,
            fastq_path=args.input_fastq,
            top_rows=top_rows,
            plots_dir=f'{output}_top_{args.plot_top_n}_clusters',
            method=args.method,
            resolution=args.resolution,
            node_weight_mode=args.node_weight_mode,
            log_prominence=args.log_prominence,
            alignment_threshold=args.alignment_threshold,
            overlap_threshold=args.overlap_threshold,
            v=v,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
