#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot summary statistics for a single UMI cluster.
Takes as input a graph, a cluster representative UMI identifier, the coverage-
connex UMIs, and the reads statistics, and plots side-by-side (i) the
histogram of the number of reads per unique raw UMI sequence in that cluster
(each unique UMI sequence is binned by the sum of reads of all nodes bearing
that sequence) and (ii) the coverage depth profile of the cluster on its
reference contig.
"""
 
 
import argparse
import os
import pickle
import sys
import time
 
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
 
import utils
 
 
def find_cluster(graph, representative, v=0):
    """Retrieve the connected component containing the given representative.
    Warns if the supplied node is not the canonical representative of its
    component (i.e. not the node with the highest n_reads).
 
    Arguments:
    graph (networkx graph) - Networkx graph object.
    representative   (str) - Name of the cluster representative UMI node.
    v                (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    sub  (networkx graph) - Subgraph of the connected component containing
                            the representative.
    """
    if representative not in graph:
        raise ValueError(
            f'Representative {representative} not found in graph nodes.'
        )
    component = nx.node_connected_component(graph, representative)
    sub = graph.subgraph(component)
 
    # Warn (don't fail) if the supplied node isn't the canonical rep -- the
    # caller may legitimately want to inspect a non-rep node's component
    canonical = max(sub.nodes, key=lambda n: sub.nodes[n]['n_reads'])
    if canonical != representative:
        utils.send_text(
            f'Warning: {representative} is not the canonical representative '
            + f'of its component (which is {canonical}). Plotting anyway.',
            v, 1, 1
        )
 
    return sub
 
 
def compute_unique_umi_reads(sub):
    """Sum reads across nodes sharing the same raw UMI sequence.
    Each node name follows the pattern '<umi_seq>_<ref>[_<cc_idx>]', so the
    UMI sequence is the first underscore-separated token. Multiple nodes may
    share a UMI sequence when the same raw UMI got split during the
    deduplication pipeline (different reference contigs, or multiple coverage
    connected components on the same contig).
 
    Arguments:
    sub  (networkx graph) - Connected (sub-)graph.
    
    Return:
        umi_reads  (dict) - Dict mapping each unique raw UMI sequence to the
                            total number of reads bearing it within the
                            cluster.
    """
    umi_reads = {}
    for node in sub:
        umi_seq = node.split('_')[0]
        umi_reads[umi_seq] = umi_reads.get(umi_seq, 0) \
            + sub.nodes[node]['n_reads']
    return umi_reads
 
 
def plot_unique_umis_histogram(umi_reads, ax, *, title_font=12,
                               label_font=10):
    """Plot the histogram of reads-per-unique-UMI for a cluster.
    Each unique UMI sequence is placed in the bin corresponding to its total
    read count in the cluster. Integer bins below 30, log-spaced bins above
    so very wide distributions stay legible.
 
    Arguments:
    umi_reads   (dict) - Dict mapping each unique raw UMI sequence to the
                        total number of reads bearing it.
    ax      (plt.Axes) - Matplotlib subplot.
    title_font (float) - Size of title.
    label_font (float) - Size of axis labels font.
    """
    counts = list(umi_reads.values())
    n_umis = len(counts)
    max_count = max(counts) if counts else 1
    total_reads = sum(counts)
 
    # Integer bins for narrow distributions, log-spaced bins for wide ones.
    # The threshold (30) is small enough that we get one bin per integer for
    # any realistic "small cluster" view, and large enough that the log
    # branch only kicks in when linear bins would crush the tail.
    if max_count <= 30:
        bins = np.arange(0.5, max_count + 1.5, 1)
    else:
        bins = np.logspace(0, np.log10(max_count + 1), 30)
        ax.set_xscale('log')
 
    ax.hist(counts, bins=bins, color='#4C9CD6', edgecolor='#333333',
            alpha=0.85)
    ax.set_xlabel('Reads per unique UMI sequence', fontsize=label_font)
    ax.set_ylabel('Number of unique UMI sequences', fontsize=label_font)
    ax.set_title(
        f'Reads per unique UMI ({n_umis} unique sequences, '
        + f'{total_reads} reads total)',
        fontsize=title_font
    )
    ax.grid(True, alpha=0.3)
 
    return 0
 
 
def compute_coverage_profile(sub, umis, reads_stats):
    """Compute the coverage depth profile of all aligned reads in the cluster.
    Unaligned reads are skipped silently.
 
    Arguments:
    sub  (networkx graph) - Connected (sub-)graph.
    umis            (dict) - Dictionnary containing the UMIs as keys and the
                             ids of the reads associated with each UMI.
    reads_stats     (dict) - Dictionnary containing the reads ids as keys and
                             alignment statistics of the reads.
    
    Return:
    positions (np.array) - Reference positions covered by the cluster reads
                           (empty if no aligned read).
    depth     (np.array) - Coverage depth at each position.
    """
    # Collect aligned read intervals
    intervals = []
    for node in sub:
        for read_id in umis[node]:
            stats = reads_stats[read_id]
            if stats['ref'] == 'unaligned':
                continue
            start = stats['ref_start']
            end = start + stats['aligned_len']
            intervals.append((start, end))
 
    if not intervals:
        return np.array([]), np.array([])
 
    # Build coverage profile by stacking reads on a numpy array. The window
    # is just [min_start, max_end), which keeps the array small enough for
    # any realistic capture-sequencing cluster.
    min_pos = min(s for s, _ in intervals)
    max_pos = max(e for _, e in intervals)
    depth = np.zeros(max_pos - min_pos, dtype=int)
    for start, end in intervals:
        depth[start - min_pos:end - min_pos] += 1
    positions = np.arange(min_pos, max_pos)
 
    return positions, depth
 
 
def plot_coverage_profile(positions, depth, ax, ref='', *, title_font=12,
                          label_font=10):
    """Plot the coverage depth profile of a cluster's reads.
 
    Arguments:
    positions (np.array) - Reference positions covered by the cluster reads.
    depth     (np.array) - Coverage depth at each position.
    ax         (plt.Axes) - Matplotlib subplot.
    ref            (str) - Reference contig name, used in axis label.
    title_font   (float) - Size of title.
    label_font   (float) - Size of axis labels font.
    """
    if len(positions) == 0:
        ax.text(
            0.5, 0.5, 'No aligned reads in cluster',
            ha='center', va='center', transform=ax.transAxes,
            fontsize=label_font
        )
        ax.set_axis_off()
        return 0
 
    ax.fill_between(positions, depth, color='#4C9CD6', alpha=0.6,
                    step='mid')
    ax.plot(positions, depth, color='#1F5F8B', linewidth=0.6)
    ax.set_xlabel(
        f'Position on {ref}' if ref else 'Reference position',
        fontsize=label_font
    )
    ax.set_ylabel('Depth', fontsize=label_font)
    breadth = int((depth > 0).sum())
    ax.set_title(
        f'Coverage profile (max depth: {int(depth.max())}, breadth: '
        + f'{breadth} bp)',
        fontsize=title_font
    )
    ax.grid(True, alpha=0.3)
    ax.set_xlim(positions.min(), positions.max())
    ax.set_ylim(bottom=0)
 
    return 0
 
 
def plot_cluster_summary(graph, umis, reads_stats, representative, out_path,
                         v=0):
    """Plot histogram of reads-per-unique-UMI and coverage profile for one
    cluster, side by side, in a single PNG file.
 
    Arguments:
    graph (networkx graph) - Networkx graph object.
    umis            (dict) - Dictionnary containing the UMIs as keys and the
                             ids of the reads associated with each UMI.
    reads_stats     (dict) - Dictionnary containing the reads ids as keys and
                             alignment statistics of the reads.
    representative   (str) - Name of the cluster representative UMI node.
    out_path         (str) - Path to output PNG file.
    v                (int) - Level of verbosity (default: 0 = muted)
    """
    if v > 0:
        t_zero = time.perf_counter()
 
    # Find cluster
    sub = find_cluster(graph, representative, v=v)
    ref = representative.split('_')[1]
    n_reads_total = sum(sub.nodes[n]['n_reads'] for n in sub)
 
    utils.send_text(
        f'Cluster contains {sub.number_of_nodes()} nodes, '
        + f'{sub.number_of_edges()} edges, {n_reads_total} reads on ref '
        + f'{ref}',
        v, 2, 1
    )
 
    # Compute unique UMI reads
    umi_reads = compute_unique_umi_reads(sub)
 
    # Compute coverage profile
    positions, depth = compute_coverage_profile(sub, umis, reads_stats)
 
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    plot_unique_umis_histogram(umi_reads, axes[0])
    plot_coverage_profile(positions, depth, axes[1], ref=ref)
 
    fig.suptitle(
        f'Cluster {representative}  -  '
        + f'{sub.number_of_nodes()} nodes, {n_reads_total} reads',
        fontsize=14
    )
    plt.tight_layout()
 
    utils.send_text('Saving cluster summary plot', v, 3, 1)
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
 
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Plotted cluster summary in {elapsed:.3f} s', v, 2, 1
        )
 
    return 0
 
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--representative', type=str, required=True,
                        help='Name of the cluster representative UMI node '
                        '(e.g. <umi_seq>_<ref>[_<cc_idx>]). Must be a node '
                        'present in the graph; ideally the canonical '
                        'representative of its connected component (the '
                        'node with the highest n_reads). A warning is '
                        'emitted if it is not, but plotting proceeds.')
    parser.add_argument('-g', '--graph_pickle', type=str, required=True,
                        help='Path to input pickle file containing the UMIs '
                        'clustering graph.')
    parser.add_argument('-u', '--cc_umis', type=str, required=True,
                        help='Path to input json file containing coverage '
                        'connex UMIs information.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        'alignment statistics.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output PNG file.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')
 
    args = parser.parse_args()
    v = args.verbose
 
    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)
 
    # Load graph
    utils.send_text(f'Loading UMIs graph', v, 1, 0)
    with open(args.graph_pickle, 'rb') as graph_pkl_file:
        graph = pickle.load(graph_pkl_file)
 
    # Load coverage connex UMIs
    utils.send_text(f'Loading coverage connex UMIs', v, 1, 0)
    cc_umis = utils.load_json(args.cc_umis)
 
    # Load reads statistics
    utils.send_text(f'Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)
 
    # Plot cluster summary
    utils.send_text(
        f'Plotting cluster summary for {args.representative}', v, 1, 0
    )
    res = plot_cluster_summary(
        graph,
        cc_umis,
        reads_stats,
        args.representative,
        args.output,
        v=v
    )
 
    return 0
 
 
if __name__ == "__main__":
    sys.exit(main())
