#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Investigate and split a single UMI cluster using community detection.
Takes an existing cluster representative, retrieves its connected component
in the UMIs graph, and re-partitions it with one of three methods: Leiden
(requires igraph+leidenalg), Louvain (uses networkx's built-in), or
density-peaks with log-prominence filtering (a topographic / persistent-
homology approach that uses node weights as elevation). Builds a SPOA draft
consensus per community, computes pairwise consensus dissimilarities, and
produces a single PNG combining (i) the cluster graph coloured by
community, (ii) the stratified coverage profile, (iii) the inter-community
consensus dissimilarity heatmap, and (for density-peaks only) (iv) the
persistence diagram showing every peak's birth and death in log10-read-count
space.

Designed to be ran on hand-picked clusters to calibrate the splitting
behaviour: a "looks fine" cluster to confirm the method does not over-split,
and an "over-merged" cluster to confirm the method does split it.
"""


import argparse
import math
import os
import pickle
import random
import sys
import time

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pysam
import rapidfuzz.distance
import spoa

import utils
import deduplicate_umis


random.seed(0)
np.random.seed(0)


def edge_similarity(data, alignment_threshold, overlap_threshold):
    """Convert edge alignment+overlap distances to a single similarity weight.
    Edges very close to threshold contribute weight ~0; edges with both
    distances near zero contribute weight ~2. Used as edge weight for Leiden
    community detection: high-weight edges keep their endpoints together,
    low-weight edges (= boundary edges = bridges) are the natural cut points.

    Arguments:
    data               (dict) - Networkx edge data dictionnary.
    alignment_threshold (int) - Alignment-distance threshold used at the
                                clustering step, used here for normalization.
    overlap_threshold (float) - Overlap-distance threshold used at the
                                clustering step, used here for normalization.
    
    Return:
    weight            (float) - Edge weight in [0, 2].
    """
    align_sim = max(
        0.0, 1.0 - data['alignment_distance'] / alignment_threshold
    )
    overlap_sim = max(
        0.0, 1.0 - data['overlap_distance'] / overlap_threshold
    )
    return align_sim + overlap_sim


def node_weight_boost(c_u, c_v, mode):
    """Multiplicative boost based on edge endpoints' node weights.
    Applied on top of edge_similarity to incorporate node-weight information
    into edge weights for modularity-based community detection (Leiden,
    Louvain). The intent is to anchor true peaks' halos into communities by
    making edges around heavy nodes harder to cut, while keeping bridge
    chains between low-count nodes weak.

    Arguments:
    c_u, c_v (int) - Read counts at edge endpoints.
    mode     (str) - One of:
                     'none'       (boost = 1, current behaviour);
                     'log_max'    (boost = log10(max(c_u, c_v) + 1), boosts
                                   any edge touching a heavy node);
                     'geomean'    (boost = sqrt(c_u * c_v), aggressive, heavy-
                                   heavy edges become very heavy);
                     'asymmetric' (boost = 1 + log10(max/min), rewards
                                   "child of a heavy parent" edges, continuous
                                   analog of umi-tools' 2N-1 rule).
    
    Return:
    boost  (float) - Multiplicative weight boost in (0, +inf).
    """
    if mode == 'none':
        return 1.0
    if mode == 'log_max':
        return math.log10(max(c_u, c_v) + 1)
    if mode == 'geomean':
        return math.sqrt(max(1, c_u) * max(1, c_v))
    if mode == 'asymmetric':
        hi = max(c_u, c_v)
        lo = max(1, min(c_u, c_v))  # avoid div-by-zero / log of 0
        return 1.0 + math.log10(hi / lo)
    raise ValueError(
        f"Unknown node_weight_mode: {mode!r}. "
        "Must be 'none', 'log_max', 'geomean', or 'asymmetric'."
    )


def run_leiden(sub, alignment_threshold=3, overlap_threshold=0.75,
               resolution=1.0, node_weight_mode='none', seed=0, v=0):
    """Run Leiden community detection on a cluster subgraph.
    Networkx graphs are converted to igraph for leidenalg, then communities
    are translated back into lists of node names.

    Arguments:
    sub      (networkx graph) - Connected (sub-)graph.
    alignment_threshold (int) - See edge_similarity (default: 3).
    overlap_threshold (float) - See edge_similarity (default: 0.75).
    resolution        (float) - Leiden RB resolution. Higher values yield more
                                (smaller) communities. Default: 1.0.
    node_weight_mode    (str) - See node_weight_boost. Default 'none'.
    seed                (int) - Random seed for Leiden initialization
                                (default: 0).
    v                   (int) - Level of verbosity (default: 0 = muted).
    
    Return:
    communities (list of list of str) - Communities sorted by total reads
                                        descending.
    """
    # Convert networkx -> igraph (leidenalg consumes igraph)
    import igraph as ig
    import leidenalg
    sorted_nodes = sorted(sub.nodes)
    nx_to_ig = {n: i for i, n in enumerate(sorted_nodes)}
    ig_to_nx = {i: n for n, i in nx_to_ig.items()}
    edges = []
    weights = []
    for u, v_, data in sub.edges(data=True):
        edges.append((nx_to_ig[u], nx_to_ig[v_]))
        base = edge_similarity(data, alignment_threshold, overlap_threshold)
        boost = node_weight_boost(
            sub.nodes[u]['n_reads'],
            sub.nodes[v_]['n_reads'],
            node_weight_mode,
        )
        weights.append(base * boost)
    ig_graph = ig.Graph(
        n=sub.number_of_nodes(), edges=edges, directed=False
    )
    ig_graph.es['weight'] = weights

    # Run Leiden using RBConfigurationVertexPartition (the resolution-aware
    # quality function; resolution=1.0 reproduces standard modularity
    # weighting with the configuration null model).
    if v > 0:
        t_zero = time.perf_counter()
    partition = leidenalg.find_partition(
        ig_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed,
    )
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Leiden produced {len(partition)} communities in {elapsed:.3f} s '
            + f'(modularity = {partition.modularity:.4f})',
            v, 3, 2
        )

    # Convert back to lists of node names, sorted by total reads desc.
    communities = []
    for comm in partition:
        communities.append([ig_to_nx[i] for i in comm])
    communities.sort(
        key=lambda c: sum(sub.nodes[n]['n_reads'] for n in c),
        reverse=True
    )
    return communities


def run_louvain(sub, alignment_threshold=3, overlap_threshold=0.75,
                resolution=1.0, node_weight_mode='none', seed=0, v=0):
    """Run Louvain community detection on a cluster subgraph.
    Uses networkx's built-in implementation. The resolution parameter has the
    following convention: 1.0 = standard modularity, higher = more / smaller
    communities.

    Arguments:
    sub      (networkx graph) - Connected (sub-)graph.
    alignment_threshold (int) - See edge_similarity (default: 3).
    overlap_threshold (float) - See edge_similarity (default: 0.75).
    resolution        (float) - Louvain resolution. Higher values yield
                                more (smaller) communities. Default: 1.0.
    node_weight_mode    (str) - See node_weight_boost. Default 'none'.
    seed                (int) - Random seed for Louvain initialization
                                (default: 0).
    v                   (int) - Level of verbosity (default: 0 = muted).
    
    Return:
    communities (list of list of str) - Communities sorted by total reads
                                        descending.
    """
    # Build a weighted copy so we don't mutate the caller's graph
    weighted = nx.Graph()
    sorted_nodes = sorted(sub.nodes)
    weighted.add_nodes_from(sorted_nodes(data=True))
    for u, v_, data in sub.edges(data=True):
        base = edge_similarity(data, alignment_threshold, overlap_threshold)
        boost = node_weight_boost(
            sub.nodes[u]['n_reads'],
            sub.nodes[v_]['n_reads'],
            node_weight_mode,
        )
        weighted.add_edge(u, v_, weight=base * boost)

    if v > 0:
        t_zero = time.perf_counter()
    communities_raw = nx.community.louvain_communities(
        weighted, weight='weight', resolution=resolution, seed=seed
    )
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        modularity = nx.community.modularity(
            weighted, communities_raw, weight='weight',
            resolution=resolution
        )
        utils.send_text(
            f'Louvain produced {len(communities_raw)} communities in '
            + f'{elapsed:.3f} s (modularity = {modularity:.4f})',
            v, 3, 2
        )

    # Convert sets -> lists, sorted by total reads desc.
    communities = [list(c) for c in communities_raw]
    communities.sort(
        key=lambda c: sum(sub.nodes[n]['n_reads'] for n in c),
        reverse=True
    )
    return communities


def run_density_peaks(sub, log_prominence=1.0, v=0):
    """Density-peak community detection with log-prominence filtering.
    Treats node read count (on a log10 scale) as topographic elevation, builds
    the persistent-homology merge tree by sweeping nodes in decreasing
    elevation, and keeps only peaks whose log-prominence (= log10 of the ratio
    peak_count / saddle_count) exceeds the threshold. Low-prominence peaks
    are absorbed into the higher peak they merged into.

    Edge weights are not used: only the graph topology and the node weights
    matter for the watershed. This is by design — the failure mode this
    method targets is precisely the one where edge weights alone cannot tell
    a real molecule's error cloud from a spurious bridge between two clouds.

    Arguments:
    sub      (networkx graph) - Connected (sub-)graph.
    log_prominence    (float) - Log10-prominence threshold. A peak survives if
                                log10(peak_count / saddle_count) >= threshold.
                                1.0 means "peak must be >=10x its saddle";
                                0.5 means ">=~3x"; 2.0 means ">=100x".
                                Default: 1.0.
    v                   (int) - Level of verbosity (default: 0 = muted).
    
    Return:
    communities (list of list of str) - Communities sorted by total reads
                                        descending.
    persistence        (dict) - Persistence data for plotting:
                                'peak_birth' (peak -> log10 count of peak),
                                'peak_death' (peak -> log10 count of saddle at
                                              which it merged, or 0 for the
                                              global max),
                                'prominence' (peak -> birth - death),
                                'survivors' (set of peaks above threshold),
                                'peak_per_community' (list of surviving peaks
                                                      in the same order as
                                                      communities),
                                'log_prominence' (threshold used).
    """
    if v > 0:
        t_zero = time.perf_counter()

    # Log-scale elevation. max(1, ...) is defensive against count-0 nodes,
    # which shouldn't exist but cost nothing to guard against.
    elev = {n: math.log10(max(1, sub.nodes[n]['n_reads'])) for n in sub}
    # Tie-break by node name for determinism: with two nodes at the same
    # elevation, this picks a stable winner-takes-peak rather than relying on
    # dict insertion order.
    nodes_by_elev = sorted(sub.nodes, key=lambda n: (-elev[n], n))

    parent = {}            # union-find parent pointers
    peak_birth = {}        # peak node -> elevation at birth
    peak_death = {}        # peak node -> elevation at death (merge into elder)
    absorbed_into = {}     # killed peak -> peak that absorbed it
    initial_root = {}      # node -> peak it first joined (used for assignment)

    def find(x):
        # Path-compressing find: walks to root, then re-points every node
        # along the path directly at the root.
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    processed = set()
    for n in nodes_by_elev:
        # Find the unique components (roots) of already-processed neighbours
        nb_roots = set()
        for m in sub.neighbors(n):
            if m in processed:
                nb_roots.add(find(m))

        if not nb_roots:
            # No higher neighbour yet -> n is a new local maximum
            parent[n] = n
            peak_birth[n] = elev[n]
            initial_root[n] = n
        else:
            # Elder rule: the peak with highest birth survives, others die at
            # elev[n] (which is by construction the lowest elevation among the
            # nodes joining n's neighbourhood and is therefore the saddle).
            roots_sorted = sorted(
                nb_roots, key=lambda r: (peak_birth[r], r), reverse=True
            )
            surviving = roots_sorted[0]
            for r in roots_sorted[1:]:
                peak_death[r] = elev[n]
                absorbed_into[r] = surviving
                parent[r] = surviving
            parent[n] = surviving
            initial_root[n] = surviving

        processed.add(n)

    # The global maximum never dies; assign conventional death at log10(1) = 0
    # so its prominence equals log10(count) and it always wins the threshold
    # gate.
    for p in list(peak_birth):
        if p not in peak_death:
            peak_death[p] = 0.0

    prominence = {
        p: peak_birth[p] - peak_death[p] for p in peak_birth
    }

    # Survivors: peaks above the log-prominence threshold. The global maximum
    # is always retained, even when the user picks an aggressive threshold,
    # so that the function always returns at least one community.
    survivors = {p for p, prom in prominence.items() if prom >= log_prominence}
    global_max = max(peak_birth, key=peak_birth.get)
    survivors.add(global_max)

    def find_surviving(peak):
        # Walk up the absorbed_into chain until we hit a surviving peak.
        # By construction (global max is always a survivor), this terminates.
        while peak not in survivors:
            peak = absorbed_into[peak]
        return peak

    # Assign each node to its surviving peak via its initial root
    comm_by_peak = {p: [] for p in survivors}
    for n in sub:
        comm_by_peak[find_surviving(initial_root[n])].append(n)

    # Sort communities by total reads desc, matching leiden/louvain convention
    items = sorted(
        comm_by_peak.items(),
        key=lambda kv: sum(sub.nodes[node]['n_reads'] for node in kv[1]),
        reverse=True,
    )
    communities = [nodes for _peak, nodes in items]
    peak_per_community = [peak for peak, _nodes in items]

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Density-peaks produced {len(communities)} communities in '
            + f'{elapsed:.3f} s (from {len(peak_birth)} candidate peaks, '
            + f'{len(survivors)} surviving log-prominence >= '
            + f'{log_prominence})',
            v, 3, 2
        )

    persistence = {
        'peak_birth': peak_birth,
        'peak_death': peak_death,
        'prominence': prominence,
        'survivors': survivors,
        'peak_per_community': peak_per_community,
        'log_prominence': log_prominence,
    }

    return communities, persistence


def load_cluster_reads(fastq_path, read_ids_needed, reads_stats, v=0):
    """Load read sequences and qualities from fastq, applying same per-read
    treatment as deduplicate_umis (reverse-complement if reverse strand,
    poly-A trimming, GA-repeat filtering). Only reads whose ids appear in
    `read_ids_needed` are kept.

    Arguments:
    fastq_path           (str) - Path to fastq file.
    read_ids_needed (iterable) - Iterable of read ids to keep.
    reads_stats         (dict) - Dictionnary containing the reads ids as keys
                                 and alignment statistics of the reads.
    v                    (int) - Level of verbosity (default: 0 = muted).
    
    Return:
    seqs                (dict) - Dict mapping read id to processed sequence.
    quals               (dict) - Dict mapping read id to processed quality
                                 string.
    """
    seqs = {}
    quals = {}
    needed = set(read_ids_needed)
    if v > 0:
        t_zero = time.perf_counter()
    with pysam.FastxFile(fastq_path) as fastq_stream:
        for read in fastq_stream:
            if read.name not in needed:
                continue
            stats = reads_stats[read.name]
            if stats['orientation'] != 'reverse':
                seq = read.sequence
                qual = read.quality
            else:
                seq = utils.rev_comp(read.sequence)
                qual = read.quality[::-1]
            if stats['ref'] != 'unaligned':
                seq, qual = deduplicate_umis.trim_polyA(seq, qual)
            seq, qual = deduplicate_umis.filter_GA(seq, qual)
            seqs[read.name] = seq
            quals[read.name] = qual
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Loaded {len(seqs)} cluster reads from fastq in '
            + f'{elapsed:.3f} s', v, 3, 2
        )
    return seqs, quals


def build_community_consensuses(communities, umis, seqs, max_seed=100, v=0):
    """Build a SPOA draft consensus per community.
    Up to max_seed reads per community are used as SPOA seeds (matching the
    pattern used in deduplicate_umis to keep SPOA fast). No racon polishing
    here — the draft is enough for a between-community comparison given how
    few base differences we are looking for.

    Arguments:
    communities (list of list of str) - Communities (lists of node names).
    umis        (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    seqs        (dict) - Dict mapping read id to processed sequence.
    max_seed     (int) - Max number of reads to seed SPOA per community
                         (default: 100).
    v            (int) - Level of verbosity (default: 0 = muted).
    
    Return:
    consensuses (dict) - Dict mapping community index to consensus sequence.
    """
    consensuses = {}
    rng = random.Random(0)
    if v > 0:
        t_zero = time.perf_counter()
    for i, comm in enumerate(communities):
        # Gather all reads belonging to community
        read_ids = []
        for node in comm:
            read_ids.extend(umis[node])
        # Subsample to max_seed for SPOA
        if len(read_ids) > max_seed:
            seed_reads = rng.sample(read_ids, max_seed)
        else:
            seed_reads = read_ids
        sequences = [seqs[r] for r in seed_reads if r in seqs]
        if not sequences:
            consensuses[i] = ''
            continue
        try:
            consensus, _msa = spoa.poa(sequences)
            consensuses[i] = consensus if isinstance(consensus, str) else ''
        except Exception as e:
            utils.send_text(
                f'SPOA error on community {i}: {e}', v, 3, 2
            )
            consensuses[i] = ''
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Built {len(consensuses)} community consensuses in '
            + f'{elapsed:.3f} s', v, 3, 2
        )
    return consensuses


def compute_dissimilarity_matrix(consensuses):
    """Pairwise Levenshtein distance between consensuses, normalized by the
    mean of the two consensus lengths. Empty / failed consensuses get
    dissimilarity 1.

    Arguments:
    consensuses (dict) - Dict mapping community index to consensus sequence.
    
    Return:
    matrix (np.array) - Symmetric dissimilarity matrix, indexed by community
                        order in `consensuses`.
    """
    n = len(consensuses)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            seq_i = consensuses[i]
            seq_j = consensuses[j]
            if not seq_i or not seq_j:
                d = 1.0
            else:
                raw = rapidfuzz.distance.Levenshtein.distance(seq_i, seq_j)
                mean_len = (len(seq_i) + len(seq_j)) / 2
                d = raw / mean_len if mean_len > 0 else 1.0
            matrix[i, j] = d
            matrix[j, i] = d
    return matrix


def restrict_communities_to_reps(communities, sub, keep_tied=False):
    """Restrict each community to its representative node(s).
    For each community, retain only the heaviest node (when keep_tied=False)
    or all nodes tied for max n_reads (when keep_tied=True). Used as a
    sharpening post-processing step: a hard cut-off at the rep node avoids
    boundary-mushiness between communities at the cost of dropping reads.

    Arguments:
    communities (list of list of str) - Communities (lists of node names).
    sub      (networkx graph) - Connected (sub-)graph (used for n_reads).
    keep_tied          (bool) - If True, keep every node tied for max
                                n_reads. If False, keep only the first such
                                node (matches Python's max() on iteration
                                order, deterministic given a deterministic
                                upstream method). Default: False.

    Return:
    restricted (list of list of str) - Restricted communities (typically
                                       one node per community, more when
                                       keep_tied and there are ties).
    """
    restricted = []
    for comm in communities:
        max_reads = max(sub.nodes[n]['n_reads'] for n in comm)
        reps = [n for n in comm if sub.nodes[n]['n_reads'] == max_reads]
        if not keep_tied:
            reps = reps[:1]
        restricted.append(reps)
    return restricted


def plot_graph_by_community(sub, communities, ax, fig, *, cmap):
    """Draw the cluster graph with nodes coloured by Leiden community.

    Arguments:
    sub  (networkx graph) - Connected (sub-)graph.
    communities (list of list of str) - Leiden communities.
    ax        (plt.Axes) - Matplotlib subplot.
    fig     (plt.Figure) - Matplotlib figure (kept for symmetry with other
                           plotting helpers; not used here).
    cmap (plt.Colormap) - Discrete colormap used to color communities.
    """
    pos = nx.spring_layout(sub, weight=None, seed=0)
    node_reads = {n: sub.nodes[n]['n_reads'] for n in sub}
    max_reads = max(node_reads.values())

    nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.45, width=0.8)

    n_comms = len(communities)
    colors = [cmap(i % cmap.N) for i in range(n_comms)]
    for i, comm in enumerate(communities):
        # Sort by reads asc so heavies end up on top of singletons
        comm_sorted = sorted(comm, key=lambda n: node_reads[n])
        sizes = [
            40 + 200 * (node_reads[n] / max_reads) for n in comm_sorted
        ]
        n_reads_comm = sum(node_reads[n] for n in comm)
        nx.draw_networkx_nodes(
            sub, pos, ax=ax,
            nodelist=comm_sorted,
            node_size=sizes,
            node_color=[colors[i]] * len(comm_sorted),
            edgecolors='#333333', linewidths=0.5,
            label=f'C{i}: {len(comm)}n / {n_reads_comm}r',
        )

    # Label the heaviest node of each community
    labels = {}
    for comm in communities:
        rep = max(comm, key=lambda n: node_reads[n])
        labels[rep] = node_reads[rep]
    nx.draw_networkx_labels(
        sub, pos, labels=labels, ax=ax, font_size=8
    )

    ax.legend(fontsize=7, loc='best')
    ax.set_title('Cluster graph coloured by community',
                 fontsize=11)
    ax.set_axis_off()


def plot_stratified_coverage(communities, umis, reads_stats, ax, *, cmap):
    """Plot coverage profile, one filled curve per community, semi-transparent
    and overlaid. The x-axis spans only the cluster's own coverage window.

    Arguments:
    communities (list of list of str) - Leiden communities.
    umis        (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    ax        (plt.Axes) - Matplotlib subplot.
    cmap (plt.Colormap) - Discrete colormap used to color communities.
    """
    n_comms = len(communities)
    colors = [cmap(i % cmap.N) for i in range(n_comms)]

    # Collect aligned intervals per community
    all_intervals = []
    per_comm_intervals = []
    for comm in communities:
        intervals = []
        for node in comm:
            for read_id in umis[node]:
                stats = reads_stats[read_id]
                if stats['ref'] == 'unaligned':
                    continue
                start = stats['ref_start']
                end = start + stats['aligned_len']
                intervals.append((start, end))
        per_comm_intervals.append(intervals)
        all_intervals.extend(intervals)

    if not all_intervals:
        ax.text(
            0.5, 0.5, 'No aligned reads in cluster',
            ha='center', va='center', transform=ax.transAxes
        )
        ax.set_axis_off()
        return

    min_pos = min(s for s, _ in all_intervals)
    max_pos = max(e for _, e in all_intervals)
    positions = np.arange(min_pos, max_pos)

    for i, intervals in enumerate(per_comm_intervals):
        if not intervals:
            continue
        depth = np.zeros(max_pos - min_pos, dtype=int)
        for start, end in intervals:
            depth[start - min_pos:end - min_pos] += 1
        n_reads_comm = len(intervals)
        ax.fill_between(
            positions, depth, color=colors[i], alpha=0.45,
            label=f'C{i} ({n_reads_comm} aligned reads)', step='mid'
        )

    ax.set_xlabel('Reference position', fontsize=10)
    ax.set_ylabel('Depth', fontsize=10)
    ax.set_title('Coverage profile stratified by community', fontsize=11)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(min_pos, max_pos)
    ax.set_yscale('log')
    ax.set_ylim(bottom=0.9)


def plot_dissimilarity_matrix(matrix, communities, sub, ax, fig):
    """Plot inter-community consensus dissimilarity heatmap with annotated
    cells. White text on green-marked cells (low dissimilarity), black on
    red (high) for readability.

    Arguments:
    matrix    (np.array) - Symmetric dissimilarity matrix.
    communities (list of list of str) - Leiden communities (used for labels).
    sub (networkx graph) - Connected (sub-)graph (used to count reads per
                            community for tick labels).
    ax        (plt.Axes) - Matplotlib subplot.
    fig     (plt.Figure) - Matplotlib figure (used for colorbar).
    """
    n = len(communities)
    if n < 2:
        ax.text(
            0.5, 0.5, 'Single community: no pairwise comparison',
            ha='center', va='center', transform=ax.transAxes,
            fontsize=10
        )
        ax.set_axis_off()
        return

    # Cap colorbar at max(0.05, observed max) so error-rate-level
    # dissimilarities (~0.01-0.02) are still distinguishable from zero
    vmax = max(0.05, float(matrix.max()))
    im = ax.imshow(matrix, cmap='RdYlGn_r', vmin=0, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    labels = []
    for i, comm in enumerate(communities):
        n_reads = sum(sub.nodes[node]['n_reads'] for node in comm)
        labels.append(f'C{i}\n{n_reads}r')
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    # Annotate cells with the dissimilarity value
    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            text_color = 'white' if val > vmax * 0.5 else 'black'
            ax.text(
                j, i, f'{val:.3f}', ha='center', va='center',
                fontsize=8, color=text_color
            )

    ax.set_title(
        'Inter-community consensus dissimilarity\n'
        + '(Levenshtein / mean consensus length)',
        fontsize=11
    )
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)


def plot_persistence_diagram(persistence, sub, ax, *, cmap):
    """Plot the persistence diagram of the density-peak filtration.
    Each peak is a (birth, death) point in log10(read count) space. Distance
    from the diagonal = log-prominence. Surviving peaks (= those above the
    log-prominence threshold) are drawn in their community colour and
    annotated with the peak's read count; absorbed peaks are drawn as small
    grey crosses. The dashed red line is the threshold.

    Arguments:
    persistence (dict) - Persistence data dict as returned by
                         run_density_peaks (keys: peak_birth, peak_death,
                         survivors, peak_per_community, log_prominence).
    sub  (networkx graph) - Connected (sub-)graph (used to read read counts).
    ax        (plt.Axes) - Matplotlib subplot.
    cmap (plt.Colormap) - Discrete colormap, same one used for the graph and
                          coverage plots so the persistence diagram colours
                          match.
    """
    peak_birth = persistence['peak_birth']
    peak_death = persistence['peak_death']
    survivors = persistence['survivors']
    peak_per_community = persistence['peak_per_community']
    threshold = persistence['log_prominence']

    # Map each surviving peak to its community colour. Communities are sorted
    # by total reads desc, so community 0 = heaviest, matching the rest of
    # the figure.
    peak_to_color = {
        peak: cmap(i % cmap.N) for i, peak in enumerate(peak_per_community)
    }

    # Determine plot bounds with a small margin
    all_b = list(peak_birth.values())
    all_d = list(peak_death.values())
    lo = min(min(all_b), min(all_d)) - 0.15
    hi = max(all_b) + 0.15

    # Diagonal (birth = death = noise)
    ax.plot(
        [lo, hi], [lo, hi], 'k--', alpha=0.35, linewidth=0.8,
        label='Diagonal (peak == saddle)'
    )

    # Threshold line: locus of points with birth - death = threshold.
    # On a diagonal-of-slope-1 plot this is a line parallel to the diagonal,
    # shifted down by `threshold` units.
    ax.plot(
        [lo + threshold, hi], [lo, hi - threshold],
        color='red', alpha=0.55, linestyle=':',
        label=f'Log-prominence threshold = {threshold}'
    )

    # Plot every peak. Survivors first into the legend-eligible scatter so
    # they sit on top of the absorbed-peak cloud.
    for peak in peak_birth:
        b, d = peak_birth[peak], peak_death[peak]
        if peak in survivors:
            color = peak_to_color.get(peak, 'gray')
            ax.scatter(
                [b], [d], c=[color], s=120, marker='o',
                edgecolors='#333333', linewidths=0.6, zorder=3
            )
            count = sub.nodes[peak]['n_reads']
            ax.annotate(
                f'{count}r',
                (b, d), textcoords='offset points', xytext=(6, 4),
                fontsize=8,
            )
        else:
            ax.scatter(
                [b], [d], c='lightgray', s=25, marker='x', linewidths=0.8,
                zorder=2
            )

    ax.set_xlabel('Birth — log10(peak read count)', fontsize=10)
    ax.set_ylabel('Death — log10(saddle read count)', fontsize=10)
    ax.set_title(
        f'Persistence diagram ({len(survivors)} surviving / '
        + f'{len(peak_birth)} total peaks)',
        fontsize=11
    )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='lower right')


def investigate_cluster(graph, umis, reads_stats, fastq_path, representative,
                        out_path, method='leiden', alignment_threshold=3,
                        overlap_threshold=0.75, resolution=1.0,
                        log_prominence=1.0, node_weight_mode='none',
                        restrict_to_reps=False, keep_tied_reps=False, v=0):
    """Re-partition the connected component containing `representative` with
    Leiden, Louvain, or density-peaks, build per-community SPOA consensuses,
    compute pairwise dissimilarity, and produce a single PNG summary plus a
    per-node TSV. When method='density_peaks', the PNG also contains a
    persistence diagram showing every peak's (birth, death) in log10 read-
    count space.

    Arguments:
    graph    (networkx graph) - Networkx graph object.
    umis               (dict) - Dictionnary containing the UMIs as keys and the
                                ids of the reads associated with each UMI.
    reads_stats        (dict) - Dictionnary containing the reads ids as keys
                                and alignment statistics of the reads.
    fastq_path          (str) - Path to fastq file with raw reads.
    representative      (str) - Name of a node belonging to the cluster of
                                interest.
    out_path            (str) - Path to output PNG file.
    method              (str) - Community detection method. Must be 'leiden'
                                (requires igraph+leidenalg) or 'louvain' (uses
                                networkx's built-in), or 'density_peaks'
                                (topographic clustering on log-scaled node
                                weights with log-prominence filtering).
                                Default: 'leiden'.
    alignment_threshold (int) - See edge_similarity (default: 3). Unused for
                                density_peaks.
    overlap_threshold (float) - See edge_similarity (default: 0.75). Unused
                                for density_peaks.
    resolution        (float) - Leiden/Louvain resolution (default: 1.0).
                                Unused for density_peaks.
    log_prominence    (float) - Log10-prominence threshold for density_peaks
                                (default: 1.0 = peaks must be >=10x their
                                saddle). Unused for leiden/louvain.
    node_weight_mode    (str) - See node_weight_boost. Used for leiden/louvain
                                only. Default 'none'.
    restrict_to_reps   (bool) - Restrict each community's consensus,
                                dissimilarity, coverage profile, and clustered-
                                UMIs JSON to just its representative node(s).
                                The graph plot and per-node TSV still show the
                                full community. Default: False.
    keep_tied_reps     (bool) - If True, all nodes tied for max n_reads
                                are kept as co-representatives (reads
                                pooled). Only meaningful when
                                restrict_to_reps is True. Default: False.
    v                   (int) - Level of verbosity (default: 0 = muted).
    """
    # Find the cluster's connected component
    if representative not in graph:
        raise ValueError(
            f'Representative {representative} not found in graph nodes.'
        )
    component = nx.node_connected_component(graph, representative)
    sub = graph.subgraph(component).copy()

    n_reads_total = sum(sub.nodes[n]['n_reads'] for n in sub)
    utils.send_text(
        f'Cluster: {sub.number_of_nodes()} nodes, '
        + f'{sub.number_of_edges()} edges, {n_reads_total} reads',
        v, 2, 1
    )

    # Run community detection
    persistence = None  # set by density_peaks; otherwise stays None
    if method == 'leiden':
        utils.send_text('Running Leiden community detection', v, 2, 1)
        communities = run_leiden(
            sub,
            alignment_threshold=alignment_threshold,
            overlap_threshold=overlap_threshold,
            resolution=resolution,
            node_weight_mode=node_weight_mode,
            v=v,
        )
    elif method == 'louvain':
        utils.send_text('Running Louvain community detection', v, 2, 1)
        communities = run_louvain(
            sub,
            alignment_threshold=alignment_threshold,
            overlap_threshold=overlap_threshold,
            resolution=resolution,
            node_weight_mode=node_weight_mode,
            v=v,
        )
    elif method == 'density_peaks':
        utils.send_text(
            'Running density-peak community detection', v, 2, 1
        )
        communities, persistence = run_density_peaks(
            sub,
            log_prominence=log_prominence,
            v=v,
        )
    else:
        raise ValueError(
            f"Unknown community detection method: {method!r}. "
            "Must be 'leiden' or 'louvain'."
        )

    # Print community summary
    utils.send_text(
        f'Found {len(communities)} communities:', v, 2, 1
    )
    for i, comm in enumerate(communities):
        n_nodes = len(comm)
        n_reads = sum(sub.nodes[node]['n_reads'] for node in comm)
        n_unique = len({node.split('_')[0] for node in comm})
        utils.send_text(
            f'  C{i}: {n_nodes} nodes, {n_unique} unique UMIs, '
            + f'{n_reads} reads',
            v, 2, 1
        )
    
    # Collapse each community to its representative node(s) when
    # restrict_to_reps is on
    if restrict_to_reps:
        effective_communities = restrict_communities_to_reps(
            communities, sub, keep_tied=keep_tied_reps
        )
        n_kept = sum(len(c) for c in effective_communities)
        utils.send_text(
            f'Restricting to representatives: {n_kept} node'
            + ('s' if n_kept != 1 else '')
            + f' kept across {len(communities)} communities '
            + ('(ties pooled)' if keep_tied_reps else '(first max each)'),
            v, 2, 1
        )
    else:
        effective_communities = communities

    # Load read sequences for the cluster
    all_read_ids = []
    for node in sub:
        all_read_ids.extend(umis[node])
    utils.send_text('Loading cluster reads from fastq', v, 2, 1)
    seqs, _quals = load_cluster_reads(
        fastq_path, all_read_ids, reads_stats, v=v
    )

    # Build community consensuses
    utils.send_text(
        'Building community consensuses with SPOA', v, 2, 1
    )
    consensuses = build_community_consensuses(
        effective_communities, umis, seqs, v=v
    )

    # Compute pairwise consensus dissimilarity
    utils.send_text(
        'Computing inter-community consensus dissimilarity', v, 2, 1
    )
    matrix = compute_dissimilarity_matrix(consensuses)

    # Plot summary. Density-peaks adds a persistence diagram as a new top
    # row spanning both columns; other methods keep the original 2-row
    # layout. The graph + dissimilarity matrix + coverage row arrangement is
    # otherwise identical.
    utils.send_text('Plotting investigation summary', v, 2, 1)
    if persistence is not None:
        fig = plt.figure(figsize=(18, 15))
        gs = fig.add_gridspec(3, 2, height_ratios=[0.9, 1.0, 0.9])
        ax_pers = fig.add_subplot(gs[0, :])
        ax_graph = fig.add_subplot(gs[1, 0])
        ax_dissim = fig.add_subplot(gs[1, 1])
        ax_cov = fig.add_subplot(gs[2, :])
    else:
        fig = plt.figure(figsize=(18, 11))
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.9])
        ax_pers = None
        ax_graph = fig.add_subplot(gs[0, 0])
        ax_dissim = fig.add_subplot(gs[0, 1])
        ax_cov = fig.add_subplot(gs[1, :])

    cmap = plt.get_cmap('tab10')
    if ax_pers is not None:
        plot_persistence_diagram(persistence, sub, ax_pers, cmap=cmap)
    plot_graph_by_community(sub, communities, ax_graph, fig, cmap=cmap)
    plot_dissimilarity_matrix(
        matrix, effective_communities, sub, ax_dissim, fig
    )
    plot_stratified_coverage(
        effective_communities, umis, reads_stats, ax_cov, cmap=cmap
    )

    # The bottom-of-suptitle tuning knob differs by method
    if method == 'density_peaks':
        knob = f'log-prominence={log_prominence}'
    else:
        knob = (f'resolution={resolution}, '
                f'node_weight_mode={node_weight_mode}')
    if restrict_to_reps:
        knob += ', restricted_to_reps'
        if keep_tied_reps:
            knob += '+tied'
    fig.suptitle(
        f'Cluster {representative}  —  {method.capitalize()} investigation\n'
        + f'{len(communities)} communit'
        + ('y' if len(communities) == 1 else 'ies')
        + f' from a single cluster of {sub.number_of_nodes()} nodes / '
        + f'{n_reads_total} reads ({knob})',
        fontsize=13
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # Companion TSV with per-node community assignment
    tsv_path = out_path.rsplit('.', 1)[0] + '_communities.tsv'
    with open(tsv_path, 'w') as tsv_out:
        tsv_out.write(
            'community\tnode\tn_reads\tn_unique_umis_in_community\t'
            + 'community_consensus_length\n'
        )
        for i, comm in enumerate(communities):
            n_unique = len({node.split('_')[0] for node in comm})
            cons_len = len(consensuses.get(i, ''))
            for node in comm:
                tsv_out.write(
                    f'{i}\t{node}\t{sub.nodes[node]["n_reads"]}\t'
                    + f'{n_unique}\t{cons_len}\n'
                )
    utils.send_text(f'Wrote community assignment to {tsv_path}', v, 2, 1)

    # Companion JSON with newly identified UMIs clusters
    json_path = out_path.rsplit('.', 1)[0] + '_clustered_umis.json'
    clustered_umis = {}
    for comm in effective_communities:
        rep = max(comm, key=lambda n: sub.nodes[n]['n_reads'])
        clustered_umis[rep] = [
            read_id for node in comm for read_id in umis[node]
        ]
    utils.save_json(json_path, clustered_umis)
    utils.send_text(
        f'Wrote clustered UMIs JSON to {json_path}', v, 2, 1
    )

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--representative', type=str, required=True,
                        help='Name of a node belonging to the cluster to '
                        "investigate. The tool retrieves that node's "
                        'connected component, re-partitions it with community '
                        'detection, and reports.')
    parser.add_argument('-g', '--graph_pickle', type=str, required=True,
                        help='Path to input pickle file containing the UMIs '
                        'clustering graph.')
    parser.add_argument('-u', '--cc_umis', type=str, required=True,
                        help='Path to input json file containing coverage '
                        'connex UMIs information.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        'alignment statistics.')
    parser.add_argument('-f', '--input_fastq', type=str, required=True,
                        help='Path to input fastq file with raw reads, used '
                        'to build SPOA community consensuses.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output PNG file. A TSV with the '
                        'per-node community assignment is written alongside '
                        '(same name with _communities.tsv suffix).')
    parser.add_argument('-at', '--alignment_threshold', type=int, default=3,
                        help='Alignment-distance threshold used to normalize '
                        'edge similarity weights for Louvain/Leiden. Should '
                        'match the value passed to deduplicate_umis to reflect'
                        ' the edge geometry actually used during clustering. '
                        'Default: 3.')
    parser.add_argument('-ot', '--overlap_threshold', type=float, default=0.75,
                        help='Overlap pass threshold (similarity, not '
                        'distance), used to normalize edge similarity '
                        'weights for Louvain/Leiden. Default: 0.75.')
    parser.add_argument('-m', '--method', type=str, default='leiden',
                        choices=['leiden', 'louvain', 'density_peaks'],
                        help='Community detection method. leiden requires '
                        'igraph and leidenalg installed; louvain uses '
                        "networkx's built-in implementation; density_peaks "
                        'is a topographic / persistent-homology approach '
                        'that uses node weights as elevation and filters '
                        'peaks by log-prominence (see --log_prominence). '
                        'Default: leiden.')
    parser.add_argument('--resolution', type=float, default=1.0,
                        help='Leiden RB / Louvain resolution parameter. Higher'
                        ' values yield more (smaller) communities. Unused when'
                        ' method is density_peaks. Default: 1.0.')
    parser.add_argument('-lp', '--log_prominence', type=float, default=1.0,
                        help='Log10-prominence threshold for density_peaks. '
                        'A peak survives if log10(peak_count / saddle_count) '
                        '>= threshold. 1.0 means "peak must be >=10x its '
                        'saddle"; 0.5 means ">=~3x"; 2.0 means ">=100x". '
                        'Unused for leiden / louvain. Default: 1.0.')
    parser.add_argument('-nwm', '--node_weight_mode', type=str,
                        default='none',
                        choices=['none', 'log_max', 'geomean', 'asymmetric'],
                        help='Multiplicative boost applied to edge weights '
                        'using endpoint node weights, when method is leiden '
                        'or louvain. none = read counts information unused; '
                        'log_max =  log10(max(c_u, c_v) + 1), boosts edges '
                        'touching any heavy node; geomean = sqrt(c_u * c_v), '
                        'aggressive; asymmetric = 1 + log10(max/min), '
                        "continuous analog of umi-tools' 2N-1 rule. Unused for"
                        ' density_peaks. Default: none.')
    parser.add_argument('-rr', '--restrict_to_reps', action='store_true',
                        help="Restrict each community's consensus, "
                        'dissimilarity, coverage profile and clustered-UMIs '
                        'JSON to just the representative node (= heaviest node'
                        ' by n_reads). The graph plot and per-node TSV still '
                        'show the full community. Default: off.')
    parser.add_argument('--keep_tied_reps', action='store_true',
                        help='When several nodes tie for the maximum '
                        'n_reads within a community, keep them all as '
                        'co-representatives (their reads are pooled in the '
                        'JSON, coverage and consensus). When off, only the '
                        'first such node is kept (deterministic). Only '
                        'meaningful with --restrict_to_reps. Default: off.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load graph
    utils.send_text('Loading UMIs graph', v, 1, 0)
    with open(args.graph_pickle, 'rb') as graph_pkl_file:
        graph = pickle.load(graph_pkl_file)

    # Load coverage connex UMIs
    utils.send_text('Loading coverage connex UMIs', v, 1, 0)
    cc_umis = utils.load_json(args.cc_umis)

    # Load reads statistics
    utils.send_text('Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)

    # Investigate
    utils.send_text(
        f'Investigating cluster around {args.representative}', v, 1, 0
    )
    res = investigate_cluster(
        graph,
        cc_umis,
        reads_stats,
        args.input_fastq,
        args.representative,
        args.output,
        method=args.method,
        alignment_threshold=args.alignment_threshold,
        overlap_threshold=args.overlap_threshold,
        resolution=args.resolution,
        log_prominence=args.log_prominence,
        node_weight_mode=args.node_weight_mode,
        restrict_to_reps=args.restrict_to_reps,
        keep_tied_reps=args.keep_tied_reps,
        v=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
