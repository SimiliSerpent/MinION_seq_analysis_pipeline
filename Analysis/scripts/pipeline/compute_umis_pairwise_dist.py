#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute pairwise distance for a set of groups of UMIs.
Two kind of distances are computed: a distance based on alignment cost of two
distinct UMIs; and a distance based on the overlap of the coverage profile of
all the reads associated with each UMI of the pair. Optionally outputs heatmaps
representing distances for every pair of UMIs.
"""


import argparse
import functools
import multiprocessing
import os
import sys
import time

import Bio.Align
import matplotlib.pyplot as plt
import numpy as np
import rapidfuzz.distance
import scipy
import skbio

import utils


# Define aligner parameters (local/global and score parameters)
MY_ALIGNER = Bio.Align.PairwiseAligner()
MY_ALIGNER.mode = 'global'
MY_ALIGNER.match_score = 0
MY_ALIGNER.mismatch_score = -2
MY_ALIGNER.open_gap_score = -2
MY_ALIGNER.extend_gap_score = -2
MY_ALIGNER.end_gap_score = -1


def pair_generator(my_set):
    """Generator of pair of items.

    Arguments:
    my_set (list) - List of any items.
    """
    for i in range(len(my_set)):
        for j in range(i + 1, len(my_set)):
            yield (my_set[i], my_set[j])


def levenshtein_dist(pair):
    """Returns Levenshtein distance

    Arguments:
    pair (2-tuple of str) - Pair of sequence.
    
    Return:
                    (int) - Edition distance between the two input sequence.
    """
    seq1, seq2 = pair
    return rapidfuzz.distance.Levenshtein.distance(
        seq1,
        seq2,
        weights=(1,1,2),
        score_cutoff=None
    ) # weights are (insertion, deletion, substitution)


def alignment_dist(pair):
    """Returns alignnment-based distance

    Arguments:
    pair (2-tuple of str) - Pair of sequence.
    
    Return:
                    (int) - Alignment distance between the two input sequence.
    """
    seq1, seq2 = pair
    return -MY_ALIGNER.score(seq1, seq2)


def compute_alignment_dist(sequences, method='Levenshtein', cpus=None, v=0):
    """Compute matrix of pairwise distances from list of sequences
    Distances are computed using the specified method.

    Arguments:
    sequences (str) - Set of sequences to compute pairwise distance from.
    method    (str) - Name of method to compute UMIs pairwise distances.
                      Must be Levenshtein or SmithWaterman (default:
                      Levenshtein).
    cpus      (str) - Number of cpus for parallelization (optional).
    v         (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    dist_matrix (np.array) - Distance matrix.
    """
    n = len(sequences)
    # Compute distance
    if v > 0:
        t_zero = time.perf_counter()
    if method == 'Levenshtein':
        scores = []
        for umi_pair in pair_generator(sequences):
            scores.append(levenshtein_dist(umi_pair))

    elif method == 'SmithWaterman':

        # If cpus have been specified, parallelize computation
        if cpus is not None:
            with multiprocessing.Pool(processes=cpus) as pool:
                scores = pool.map(alignment_dist, pair_generator(sequences))

        else:
            scores = []
            for umi_pair in pair_generator(sequences):
                scores.append(alignment_dist(umi_pair))
    
    else:
        raise ValueError(
            'compute_dist_matrix: method must be either Levenshtein or '
            + f'SmithWaterman, not {method}'
        )
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Distances computed in {elapsed:.3f} s for {n} '
                        + 'sequences', v, 2, 1)
    
    # Populate 2D array
    if v > 0:
        t_zero = time.perf_counter()
    dist_matrix = np.zeros([n, n], dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            dist_matrix[i, j] = scores[n*i+j-1-i*(i+3)//2]
            dist_matrix[j, i] = scores[n*i+j-1-i*(i+3)//2]
    # dist_matrix = np.ma.masked_where(dist_matrix == -1, dist_matrix)
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Created distance matrix in {elapsed:.3f} s', v, 2, 1)

    return dist_matrix


def build_umi_cov_events(stats, umis, umi):
    """Builds a list of coverage event (inc / dec) using UMI reads information

    Arguments:
    stats      (dict) - Dictionnary containing the reads ids as keys and
                        alignment statistics of the reads.
    umis       (dict) - Dictionnary containing the UMIs as keys and the
                        ids of the reads associated with each UMI.
    umi         (str) - UMI "name".
    
    Return:
    stats      (dict) - Dictionnary containing the reads ids as keys and
                        alignment statistics of the reads.
    size        (int) - Number of bases in that UMI reads.
    """
    # Retrieve reads information
    reads_info = [
        (
            stats[read_id]['ref_start'],
            stats[read_id]['aligned_len']
        ) for read_id in umis[umi]
    ]
    # Build events list describing coverage profile for that UMI's reads
    cov_events = []
    size = 0
    for start_pos, read_len in reads_info:
        cov_events.append((start_pos, umi, 1))                  # read start
        cov_events.append((start_pos + read_len - 1, umi, -1))  # read end
        size += read_len
    return cov_events, size


def get_overlap(cov_events_pair):
    """Returns coverage overlap of two sets of reads
    Coverage overlap is defined as the ratio of the intersection of the two
    coverage profiles (for each list of reads) over the minimum amount of bases
    in both list of reads.

    Arguments:
    cov_events_pair (2-tuple) - Pair of (list of events, size), corresponding
                                to a coverage increase or decrease for the
                                events, and the number of bases associated with
                                one UMI. One tuple for each UMI of the pair.
    
    Return:
                      (float) - Overlap distance between the two read
                                sets.
    """
    (cov_events_A, size_A), (cov_events_B, size_B) = cov_events_pair
    umi_A, umi_B = cov_events_A[0][1], cov_events_B[0][1]
    # Merge the list of events
    cov_events = cov_events_A + cov_events_B
    # Sort events to scan reference 5' -> 3'
    cov_events.sort()

    # Loop through events and compute intersection size
    cov_A = cov_B = 0
    prev = None
    intersection = 0
    for pos, umi, delta in cov_events:
        if prev is not None and pos > prev:
            intersection += (pos - prev) * min(cov_A, cov_B)

        if umi == umi_A:
            cov_A += delta
        else:
            cov_B += delta

        prev = pos

    return 1 - (intersection / min(size_A, size_B))


def compute_overlap_dist(umi_names, umis, stats, align_mat=None,
                         align_threshold=3, cpus=None, v=0):
    """Compute matrix of pairwise coverage-based distances from sets of reads
    Distances are computed using coverage overlap between two read sets.
    If an alignment-distance matrix is provided, only pairs of UMIs that share
    the same reference contig AND whose alignment distance is below
    align_threshold are evaluated; other pairs get the max overlap distance
    (1.0). This skips computing distances for pairs that could never end up
    merged by the clustering step anyway.

    Arguments:
    umi_names (list of str) - List of UMI names (used to preserve UMIs order).
    umis        (dict) - Dictionnary containing the UMIs as keys and the
                         ids of the reads associated with each UMI.
    stats       (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    align_mat (np.array) - Alignment distance matrix. If specified, only the
                         pairs of UMIs closer than the align_threshold will be
                         considered for overlap computing. Others will be set
                         to overlap-distance 1. Defaults to None, i.e. all
                         pairs will be evaluated.
    align_threshold (float) - Threshold used together with align_mat. Defaults
                         to 3.
    cpus         (str) - Number of cpus for parallelization (optional).
    v            (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    dist_matrix (np.array) - Distance matrix.
    """
    n = len(umi_names)
    # Compute coverage profiles for all UMI's reads set
    if v > 0:
        t_zero = time.perf_counter()
    build_umi_cov_events_partial = functools.partial(
        build_umi_cov_events,
        stats,
        umis
    )
    # If cpus have been specified, parallelize computation
    if cpus is not None:
        with multiprocessing.Pool(processes=cpus) as pool:
            all_cov_events = pool.map(build_umi_cov_events_partial, umi_names)
    else:
        all_cov_events = []
        for umi_name in umi_names:
            all_cov_events.append(build_umi_cov_events_partial(umi_name))
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Built coverage profiles in {elapsed:.3f} s for {n} '
                        + 'read sets', v, 2, 1)
    
    # Build the candidate pair list. If align_mat is provided, restrict to
    # pairs that could conceivably be merged later (same reference contig AND
    # alignment distance below threshold). Otherwise, every pair is a
    # candidate (legacy behaviour, useful for full-matrix statistics).
    if v > 0:
        t_zero = time.perf_counter()
    umi_refs = [umi.split('_')[1] for umi in umi_names]
    if align_mat is None:
        candidate_pairs = [
            (i, j) for i in range(n) for j in range(i + 1, n)
            if umi_refs[i] == umi_refs[j]
        ]
    else:
        candidate_pairs = [
            (i, j) for i in range(n) for j in range(i + 1, n)
            if umi_refs[i] == umi_refs[j]
            and align_mat[i, j] <= align_threshold
        ]
    n_pairs = len(candidate_pairs)
    n_total = n * (n - 1) // 2
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        pct = 100 * n_pairs / max(1, n_total)
        utils.send_text(
            f'Selected {n_pairs} candidate pairs out of {n_total} '
            f'({pct:.2f}%) in {elapsed:.3f} s',
            v, 2, 1
        )

    # Compute overlap distance for each candidate pair
    if v > 0:
        t_zero = time.perf_counter()
    pair_inputs = [
        (all_cov_events[i], all_cov_events[j]) for i, j in candidate_pairs
    ]
    # If cpus have been specified, parallelize computation
    if cpus is not None:
        with multiprocessing.Pool(processes=cpus) as pool:
            overlaps = pool.map(get_overlap, pair_inputs)
    else:
        overlaps = [get_overlap(pair) for pair in pair_inputs]
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Distances computed in {elapsed:.3f} s for {n} read '
                        + 'sets', v, 2, 1)
    
    # Populate 2D array. Pairs filtered out get the max overlap distance
    # (1.0), so they will never pass the overlap gate during clustering.
    if v > 0:
        t_zero = time.perf_counter()
    dist_matrix = np.ones([n, n], dtype=float)
    np.fill_diagonal(dist_matrix, 0.0)
    for (i, j), overlap in zip(candidate_pairs, overlaps):
        dist_matrix[i, j] = overlap
        dist_matrix[j, i] = overlap
    # dist_matrix = np.ma.masked_where(dist_matrix == -1, dist_matrix)
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Created distance matrix in {elapsed:.3f} s', v, 2, 1)
    
    return dist_matrix


def plot_distance_matrix(dist_matrix, png_path, method='avg-linkage', v=0):
    """Plot input distance matrix as heatmap

    Arguments:
    dist_matrix (np.array) - Distance matrix.
    png_path         (str) - Path to heatmap showing distance matrix.
    method           (str) - Method to perform reordering of the distance
                             matrix indices for visualization purposes. Must be
                             one of avg-linkage, neighbour-joining. The former
                             is much faster but less exact than the latter.
    v                (int) - Level of verbosity (default: 0 = muted).
    """
    utils.send_text('Initialize distance matrix', v, 2, 1)
    n = dist_matrix.shape[0]

    # Skip plotting for matrix of size 1 by 1
    if n < 2:
        utils.send_text(
            f'plot_distance_matrix: skipping {png_path} (n={n} < 2, '
            + 'nothing to cluster).',
            v, 1, 1
        )
        return 0

    if method == 'neighbour-joining':
        # Perform Neighbor-joining on UMIs
        if v > 0:
            t_zero = time.perf_counter()
        indices = [str(i) for i in range(n)]
        skbio_dist_mat = skbio.DistanceMatrix(dist_matrix, indices)

        utils.send_text('Perform neighbour-joining', v, 2, 1)
        distance_tree = skbio.tree.nj(skbio_dist_mat)
        if v > 0:
            elapsed = time.perf_counter() - t_zero
            utils.send_text(f'Performed neighbour-joining in {elapsed:.3f} s',
                            v, 2, 1)

        # Extract clusters by distance threshold
        # groups = []
        # threshold = 5
        # for node in distance_tree.non_tips():
        #     tips = [t.name for t in node.tips()]
        #     if node.length is not None and node.length <= threshold:
        #         groups.append(tips)

        # Reorder distance matrix to cluster similar UMIs
        utils.send_text('Re-order indices', v, 2, 1)
        leaf_order = [int(tip.name) for tip in distance_tree.tips()]
    
    elif method == 'avg-linkage':
        vectorized_distances = scipy.spatial.distance.squareform(
            dist_matrix,
            checks=False
        )
        linkage_matrix = scipy.cluster.hierarchy.linkage(
            vectorized_distances,
            method='average'
        )
        leaf_order = scipy.cluster.hierarchy.leaves_list(linkage_matrix)
    
    else:
        raise ValueError('The specified clustering method is not allowed:',
                         method)

    dist_matrix_reordered = dist_matrix[np.ix_(leaf_order, leaf_order)]

    # Plot heatmap
    utils.send_text('Plot distance heatmap', v, 2, 1)
    plt.clf() # Clear existing figure
    fig = plt.figure(
        figsize=(20, 15)
    )
    cmap = plt.get_cmap('RdYlGn').copy()
    cmap.set_bad(color='black')
    ax = fig.add_subplot(1, 1, 1)
    mesh = ax.pcolormesh(
        dist_matrix,
        cmap=cmap,
        shading='auto'
    )
    fig.supxlabel('Read length')
    fig.supylabel('Read quality (Phred)')
    fig.colorbar(
        mesh,
        label='Distance',
        ax=ax,
        location='right',
        fraction = 0.06,
        pad = 0.03,
        aspect=6
    )
    fig.suptitle('UMIs distance matrix')

    # Saving figure
    utils.send_text('Saving distance matrix', v, 2, 1)
    plt.savefig(png_path, dpi=500, bbox_inches="tight")

    return 0


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
    parser.add_argument('-M', '--dist_method', type=str, default='Levenshtein',
                        help='Name of method to compute UMIs pairwise '
                        'distances. Must be Levenshtein or SmithWaterman. '
                        'Default: Levenshtein.')
    parser.add_argument('-at', '--alignment_threshold', type=int, default=3,
                        help='Two clusters of reads associated with two '
                        'distinct UMI sequences must have their UMIs distant '
                        'from at most this threshold (alignment-based '
                        'distance) for their overlap distance to be '
                        'considered. Using Levenshtein distance, a '
                        'substitution costs 2 and a translation costs 2. '
                        'Default value is 3, meaning that two clusters must '
                        'have their UMIs differ from at most 1 translation or '
                        '1 substitution for their overlap to be considered.')
    parser.add_argument('--full_overlap_matrix', action='store_true',
                        help='Ignores the `alignment_threshold` parameter and '
                        'compute the full pairwise overlap distance matrix. '
                        'Much slower. Default: False (only UMI pairs with '
                        'short enough alignment distances are evaluated, the '
                        'rest are set to overlap distance 1.0).')
    parser.add_argument('-c', '--cores', type=int, default=None,
                        help='Number of CPUs available. If specified, distance'
                        ' computation will be parralelized only if '
                        'SmithWaterman is chosen, because overhead is too '
                        'high for Levenshtein distance computations.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output files do not start with an underscore
    output = args.output_prefix
    if output.endswith('/'):
        output += 'dist_out'
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

    # Build UMI sequences alignment-based distance matrix (and plot it)
    utils.send_text(f'Computing UMIs alignment distance matrix', v, 1, 0)
    alignment_dist_matrix = compute_alignment_dist(
        sequences=[umi_name.split('_')[0] for umi_name in cc_umis_list],
        method=args.dist_method,
        cpus=args.cores,
        v=v
    )
    # Save alignment-based distance matrix
    utils.send_text(f'Saving alignment distance matrix', v, 1, 0)
    np.savez_compressed(
        f'{output}_dist_mat_alignment.npz',
        align_dist_mat=alignment_dist_matrix
    )
    # Plot alignment-based distance matrix as heatmap
    utils.send_text(f'Plotting UMIs alignment distance heatmap', v, 1, 0)
    res = plot_distance_matrix(
        alignment_dist_matrix,
        f'{output}_alignment_heatmap.png',
        v=v
    )

    # Build UMI reference coverage overlap-based distance matrix (and plot it)
    utils.send_text(f'Computing UMIs overlap distance matrix', v, 1, 0)
    overlap_dist_matrix = compute_overlap_dist(
        umi_names=cc_umis_list,
        umis=cc_umis,
        stats=reads_stats,
        align_mat=None if args.full_overlap_matrix else alignment_dist_matrix,
        align_threshold=args.alignment_threshold,
        cpus=args.cores,
        v=v
    )
    # Save overlap-based distance matrix
    utils.send_text(f'Saving overlap distance matrix', v, 1, 0)
    np.savez_compressed(
        f'{output}_dist_mat_overlap.npz',
        overlap_dist_mat=overlap_dist_matrix
    )
    # Plot overlap-based distance matrix as heatmap
    utils.send_text(f'Plotting UMIs overlap distance heatmap', v, 1, 0)
    res = plot_distance_matrix(
        overlap_dist_matrix,
        f'{output}_overlap_heatmap.png',
        v=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
