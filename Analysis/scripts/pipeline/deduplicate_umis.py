#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Homemade UMIs deduplication tool.
From a BAM file and a UMI sequence, perform UMI deduplication. Optionally
outputs a range of statistics, and optionally polishes the resulting sequences
using consensus calling.
"""


import argparse
import contextlib
import itertools
import multiprocessing
import os
import random
import subprocess
import sys
import tempfile
import time

import Bio.Align
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pysam
import rapidfuzz.distance
import skbio
import spoa

import utils


# Set to True if comparing with umi-tools to ensure identical bam input
TEST_COMPARE_UMITOOLS = True

# CIGAR operations used for CIGAR updating
M, I, D, N, S, EQ, X = 0, 1, 2, 3, 4, 7, 8
QUERY_OPS = {M, I, S, EQ, X}
REF_OPS   = {M, D, N, EQ, X}

# Define aligner parameters (local/global and score parameters)
MY_ALIGNER = Bio.Align.PairwiseAligner()
MY_ALIGNER.mode = 'global'
MY_ALIGNER.match_score = 0
MY_ALIGNER.mismatch_score = -2
MY_ALIGNER.open_gap_score = -2
MY_ALIGNER.extend_gap_score = -2
MY_ALIGNER.end_gap_score = -1


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
                        + f'{elapsed:.1f} s.',
                        v, 1, 1
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
            + f'{elapsed:.1f} s.',
            v, 1, 1
        )
            
    return reads_stats


def update_umis(old_umis, reads_stats):
    """Update UMIs with only reads in the provided statistics

    Arguments:
    old_umis    (dict) - Dictionnary containing the UMIs as keys and and the
                         ids of the reads associated with each UMI.
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    
    Return:
    umis        (dict) - Dictionnary containing the UMIs as keys and and the
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


def split_connexe_coverage(umis, reads_stats, v=0):
    """Split UMIs according to associated reads coverage connected components
    Reads associated with a given UMIs are placed in a graph. Edges are drawn
    between reads if and only if they overlap (once aligned on the reference).
    If the graph has multiple connected components, the UMI is split in
    several sub read-lists. In practice, use a 1-D line rather than a graph for
    faster execution

    Arguments:
    umis        (dict) - Dictionnary containing the UMIs as keys and and the
                         ids of the reads associated with each UMI.
    reads_stats (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    v            (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    by_cc_umis  (dict) - Updated dictionnary containing the UMIs.
    """
    n = len(umis)
    # Initiate new UMIs dictionnary
    by_ref_umis = {}
    by_cc_umis = {}

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
            + f'{elapsed:.1f} s.',
            v, 1, 1
        )
        t_zero = time.perf_counter()

    # Split UMIs on connected components
    for umi, read_ids in by_ref_umis.items():
        if umi.endswith('unaligned'):
            by_cc_umis[umi] = read_ids
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
            by_cc_umis[umi] = read_ids
        else:
            for i, component in enumerate(components):
                by_cc_umis[f'{umi}_{i+1}'] = component
        # TODO: Remove
        # if '9173c2f4-e258-4693-96b9-9c2caba906fe_1' in read_ids \
        #     or 'ccef40b9-ff6b-41ae-8351-d09b6efe0a16_1' in read_ids:
        #     print(components)
        #     print(intervals)
        #     print(umi)
    by_cc_n = len(by_cc_umis)
    
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Split {by_ref_n} ref-aware UMIs in {by_cc_n} coverage-connected '
            + f'UMIs in {elapsed:.1f} s.',
            v, 1, 1
        )
    
    return by_cc_umis


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

        # If cpus have been specified, parrellelize computation
        if cpus is not None:
            with multiprocessing.Pool(processes=cpus) as p:
                scores = p.map(alignment_dist, pair_generator(sequences))

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
        utils.send_text(f'Distances computed in {elapsed} s for {n} sequences',
                        v, 2, 1)
    
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
        utils.send_text(f'Created distance matrix in {elapsed} s', v, 2, 1)

    return dist_matrix


def get_overlap(read_list_pair):
    """Returns coverage overlap of two sets of reads
    Coverage overlap is defined as the ratio of the intersection of the two
    coverage profiles (for each list of reads) over the minimum amount of bases
    in both list of reads.

    Arguments:
    pair (2-tuple of list) - Pair of list of reads. Each read is defined as
                             (start, length) 2-tuple.
    
    Return:
                   (float) - Overlap between the two read sets.
    """
    set_A, set_B = read_list_pair
    # Define a list of events, corresponding to a coverage increase or decrease
    # for either list A or list B
    events = []

    # Fill events by looping through each set's reads
    for start_pos, read_len in set_A:
        events.append((start_pos, 0, 1))        # start A
        events.append((start_pos + read_len - 1, 0, -1))   # end A

    for start_pos, read_len in set_B:
        events.append((start_pos, 1, 1))        # start B
        events.append((start_pos + read_len - 1, 1, -1))   # end B

    # Sort events to scan reference 5' -> 3'
    events.sort()

    # Loop through events and compute intersection size
    covA = covB = 0
    prev = None
    intersection = 0
    for pos, set_id, delta in events:
        if prev is not None and pos > prev:
            intersection += (pos - prev) * min(covA, covB)

        if set_id == 0:
            covA += delta
        else:
            covB += delta

        prev = pos

    sizeA = sum([read_len for start_pos, read_len in set_A])
    sizeB = sum([read_len for start_pos, read_len in set_B])

    return 1 - (intersection / min(sizeA, sizeB))


def compute_overlap_dist(umi_names, umis, stats, cpus=None, v=0):
    """Compute matrix of pairwise coverage-based distances from sets of reads
    Distances are computed using coverage overlap between two read sets.

    Arguments:
    umi_names (list of str) - List of UMI names (used to preserve UMIs order).
    umis        (dict) - Dictionnary containing the UMIs as keys and and the
                         ids of the reads associated with each UMI.
    stats       (dict) - Dictionnary containing the reads ids as keys and
                         alignment statistics of the reads.
    cpus         (str) - Number of cpus for parallelization (optional).
    v            (int) - Level of verbosity (default: 0 = muted)
    
    Return:
    dist_matrix (np.array) - Distance matrix.
    """
    n = len(umi_names)
    # Build list of list of reads to provide pair generator with
    reads_lists = [
        [
            (
                stats[read_id]['ref_start'],
                stats[read_id]['aligned_len']
            ) for read_id in umis[umi]
        ] for umi in umi_names
    ]
    
    # Compute distance
    if v > 0:
        t_zero = time.perf_counter()
    # If cpus have been specified, parrellelize computation
    if cpus is not None:
        with multiprocessing.Pool(processes=cpus) as p:
            overlaps = p.map(get_overlap, pair_generator(reads_lists))

    else:
        overlaps = []
        for reads_list_pair in pair_generator(reads_lists):
            overlaps.append(get_overlap(reads_list_pair))

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Distances computed in {elapsed} s for {n} read sets',
                        v, 2, 1)
    
    # Populate 2D array
    if v > 0:
        t_zero = time.perf_counter()
    dist_matrix = np.zeros([n, n], dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            dist_matrix[i, j] = overlaps[n*i+j-1-i*(i+3)//2]
            dist_matrix[j, i] = overlaps[n*i+j-1-i*(i+3)//2]
    # dist_matrix = np.ma.masked_where(dist_matrix == -1, dist_matrix)
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Created distance matrix in {elapsed} s', v, 2, 1)
    
    return dist_matrix


def plot_distance_matrix(dist_matrix, png_path, v=0):
    """Plot input distance matrix as heatmap

    Arguments:
    dist_matrix (np.array) - Distance matrix.
    png_path         (str) - Path to heatmap showing distance matrix.
    v                (int) - Level of verbosity (default: 0 = muted).
    """
    n = dist_matrix.shape[0]
    # Perform Neighbor-joining on UMIs
    if v > 0:
        t_zero = time.perf_counter()
    indices = [str(i) for i in range(n)]
    skbio_dist_mat = skbio.DistanceMatrix(dist_matrix, indices)
    distance_tree = skbio.tree.nj(skbio_dist_mat)
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(f'Performed neighbour-joining in {elapsed} s',
                        v, 2, 1)

    # Extract clusters by distance threshold
    # groups = []
    # threshold = 5
    # for node in distance_tree.non_tips():
    #     tips = [t.name for t in node.tips()]
    #     if node.length is not None and node.length <= threshold:
    #         groups.append(tips)

    # Reorder distance matrix to cluster similar UMIs
    leaf_order = [int(tip.name) for tip in distance_tree.tips()]
    dist_matrix_reordered = dist_matrix[np.ix_(leaf_order, leaf_order)]

    # Plot heatmap
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
    utils.send_text('Saving distance matrix', v, 3, 1)
    plt.savefig(png_path, dpi=500, bbox_inches="tight")

    return 0


def cluster_umis(umi_names, umis, alignment_mat, overlap_mat,
                 alignment_threshold=3, overlap_threshold=0.7,
                 umi_stats_path=None, reads_stats=None, groups_path=None, v=0):
    """Cluster UMIs based on sequence and coverage proximity

    Arguments:
    umi_names  (list of str) - List of UMI names (used to preserve UMIs order).
    umis               (str) - Dictionnary containing the UMIs as keys and and
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
    """
    # Initialize graph
    n = len(umis)
    umis_graph = nx.Graph()
    umis_graph.add_nodes_from(umi_names)
    if v > 0:
        t_zero = time.perf_counter()

    # Connect UMIs based on pairwise distances and thresholds
    # if and only if they are aligned to the same reference contig
    # (or both unaligned)
    umi_refs = {umi: umi.split('_')[1] for umi in umi_names}
    for i in range(n):
        for j in range(i + 1, n):
            if umi_refs[umi_names[i]] == umi_refs[umi_names[j]]:
                if alignment_mat[i, j] <= alignment_threshold:
                    if overlap_mat[i, j] <= (1 - overlap_threshold):
                        umis_graph.add_edge(umi_names[i], umi_names[j])
    # Display processing time
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Connected UMIs using distance in {elapsed:.1f} s.',
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

    return clustered_umis


def trim_polyA(seq, qual):
    """Trim 3' trailing A

    Arguments:
    seq  (str) - Sequence to trim.
    qual (arr) - Associated quality string.
    
    Return:
    seq  (str) - Trimmed sequence.
    qual (arr) - Updated quality string.
    """
    # Trim sequence
    right_trimmed = 0
    length = len(seq)
    while seq[-5:].count('A') >= 3:
        if seq.endswith('A'):
            seq = seq.rstrip('A')
        else:
            seq = seq[:-1]
    right_trimmed = length - len(seq)
    # Update quality array
    if right_trimmed > 0:
        qual = qual[:-right_trimmed]
    return seq, qual


def filter_GA(seq, qual, threshold=10):
    """Filter GA(/AG) repeats longer than specified length

    Arguments:
    seq  (str) - Sequence to filter.
    qual (arr) - Associated quality string.
    
    Return:
    seq  (str) - Filtered sequence.
    qual (arr) - Updated quality string.
    """
    repeats = []
    for ind, char in enumerate(seq):
        if ind < len(seq)-1:
        # Test for GA/AG
            if seq[ind:ind+2] == 'GA' or seq[ind:ind+2] == 'AG':
                # Test if still in previously detected repeat
                if repeats and ind + 1 == repeats[-1][0] + repeats[-1][1]:
                    # update current repeat...
                    repeats[-1][1] += 1
                # ...else, start new repeat
                else:
                    repeats.append([ind, 2])
    repeats = [rep for rep in repeats if rep[1] >= threshold]
    for rep in repeats[::-1]:
        seq = seq[:rep[0]] + seq[sum(rep):]
        qual = qual[:rep[0]] + qual[sum(rep):]
    return seq, qual


def deduplicate(umis, reads_stats, fastq_path, out_prefix, cpus=None, v=0):
    """Deduplicate reads from UMI's clusters and write out fasta file

    Arguments:
    umis      (dict) - Dictionnary containing the UMIs as keys and and the ids
                       of the reads associated with each UMI.
    stats     (dict) - Dictionnary containing the reads ids as keys and
                       alignment statistics of the reads.
    fastq_path (str) - Path to the input fastq file with raw sequences.
    out_prefix (str) - Prefix for output files.
    cpus       (str) - Number of cpus for parallelization (optional).
    v          (int) - Level of verbosity (default: 0 = muted).
    """
    # Gather reads full sequences and quality strings from fastq
    with pysam.FastxFile(fastq_path) as fastq_stream:
        for read in fastq_stream:
            # Only retain reads that have their ids potentially used in
            # clusters
            if read.name in reads_stats.keys():
                if reads_stats[read.name]['orientation'] != 'reverse':
                    seq = read.sequence
                    qual = read.quality
                else:
                    seq = utils.rev_comp(read.sequence)
                    qual = read.quality[::-1]
                # Trim poly-A
                if reads_stats[read.name]['ref'] != 'unaligned':
                    seq, qual = trim_polyA(seq, qual)
                # Filter GA repeats
                seq, qual = filter_GA(seq, qual)
                # Register read sequence and quality
                reads_stats[read.name]['seq'] = seq
                reads_stats[read.name]['qual'] = qual

    # Define thread variable
    threads = str(cpus) if cpus is not None else '1'

    # Loop through reads clusters to compute consensus
    consensus = {}
    total = len(umis)
    tenth = total // 10
    error_counts = 0
    error_log = f'{out_prefix}_deduplication_error_log.txt'
    for umi, read_ids in umis.items():

        # Display progression
        partial = len(consensus)
        if partial % tenth == 0:
            utils.send_text(f'Deduplicated {partial} out of {total} UMIs.',
                            v, 1, 1)

        # If there is only one sequence in the cluster, job's done! :)
        if len(read_ids) == 1:
            consensus[umi] = reads_stats[read_ids[0]]['seq']
            continue

        # Compute draft alignment to then align reads to it
        # (if there are too many reads, use a subset for the draft reference)
        seed_reads = random.sample(read_ids, min(100, len(read_ids)))
        sequences = [reads_stats[read_id]['seq'] for read_id in seed_reads]
        draft_consensus, msa = spoa.poa(sequences)
        # Ensure there is only one sequence in draft result
        if not isinstance(draft_consensus, str):
            raise ValueError('Something went wrong with pyspoa consensus '
                             + f'building from {sequences}.')

        # Produce consensus using a traditional minimap2/racon pipeline
        # (iterate twice)
        for iteration in range(2):

            # Use temporary directory to store intermediate files
            with tempfile.TemporaryDirectory() as tmp:

                # Define file paths
                reads_fq = f'{tmp}/reads.fq'
                draft_fa = f'{tmp}/draft.fa'
                paf = f'{tmp}/aln.paf'
                polished = f'{tmp}/consensus.fa'

                # TODO: Remove
                # reads_fq = f'TEMP/reads.fq'
                # draft_fa = f'TEMP/draft.fa'
                # paf = f'TEMP/aln.paf'
                # polished = f'TEMP/consensus.fa'


                # Write reads to temporary fastq file
                with open(reads_fq, 'w') as tmp_fq:
                    for read_id in read_ids:
                        seq = reads_stats[read_id]['seq']
                        qual = reads_stats[read_id]['qual']
                        res = tmp_fq.write(f'@{read_id}\n{seq}\n+\n{qual}\n')

                # Write draft consensus to temporary fasta file
                with open(draft_fa, 'w') as tmp_fa:
                    res = tmp_fa.write(f'>draft_ref\n{draft_consensus}\n')

                # Align reads to draft reference
                try:
                    subprocess.run(
                        ['minimap2', '-x', 'map-ont', '-t', threads, draft_fa,
                        reads_fq],
                        stdout=open(paf, 'w'),
                        stderr=open(f'{out_prefix}_minimap2_log.txt', 'a'),
                        check=True
                    )
                except:
                    erroneous_reads = str(read_ids)
                    with open(error_log, 'a') as err:
                        err.write(
                            '#####################################\n'
                            + f'{umi} had minimap2 error. Reads:\n'
                            + f'{erroneous_reads}\n'
                        )
                    error_counts += 1
                    continue

                # Polish alignment to generate consensus
                try:
                    subprocess.run(
                        ['racon', '-t', threads, reads_fq, paf, draft_fa],
                        stdout=open(polished, 'w'),
                        stderr=open(f'{out_prefix}_racon_log.txt', 'a'),
                        check=True
                    )
                except:
                    erroneous_reads = str(read_ids)
                    with open(error_log, 'a') as err:
                        err.write(
                            '#####################################\n'
                            + f'{umi} had racon error. Reads:\n'
                            + f'{erroneous_reads}\n'
                        )
                    error_counts += 1
                    continue

                # Retrieve consensus sequence
                with open(polished) as tmp_consensus_file:
                    lines = tmp_consensus_file.read().splitlines()
                    draft_consensus = ''.join(lines[1:])
        
        consensus[umi] = draft_consensus
    
    # Write consensus sequences to fastq file
    utils.send_text(f'Deduplicated {total} UMIs with {error_counts} errors.',
                    v, 1, 1)
    utils.send_text(f'Writing consensus obtained to fasta file.', v, 1, 1)
    with open(f'{out_prefix}_dedup.fasta', 'w') as final_fasta:
        for umi, read_ids in umis.items():
            nb_reads = len(read_ids)
            sequence = consensus[umi]
            final_fasta.write(f'>{umi}-{nb_reads}\n{sequence}\n')
    
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input_fastq', type=str, required=True,
                        help='Path to input fastq file with raw reads.')
    parser.add_argument('-b', '--input_bam', type=str, required=True,
                        help='Path to input bam file with raw reads aligned to'
                        + ' the reference.')
    parser.add_argument('-o', '--output_prefix', type=str, default=None,
                        help='Prefix for output files. Default is the path to '
                        + 'the input fastq without the fastq extension. Note '
                        + 'that this can be used to redirect outputs to a '
                        + 'given directory. Output files are: '
                        + '_reads_stats.json, _alignment_heatmap.png, '
                        + '_overlap_heatmap.png, _umis_stats.tsv, _groups.tsv,'
                        ' (a tsv file similar to umi-tools group tsv output), '
                        + '_minimap2_log.txt, _racon_log.txt, '
                        + '_deduplication_error_log.txt')
    parser.add_argument('-r', '--raw_umis', type=str, required=True,
                        help='Path to input json file containing raw UMIs '
                        + 'information.')
    parser.add_argument('-s', '--reads_stats', type=str, default=None,
                        help='Path to input json file containing reads '
                        + 'statistics to skip the bam file parsing step.')
    parser.add_argument('-m', '--dist_method', type=str, default='Levenshtein',
                        help='Name of method to compute UMIs pairwise '
                        + 'distances. Must be Levenshtein or SmithWaterman.')
    parser.add_argument('-c', '--cores', type=int, default=None,
                        help='Number of CPUs available. If specified, distance'
                        + ' computation will be parralelized only if '
                        + 'SmithWaterman is chosen, because overhead is too '
                        + 'high for Levenshtein distance computations.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Define output path
    if args.output_prefix is not None:
        output = args.output_prefix
    else:
        output = args.input_fastq.replace('.fastq', '').replace('.fq', '')
    # Ensure directory exists
    output_dir = os.path.split(output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load raw UMIs
    utils.send_text(f'Loading raw UMIs statistics', v, 1, 0)
    umis = utils.load_json(args.raw_umis)

    # Skip bam parsing if information is already available
    if args.reads_stats is not None:
        utils.send_text(f'Loading reads statistics', v, 1, 0)
        reads_stats = utils.load_json(args.reads_stats)
    
    else:
        # Retrieve reads statistics
        utils.send_text(f'Retrieving reads statistics from {args.input_bam}',
                        v, 1, 0)
        reads_stats = retrieve_reads_stats(args.input_bam, v)
        # Save UMI and reads data json for future time save
        utils.send_text('Saving reads statistics information', v, 1, 0)
        utils.save_json(f'{output}_reads_stats.json', reads_stats)

    # Update UMIs to filter out reads that aren't in the statistics
    utils.send_text(f'Filtering out =unmapped reads from UMIs', v, 1, 0)
    umis = update_umis(umis, reads_stats)
    
    # Divide UMI-associated read groups that have distinct connexe components
    # in their coverage profile
    utils.send_text(f'Splitting umis on coverage connexe components', v, 1, 0)
    umis = split_connexe_coverage(umis, reads_stats, v)
    connex_umis_list = list(umis.keys())
    connex_umis_list.sort()

    # Build UMI sequences alignment-based distance matrix (and plot it)
    utils.send_text(f'Computing UMIs alignment distance matrix', v, 1, 0)
    alignment_dist_matrix = compute_alignment_dist(
        sequences=[umi_name.split('_')[0] for umi_name in connex_umis_list],
        method=args.dist_method,
        cpus=args.cores,
        v=v
    )
    res = plot_distance_matrix(
        alignment_dist_matrix,
        f'{output}_alignment_heatmap.png',
        v
    )

    # Build UMI reference coverage overlap-based distance matrix (and plot it)
    utils.send_text(f'Computing UMIs overlap distance matrix', v, 1, 0)
    overlap_dist_matrix = compute_overlap_dist(
        umi_names=connex_umis_list,
        umis=umis,
        stats=reads_stats,
        cpus=args.cores,
        v=v
    )
    res = plot_distance_matrix(
        overlap_dist_matrix,
        f'{output}_overlap_heatmap.png',
        v
    )

    # Cluster UMIs
    utils.send_text(f'Cluster UMIs using computed distances', v, 1, 0)
    umis = cluster_umis(
        umi_names=connex_umis_list,
        umis=umis,
        alignment_mat=alignment_dist_matrix,
        overlap_mat=overlap_dist_matrix,
        umi_stats_path=f'{output}_umis_stats.tsv',
        reads_stats=reads_stats,
        groups_path=f'{output}_groups.tsv',
        v=v
    )

    # Deduplicate reads
    utils.send_text(f'Deduplicating reads (calling consensus per UMI)',
                    v, 1, 0)
    res = deduplicate(
        umis,
        reads_stats,
        args.input_fastq,
        output,
        args.cores,
        v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
