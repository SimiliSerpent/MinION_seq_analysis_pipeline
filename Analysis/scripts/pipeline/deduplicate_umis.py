#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Deduplicate reads assigned to a single original molecule.
Retrieve reads associated to closely related UMIs (in terms of sequence
similarity and mapping proximity) and thus considered originating from the same
original molecule and call the consensus of all thoses reads.
"""


import argparse
import functools
import glob
import multiprocessing
import os
import random
import subprocess
import sys
import tempfile
import time

import pysam
import spoa

import utils

# Module-level globals populated by the worker initializer. Workers access
# read sequences and quality strings through these to avoid pickling the
# whole dicts per task. On Linux, fork makes this near-free; on spawn-based
# platforms, dicts are pickled once per worker at pool startup.
WORKER_SEQS = None
WORKER_QUALS = None


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


def init_consensus_worker(seqs, quals):
    """Pool initializer: share sequence/quality dicts read-only with workers.
    """
    global WORKER_SEQS, WORKER_QUALS
    WORKER_SEQS = seqs
    WORKER_QUALS = quals


def call_consensus_for_cluster(umi_and_ids, polish_iterations,
                               polish_min_reads, log_prefix, threads=1):
    """Build a consensus for one UMI cluster.
    Runs SPOA, then optionally one or more rounds of minimap2 + racon polish.

    Arguments:
    umi_and_ids     (tuple) - (umi, read_ids): umi = Name of the UMI associated
                              with the reads cluster used for consensus
                              calling; read_ids =  List of ids of reads used
                              for consensus calling.
    polish_iterations (int) - Number of minimap2+racon polish iterations to run
                              on top of the SPOA draft.
    polish_min_reads  (int) - Clusters with strictly fewer reads than this skip
                              polish entirely - the SPOA draft is returned
                              as-is.
    log_prefix        (str) - Prefix for output log files.
    threads           (int) - Threads passed to minimap2 and racon.

    Return:
    umi           (str) - Name of the UMI associated with the reads cluster
                          used for consensus calling.
    consensus_seq (str) - Consensus sequence called. None only if SPOA failed;
                          partial polish failures still return the best draft
                          obtained so far.
    error_msg     (str) - Error log message. None on full success.
    """
    umi, read_ids = umi_and_ids
    seqs, quals = WORKER_SEQS, WORKER_QUALS

    # SPOA draft from up to 100 reads (the rest of the cluster is only
    # used in polishing if any)
    if len(read_ids) <= 100:
        seed_reads = read_ids
    else:
        seed_reads = random.sample(read_ids, 100)
    sequences = [seqs[read_id] for read_id in seed_reads]
    try:
        consensus_seq, _msa = spoa.poa(sequences)
    except Exception as e:
        return (umi, None, f'SPOA error: {e}')
    if not isinstance(consensus_seq, str):
        return (umi, None, 'SPOA returned non-string consensus')

    # Skip polish for small clusters - SPOA alone is good enough on a handful
    # of ONT reads, and the minimap2+racon fork overhead dominates the actual
    # work below this threshold
    if polish_iterations == 0 or len(read_ids) < polish_min_reads:
        return (umi, consensus_seq, None)

    # Per-PID aligner logs avoid contention between workers; they get
    # concatenated into the canonical log files in the parent at the end
    pid = os.getpid()
    minimap2_log = f'{log_prefix}_minimap2_log.{pid}.txt'
    racon_log = f'{log_prefix}_racon_log.{pid}.txt'

    error_msg = None
    with tempfile.TemporaryDirectory() as tmp:
        reads_fq = f'{tmp}/reads.fq'
        draft_fa = f'{tmp}/draft.fa'
        paf = f'{tmp}/aln.paf'
        polished = f'{tmp}/consensus.fa'

        # Reads don't change across polish iterations - write them once
        with open(reads_fq, 'w') as tmp_fq:
            for read_id in read_ids:
                tmp_fq.write(
                    f'@{read_id}\n{seqs[read_id]}\n+\n{quals[read_id]}\n'
                )

        for iteration in range(polish_iterations):
            with open(draft_fa, 'w') as tmp_fa:
                tmp_fa.write(f'>draft_ref\n{consensus_seq}\n')

            # Align reads to draft
            try:
                with open(paf, 'w') as paf_out, \
                     open(minimap2_log, 'a') as mm2_err:
                    subprocess.run(
                        ['minimap2', '-x', 'map-ont', '-t', str(threads),
                         draft_fa, reads_fq],
                        stdout=paf_out, stderr=mm2_err, check=True
                    )
            except Exception as e:
                error_msg = f'minimap2 error at iter {iteration}: {e}'
                break

            # Polish alignment to generate consensus
            try:
                with open(polished, 'w') as pol_out, \
                     open(racon_log, 'a') as rc_err:
                    subprocess.run(
                        ['racon', '-t', str(threads), reads_fq, paf, draft_fa],
                        stdout=pol_out, stderr=rc_err, check=True
                    )
            except Exception as e:
                error_msg = f'racon error at iter {iteration}: {e}'
                break

            # Retrieve polished consensus
            with open(polished) as pol_in:
                lines = pol_in.read().splitlines()
                consensus_seq = ''.join(lines[1:])

    return umi, consensus_seq, error_msg


def deduplicate(umis, reads_stats, fastq_path, out_prefix, cpus=None,
                polish_iter=1, polish_min=6, polish_big_min=200, v=0):
    """Deduplicate reads from UMI's clusters and write out fasta file

    Arguments:
    umis       (dict) - Dictionnary containing the UMIs as keys and the ids
                       of the reads associated with each UMI.
    stats      (dict) - Dictionnary containing the reads ids as keys and
                       alignment statistics of the reads.
    fastq_path  (str) - Path to the input fastq file with raw sequences.
    out_prefix  (str) - Prefix for output files.
    cpus        (str) - Number of cpus for parallelization (optional).
    polish_iter (int) - Number of minimap2+racon polish iterations to run on
                        top of the SPOA draft (default: 1).
    polish_min  (int) - Clusters with strictly fewer reads than this skip
                        polish entirely - the SPOA draft is returned as-is
                        (default: 6).
    polish_big_min (int) - Cluster-size cutoff between "small" and "big"
                           multi-read clusters (default: 200). Small clusters
                           are polished in a parallel pool with 1 thread per
                           worker; big clusters are polished serially with
                           `cpus` threads each.
    v           (int) - Level of verbosity (default: 0 = muted).
    """
    # Gather reads full sequences and quality strings from fastq, in dedicated
    # dicts so the per-worker memory footprint stays small (workers only need
    # seqs/quals, not the alignment metadata in reads_stats)
    seqs = {}
    quals = {}
    if v > 0:
        t_zero = time.perf_counter()
    with pysam.FastxFile(fastq_path) as fastq_stream:
        for read in fastq_stream:
            # Only retain reads that have their ids potentially used in
            # clusters
            if read.name not in reads_stats:
                continue
            stats = reads_stats[read.name]
            if stats['orientation'] != 'reverse':
                seq = read.sequence
                qual = read.quality
            else:
                seq = utils.rev_comp(read.sequence)
                qual = read.quality[::-1]
            # Trim poly-A
            if stats['ref'] != 'unaligned':
                seq, qual = trim_polyA(seq, qual)
            # Filter GA repeats
            seq, qual = filter_GA(seq, qual)
            # Register read sequence and quality
            seqs[read.name] = seq
            quals[read.name] = qual
    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Retrieved reads properties from fastq in {elapsed:.3f} s.',
            v, 2, 1
        )
    
    # Compute trivial single-read consensuses up front (no IPC overhead) and
    # build work items for every multi-read cluster for later multi-threaded
    # consensus calling
    consensus = {}
    work_items = []
    for umi, read_ids in umis.items():
        if len(read_ids) == 1:
            consensus[umi] = seqs[read_ids[0]]
        else:
            work_items.append((umi, read_ids))
            #     (umi, read_ids, polish_iterations, polish_min_reads,
            #      out_prefix)
            # )
    
    # Define partial function for the consensus calling function to only change
    # arguments specific to one reads cluster
    call_consensus_partial = functools.partial(
        call_consensus_for_cluster,
        polish_iterations=polish_iter,
        polish_min_reads=polish_min,
        log_prefix=out_prefix,
    )

    # Split multi-read clusters into "big" and "small". Big clusters are
    # processed serially in the main process with full CPU budget per cluster
    # because for a tail-dominated workload the single biggest cluster gates
    # wall time, and racon scales better with threads than with worker count.
    # Small clusters are processed in a worker pool with 1 thread each.
    big_items = [item for item in work_items
                 if len(item[1]) >= polish_big_min]
    small_items = [item for item in work_items
                   if len(item[1]) < polish_big_min]
    big_items.sort(key=lambda item: len(item[1]), reverse=True)

    total = len(umis)
    n_workers = cpus if cpus is not None else 1
    utils.send_text(
        f'{total} UMIs: {len(consensus)} singletons, {len(small_items)} small '
        + f'multi-read clusters (parallel, 1 thread each), {len(big_items)} '
        + f'large clusters with >= {polish_big_min} reads (serial, {n_workers}'
        + ' threads each).',
        v, 2, 1
    )

    error_counts = 0
    error_log = f'{out_prefix}_deduplication_error_log.txt'

    # Define helper function to process outputs of consensus calling
    def handle_result(idx, total_for_label, label, umi, cons, err):
        nonlocal error_counts
        if err is not None:
            error_counts += 1
            with open(error_log, 'a') as err_out:
                err_out.write(
                    '#####################################\n'
                    + f'{umi}: {err}\n'
                )
        if cons is not None:
            consensus[umi] = cons
        report_step = max(1, total_for_label // 20)
        if (idx + 1) % report_step == 0:
            utils.send_text(
                f'Processed {idx+1} / {total_for_label} {label} clusters.',
                v, 2, 1
            )
    
    # Big clusters: serial, all cores per cluster
    if big_items:
        if v > 0:
            t_zero = time.perf_counter()
        init_consensus_worker(seqs, quals)
        # Define partial function for the consensus calling function to only
        # change arguments specific to one reads cluster
        big_partial = functools.partial(
            call_consensus_for_cluster,
            polish_iterations=polish_iter,
            polish_min_reads=polish_min,
            log_prefix=out_prefix,
            threads=n_workers,
        )
        for i, item in enumerate(big_items):
            umi, cons, err = big_partial(item)
            handle_result(i, len(big_items), 'large', umi, cons, err)
        if v > 0:
            elapsed = time.perf_counter() - t_zero
            utils.send_text(
                f'Polished {len(big_items)} large clusters in '
                + f'{elapsed:.3f} s.',
                v, 2, 1
            )
    
    # Small clusters: parallel pool, 1 thread each
    if small_items:
        if v > 0:
            t_zero = time.perf_counter()
        chunksize = max(1, len(small_items) // (4 * n_workers))
        # Again, partial function with single thread
        small_partial = functools.partial(
            call_consensus_for_cluster,
            polish_iterations=polish_iter,
            polish_min_reads=polish_min,
            log_prefix=out_prefix,
            threads=1,
        )
        if cpus is not None and cpus > 1:
            with multiprocessing.Pool(
                cpus,
                initializer=init_consensus_worker,
                initargs=(seqs, quals)
            ) as pool:
                for i, (umi, cons, err) in enumerate(pool.imap_unordered(
                    small_partial, small_items, chunksize=chunksize
                )):
                    handle_result(i, len(small_items), 'small',
                                  umi, cons, err)
        else:
            init_consensus_worker(seqs, quals)
            for i, item in enumerate(small_items):
                umi, cons, err = small_partial(item)
                handle_result(i, len(small_items), 'small', umi, cons, err)
        if v > 0:
            elapsed = time.perf_counter() - t_zero
            utils.send_text(
                f'Polished {len(small_items)} small clusters in '
                + f'{elapsed:.3f} s.',
                v, 2, 1
            )

    if big_items or small_items:
        # Concatenate per-PID aligner logs into the canonical log files,
        # then clean up
        for tool in ('minimap2', 'racon'):
            main_log = f'{out_prefix}_{tool}_log.txt'
            per_pid = sorted(glob.glob(f'{out_prefix}_{tool}_log.*.txt'))
            if not per_pid:
                continue
            with open(main_log, 'a') as main_out:
                for fname in per_pid:
                    with open(fname) as f_in:
                        main_out.write(f_in.read())
                    os.remove(fname)
    
    # Write consensus sequences to fastq file
    utils.send_text(f'Writing consensus obtained to fasta file.', v, 2, 1)
    with open(f'{out_prefix}_dedup.fasta', 'w') as final_fasta:
        for umi, read_ids in umis.items():
            if umi not in consensus:
                continue  # SPOA failed on this cluster, nothing to write
            nb_reads = len(read_ids)
            sequence = consensus[umi]
            final_fasta.write(f'>{umi}-{nb_reads}\n{sequence}\n')
    
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input_fastq', type=str, required=True,
                        help='Path to input fastq file with raw reads.')
    parser.add_argument('-o', '--output_prefix', type=str, required=True,
                        help='Prefix for output files.')
    parser.add_argument('-u', '--clustered_umis', type=str, required=True,
                        help='Path to input json file containing clustered '
                        'UMIs information.')
    parser.add_argument('-s', '--reads_stats', type=str, required=True,
                        help='Path to input json file containing reads '
                        'alignment statistics.')
    parser.add_argument('-i', '--polish_iterations', type=int, default=1,
                        help='Number of minimap2+racon polish iterations to '
                        'run on top of the SPOA draft consensus. Set to 0 to '
                        'skip polish entirely and return the SPOA draft, '
                        'which is much faster but slightly lower quality. '
                        + 'Default: 1.')
    parser.add_argument('-m', '--polish_min_reads', type=int, default=6,
                        help='Clusters with strictly fewer reads than this '
                        + 'skip the minimap2+racon polish step entirely - the '
                        + 'SPOA draft is returned as-is. SPOA alone is good '
                        + 'enough on a handful of ONT reads, and avoiding the '
                        'fork overhead of minimap2/racon on tiny clusters '
                        'is a major speedup. Default: 6.')
    parser.add_argument('-c', '--cores', type=int, default=None,
                        help='Number of CPUs available.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output files do not start with an underscore
    output = args.output_prefix
    if output.endswith('/'):
        output += 'dedup_out'
    # Ensure output directory exists
    output_dir = os.path.split(output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load reads statistics
    utils.send_text(f'Loading reads statistics', v, 1, 0)
    reads_stats = utils.load_json(args.reads_stats)

    # Load clustered UMIs
    utils.send_text(f'Loading clustered UMIs', v, 1, 0)
    clustered_umis =  utils.load_json(args.clustered_umis)

    # Deduplicate reads
    utils.send_text(f'Deduplicating reads (calling consensus per UMI)',
                    v, 1, 0)
    res = deduplicate(
        clustered_umis,
        reads_stats,
        args.input_fastq,
        output,
        args.cores,
        polish_iter=args.polish_iterations,
        polish_min=args.polish_min_reads,
        v=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
