#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Remove poly-A tails and plot trimmed tails lengths ditribution.
Takes as input the fastq file to be trimmed and output the fasta without poly-
A tail. Looks for poly-A tail at the end of reads, and for poly-T tails at the
start of reads. Consider single non-A nucleotides inside a poly-A tail to be a
sequencing/polymerization error. The tail must start with at least 'AAA'. Also
plots a histogram of the sizes of the trimmed extremities.
"""


import os
import sys

import argparse
import matplotlib.patches as pcs
import matplotlib.pyplot as plt

import utils # import local utility functions


def trim_tails(input, output, verbosity=0):
    """Trim the 3'-polyA / 5'-polyT tails.

    Positionnal arguments:
    input     (str) - path to the input fastq to trim reads from.
    output    (str) - path to the output fastq where to write trimmed reads to.
    verbosity (int) - level of verbosity.

    Return:
    trimmed_tail_lengths (list of int) - list of size of trimmed extremities
        (one read can have its two ends trimmed, if at least 3 As are in its 5
        3'-most bases AND at least 3 Ts are in its 5 5'-most bases).
    """

    # Open input fastq
    with open(input, 'r') as fastq_in:

        with open(output, 'w') as fastq_out:

            # Use counter for statistics and booleans for localization
            seq_line = False
            score_line = False
            trimmed_size = 0
            nb_trimmed_polyA = 0
            nb_trimmed_polyT = 0
            nb_trimmed_bases = 0
            trimmed_tail_lengths = []

            # Loop through input lines
            line = fastq_in.readline()
            while line:

                # Copy header and detect coming sequence
                if not score_line and line.startswith('@'):
                    seq_name = line
                    seq_line = True
                
                # Trimm and detect coming sequencing score
                elif seq_line:
                    seq = line.strip()
                    seq_size = len(seq)
                    right_trimmed = 0
                    left_trimmed = 0

                    # Trimm poly-A
                    trimmed_polyA = False
                    while seq[-5:].count('A') >= 3:

                        # Remove trailing A-s
                        if seq.endswith('A'):
                            trimmed_polyA = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq.rstrip('A') # Remove A trailing bases
                            right_trimmed += seq_size - len(seq) # Count the number of trimmed bases
                            seq_size = len(seq) # Refresh the size marker
                        
                        # Remove what is thought to be a falsely copied/sequenced A
                        else:
                            trimmed_polyA = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq[:-1] # Remove the last base
                            right_trimmed += 1 # Count the number of trimmed bases
                            seq_size -= 1 # Refresh the size marker

                    # Trimm poly-T
                    trimmed_polyT = False
                    while seq[:5].count('T') >= 3:

                        # Remove trailing T-s
                        if seq.startswith('T'):
                            trimmed_polyT = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq.lstrip('T') # Remove A trailing bases
                            left_trimmed += seq_size - len(seq) # Count the number of trimmed bases
                            seq_size = len(seq) # Refresh the size marker
                        
                        # Remove what is thought to be a falsely copied/sequenced T
                        else:
                            trimmed_polyT = True # Set boolean to True to later increment the trimmed poly-A counter
                            seq = seq[1:] # Remove the first base
                            left_trimmed += 1 # Count the number of trimmed bases
                            seq_size -= 1 # Refresh the size marker
                    
                    # Increment counters
                    if trimmed_polyA: nb_trimmed_polyA += 1
                    if trimmed_polyT: nb_trimmed_polyT += 1
                    nb_trimmed_bases += left_trimmed + right_trimmed
                    trimmed_tail_lengths.append(left_trimmed)
                    trimmed_tail_lengths.append(right_trimmed)

                    # Copy sequence if not empty
                    if seq != '':
                        fastq_out.write(seq_name)
                        fastq_out.write(seq + '\n+\n')
                        score_line = True
                        # Skip '+' line
                        line = fastq_in.readline()
                    
                    seq_line = False
                
                # Trim sequencing score accordingly
                elif score_line:
                    score = line.strip()[left_trimmed:]
                    if right_trimmed > 0: score = score[:-right_trimmed]
                    fastq_out.write(score + '\n')
                    score_line = False

                line = fastq_in.readline()
        
        # Output trimming statistics
        if verbosity > 0:
            print('Trimmed', nb_trimmed_polyA, "polyA ending tail.")
            print('Trimmed', nb_trimmed_polyT, "polyT ending tail.")
            print('Total trimmed bases:', nb_trimmed_bases)
    
    trimmed_tail_lengths.sort()
    
    return(trimmed_tail_lengths)


def plot_tail_histo(sizes, title, output_file, verbosity=0):
    """Plot the histogram of trimmed tail lengths.

    Positionnal arguments:
    sizes (list of int) - list of all tail sizes
    title         (str) - plot title
    output_file   (str) - path to output file
    log_scale    (bool) - use True for logarithmic scale on y axis
    verbosity     (int) - level of verbosity
    """

    plt.close() # ensure no figure is already open

    # Build figure
    bins = [1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121, 131, 141, 
        151, 201, 301, 401, 501]
    labels = [
        str(bins[i]) + '-' + str(bins[i+1]-1) for i in range(len(bins) - 1)
    ] + [str(bins[-1]) + '+']

    # Count the number of sizes in each bin
    counts = [0 for bin in bins]
    size_index = 0
    for bin_index in range(len(bins) - 1):
        while size_index < len(sizes) and sizes[size_index] < bins[bin_index + 1]:
            counts[bin_index] += 1
            size_index += 1
    if size_index < len(sizes):
        counts[-1] += len(sizes) - size_index

    fig, ax = plt.subplots() # initialize figure
    #dummy_var = fig.set_size_inches(10, 5) # (width, height) set figure size

    # plot tails sizes histogram
    dummy_var = ax.bar(labels, counts)
    
    ax.set_ylim([1, 1000000])
    ax.set_yscale('log')
    ax.set_ylabel('Number of trimmed tails')
    ax.set_xlabel('Trimmed tail lengths')
    ax.set_xticks(
        range(len(labels)),
        labels=labels,
        rotation=60,
        ha='right',
        fontdict={'fontsize': 7}
    )
    ax.set_title(title)

    # Save plot
    utils.send_text('Saving plot', verbosity, 2, 1)
    plt.tight_layout()
    plt.savefig(output_file, dpi=500)

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_fastq', type=str, required=True,
                        help='Path to the input fastq file.')
    parser.add_argument('-o', '--out_fastq', type=str,
                        help='Path to the output fastq file (default: <in_name'
                        + '>_no_tail.fastq)')
    parser.add_argument('-p', '--histo_png_out', type=str,
                        help='Path to the output histogram png file (default: '
                        + '<in_name>_no_tail_trimmed_tail_lengths.png)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Open output fastq to write to
    if args.out_fastq is None:
        out_path = args.in_fastq.rsplit('.fastq', 1)[0] + '_no_tail.fastq'
    else:
        out_path = args.out_fastq
    
    # Trim extremities and retrieve sizes of trimmed ends
    trimmed_tail_lengths = trim_tails(args.in_fastq, out_path, v)

    # Plot histogram of trimmed tail lengths
    if args.histo_png_out is None:
        histo_out = out_path[:-6] + '_trimmed_tail_lengths.png'
    else:
        histo_out = args.histo_png_out
    title = out_path.split('/')[-1][:-15] \
            + ' - ' \
            + str(len(trimmed_tail_lengths)) \
            + ' trimmed tails, ' \
            + str(sum(trimmed_tail_lengths)) \
            + ' trimmed bases'
    
    return(plot_tail_histo(trimmed_tail_lengths, title, histo_out, v))


if __name__ == "__main__":
    sys.exit(main())
