#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Remove GA-repeats and plot removed subsequences lengths ditribution.
Takes as input the fastq/fasta file with the read to be filtered and output the
fastq/fasta without GA-repeats. Filtered repeats are longer than a 
ser-specified threshold, and the user can chose to remove only the longest one.
Also plots a histogram of the sizes of the filtered GA-repeats.
"""


import os
import sys

import argparse
import matplotlib.patches as pcs
import matplotlib.pyplot as plt

import utils # import local utility functions


def filter_GA(input, output, threshold=10, longest=False, verbosity=0):
    """Filter GA repeats.

    Positionnal arguments:
    input     (str) - path to the input fastq/fasta to filter reads from.
    output    (str) - path to the output fastq/fasta where to write filtered
                      reads to
    threhsold (int) - minimum size of a repeat for it to get filtered (default
                      to 10)
    longest  (bool) - set to True to only filter the longest GA repeat in each
                      read (default to False)
    verbosity (int) - level of verbosity (default to 0 = mute)

    Return:
    GA_repeats_lengths (list of int) - list of size of filtered repeats
    """

    # Open input fastq/fasta
    with open(input, 'r') as reads_in:

        with open(output, 'w') as reads_out:

            # Use counter for statistics and booleans for localization
            seq_line = False
            score_line = False
            trimmed_size = 0
            nb_trimmed_repeats = 0
            nb_trimmed_bases = 0
            GA_repeats_lengths = []

            # Loop through input lines
            line = reads_in.readline()

            # Decide whether we deal with fastq or fasta
            if line[0] == '@':
                read_format = 'fastq'
                starter = '@'
            elif line[0] == '>':
                read_format = 'fasta'
                starter = '>'
            else:
                raise ValueError(input, 'is no valid fastq/fasta format.')

            while line:

                # Copy header and detect coming sequence
                if not score_line and line.startswith(starter):
                    seq_name = line
                    seq_line = True
                
                # Trimm and detect coming sequencing score
                elif seq_line:
                    seq = line.strip()
                    right_trimmed = 0
                    left_trimmed = 0

                    # Look for GA repeats (store list of repeats start and len)
                    repeats = []
                    longest_repeat = []
                    for ind, char in enumerate(seq):
                        # Test for GA
                        if ind < len(seq)-1 and seq[ind:ind+2] == 'GA':
                            # Test if still in previously detected repeat
                            if repeats and ind == repeats[-1][0] + repeats[-1][1]:
                                # update current repeat...
                                repeats[-1][1] += 2
                            # ...else, start new repeat
                            else:
                                repeats.append([ind, 2])
                            # Update longest repeat
                            if longest_repeat and repeats[-1][1] > longest_repeat[1]:
                                longest_repeat = repeats[-1]

                    # Limit to longest repeat if requested
                    if longest:
                        repeats = longest_repeat

                    # Filter repeats based on their size
                    repeats = [rep for rep in repeats if rep[1] >= threshold]

                    # Filter the read based on the repeats
                    for rep in repeats[::-1]:
                        seq = seq[:rep[0]] + seq[rep[0] + rep[1]:]
                    
                    # Increment counters
                    if repeats:
                        nb_trimmed_repeats += len(repeats)
                        nb_trimmed_bases += sum([rep[1] for rep in repeats])
                        GA_repeats_lengths += [rep[1] for rep in repeats]

                    # Copy sequence if not empty
                    if seq != '':
                        reads_out.write(seq_name)
                        reads_out.write(seq + '\n')
                        if read_format == 'fastq':
                            reads_out.write('+\n') # write line with plus sign
                            score_line = True
                            line = reads_in.readline() # skip '+' line
                    
                    seq_line = False
                
                # Filter sequencing score accordingly if fastq format
                elif score_line:
                    score = line.strip()
                    for rep in repeats[::-1]:
                        score = score[:rep[0]] + score[rep[0] + rep[1]:]
                    reads_out.write(score + '\n')
                    score_line = False

                line = reads_in.readline()
        
        # Output GA repeats filtering statistics
        if verbosity > 0:
            print('Filtered', nb_trimmed_repeats, "GA repeats.")
            print('Total filtered bases:', nb_trimmed_bases)
    
    GA_repeats_lengths.sort()
    
    return(GA_repeats_lengths)


def plot_repeats_histo(sizes, title, output_file, verbosity=0):
    """Plot the histogram of filtered GA repeats lengths.

    Positionnal arguments:
    sizes (list of int) - list of all repeats sizes
    title         (str) - plot title
    output_file   (str) - path to output file
    log_scale    (bool) - use True for logarithmic scale on y axis
    verbosity     (int) - level of verbosity (default to 0 = mute)
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

    # plot repeats sizes histogram
    dummy_var = ax.bar(labels, counts)
    
    ax.set_ylim([1, 1000000])
    ax.set_yscale('log')
    ax.set_ylabel('Number of filtered repeats')
    ax.set_xlabel('Filtered GA repeats lengths')
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
    parser.add_argument('-i', '--in_reads', type=str, required=True,
                        help='Path to the input fastq/fastq file.')
    parser.add_argument('-o', '--out_reads', type=str,
                        help='Path to the output fastq/fastq file (default: '
                        + '<input_name>_no_GA with same extension suffix)')
    parser.add_argument('-t', '--threshold', type=int, default=10,
                        help='Min size to filter a GA repeat (default: 10)')
    parser.add_argument('-l', '--longest', action='store_true',
                        help='Use to only filter the longest repeat in each '
                        + 'read')
    parser.add_argument('-p', '--histo_png_out', type=str,
                        help='Path to the output histogram png file (default: '
                        + '<input_name>_no_GA_filtered_lengths.png)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Open output fastq to write to
    if args.out_reads is None:
        ext = args.in_reads.rsplit('.', 1)[1]
        out_path = args.in_reads.rsplit('.', 1)[0] + '_no_GA.' + ext
    else:
        out_path = args.out_reads
    
    # Filter GA repeats and retrieve sizes of filtered subsequences
    GA_repeats_lengths = filter_GA(
        args.in_reads, out_path, args.threshold, args.longest, v
    )

    # Plot histogram of filtered GA repeats lengths
    if args.histo_png_out is None:
        histo_out = out_path.rsplit('.', 1)[0] + '_filtered_lengths.png'
    else:
        histo_out = args.histo_png_out
    title = out_path.split('/')[-1].rsplit('_no_GA', 1)[0] \
            + ' - ' \
            + str(len(GA_repeats_lengths)) \
            + ' filtered repeats, ' \
            + str(sum(GA_repeats_lengths)) \
            + ' filtered bases'
    
    return(plot_repeats_histo(GA_repeats_lengths, title, histo_out, v))


if __name__ == "__main__":
    sys.exit(main())
