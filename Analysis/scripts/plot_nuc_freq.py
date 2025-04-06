#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot nucleotides frequencies in raw and trimmed sequencies.
Takes as input the directory with nucleotides frequencies (as output by the
compute_nuc_freq.py script), and outputs a png graphic showing the evolution
of nucleotides frequencies before and after trimming for all found samples.
"""


import os
import sys

import argparse
import json
import matplotlib.patches as pcs
import matplotlib.pyplot as plt
import numpy as np

import utils # import local utility functions


MY_NUC_COLORS = {
    'A': 'tomato',
    'T': 'springgreen',
    'G': 'gold',
    'C': 'mediumslateblue'
}


def get_samples_freq(nuc_freq_dir):
    """Get the list of samples and their nucleotides frequencies.

    Positionnal arguments:
    nuc_freq_dir (str) - path to directory containing the 'trimmed_fastq' sub-
        directory with the nucleotides frequencies '..._freq.txt' files.

    Return:
    samples_freq_dict (dict) - dictionnary containing the samples frequencies
        in nucleotides.
    """

    # Initialize dictionnary with samples
    samples_freq_dict = {}

    # Loop through files and parse the names to find samples and targets
    for file_path in utils.find_file(nuc_freq_dir + '/trimmed_fastq', '_freq.txt'):
        if file_path.endswith('_freq.txt'):

            # Parse the name to retrieve sample name
            sample = file_path.split('/')[-1].split('_freq.')[0]

            # Open file and retrieve the number of mapped reads
            with open(nuc_freq_dir+'/trimmed_fastq/'+file_path, 'r') as file:
                samples_freq_dict[sample] = {
                    line[0]: float(line.split('\t')[1].split(' ')[0]) for line in file.readlines()
                }

    return samples_freq_dict


def plot_frequencies(samples, frequencies, bar_width, bar_colors,
    legend_elements, title, output_file, y_label='Nucleotide frequency (%)',
    verbosity=0):
    """Returns a dict containing the freq. of each nuc. for each sample.

    Positionnal arguments:
    samples          (list) - sorted list of all samples
    frequencies      (dict) - dictionnary containing the freq. of nuc. per
        sample
    bar_width       (float) - width of bars in the figure
    bar_colors       (dict) - dictionnary containing one color per nucleotide
    legend_elements  (list) - list of the patch of colors and the name of each
        specie
    title             (str) - plot title
    output_file       (str) - path to output file
    log_scale        (bool) - use True for logarithmic scale on y axis
    y_label           (str) - y axis name
    verbosity         (int) - level of verbosity
    """

    plt.close() # ensure no figure is already open
    fig, ax = plt.subplots() # initialize figure
    dummy_var = fig.set_size_inches(10, 5) # (width, height) set figure size

    # plot number of reads per species
    multiplier = 0
    for nuc in 'ATGC':
        offset = multiplier * bar_width
        dummy_var = ax.bar(
            np.arange(len(samples)) * ((len('ATGC') + 2.5) * bar_width) + offset,
            [frequencies[sample][nuc] for sample in samples],
            bar_width,
            color=bar_colors[nuc]
        )
        multiplier += 1

    # ax.set_ylim([1, 1000000])
    ax.set_ylabel(y_label)
    # ax.set_xlabel('Read size')
    ax.set_xticks(
        np.arange(len(samples)) * ((len('ATGC') + 2.5) * bar_width) + (((len('ATGC') - 1) * bar_width) / 2),
        labels=samples,
        rotation=60,
        ha='right'
    ) # position labels along x axis
    ax.set_title(title)
    ax.legend(handles=legend_elements, framealpha=1)
    ax.set_axisbelow(True) # prevent grid from being printed ahead of figure
    ax.yaxis.grid(color='lightgray') # plot the horizontal grey lines

    # Save plot
    utils.send_text('Saving plot', verbosity, 2, 2)
    plt.tight_layout()
    plt.savefig(output_file,
                dpi=500)

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', type=str, required=True,
                        help='Path to the input directory containing the '
                        + 'nucleotide frequency per sample (files must be '
                        + 'located in sub-directory \"trimmed_fastq\")')
    parser.add_argument('-o', '--output', type=str,
                        help='Path to the output directory (and prefix) '
                        '(default: same directory)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Retrieve the samples composition in nucleotides
    samples_freq_dict = get_samples_freq(args.input_dir.rstrip('/'))

    # Retrieve sorted list of all samples
    samples = [key for key in samples_freq_dict.keys()]
    samples.sort()

    # Determine output file location
    output_path = args.input_dir.rstrip('/') + '/'
    if args.output is not None:
        output_path = args.output

    # Build figures
    bar_width = 0.8 # set width or bars
    legend_elements = [
        pcs.Patch(facecolor=MY_NUC_COLORS[nuc], label=nuc) for nuc in 'ATGC'
    ] # set patch of colors to be used in the figure legend

    # Plot samples composition in terms of nucleotides
    title = 'Nucleotide frequencies'
    dummy_var = plot_frequencies(
        samples,
        samples_freq_dict,
        bar_width,
        MY_NUC_COLORS,
        legend_elements,
        title,
        output_file=output_path + '_nuc_freq.png',
        verbosity=v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
