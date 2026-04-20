#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot statistics of adapters, tails and repeats cleaning.
Takes as input a directory and recursively looks for statistics of adapters
trimming, tail cleaning and repeats filtering. Load them and plots the
statistics for all the samples identified.
"""


import os
import sys

import argparse
import json
import matplotlib.patches as pcs
import matplotlib.pyplot as plt
import numpy as np
import re

import utils # import local utility functions

# Define the different cleaning categories
C_TYPES = [
    'seq_adapt',
    'twist_outer',
    'twist_primers',
    'ont_primer_tails',
    'ont_unknown_seq',
    'tso',
    'ONT_barcode',
    'ONT_sequences',
    'polyA/T_tails',
    'GA_repeats'
]

# Define the color for each cleaning type and statistic
PLOT_COLORS = {
    'seq_adapt': {
        'title_background': 'papayawhip',
        'start_trim_seq': 'chocolate',
        'start_trim_bases': 'orangered',
        'end_trim_seq': 'lightsalmon',
        'end_trim_bases': 'peachpuff',
        'split': 'firebrick'
    },
    'twist_outer': {
        'title_background': 'thistle',
        'start_trim_seq': 'purple',
        'start_trim_bases': 'darkorchid',
        'end_trim_seq': 'orchid',
        'end_trim_bases': 'violet',
        'split': 'mediumvioletred'
    },
    'twist_primers': {
        'title_background': 'thistle',
        'start_trim_seq': 'purple',
        'start_trim_bases': 'darkorchid',
        'end_trim_seq': 'orchid',
        'end_trim_bases': 'violet',
        'split': 'mediumvioletred'
    },
    'ONT_barcode': {
        'title_background': 'papayawhip',
        'start_trim_seq': 'chocolate',
        'start_trim_bases': 'orangered',
        'end_trim_seq': 'lightsalmon',
        'end_trim_bases': 'peachpuff',
        'split': 'firebrick'
    },
    'ont_primer_tails': {
        'title_background': 'papayawhip',
        'start_trim_seq': 'chocolate',
        'start_trim_bases': 'orangered',
        'end_trim_seq': 'lightsalmon',
        'end_trim_bases': 'peachpuff',
        'split': 'firebrick'
    },
    'ont_unknown_seq': {
        'title_background': 'papayawhip',
        'start_trim_seq': 'chocolate',
        'start_trim_bases': 'orangered',
        'end_trim_seq': 'lightsalmon',
        'end_trim_bases': 'peachpuff',
        'split': 'firebrick'
    },
    'ONT_sequences': {
        'title_background': 'papayawhip',
        'start_trim_seq': 'chocolate',
        'start_trim_bases': 'orangered',
        'end_trim_seq': 'lightsalmon',
        'end_trim_bases': 'peachpuff',
        'split': 'firebrick'
    },
    'tso': {
        'title_background': 'palegreen',
        'start_trim_seq': 'darkolivegreen',
        'start_trim_bases': 'forestgreen',
        'end_trim_seq': 'xkcd:irish green',
        'end_trim_bases': 'xkcd:seaweed'
    },
    'polyA/T_tails': {
        'title_background': 'lavender',
        'start_trim_seq': 'mediumblue',
        'start_trim_bases': 'royalblue',
        'end_trim_seq': 'skyblue',
        'end_trim_bases': 'lightskyblue'
    },
    'GA_repeats': {
        'title_background': 'xkcd:ecru',
        'sequences': 'gold',
        'bases': 'yellow'
    }
}

# Define the labels for each cleaning type
C_LABELS = {
    'seq_adapt': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'twist_outer': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'twist_primers': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'ONT_barcode': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'ont_primer_tails': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'ont_unknown_seq': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'ONT_sequences': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases',
        'split'
    ],
    'tso': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases'
    ],
    'polyA/T_tails': [
        'start_trim_seq',
        'start_trim_bases',
        'end_trim_seq',
        'end_trim_bases'
    ],
    'GA_repeats': [
        'sequences',
        'bases'
    ]
}

# Define the text to plot for each cleaning step
C_CLEAN_NAMES = {
    'seq_adapt': 'Sequencing adapters',
    'twist_outer': 'Twist outer univ. sequences',
    'twist_primers': 'Twist\nprimers',
    'ont_primer_tails': 'ONT primer tails',
    'ont_unknown_seq': 'ONT unknown sequences',
    'tso': 'Takara TSO',
    'ONT_barcode': 'ONT barcodes',
    'ONT_sequences': 'ONT sequences',
    'polyA/T_tails': 'polyA/T tails',
    'GA_repeats': 'GA repeats'
}

FONTSIZE = 7


def get_cleaning_stats(statistics_dir):
    """Get the statistics of reads cleaning from all the statistics files.

    Positionnal arguments:
    statistics_dir (str) - path to directory with reads cleaning statistics

    Return:
    clean_stats   (dict) - dictionnary with the reads cleaning statistics
    """

    # Find the paths to the statistics files
    clean_paths = utils.find_path(statistics_dir, 'porechopping_stats.json')
    clean_paths += utils.find_path(statistics_dir, '_polyA_trimming_stats.json')
    clean_paths += utils.find_path(statistics_dir, '_GA_filtering_stats.json')
    
    ### Store the statistics according to their experiments

    # Initialize statistics dictionnary
    clean_stats = {}

    # Loop through paths
    for c_path in clean_paths:

        # Identify the experiment (if none, use the "empty string" experiment)
        no_exp_found = True
        exp_pattern = re.compile(r'^\d{8}_[A-Z0-9]+_\d+$') # match an exp name
        # Find the exp name in the path
        for dir in c_path.split('/'):
            if bool(exp_pattern.match(dir)):
                exp_name = dir
                no_exp_found = False
        if no_exp_found:
            exp_name = ''

        # Initialize the experiment sub-dictionnary if not done before
        if not exp_name in clean_stats.keys():
            clean_stats[exp_name] = {}
        
        # Load the cleaning statistics
        clean_stat = utils.load_json(c_path)
        for barcode in clean_stat.keys():
            # Initialize barcode if not done before
            if not barcode in clean_stats[exp_name].keys():
                clean_stats[exp_name][barcode] = {}
            # Add the cleaning stats to the global dictionnary
            for adapter in clean_stat[barcode].keys():
                new_stat = clean_stat[barcode][adapter]
                clean_stats[exp_name][barcode][adapter] = new_stat
            
    return(clean_stats)


def plot_cleaning_stats(clean_stats, output_file, sample_names=None,
    c_types=None, bar_width=0.8, total=False, log_scale=False, verbosity=0):
    """Plost the statistics of reads cleaning.

    Positionnal arguments:
    clean_stats    (dict) - dictionnary containing the cleaning statistics for
        each barcode of each found experiment
    output_file     (str) - path to output file
    sample_names   (dict) - dictionnary containing the name of each barcode for
        any found experiment - default is None - if specified, only the named
        barcode will appear
    c_types (list of str) - listo of c_types to use - if none, all c_types are
        used
    bar_width     (float) - width of bars in the figure
    output_file     (str) - path to output file
    total          (bool) - plot total reads as well
    log_scale      (bool) - use True for logarithmic scale on y axis
    verbosity       (int) - level of verbosity
    """

    # Use all cleaning types if none are specified
    if c_types is None:
        c_types = C_TYPES

    # Retrieve the experiments names
    if sample_names is not None:
        # If sample names have been provided, only use the experiments for
        # which it is the case, if they are available in the gathered cleaning
        # statistics json
        exps = [
            exp for exp in sample_names.keys() if exp in clean_stats.keys()
        ]
    else:
        exps = [exp for exp in clean_stats.keys()]
    exps.sort()

    # Retrieve the barcodes
    barcodes = {}
    if sample_names is not None:
        # Same as above in the case sample names have been provided
        for exp in exps:
            bcs = [
                bc for bc in sample_names[exp].keys() if bc in clean_stats[exp].keys()
            ]
            bcs.sort()
            barcodes[exp] = bcs
    else:
        for exp in exps:
            bcs = [bc for bc in clean_stats[exp].keys()]
            bcs.sort()
            barcodes[exp] = bcs
    
    # Restrict cleaning types to those actually found in the cleaning stats
    available_c_types = []
    for exp in exps:
        for bc in barcodes[exp]:
            available_c_types += [c for c in clean_stats[exp][bc].keys()]
    c_types = [c for c in c_types if c in available_c_types]
    
    # Define the indices to plot all cleaning stats
    c_indices = [
        [i+1 for i in range(len(C_LABELS[c]))] for c in c_types
    ]
    # Increase the indices for each cleaning type given the other clean types
    # already plotted on the left
    for c_i in range(1, len(c_indices)):
        c_indices[c_i] = [
            c_indices[c_i - 1][-1] + 2 + i for i in c_indices[c_i]
        ]
    # Concatenate the indices of all the cleaning types
    bar_x = []
    for c_ind in c_indices:
        bar_x += c_ind

    # Define the x labels for all the plots
    bar_labels = []
    for c_type in c_types:
        bar_labels += C_LABELS[c_type]
    
    # Define the x colors for all the plots
    bar_col = []
    for c_type in c_types:
        for label in C_LABELS[c_type]:
            bar_col.append(PLOT_COLORS[c_type][label])
    
    # Split the previously defined list between sequences and bases, since 
    # these two will not have the same scale (and thus we have to plot to 
    # different bar figures for each subplot)
    bar_x_seq, bar_x_base = [], []
    bar_lab_seq, bar_lab_base = [], []
    bar_col_seq, bar_col_base = [], []
    for lab_ind, label in enumerate(bar_labels):
        if 'bases' in label:
            bar_x_base.append(bar_x[lab_ind])
            bar_lab_base.append(bar_labels[lab_ind])
            bar_col_base.append(bar_col[lab_ind])
        else:
            bar_x_seq.append(bar_x[lab_ind])
            bar_lab_seq.append(bar_labels[lab_ind])
            bar_col_seq.append(bar_col[lab_ind])

    # Create the figure
    plt.close() # ensure no figure is already open
    fig, axs = plt.subplots(
        len(exps),
        max([len(barcodes[exp]) for exp in exps]),
        squeeze=False
    ) # initialize figure
    # (squeeze set to False allows to always generate a 2D-array for the
    # subplots, thus allowing us to deal with the subplots in the same way no
    # matter if there is only one experiment or only one barcode (or both))
    set_size = fig.set_size_inches(
        (1+max(bar_x)/3)*max([len(barcodes[exp]) for exp in exps]),
        4*len(exps)
    ) # (width, height) set figure size

    # Retrieve the maximum height of all the barcodes seq/base statistics
    max_seq_amount, max_base_amount = 0, 0
    for exp in exps:
        for bc in barcodes[exp]:
            for c_type in clean_stats[exp][bc].keys():
                if c_type in c_types:
                    for label in clean_stats[exp][bc][c_type].keys():
                        if 'bases' in label \
                        and clean_stats[exp][bc][c_type][label] > max_base_amount:
                            max_base_amount = clean_stats[exp][bc][c_type][label]
                        elif not 'bases' in label and label != 'input_seq' \
                        and clean_stats[exp][bc][c_type][label] > max_seq_amount:
                            max_seq_amount = clean_stats[exp][bc][c_type][label]

    # Plot for each experiment
    for exp_ind, exp in enumerate(exps):

        # Plot for each barcode:
        for bc_ind, bc_name in enumerate(barcodes[exp]):

            # Retrieve the statistics for this barcode
            bar_y = []
            for c_type in c_types:
                if c_type in clean_stats[exp][bc_name].keys():
                    for label in C_LABELS[c_type]:
                        # Test if the splitting parameter exists in the log
                        # file - else, assume it is 0 (this happens e.g. in
                        # porechop if no adapter is found, then splitting
                        # is not performed and thus the splitting statistics
                        # cannot be parsed from the log file)
                        if label in clean_stats[exp][bc_name][c_type].keys():
                            bar_y.append(
                                clean_stats[exp][bc_name][c_type][label]
                            )
                        else:
                            bar_y.append(0)
                else:
                    for label in C_LABELS[c_type]:
                        bar_y.append(0)
            
            # Split the bar heights between sequence and base scales
            bar_y_seq, bar_y_base = [], []
            for lab_ind, label in enumerate(bar_labels):
                if 'bases' in label:
                    bar_y_base.append(bar_y[lab_ind])
                else:
                    bar_y_seq.append(bar_y[lab_ind])
            
            # Plot the barplot for this barcode (read scale)
            plotting_seq_barplot = axs[exp_ind, bc_ind].bar(
                bar_x_seq,
                bar_y_seq,
                bar_width,
                color=bar_col_seq,
                label=bar_lab_seq
            )
            axs[exp_ind, bc_ind].set_ylim(top=max_seq_amount) # set the right
            # y limit

            # Set the ticks and labels along the x axis
            axs[exp_ind, bc_ind].set_xticks(
                bar_x,
                labels=bar_labels,
                rotation=60,
                ha='right',
                fontdict={'fontsize': FONTSIZE}
            ) # position labels along x axis

            # Add ax titles
            if bc_ind == 0:
                axs[exp_ind, bc_ind].set_ylabel(exp)
            if sample_names is not None:
                x_label = sample_names[exp][bc_name]
            else:
                x_label = bc_name
            axs[exp_ind, bc_ind].set_xlabel(x_label.replace('\\n', '\n'))

            # Instantiate another ax for the base scale
            ax_base = axs[exp_ind, bc_ind].twinx()
            plotting_bases_barplot = ax_base.bar(
                bar_x_base,
                bar_y_base,
                bar_width,
                color=bar_col_base,
                label=bar_lab_base
            )
            ax_base.set_ylim(top=max_base_amount) # set the right y limit

            # Add some boxes with the cleaning steps above the bars
            for c_i, c_type in enumerate(c_types):
                # Retrieve the indices of that cleaning type
                c_bars_all_i = c_indices[c_i]
                # Retrieve the bars corresponding to that cleaning type
                c_seq_bars_i = [i for i in range(len(bar_x_seq)) if bar_x_seq[i] in c_bars_all_i]
                c_base_bars_i = [i for i in range(len(bar_x_base)) if bar_x_base[i] in c_bars_all_i]
                c_bars = [
                    plotting_seq_barplot[i] for i in c_seq_bars_i
                ] + [
                    plotting_bases_barplot[i] for i in c_base_bars_i
                ]
                # Get the average x position
                c_bars_all_x = [bar.get_x() + bar_width/2 for bar in c_bars]
                c_box_x = sum(c_bars_all_x) / len(c_bars_all_x)
                # Get the maximum height of the bars
                # ...and convert them in data coordinate of the second=base ax
                c_heights = [
                    plotting_seq_barplot[i].get_height() \
                    / axs[exp_ind, bc_ind].get_ylim()[1] \
                    * ax_base.get_ylim()[1] \
                    for i in c_seq_bars_i
                ] + [
                    plotting_bases_barplot[i].get_height() \
                    for i in c_base_bars_i
                ]
                # ...then add a lil something to put the box slightly over the
                # highest bar
                c_box_y = max(c_heights) + 0.05 * ax_base.get_ylim()[1]
                # Define the text to be shown
                if c_type in clean_stats[exp][bc_name].keys():
                    in_seq = clean_stats[exp][bc_name][c_type]['input_seq']
                else:
                    in_seq = 0
                # Add the text
                c_box = ax_base.text(
                    c_box_x,
                    c_box_y,
                    C_CLEAN_NAMES[c_type] + '\ninput seq: ' + str(in_seq),
                    ha="center", va="bottom",
                    fontsize=FONTSIZE, color="black",
                    bbox=dict(
                        boxstyle="square,pad=0.5",
                        edgecolor="black",
                        facecolor=PLOT_COLORS[c_type]['title_background'],
                        alpha=1
                    )
                )
            
            # Set scale to logarithmic if specified
            # (someth goes wrong with the boxes with log scale but idk why and
            # tbh idc since i am gonna use linear scale)
            if log_scale:
                axs[exp_ind, bc_ind].set_yscale('log')
                ax_base.set_yscale('log') # set log scale on y axis
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=500)

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', type=str, required=True,
                        help='Path to the input directory. This must contain '
                        + 'the cleaning statistics of the samples.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output png file.')
    parser.add_argument('-s', '--samples', type=str,
                        help='Path to the directory containing text files with'
                        + ' the sample names corresponding to each barcode '
                        + 'used at barcoding step, each in the form: '
                        + '<barcode_nb>=<sample_name> where <barcode_nb> '
                        + 'should be an integer. Text files are searched for '
                        + 'recursively. They should be named after their '
                        + 'experiment name as follow: <exp_name>.barcodes')
    parser.add_argument('-c', '--c_types', type=str, nargs='+',
                        help='Coma-separated list of cleaning types to use - '
                        + 'this must belong to the following: seq_adapt, '
                        + 'twist_outer, twist_primers, ont_primer_tails, '
                        + 'ont_unknown_seq, tso, polyA/T_tails, GA_repeats')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Check input directory
    utils.check_dir(args.input_dir)

    # Retrieve the reads cleaning statistics
    text = 'Looking for filtering/trimming/cleaning statistics in ' \
           + args.input_dir
    utils.send_text(text, v, 1, 0)
    clean_stats = get_cleaning_stats(args.input_dir)

    utils.send_text('Found cleaning stats are:', v, 2, 0)
    if v >= 2: utils.print_json(clean_stats)

    if args.samples is not None:
        # Retrieve the barcodes names
        utils.send_text('Retrieving the name chosen for each barcode', v, 1, 0)
        sample_names = utils.get_sample_names(args.samples)
        utils.send_text('Found sample names are:', v, 2, 0)
        if v >= 2: utils.print_json(sample_names)
    else:
        sample_names = None

    # Plot the statistics
    utils.send_text(
        'Plotting the trimming statistics (this can take some time...)',
        v, 1, 0
    )
    plotting_res = plot_cleaning_stats(
        clean_stats, args.output, sample_names, args.c_types
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
