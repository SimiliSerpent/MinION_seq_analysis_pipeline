#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot the evolution of read qualities versus read lengths.
Takes as input a '.pkl' file with all the information in the form of a
pandas data frame, i.e. reads lengths, the reads nucleotides frequencies, the
reads species it was assigned to, and the reads bases quality (Phred score),
and plots one scatter plot for each barcode (lines) and trimming condition
(columns).
"""


import os
import sys

import argparse
import matplotlib.pyplot as plt
import pandas as pd

import utils # import local utility functions


MY_MARKERS = ['.', 'o', '^', 's', 'P', '*', 'X', 'd', '2']
MY_COLORS = [
    'cornflowerblue', 'cornflowerblue', 'orangered', 'orange', 'yellow',
    'springgreen', 'aqua', 'magenta', 'indigo'
]

# Define the label to plot for each cleaning step
C_CLEAN_NAMES = {
    '1_seq_adapt': 'Sequencing\nadapters',
    '2_twist_outer': 'Twist outer\nuniv. sequences',
    '3_twist_primers': 'Twist\nprimers',
    '2_ont_primer_tails': 'ONT primer\ntails',
    '3_ont_unknown_seq': 'ONT unknown\nsequences',
    '4_tso': 'Takara\nTSO',
    'ONT_barcode': 'ONT barcodes',
    'ONT_sequences': 'ONT sequences',
    '5_polyA': 'polyA/T\ntails',
    '6_GA': 'GA repeats'
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pickle_dir', type=str, required=True,
                        help='Path to the directory containing the input '
                        + 'pickle files. The path should be at least of '
                        + 'length 2, with the parent of the pickle directory '
                        + 'being named after the experiment')
    parser.add_argument('-o', '--output', type=str,
                        help="Path to the output '.png' file (default: same "
                        + "directory as pickle files, name 'size_quality_evolu"
                        + "tion.png').")
    parser.add_argument('-s', '--samples', type=str,
                        help='Path to the directory containing text files with'
                        + ' the sample names corresponding to each barcode '
                        + 'used at barcoding step, each in the form: '
                        + '<barcode_nb>=<sample_name> where <barcode_nb> '
                        + 'should be an integer. Text files are searched for '
                        + 'recursively. They should be named after their '
                        + 'experiment name as follow: <exp_name>.barcodes')
    parser.add_argument('-n', '--sampling_flat', type=int,
                        help='Plot using a subsample of size n of each '
                        + 'barcode.')
    parser.add_argument('-f', '--sampling_frac', type=float,
                        help="Plot using a fraction f of each barcode.")
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Determine output directory
    if args.output is None:
        args.output = args.pickle_dir.rstrip('/') + '/size_quality_evolution.png'

    # Print script starting informations
    command = 'plot_sequence_stats.py' \
              + ' -p ' + args.pickle_dir \
              + ' -o ' + args.output \
              + ' -v ' + str(v)
    utils.send_text('plot_sequence_stats.py: Running with following '
                    + 'command:', v, 3, 0)
    utils.send_text(command, v, 3, 0)

    # Retrieve the barcodes (cutting the size of '_untrimmed.pkl' = 14
    # letters from the end of the path)
    utils.send_text('compute_sequence_stats.py: Getting barcodes', v, 1, 1)
    bc_unt_pkl_paths = utils.find_path(args.pickle_dir, '_untrimmed.pkl')
    barcodes = [r_path.split('/')[-1][:-14] for r_path in bc_unt_pkl_paths]
    barcodes.sort()
    utils.send_text('compute_sequence_stats.py: Found barcodes are: ' \
                    +  ', '.join(barcodes), v, 2, 2)

    # Retrieve the different trimming states
    utils.send_text(
        'compute_sequence_stats.py: Getting trimming states', v, 1, 1
    )
    all_trim_states = utils.find_file(
        args.pickle_dir, barcodes[0] + '*.pkl'
    )
    # Remove the name of the barcode from the beginning of the file name (+ 1
    # for the '_')...
    # ...and the '.pkl' suffix (size=4) from the end of the file name
    all_trim_states = [
        name[len(barcodes[0])+1 : -4] for name in all_trim_states
    ]
    all_trim_states.remove('untrimmed')
    all_trim_states.sort() # sort...
    all_trim_states = ['untrimmed'] + all_trim_states # ...& add untrimmed 1st
    utils.send_text('compute_sequence_stats.py: Found trimming states are: ' \
                    +  ', '.join(all_trim_states), v, 2, 2)

    # Retrieve the name of the experiment
    exp_name = args.pickle_dir.rstrip('/').split('/')[-3]
    
    if args.samples is not None:
        # Retrieve the barcodes names
        utils.send_text('Retrieving the name chosen for each barcode', v, 1, 0)
        sample_names = utils.get_sample_names(args.samples)[exp_name]
        utils.send_text('Found sample names are:', v, 2, 0)
        if v >= 2: utils.print_json(sample_names)
        # Filter the barcodes depending on whether they are named
        barcodes = [bc for bc in barcodes if bc in sample_names.keys()]
    else:
        sample_names = None
    

    # Loop through barcodes and trimming states and load data frames to store
    # extreme values and found species
    utils.send_text('compute_sequence_stats.py: Retrieving extreme x & y values',
                    v, 1, 1)
    min_len = 10000
    max_len = 0
    min_qual = 200
    max_qual = 0
    species = []
    for barcode in barcodes:
        for trim_state in all_trim_states:
            
            # Load data frame
            stored_df_path = args.pickle_dir.rstrip('/') +  '/' + barcode + '_' \
                             + trim_state + '.pkl'
            qual_df = pd.read_pickle(stored_df_path)

            # Update extrema
            if len(qual_df['size']) > 0:
                min_len = min(min_len, min(qual_df['size']))
                max_len = max(max_len, max(qual_df['size']))
                min_qual = min(min_qual, min(qual_df['qual']))
                max_qual = max(max_qual, max(qual_df['qual']))
            # Update found species
            species += list(
                qual_df['species'][~qual_df['species'].isin(species)].unique()
            )
    # Place 'Trimmed' and 'Unmapped' at the front of the species list, because
    # we want to plot them in a given way later
    species.sort()
    if 'Unmapped' in species:
        species.insert(0, species.pop(species.index('Unmapped')))
    if 'Trimmed' in species:
        species.insert(0, species.pop(species.index('Trimmed')))

    # Initialize plot
    utils.send_text('compute_sequence_stats.py: Starting to plot...', v, 1, 1)
    plt.clf()
    nb_barcodes = len(barcodes)
    nb_trim_states = len(all_trim_states)
    fig, axs = plt.subplots(
        nb_barcodes,
        nb_trim_states,
        figsize=(20, 15),
        constrained_layout=True
    ) # sharex=True, sharey=True
    axs = axs.ravel() # flattens the multidimensional array of axs
    
    # Loop through barcodes and trimming states, load data frames and plot
    # scatter plots
    for bc_ind, bc in enumerate(barcodes):

        for trim_state_ind, trim_state in enumerate(all_trim_states):

            utils.send_text('compute_sequence_stats.py: Plotting '
                            + bc + ' - ' + trim_state,
                            v, 2, 2)

            # Load data frame
            stored_df_path = args.pickle_dir.rstrip('/') + '/' + bc + '_' \
                             + trim_state + '.pkl'
            qual_df = pd.read_pickle(stored_df_path)

            # Sub-sample if requested
            if args.sampling_flat is not None:
                tot_reads = len(qual_df.index)
                sampling_size = min(tot_reads, args.sampling_flat)
                inds = list(qual_df.sample(sampling_size).index)
                qual_df = qual_df[qual_df.index.isin(inds)]
                if args.sampling_frac is not None:
                    print('Both flat and proportional subsampling were used.' \
                        + ' Flat subsampling has been priorized')
            elif args.sampling_frac is not None:
                sampling_size = int(args.sampling_frac * qual_df.size)
                inds = list(qual_df.sample(sampling_size).index)
                qual_df = qual_df[qual_df.index.isin(inds)]

            # Plot the quality versus read length
            ax_id = bc_ind * nb_trim_states + trim_state_ind
            axs[ax_id].grid(True, alpha=0.5)

            # Plot reads belonging to each species with different colors
            for sp_ind, specie in enumerate(species):
                utils.send_text('compute_sequence_stats.py: Plotting '
                                + bc + ' - ' + trim_state
                                + ' (' + specie + ' reads)',
                                v, 3, 3)
                dummy = axs[ax_id].scatter(
                    x=qual_df['size'][qual_df['species'] == specie].values,
                    y=qual_df['qual'][qual_df['species'] == specie].values,
                    s=2,
                    c=MY_COLORS[sp_ind],
                    marker=MY_MARKERS[sp_ind]
                )

            # Configure plot
            axs[ax_id].set_xlim(min_len, max_len)
            axs[ax_id].set_ylim(min_qual, max_qual)
            # Remove x ticks on all but last row
            if bc_ind != (nb_barcodes-1):
                axs[ax_id].set_xticks([])
                axs[ax_id].tick_params(bottom=False, labelbottom=False)
            # Remove y ticks on all but first column
            if trim_state_ind != 0:
                axs[ax_id].set_yticks([])
                axs[ax_id].tick_params(left=False, labelleft=False)
            # Add x labels on the last row
            if bc_ind == (nb_barcodes-1):
                if all_trim_states[trim_state_ind] in C_CLEAN_NAMES.keys():
                    clean_label = C_CLEAN_NAMES[all_trim_states[trim_state_ind]]
                else :
                    clean_label = all_trim_states[trim_state_ind]
                axs[ax_id].set_xlabel(clean_label)
            # Add y labels on the first column
            if trim_state_ind == 0:
                if sample_names:
                    axs[ax_id].set_ylabel(
                        sample_names[barcodes[bc_ind]].replace('\\n', '\n')
                    )
                else:
                    axs[ax_id].set_ylabel(barcodes[bc_ind])
            # Set log scales if data not empty
            if len(qual_df['size']) > 0:
                axs[ax_id].set_xscale('log')
                axs[ax_id].set_yscale('log')

    # Set labels, legend and title
    fig.supxlabel('Read length')
    fig.supylabel('Read quality (Phred)')
    fig.legend(species, loc='outside right')
    # (Explanation:
    #    Here we are using a discouraged matplotlib functionnality of the
    # figure 'legend' method. By passing a list as argument, matplotlib
    # associates one label with one artist, artists basically being the scatter
    # plots defined above. Here, we have plotted nb_barcodes x nb_trim_states x
    # nb_species scatter plots, and the list we pass to the legend method
    # is of size nb_species. The list of artists will be truncated, and our
    # labels will be associated with the nb_species first artists.
    # Theoretically, it should map fine. However, the ordering of the artists
    # and of our labels list is the only reason why it works! Hence why it is
    # discouraged by matplotlib.
    # )
    fig.suptitle(
        exp_name + ' - Evolution of read size & quality during trimming'
    )

    # Force matplotlib to draw the figure, which will force it into calculating
    # the layout, thus applying the `contraint_layout` options. Since this opt-
    # ion re-defines the ticks, we can then remove the undesired ticks, before
    # saving.
    # fig.canvas.draw()

    # Saving figure
    utils.send_text('compute_sequence_stats.py: Saving fig', v, 3, 2)
    # plt.tight_layout()
    plt.savefig(args.output, dpi=500, bbox_inches="tight")

    return 0


if __name__ == "__main__":
    sys.exit(main())
