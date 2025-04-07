#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot the evolution of read qualities versus read lengths.
Takes as input a '.pkl' file with all the information in the form of a
pandas data frame, i.e. reads lengths, the reads nucleotides frequencies, the
reads species it was assigned to, and the reads bases quality (Phred score),
and plots one scatter plot for each sample (lines) and trimming condition
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pickle_dir', type=str, required=True,
                        help='Path to the directory containing the input pickle'
                        + 'files.')
    parser.add_argument('-o', '--output', type=str,
                        help="Path to the output '.png' file (default: same "
                        + "directory as pickle files, name 'size_quality_evolu"
                        + "tion.png').")
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

    # Retrieve the samples names (cutting the size of '_untrimmed.pkl' = 14
    # letters from the end of the path)
    utils.send_text('compute_sequence_stats.py: Getting sample names', v, 1, 1)
    sample_paths = utils.find_path(args.pickle_dir, '_untrimmed.pkl')
    sample_names = [r_path.split('/')[-1][:-14] for r_path in sample_paths]
    sample_names.sort()
    utils.send_text('compute_sequence_stats.py: Found samples are: ' \
                    +  ', '.join(sample_names), v, 2, 2)

    # Retrieve the different trimming states
    utils.send_text(
        'compute_sequence_stats.py: Getting trimming states', v, 1, 1
    )
    all_trim_states = utils.find_file(
        args.pickle_dir, sample_names[0] + '*.pkl'
    )
    # Remove the name of the sample from the beginning of the file name (+ 1
    # for the '_')...
    # ...and the '.pkl' suffix (size=4) from the end of the file name
    all_trim_states = [
        name[len(sample_names[0])+1 : -4] for name in all_trim_states
    ]
    all_trim_states.remove('untrimmed')
    all_trim_states.sort() # sort...
    all_trim_states = ['untrimmed'] + all_trim_states # ...& add untrimmed 1st
    utils.send_text('compute_sequence_stats.py: Found trimming states are: ' \
                    +  ', '.join(all_trim_states), v, 2, 2)

    # Loop through samples and trimming states and load data frames to store
    # extreme values and found species
    utils.send_text('compute_sequence_stats.py: Retrieving extreme x & y values',
                    v, 1, 1)
    min_len = 10000
    max_len = 0
    min_qual = 200
    max_qual = 0
    species = []
    for sample in sample_names:
        for trim_state in all_trim_states:
            
            # Load data frame
            stored_df_path = args.pickle_dir.rstrip('/') +  '/' + sample + '_' \
                             + trim_state + '.pkl'
            qual_df = pd.read_pickle(stored_df_path)

            # Update extrema
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
    nb_samples = len(sample_names)
    nb_trim_states = len(all_trim_states)
    fig, axs = plt.subplots(
        nb_samples,
        nb_trim_states,
        figsize=(20, 15),
        constrained_layout=True
    ) # sharex=True, sharey=True
    axs = axs.ravel() # flattens the multidimensional array of axs
    
    # Loop through samples and trimming states, load data frames and plot
    # scatter plots
    for sample_ind, sample in enumerate(sample_names):

        sampled = False # TODO: remove line

        for trim_state_ind, trim_state in enumerate(all_trim_states):

            utils.send_text('compute_sequence_stats.py: Plotting '
                            + sample + ' - ' + trim_state,
                            v, 2, 2)

            # Load data frame
            stored_df_path = args.pickle_dir.rstrip('/') + '/' + sample + '_' \
                             + trim_state + '.pkl'
            qual_df = pd.read_pickle(stored_df_path)

            # if not sampled: # TODO: remove if statement
            #     inds = list(qual_df.sample(1500).index)
            #     sampled = True
            # qual_df = qual_df[qual_df.index.isin(inds)]

            # Plot the quality versus read length
            ax_id = sample_ind * nb_trim_states + trim_state_ind
            axs[ax_id].grid(True, alpha=0.5)

            # Plot reads belonging to each species with different colors
            for sp_ind, specie in enumerate(species):
                utils.send_text('compute_sequence_stats.py: Plotting '
                                + sample + ' - ' + trim_state
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
            if sample_ind != (nb_samples-1):
                axs[ax_id].set_xticks([])
                axs[ax_id].tick_params(bottom=False, labelbottom=False)
            # Remove y ticks on all but first column
            if trim_state_ind != 0:
                axs[ax_id].set_yticks([])
                axs[ax_id].tick_params(left=False, labelleft=False)
            # Add x labels on the last row
            if sample_ind == (nb_samples-1):
                axs[ax_id].set_xlabel(all_trim_states[trim_state_ind])
            # Add y labels on the first column
            if trim_state_ind == 0:
                axs[ax_id].set_ylabel(sample_names[sample_ind])
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
    # plots defined above. Here, we have pltoted nb_samples x nb_trim_states x
    # nb_species scatter plots, and the list we pass to the legend method
    # is of size nb_species. The list of artists will be truncated, and our
    # labels will be associated with the nb_species first artists.
    # Theoretically, it should map fine. However, the ordering of the artists
    # and of our labels list is the only reason why it works! Hence why it is
    # discouraged by matplotlib.
    # )
    exp_name = args.pickle_dir.rstrip('/').split('/')[-2]
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
