#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot the heatmap of read qualities versus read lengths for each target.
Takes as input a '.pkl' file with all the information in the form of a
pandas data frame, i.e. reads lengths, the reads nucleotides frequencies, the
reads species it was assigned to, and the reads bases quality (Phred score),
and plots one heatmap for each barcode (lines) and target / all / unclassified
(columns).
"""


import os
import sys

import argparse
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utils # import local utility functions


# Define number of bins for quality axis
QBIN = 250
# Define number of bins for size axis
SBIN = 250
# Define the label to plot for each target
C_CLEAN_NAMES = {
    'All': 'All',
    'Unmapped': 'Unmapped',
    'human': 'Human',
    'mouse': 'Mouse',
    'SARSCoV2': 'SARS-CoV-2'
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pickle_dir', type=str, required=True,
                        help='Path to the directory containing the input '
                        + 'pickle files. The path should be at least of '
                        + 'length 2, with the parent of the pickle directory '
                        + 'being named after the experiment. Pickle files '
                        + 'should be named as follow: <sample>_6_GA.pkl')
    parser.add_argument('-o', '--output', type=str,
                        help="Path to the output '.png' file (default: same "
                        + "directory as pickle files, name 'size_quality_heat"
                        + "map.png').")
    parser.add_argument('-s', '--samples', type=str,
                        help='Path to the directory containing text files with'
                        + ' the sample names corresponding to each barcode '
                        + 'used at barcoding step, each in the form: '
                        + '<barcode_nb>=<sample_name> where <barcode_nb> '
                        + 'should be an integer. Text files are searched for '
                        + 'recursively. They should be named after their '
                        + 'experiment name as follow: <exp_name>.barcodes')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Determine output directory
    if args.output is None:
        args.output = args.pickle_dir.rstrip('/') + '/size_quality_heatmap.png'

    # Print script starting informations
    command = 'plot_sequence_stats_heatmap.py' \
              + f' -p {args.pickle_dir} -o {args.output} -s {args.samples}' \
              + f' -v {v}'
    utils.send_text('plot_sequence_stats_heatmap.py: Running with following '
                    + 'command:', v, 3, 0)
    utils.send_text(command, v, 3, 0)

    # Retrieve the barcodes (cutting the size of '_6_GA.pkl' = 9
    # letters from the end of the path)
    utils.send_text('plot_sequence_stats_heatmap.py: Getting samples', v, 1, 1)
    bc_unt_pkl_paths = utils.find_path(args.pickle_dir, '_6_GA.pkl')
    barcodes = [r_path.split('/')[-1][:-9] for r_path in bc_unt_pkl_paths]
    barcodes.sort()
    utils.send_text('plot_sequence_stats_heatmap.py: Found barcodes are: ' \
                    +  ', '.join(barcodes), v, 2, 2)

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
    

    # Loop through barcodes and load data frames to store extreme values and
    # found species
    utils.send_text('plot_sequence_stats_heatmap.py: Computing extreme x/y',
                    v, 1, 1)

    min_len = 10000
    max_len = 0
    min_qual = 200
    max_qual = 0
    species = []

    for barcode in barcodes:
            
        # Load data frame
        stored_df_path = f"{args.pickle_dir.rstrip('/')}/{barcode}_6_GA.pkl"
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
    
    utils.send_text(
        f'Lengths: [{min_len}; {max_len}] - Quality: [{min_qual}; {max_qual}]',
        v, 3, 3
    )

    # Place 'Trimmed' and 'Unmapped' at the front of the species list, because
    # we want to plot them in a given way later
    species.sort()
    if 'Unmapped' in species:
        species.insert(0, species.pop(species.index('Unmapped')))
    species.insert(0, 'All')

    # Loop through barcodes, load data frames and compute heatmaps data
    hm_data = {
        bc: {
            target: [] for target in species
        } for bc in barcodes
    }
    max_bin_val = 0
    q_bins = np.linspace(min_qual, max_qual, num=QBIN+1)
    s_bins = np.exp(
        np.linspace(
            np.log(min_len),
            np.log(max_len),
            num=SBIN+1
        )
    )
    # q_mesh, s_mesh = np.meshgrid(q_bins, s_bins)
    # q_bins = np.append(q_bins, float('+inf'))
    # s_bins = np.append(s_bins, float('+inf'))
    for bc_ind, bc in enumerate(barcodes):

        # Load data frame
        stored_df_path = f"{args.pickle_dir.rstrip('/')}/{bc}_6_GA.pkl"
        qual_df = pd.read_pickle(stored_df_path)

        utils.send_text(
            f'plot_sequence_stats_heatmap.py: Prepare data for plotting'
            + f' {bc} - All targets',
            v, 2, 2
        )

        # Build heatmap images
        hm_data[bc]['All'], q_edges_out, s_edges_out = np.histogram2d(
            qual_df['qual'].to_numpy(),
            qual_df['size'].to_numpy(),
            bins=[q_bins, s_bins]
        )

        # Mask 0 data
        hm_data[bc]['All'] = np.ma.masked_where(
            hm_data[bc]['All'] == 0, 
            hm_data[bc]['All']
        )

        for target_ind, target in enumerate(species[1:]):

            utils.send_text(
                f'plot_sequence_stats_heatmap.py: Prepare data for plotting'
                + f' {bc} - {target}',
                v, 2, 2
            )

            # Build heatmap images
            hm_data[bc][target], q_edges_out, s_edges_out = np.histogram2d(
                qual_df['qual'][qual_df['species'] == target].to_numpy(),
                qual_df['size'][qual_df['species'] == target].to_numpy(),
                bins=[q_bins, s_bins]
            )

            # Mask 0 data
            hm_data[bc][target] = np.ma.masked_where(
                hm_data[bc][target] == 0, 
                hm_data[bc][target]
            )

    # Initialize plot
    utils.send_text(
        'plot_sequence_stats_heatmap.py: Starting to plot...', v, 1, 1
    )
    plt.clf() # Clear existing figure
    nb_barcodes = len(barcodes)
    nb_targets = len(species)
    fig, axs = plt.subplots(
        nb_barcodes,
        nb_targets,
        figsize=(20, 15),
        constrained_layout=True
    ) # sharex=True, sharey=True
    axs = axs.ravel() # flattens the multidimensional array of axs

    # Modify color map to display 0-counts bins as white
    cmap = plt.get_cmap('RdYlGn').copy()
    cmap.set_bad(color='white')

    # Compute a matlplotlib normalization to get common colormap legend
    all_hm = np.stack(
        [hm_data[bc][target] for (bc, target) in list(
            itertools.product(barcodes, species)
        )]
    ) # all heatmap data arrays
    norm = plt.Normalize(vmin=0, vmax=all_hm.max())

    # Determine ax ids where to plot x-axis labels
    # ...because if the plot has no data, it is almost impossible to plot the
    # x-ticks nicely in a log-way - and that can happen for the last plot of
    # any column!
    # Add labels on last non-empty subplot for each column
    x_lab_inds = []
    for sp_ind, specie in enumerate(species):
        x_lab_inds.append([
            bc_i * nb_targets + sp_ind for bc_i, bc in enumerate(barcodes) \
            if hm_data[bc][specie].sum() > 0
        ][-1])
    
    # Loop through barcodes, load data frames and plot heatmaps
    for bc_ind, bc in enumerate(barcodes):

        for target_ind, target in enumerate(species):

            utils.send_text(
                f'plot_sequence_stats_heatmap.py: Plotting {bc} - {target}',
                v, 2, 2
            )

            # Plot the quality versus read length
            ax_id = bc_ind * nb_targets + target_ind
            if hm_data[bc][target].sum() > 0:
                axs[ax_id].grid(True, axis='both', alpha=0.5)

            mesh = axs[ax_id].pcolormesh(
                s_bins,
                q_bins,
                hm_data[bc][target],
                cmap=cmap,
                norm=norm,
                shading='auto'
            )

            # Configure plot
            # axs[ax_id].set_xlim(min_len, max_len)
            # axs[ax_id].set_ylim(min_qual, max_qual)
            # Remove x tick labels on all but last row
            # if bc_ind != (nb_barcodes-1):
            if not (ax_id in x_lab_inds):
                axs[ax_id].tick_params(labelbottom=False)
            # Remove y tick labels on all but first column
            if target_ind != 0:
                axs[ax_id].tick_params(labelleft=False)
            # Add x labels on the last row
            if bc_ind == (nb_barcodes-1):
                if target in C_CLEAN_NAMES.keys():
                    clean_label = C_CLEAN_NAMES[target]
                else :
                    clean_label = target
                axs[ax_id].set_xlabel(clean_label)
            # Add y labels on the first column
            if target_ind == 0:
                if sample_names:
                    axs[ax_id].set_ylabel(
                        sample_names[bc].replace('\\n', '\n'),
                        rotation=60,
                        labelpad=20
                    )
                else:
                    axs[ax_id].set_ylabel(bc)
            # Set log scales if data not empty
            if hm_data[bc][target].sum() > 0:
                axs[ax_id].set_xscale('log')
                #axs[ax_id].set_yscale('log')
            else:
                # Remove x ticks
                axs[ax_id].tick_params(bottom=False, labelbottom=False)

    # Set labels, legend and title
    fig.supxlabel('Read length')
    fig.supylabel('Read quality (Phred)')
    fig.colorbar(
        mesh,
        label='Read count',
        ax=axs.tolist(),
        location='right',
        fraction = 0.02,
        pad = 0.01,
        aspect=6
    )
    fig.suptitle(
        exp_name + ' - Read size vs quality after trimming'
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
