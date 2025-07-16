#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot composition of mixture.
Takes as input the mapping statistics of several samples agaisnt the different
targets in the form given by the `samtools stats` command, with the option `-t
chromosomes_sizes_file.txt`, before running `cat ... | grep ^RL | cut -f 2- >
...`, and plot the samples composition in the form of a bar plot.
"""


import os
import sys

import argparse
import json
import matplotlib.patches as pcs
import matplotlib.pyplot as plt
import numpy as np

import utils # import local utility functions


MY_COLORS = [
    'lime', 'mediumseagreen', 'aquamarine', 'cyan', 'deepskyblue', 'steelblue'
]


def get_samples_composition(samstats_summary_dir):
    """Get the composition of every samples after mapping and samtools stats.

    Positionnal arguments:
    samstats_summary_dir (str) - path to directory with output of the `cat ... |
        grep ^RL | cut -f 2- > ...` command applied to samtools stats output.
        The file names must comply to the following naming convention:
        '<target>/<sample>.<target>.summary'.

    Return:
    samples_dict (dict) - dictionnary containing the samples composition in
        targets.
    """

    # Initialize dictionnary with samples
    samples_dict = {}

    # Loop through files and parse the names to find cross-species alignemnt stats
    for file_path in utils.find_file(samstats_summary_dir, '_snakeref.samstats'):
        if file_path.endswith('_snakeref.samstats'):

            # Parse the name to retrieve sample
            sample = file_path.split('/')[-1].split('_map2_')[0]

            # Open file and retrieve the number of total reads
            with open(samstats_summary_dir.rstrip('/')+'/'+file_path, 'r') as file:
                line = file.readline()
                while line:
                    if line.startswith('SN\traw total sequences:'):
                        total_reads = int(line.strip().split('\t')[2])
                        break
                    line = file.readline()

            # Update the samples dictionnary accordingly
            samples_dict[sample] = {'all': total_reads}

    # Loop through files and parse the names to find samples and targets
    for file_path in utils.find_file(samstats_summary_dir, '.summary'):
        if file_path.endswith('.summary'):

            # Parse the name to retrieve sample and target
            sample, target, dummy = file_path.split('/')[-1].split('.')

            # Open file and retrieve the number of mapped reads
            with open(samstats_summary_dir.rstrip('/')+'/'+target+'/'+file_path, 'r') as file:
                mapped_reads = int([
                    line.strip().split('\t')[1] for line in file.readlines() if line.startswith('reads mapped:')
                ][0])

            # Update the samples dictionnary accordingly
            if not sample in samples_dict.keys():
                samples_dict[sample] = {target: mapped_reads}
            else:
                samples_dict[sample][target] = mapped_reads

    return samples_dict


def get_list_of_samples_and_species(samples_dict):
    """Returns the sorted list of all samples and all species.

    Positionnal arguments:
    samples_dict (dict) - dictionnary containing the samples composition in
        targets

    Return:
    samples (list) - sorted list of all samples
    species (list) - sorted list of all species
    """

    samples = [sample for sample in samples_dict.keys()]
    samples.sort()
    species = []
    for sample in samples:
        species += [specie for specie in samples_dict[sample].keys()]
    species = list(set(species))
    species.sort()

    return samples, species


def get_species_repartition(samples, species, samples_dict):
    """Returns a dict containing the number of reads per sample for each target.

    Positionnal arguments:
    samples      (list) - sorted list of all samples
    species      (list) - sorted list of all species
    samples_dict (dict) - dictionnary containing the samples composition in
        targets

    Return:
    composition (dict) - dictionnary containing the number of reads per sample
        for each target
    """

    composition = {specie: [] for specie in species}
    for specie in species:
        for sample in samples:
            if specie in samples_dict[sample].keys():
                composition[specie].append(samples_dict[sample][specie])
            else:
                composition[specie].append(0)

    return composition


def plot_composition(samples, species, composition, bar_width, bar_colors,
    legend_elements, title, output_file, total=False, log_scale=False,
    y_label='Read count', verbosity=0):
    """Returns a dict containing the number of reads per sample for each target.

    Positionnal arguments:
    samples          (list) - sorted list of all samples
    species          (list) - sorted list of all species
    composition      (dict) - dictionnary containing the number of reads per
        sample for each target
    bar_width       (float) - width of bars in the figure
    bar_colors       (dict) - dictionnary containing one color per specie
    legend_elements  (list) - list of the patch of colors and the name of each
        specie
    title             (str) - plot title
    output_file       (str) - path to output file
    total            (bool) - plot total reads as well
    log_scale        (bool) - use True for logarithmic scale on y axis
    y_label           (str) - y axis name
    verbosity         (int) - level of verbosity
    """

    plt.close() # ensure no figure is already open
    fig, ax = plt.subplots() # initialize figure
    dummy_var = fig.set_size_inches(10, 5) # (width, height) set figure size

    if total:
        # plot total number of reads per sample
        dummy_var = ax.bar(
            np.arange(len(samples)) * ((len(species) + 2.5) * bar_width) + (len(species[:-1]) * bar_width / 2),
            composition['all'],
            bar_width * (len(species) + 1),
            color='gainsboro'
        )

    # plot number of reads per species
    multiplier = 0
    for specie in species:
        offset = multiplier * bar_width
        dummy_var = ax.bar(
            np.arange(len(samples)) * ((len(species) + 2.5) * bar_width) + offset,
            composition[specie],
            bar_width,
            color=bar_colors[specie]
        )
        multiplier += 1

    # ax.set_ylim([1, 1000000])
    if log_scale:
        ax.set_yscale('log') # set log scale on y axis
    ax.set_ylabel(y_label)
    # ax.set_xlabel('Read size')
    ax.set_xticks(
        np.arange(len(samples)) * ((len(species) + 2.5) * bar_width) + (((len(species) - 1) * bar_width) / 2),
        labels=samples,
        rotation=60,
        ha='right',
        fontdict={'fontsize': 7}
    ) # position labels along x axis
    ax.set_title(title)
    ax.legend(handles=legend_elements, framealpha=0.6)
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
                        help='Path to the input directory. /!\\ Must contain the'
                        + 'mapping statistics summary of all the samples '
                        + 'against all the targets. The files must respect the '
                        + 'following naming convention: '
                        + "'<target>/<sample>.<target>.summary'.")
    parser.add_argument('-o', '--output', type=str,
                        help='Path to the output directory (and prefix) '
                        '(default: same directory)')
    parser.add_argument('-s', '--samples', type=str,
                        help='Coma-separated list of barcode/samples used at '
                        + 'barcoding step, each in the form: '
                        + '<barcode_nb>=<sample_name> where <barcode_nb> '
                        + 'should be an integer')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Retrieve the samples and targets
    samples_dict = get_samples_composition(args.input_dir)

    # Retrieve sorted list of all samples and all species
    samples, species = get_list_of_samples_and_species(samples_dict)

    # Use only the specified samples
    if args.samples is not None and len(args.samples) > 0:
        samples_names = {}
        try:
            for bc_samp_str in args.samples.split(','):
                bc_samp = bc_samp_str.split('=')
                samples_names[
                    'barcode' + ('0' + bc_samp[0])[-2:]
                ] = bc_samp[1]
        except ValueError:
            print('String, coma-separated list of barcode/sample couples, '
                    + 'each in the form: <barcode_nb>=<sample_name> where '
                    + '<barcode_nb> should be an integer')
        samples = [
            sample for sample in samples if sample in samples_names.keys()
        ]

    # Prepare data for plotting
    composition = get_species_repartition(samples, species, samples_dict)
    norm_comp = {}
    for sp in composition.keys():
        if max(composition[sp]) > 0:
            norm_comp[sp] = [
                float(count)/sum(composition[sp]) for count in composition[sp]
            ]
        else:
            norm_comp[sp] = composition[sp]
    
    # Remove the 'all' artificial species
    species.remove('all')

    # Compute correlation matrix
    # species_count = np.array([composition[specie][:-1] for specie in species])
    # print(species_count)
    # print(species)
    # np.set_printoptions(precision=2)
    # print(np.corrcoef(species_count))

    # fig, ax = plt.subplots()
    # ax.scatter(composition['mouse'], composition['human'])
    # for i, sample in enumerate(samples):
    #     ax.annotate(sample, (composition['mouse'][i], composition['human'][i]))
    # plt.show()

    # Determine output file location
    output_path = args.input_dir.rstrip('/') + '/'
    if args.output is not None:
        output_path = args.output

    # Build figures
    bar_width = 0.8 # set width or bars
    bar_colors = {
        species[sp_ind]: MY_COLORS[sp_ind] for sp_ind in range(len(species))
    } # set colors of species
    legend_elements = [
        pcs.Patch(facecolor=bar_colors[sp], label=sp) for sp in species
    ] # set patch of colors to be used in the figure legend

    if args.samples is not None and len(args.samples) > 0:
        # Use samples names as labels for x axis
        samples = [samples_names[sample] for sample in samples]
        # Change plot file name accordingly
        output_path += '_samples'
    else:
        output_path += '_barcodes'

    # Plot samples composition in terms of targets
    title = 'Outcome of mapping against ' + str(len(species)) + ' species'
    dummy_var = plot_composition(
        samples,
        species,
        composition,
        bar_width,
        bar_colors,
        legend_elements,
        title,
        output_file=output_path + '_composition.png',
        log_scale=True,
        verbosity=v
    )

    # Plot normalized samples composition in terms of targets
    title = 'Outcome of mapping against ' + str(len(species)) + ' species - $\\it{normalized}$'
    dummy_var = plot_composition(
        samples,
        species,
        norm_comp,
        bar_width,
        bar_colors,
        legend_elements,
        title,
        output_file = output_path + '_normalized_composition.png',
        y_label='Proportion of species reads',
        verbosity = v
    )

    legend_elements = [
        pcs.Patch(facecolor='gainsboro', label='Total')
    ] + legend_elements # add the total reads to the figure legend

    # Plot samples composition in terms of targets (with total reads)
    title = 'Outcome of mapping against ' + str(len(species)) + ' species'
    dummy_var = plot_composition(
        samples,
        species,
        composition,
        bar_width,
        bar_colors,
        legend_elements,
        title,
        output_file=output_path + '_composition_with_total.png',
        total=True,
        log_scale=True,
        verbosity=v
    )

    # Plot normalized samples composition in terms of targets (with total reads)
    title = 'Outcome of mapping against ' + str(len(species)) + ' species - $\\it{normalized}$'
    dummy_var = plot_composition(
        samples,
        species,
        norm_comp,
        bar_width,
        bar_colors,
        legend_elements,
        title,
        output_file = output_path + '_normalized_composition_with_total.png',
        total=True,
        y_label='Proportion of species reads',
        verbosity = v
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
