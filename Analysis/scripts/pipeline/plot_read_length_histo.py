#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot read lengths distribution.
Takes as input the read length distribution in the form given by the `samtools
stats` command, before running `cat ... | grep ^RL | cut -f 2- > ...`, and plot
it in the form of a histogram.
"""


import os
import sys

import argparse
import matplotlib.patches as pcs
import matplotlib.pyplot as plt
import pandas as pd

import utils # import local utility functions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--length_distrib', type=str, required=True,
                        help='Text file with tab-separated length;nb_of_reads')
    parser.add_argument('-o', '--output', type=str,
                        help='Path to the output file (default: same directory'
                        + ', with _histo.png suffix)')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Retrieve the counts in a data frame
    distrib_df = pd.read_csv(
        args.length_distrib,
        sep='\t',
        names=['size', 'count']
    )

    # Build figure
    bins = [1, 11, 21, 31, 41, 51, 101, 201, 301, 401, 501, 1001, 2001, 3001,
        4001, 5001, 6001, 7001, 8001, 9001, 10001, 15001, 20001, 50001, 100001]
    # bins = [bin for bin in bins if bin <= distrib_df['size'].max()]
    labels = [
        str(bins[i]) + '-' + str(bins[i+1]-1) for i in range(len(bins) - 1)
    ] + [str(bins[-1]) + '+']
    bar_colors = [
        'paleturquoise', 'paleturquoise', 'paleturquoise', 'paleturquoise',
        'paleturquoise', 'darkturquoise', 'lightseagreen', 'lightseagreen',
        'lightseagreen', 'lightseagreen', 'darkcyan', 'teal', 'teal', 'teal',
        'teal', 'teal', 'teal', 'teal', 'teal', 'teal', 'darkslategray',
        'darkslategray', 'xkcd:charcoal', 'black', 'black'
    ]
    legend_elements = [
        pcs.Patch(facecolor='paleturquoise', label='10'),
        pcs.Patch(facecolor='darkturquoise', label='50'),
        pcs.Patch(facecolor='lightseagreen', label='100'),
        pcs.Patch(facecolor='darkcyan', label='500'),
        pcs.Patch(facecolor='teal', label='1000'),
        pcs.Patch(facecolor='darkslategray', label='5000'),
        pcs.Patch(facecolor='xkcd:charcoal', label='30000'),
        pcs.Patch(facecolor='black', label='50000+')
    ]

    counts = [0 for bin in bins]
    size_index = 0
    for bin_index in range(len(bins) - 1):
        while size_index < len(distrib_df['size']) and distrib_df['size'][size_index] < bins[bin_index + 1]:
            counts[bin_index] += distrib_df['count'][size_index]
            size_index += 1
    if size_index < len(distrib_df['size']):
        for count in distrib_df['count'][size_index:]:
            counts[-1] += count

    # print(bins)
    # print(counts)

    plt.close()
    fig, ax = plt.subplots()
    dummy_var = ax.bar(labels, counts, color=bar_colors)

    ax.set_ylim([1, 1000000])
    ax.set_yscale('log')
    ax.set_ylabel('Read count')
    ax.set_xlabel('Read size')
    ax.set_xticks(
        range(len(labels)),
        labels=labels,
        rotation=60,
        ha='right',
        fontdict={'fontsize': 7}
    )
    ax.set_title(args.length_distrib.split('/')[-1])
    ax.legend(
        handles=legend_elements,
        title='Bin size',
        labelspacing=0.3,
        fontsize=9,
        markerscale=0.7
    )

    # Determine output file path and name
    output_path = args.length_distrib + '_histo.png'
    if args.output is not None:
        output_path = args.output

    # Save plot
    utils.send_text('Saving plot', v, 2, 1)
    plt.tight_layout()
    plt.savefig(output_path,
                dpi=500)

    return 0


if __name__ == "__main__":
    sys.exit(main())
