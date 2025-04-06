#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute reads lengths, nuc. freq., species and quality and store in file.
Takes as input a fastq file and the sam file of its alignment versus the desired
reference, and outputs a '.pkl' file with all the information in the form of a
pandas data frame, i.e. reads lengths, the reads nucleotides frequencies, the
reads species it was assigned to, and the reads bases quality (Phred score).
"""


import os
import sys

import argparse
import pandas as pd

import utils # import local utility functions


SANGER_SCORE_OFFSET = 33
INVALID_CHAR_CODE = 200


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--exp_dir', type=str, required=True,
                        help='Path to the directory containing the experiment,'
                        + " i.e. an 'alignments' directory with the alignment "
                        + "sam files; a 'raw_fastq' directory containing the "
                        + " raw fastq files; and a 'trimmed_fastq' directory "
                        + "containing the porechopped fastq files.")
    parser.add_argument('-r', '--ref_regions_dir', type=str, required=True,
                        help='Path to the directory containing the txt files '
                        + "with each species' ROI.")
    parser.add_argument('-o', '--output', type=str,
                        help="Path to the output directory with the '.pkl' "
                        + 'files with the pandas data frames for all samples '
                        + 'and trimming state inside (default: same directory '
                        + 'as experiment).')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Determine output directory
    if args.output is None:
        args.output = args.exp_dir

    # Print script starting informations
    command = 'compute_sequence_stats.py' \
              + ' -e ' + args.exp_dir \
              + ' -r ' + args.ref_regions_dir \
              + ' -o ' + args.output \
              + ' -v ' + str(v)
    utils.send_text('compute_sequence_stats.py: Running with following '
                    + 'command:', v, 3, 0)
    utils.send_text(command, v, 3, 0)

    # Retrieve the region names for all species (cutting the size of
    # '_regions.txt' = 12 letters from the end of the path)
    utils.send_text('compute_sequence_stats.py: Getting species region names',
                    v, 1, 1)
    species_paths = utils.find_path(args.ref_regions_dir, '_regions.txt')
    species_names = [r_path.split('/')[-1][:-12] for r_path in species_paths]

    # Build a dictionnary associating its species to each region
    utils.send_text('compute_sequence_stats.py: Building species region dict',
                    v, 1, 2)
    regions_dict = {}
    for species in species_names:
        fname = args.ref_regions_dir.rstrip('/') + '/' + species + '_regions.txt'
        with open(fname, 'r') as region_file:
            line = region_file.readline()
            while line:
                regions_dict[line.strip()] = species
                line = region_file.readline()

    # Retrieve the samples names (cutting the size of '.fastq' = 6 letters
    # from the end of the path)
    utils.send_text('compute_sequence_stats.py: Getting sample names', v, 1, 1)
    sample_paths = utils.find_path(
        args.exp_dir.rstrip('/') + '/raw_fastq/',
        '.fastq'
    )
    sample_names = [r_path.split('/')[-1][:-6] for r_path in sample_paths]
    sample_names.sort()
    utils.send_text('compute_sequence_stats.py: Found samples are: ' \
                    +  ', '.join(sample_names), v, 2, 2)

    # Retrieve the different trimming states
    utils.send_text(
        'compute_sequence_stats.py: Getting trimming states', v, 1, 1
    )
    all_trim_states = utils.find_file(
        args.exp_dir.rstrip('/') + '/trimmed_fastq/',
        sample_names[0] + '_porechopped_'
    )
    # Remove the name of the sample from the beginning of the file name + the
    # length of the '_porechopped_' string (size=13)...
    # ...and the '.fastq' suffix (size=6) from the end of the file name
    all_trim_states = ['untrimmed'] + [
        name[len(sample_names[0])+13 : -6] for name in all_trim_states
    ] # and add 'untrimmed' to trimming status for untrimmed sequences
    utils.send_text('compute_sequence_stats.py: Found trimming states are: ' \
                    +  ', '.join(all_trim_states), v, 2, 2)
    
    # Define a mapping allowing ASCII characters translation in Phred scores
    q_mapping = bytes(
        (
            letter - SANGER_SCORE_OFFSET
            if SANGER_SCORE_OFFSET <= letter < 94 + SANGER_SCORE_OFFSET
            else INVALID_CHAR_CODE
        )
        for letter in range(256)
    )

    # Define a list of suffix resulting of read splitting
    split_suffixes = tuple('_' + str(i) for i in range(10))

    # Loop through samples and trimming states and build data frames
    for sample in sample_names:
        utils.send_text(
            'compute_sequence_stats.py: Analyzing ' + sample, v, 1, 1
        )

        # Parse sam file to retrieve the reads alignment status
        sam = args.exp_dir.rstrip('/') \
              + '/alignments/' \
              + sample \
              + '_map2_snakeref.sam'
        utils.send_text(
            'compute_sequence_stats.py: Parsing sam ' + sam, v, 1, 2
        )

        alignment_dict = {} # build a dictionnary to store read hits

        with open(sam, 'r') as sam_file:

            # Skip sam header
            line = sam_file.readline()
            while line.startswith('@'):
                line = sam_file.readline()
            
            # Loop through all reads
            while line:

                # Parse interesting fields (we'll get the sequence and qual
                # scores in the fastq files later on)
                read_id, flag, target = line.split()[:3]
                flag = int(flag)

                # Check if read is aligned to a reference
                # (check that using bitwise operations)
                # (see: https://stackoverflow.com/questions/28590639/python-check-if-bit-in-sequence-is-true-or-false)
                if (flag & 1<<2) == 0:

                    # Keep alignment only if read is primary aligned, i.e. that
                    #   - secondary alignment bit is unset
                    #   - supplementary alignment bit is unset
                    if ((flag & 1<<8) == 0) and ((flag & 1<<11) == 0):

                        # Retrieve the target the read was primarily aligned to
                        # ...and store in dict
                        alignment_dict[read_id] = regions_dict[target]

                        # Do the same for split reads
                        if read_id.endswith(split_suffixes):
                            # Only if assignment has not already been done for
                            # the read (we keep only the first aligned seq res-
                            # ulting from the split)
                            raw_read = read_id[:-2]
                            if not raw_read in alignment_dict.keys():
                                alignment_dict[raw_read] = regions_dict[target]
                
                # If read is unmapped, specify it in the dictionnary
                else:
                    alignment_dict[read_id] = 'Unmapped'

                line = sam_file.readline()

        # Loop through trim states and build each data frame
        for trim_state in all_trim_states:

            # Initialize dataframe, using categorical and numeric dtypes where
            # possible
            utils.send_text('compute_sequence_stats.py: Building data frame '
                            + 'for trim state ' + trim_state, v, 1, 2)
            species_cat = pd.CategoricalDtype(
                categories=species_names+['Trimmed', 'Unmapped']
            )
            qual_df = pd.DataFrame(columns=[
                'species', 'size', 'qual'
            ]).astype({
                'species': species_cat,
                'size': 'int64',
                'qual': 'float64'
            })

            # Determine path to fastq file
            if trim_state == 'untrimmed':
                fname = args.exp_dir.rstrip('/') \
                        + '/raw_fastq/' \
                        + sample \
                        + '.fastq'
            else:
                fname = args.exp_dir.rstrip('/') \
                        + '/trimmed_fastq/' \
                        + sample \
                        + '_porechopped_' \
                        + trim_state \
                        + '.fastq'
            
            utils.send_text('compute_sequence_stats.py: Parsing '
                            + fname.split('/')[-1], v, 2, 2)

            # Store reads statistics inside a dictionnary
            reads_stats = {}

            # Parse fastq file
            with open(fname, 'r') as fastq_file:

                line = fastq_file.readline()
                while line:

                    # Retrieve read id
                    read_id = line.strip().split()[0][1:]

                    line = fastq_file.readline() # skip header
                    line = fastq_file.readline() # skip sequence
                    line = fastq_file.readline() # skip '+' line

                    # Read ASCII Phred score and convert to list of int
                    score = line.strip()
                    phred = list(score.encode().translate(q_mapping))

                    # Add read to dictionnary
                    reads_stats[read_id] = [
                        alignment_dict.get(read_id, 'Trimmed'),
                        len(score),
                        sum(phred)/len(phred)
                    ]
                    
                    # Go to next read
                    line = fastq_file.readline()
            
            qual_df = pd.DataFrame.from_dict(
                reads_stats,
                orient='index',
                columns=['species', 'size', 'qual']
            ).astype({
                'species': species_cat,
                'size': 'int64',
                'qual': 'float64'
            })
            
            # Save data frame
            stored_df_path = args.output.rstrip('/') + '/' + sample + '_' \
                             + trim_state + '.pkl'
            qual_df.to_pickle(stored_df_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
