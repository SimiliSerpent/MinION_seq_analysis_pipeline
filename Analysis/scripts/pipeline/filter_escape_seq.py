#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Filter out the escape sequences from text files.
Takes as input a txt file and remove all the escape sequences. Write the output
to the specified output txt, or in place (default).
"""


import sys

import argparse
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, required=True,
                        help='Path to the txt file to filter.')
    parser.add_argument('-o', '--output_path', type=str,
                        help='Path to the filtered output txt file (default: '
                        + 'same as input = edit file in place).')

    args = parser.parse_args()

    # Read input file
    with open(args.input_path, 'r') as file:
        lines = file.readlines()
    
    # Filter escape strings
    new_lines = [
        re.sub(r'\r', r'\n', re.sub(
            r'\x1B[@-_][0-?]*[ -/]*[@-~]', '', str(line))
        ) for line in lines
    ]
    new_lines = [
        line for line in new_lines if not (
            ' / ' in line and ' (' in line and '%)\n' in line
        )
    ]
    
    # Write output file
    out_path = args.input_path
    if args.output_path is not None:
        out_path = args.output_path
    with open(out_path, 'w') as file:
        file.write(''.join(new_lines))

    return 0


if __name__ == "__main__":
    sys.exit(main())
