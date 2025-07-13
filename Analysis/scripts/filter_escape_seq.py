#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Filter out the escape sequences from text files.
Takes as input a multifasta file, loop through its lines and count the size of
each contig. Write it in a txt file with following formatting:
contig_name\tstart\tstop
"""


import sys

import argparse
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, required=True,
                        help='Path to the txt file to filter.')

    args = parser.parse_args()

    with open(args.file, 'r') as file:
        lines = file.readlines()
    
    # for line in lines[:100]:
    #     print('O_O' + re.sub(r'\x1B[@-_][0-?]*[ -/]*[@-~]', '', line))
    
    # return 0
    
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
    
    with open(args.file, 'w') as file:
        file.write(''.join(new_lines))

    return 0


if __name__ == "__main__":
    sys.exit(main())
