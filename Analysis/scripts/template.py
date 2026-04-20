#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This does that.              # short description (one-line)
With a bit more details,
you can explain here that this  # long description (multi-line)
does that in one very odd way,
taking this specific input and
outputing this uncommon output
with this significant special
feature.
# 34567891123456789212345678931234567894123456789512345678961234567897123456789
"""


import os                       # standard lib import (one per line!)
import sys

import argparse                 # 3rd party lib import
import numpy as np

import utils                    # local modules import


C_TYPE = 'donkey'               # define variables
FONTSIZE = 7


def my_func(arg1, /, arg2, *, arg3):        # define local functions
    """This does good.          # short description (one-line)
    And I tell you why in       # long description (multi-line)
    these very lines.

    Positionnal arguments:
    arg1 (str) - this means that

    Keyword arguments:
    arg2 (type) - arg desc
    arg3 (type) - arg desc
    

    Return:
    a    (type) - result desc
    """
    # Keep code commented as much as possible
    a = arg1 + 'does something'
            
    return(a)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Arg desc.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Arg desc.')
    parser.add_argument('-s', '--something', type=str, nargs='+',
                        help='Arg desc.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Do wonders
    b = some_lib.some_func(args.i)

    utils.send_text('Display something', v, 2, 0)

    return 0


if __name__ == "__main__":
    sys.exit(main())
